/**
 *  @file MultiSpeciesThermo.cpp
 *  Declarations for a thermodynamic property manager for multiple species
 *  in a phase (see \ref spthermo and
 * \link Cantera::MultiSpeciesThermo MultiSpeciesThermo\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/MultiSpeciesThermo.h"
#include "cantera/thermo/SpeciesThermoFactory.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/utilities.h"
#include "cantera/base/ctexceptions.h"

namespace Cantera
{
MultiSpeciesThermo::MultiSpeciesThermo() :
    m_tlow_max(0.0),
    m_thigh_min(1.0E30),
    m_p0(OneAtm)
{
}

void MultiSpeciesThermo::install_STIT(size_t index,
                                        shared_ptr<SpeciesThermoInterpType> stit_ptr)
{
    if (!stit_ptr) {
        throw CanteraError("MultiSpeciesThermo::install_STIT",
                           "null pointer");
    }
    AssertThrowMsg(m_speciesLoc.find(index) == m_speciesLoc.end(),
            "MultiSpeciesThermo::install_STIT",
            "Index position isn't null, duplication of assignment: {}", index);
    int type = stit_ptr->reportType();
    m_speciesLoc[index] = {type, m_sp[type].size()};
    m_sp[type].emplace_back(index, stit_ptr);
    if (m_sp[type].size() == 1) {
        m_tpoly[type].resize(stit_ptr->temperaturePolySize());
    }

    // Calculate max and min T
    m_tlow_max = std::max(stit_ptr->minTemp(), m_tlow_max);
    m_thigh_min = std::min(stit_ptr->maxTemp(), m_thigh_min);
    markInstalled(index);
}

void MultiSpeciesThermo::modifySpecies(size_t index,
                                         shared_ptr<SpeciesThermoInterpType> spthermo)
{
    if (!spthermo) {
        throw CanteraError("MultiSpeciesThermo::modifySpecies",
                           "null pointer");
    }
    if (m_speciesLoc.find(index) == m_speciesLoc.end()) {
        throw CanteraError("MultiSpeciesThermo::modifySpecies",
                           "Species with this index not previously added: {}",
                           index);
    }
    int type = spthermo->reportType();
    if (m_speciesLoc[index].first != type) {
        throw CanteraError("MultiSpeciesThermo::modifySpecies",
                           "Type of parameterization changed: {} != {}", type,
                           m_speciesLoc[index].first);
    }
    if (spthermo->minTemp() > m_tlow_max) {
        throw CanteraError("MultiSpeciesThermo::modifySpecies",
            "Cannot increase minimum temperature for phase from {} to {}",
            m_tlow_max, spthermo->minTemp());
    }
    if (spthermo->maxTemp() < m_thigh_min) {
        throw CanteraError("MultiSpeciesThermo::modifySpecies",
            "Cannot increase minimum temperature for phase from {} to {}",
            m_thigh_min, spthermo->maxTemp());
    }

    m_sp[type][m_speciesLoc[index].second] = {index, spthermo};
}

void MultiSpeciesThermo::update_single(size_t k, double t, double* cp_R,
                                       double* h_RT, double* s_R) const
{
    const SpeciesThermoInterpType* sp_ptr = provideSTIT(k);
    if (sp_ptr) {
        sp_ptr->updatePropertiesTemp(t, cp_R, h_RT, s_R);
    }
}

void MultiSpeciesThermo::update(doublereal t, doublereal* cp_R,
                                  doublereal* h_RT, doublereal* s_R) const
{
    auto iter = m_sp.begin();
    auto jter = m_tpoly.begin();
    for (; iter != m_sp.end(); iter++, jter++) {
        const std::vector<index_STIT>& species = iter->second;
        double* tpoly = &jter->second[0];
        species[0].second->updateTemperaturePoly(t, tpoly);
        for (size_t k = 0; k < species.size(); k++) {
            size_t i = species[k].first;
            species[k].second->updateProperties(tpoly, cp_R+i, h_RT+i, s_R+i);
        }
    }
}

int MultiSpeciesThermo::reportType(size_t index) const
{
    const SpeciesThermoInterpType* sp = provideSTIT(index);
    if (sp) {
        return sp->reportType();
    }
    return -1;
}

void MultiSpeciesThermo::reportParams(size_t index, int& type,
        doublereal* const c, doublereal& minTemp_, doublereal& maxTemp_,
        doublereal& refPressure_) const
{
    const SpeciesThermoInterpType* sp = provideSTIT(index);
    size_t n;
    if (sp) {
        sp->reportParameters(n, type, minTemp_, maxTemp_,
                             refPressure_, c);
    } else {
        type = -1;
    }
}

doublereal MultiSpeciesThermo::minTemp(size_t k) const
{
    if (k != npos) {
        const SpeciesThermoInterpType* sp = provideSTIT(k);
        if (sp) {
            return sp->minTemp();
        }
    }
    return m_tlow_max;
}

doublereal MultiSpeciesThermo::maxTemp(size_t k) const
{
    if (k != npos) {
        const SpeciesThermoInterpType* sp = provideSTIT(k);
        if (sp) {
            return sp->maxTemp();
        }
    }
    return m_thigh_min;
}

doublereal MultiSpeciesThermo::refPressure(size_t k) const
{
    if (k != npos) {
        const SpeciesThermoInterpType* sp = provideSTIT(k);
        if (sp) {
            return sp->refPressure();
        }
    }
    return m_p0;
}

SpeciesThermoInterpType* MultiSpeciesThermo::provideSTIT(size_t k)
{
    try {
        const std::pair<int, size_t>& loc = m_speciesLoc.at(k);
        return m_sp.at(loc.first)[loc.second].second.get();
    } catch (std::out_of_range&) {
        return 0;
    }
}

const SpeciesThermoInterpType* MultiSpeciesThermo::provideSTIT(size_t k) const
{
    try {
        const std::pair<int, size_t>& loc = m_speciesLoc.at(k);
        return m_sp.at(loc.first)[loc.second].second.get();
    } catch (std::out_of_range&) {
        return 0;
    }
}

doublereal MultiSpeciesThermo::reportOneHf298(const size_t k) const
{
    const SpeciesThermoInterpType* sp_ptr = provideSTIT(k);
    doublereal h = -1.0;
    if (sp_ptr) {
        h = sp_ptr->reportHf298(0);
    }
    return h;
}

void MultiSpeciesThermo::modifyOneHf298(const size_t k, const doublereal Hf298New)
{
    SpeciesThermoInterpType* sp_ptr = provideSTIT(k);
    if (sp_ptr) {
        sp_ptr->modifyOneHf298(k, Hf298New);
    }
}

void MultiSpeciesThermo::resetHf298(const size_t k)
{
    SpeciesThermoInterpType* sp_ptr = provideSTIT(k);
    if (sp_ptr) {
        sp_ptr->resetHf298();
    }
}

bool MultiSpeciesThermo::ready(size_t nSpecies) {
    if (m_installed.size() < nSpecies) {
        return false;
    }
    for (size_t k = 0; k < nSpecies; k++) {
        if (!m_installed[k]) {
            return false;
        }
    }
    return true;
}

void MultiSpeciesThermo::markInstalled(size_t k) {
    if (k >= m_installed.size()) {
        m_installed.resize(k+1, false);
    }
    m_installed[k] = true;
}

}
