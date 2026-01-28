/**
 *  @file MultiSpeciesThermo.cpp
 *  Declarations for a thermodynamic property manager for multiple species
 *  in a phase (see @ref spthermo and
 * @link Cantera::MultiSpeciesThermo MultiSpeciesThermo@endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/MultiSpeciesThermo.h"
#include "cantera/thermo/SpeciesThermoFactory.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/utilities.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/global.h"

namespace Cantera
{

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
    if (m_p0 == 0) {
        // First species added; use this to set the reference pressure
        m_p0 = stit_ptr->refPressure();
    } else if (fabs(m_p0 - stit_ptr->refPressure()) > 1e-6) {
        throw CanteraError("MultiSpeciesThermo::install_STIT",
            "Cannot add species {} with reference pressure {}.\n"
            "Inconsistent with previously-added species with reference pressure {}.",
            index, stit_ptr->refPressure(), m_p0);
    }
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

void MultiSpeciesThermo::update_single(size_t k, double t, double& cp_R,
                                       double& h_RT, double& s_R) const
{
    const SpeciesThermoInterpType* sp_ptr = provideSTIT(k);
    if (sp_ptr) {
        sp_ptr->updatePropertiesTemp(t, cp_R, h_RT, s_R);
    }
}

void MultiSpeciesThermo::update(double t, span<double> cp_R, span<double> h_RT,
                                span<double> s_R) const
{
    auto iter = m_sp.begin();
    auto jter = m_tpoly.begin();
    for (; iter != m_sp.end(); iter++, jter++) {
        const vector<index_STIT>& species = iter->second;
        span<double> tpoly(jter->second);
        species[0].second->updateTemperaturePoly(t, tpoly);
        for (auto& [i, spthermo] : species) {
            spthermo->updateProperties(tpoly, cp_R[i], h_RT[i], s_R[i]);
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

void MultiSpeciesThermo::reportParams(size_t index, int& type, span<double> c,
        double& minTemp_, double& maxTemp_, double& refPressure_) const
{
    const SpeciesThermoInterpType* sp = provideSTIT(index);
    size_t n;
    if (sp) {
        sp->reportParameters(n, type, minTemp_, maxTemp_,
                             refPressure_, c.first(sp->nCoeffs()));
    } else {
        type = -1;
    }
}

double MultiSpeciesThermo::minTemp(size_t k) const
{
    if (k != npos) {
        const SpeciesThermoInterpType* sp = provideSTIT(k);
        if (sp) {
            return sp->minTemp();
        }
    }
    return m_tlow_max;
}

double MultiSpeciesThermo::maxTemp(size_t k) const
{
    if (k != npos) {
        const SpeciesThermoInterpType* sp = provideSTIT(k);
        if (sp) {
            return sp->maxTemp();
        }
    }
    return m_thigh_min;
}

double MultiSpeciesThermo::refPressure() const
{
    return m_p0;
}

SpeciesThermoInterpType* MultiSpeciesThermo::provideSTIT(size_t k)
{
    try {
        auto& [iParam, jSpecies] = m_speciesLoc.at(k);
        return m_sp.at(iParam)[jSpecies].second.get();
    } catch (std::out_of_range&) {
        return 0;
    }
}

const SpeciesThermoInterpType* MultiSpeciesThermo::provideSTIT(size_t k) const
{
    try {
        auto& [iParam, jSpecies] = m_speciesLoc.at(k);
        return m_sp.at(iParam)[jSpecies].second.get();
    } catch (std::out_of_range&) {
        return 0;
    }
}

double MultiSpeciesThermo::reportOneHf298(const size_t k) const
{
    const SpeciesThermoInterpType* sp_ptr = provideSTIT(k);
    double h = -1.0;
    if (sp_ptr) {
        h = sp_ptr->reportHf298();
    }
    return h;
}

void MultiSpeciesThermo::modifyOneHf298(const size_t k, const double Hf298New)
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
