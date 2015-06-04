/**
 *  @file GeneralSpeciesThermo.cpp
 *  Declarations for a completely general species thermodynamic property
 *  manager for a phase (see \ref spthermo and
 * \link Cantera::GeneralSpeciesThermo GeneralSpeciesThermo\endlink).
 */
// Copyright 2001-2004  California Institute of Technology

#include "cantera/thermo/GeneralSpeciesThermo.h"
#include "cantera/thermo/SpeciesThermoFactory.h"

namespace Cantera
{
GeneralSpeciesThermo::GeneralSpeciesThermo() :
    m_tlow_max(0.0),
    m_thigh_min(1.0E30),
    m_p0(OneAtm)
{
}

GeneralSpeciesThermo::GeneralSpeciesThermo(const GeneralSpeciesThermo& b) :
    SpeciesThermo(b),
    m_tpoly(b.m_tpoly),
    m_speciesLoc(b.m_speciesLoc),
    m_tlow_max(b.m_tlow_max),
    m_thigh_min(b.m_thigh_min),
    m_p0(b.m_p0)
{
    m_sp.clear();
    // Copy SpeciesThermoInterpTypes from 'b'
    for (STIT_map::const_iterator iter = b.m_sp.begin();
         iter != b.m_sp.end();
         iter++) {
        for (size_t k = 0; k < iter->second.size(); k++) {
            size_t i = iter->second[k].first;
            shared_ptr<SpeciesThermoInterpType> spec(
                iter->second[k].second->duplMyselfAsSpeciesThermoInterpType());
            m_sp[iter->first].push_back(std::make_pair(i, spec));
        }
    }
}

GeneralSpeciesThermo&
GeneralSpeciesThermo::operator=(const GeneralSpeciesThermo& b)
{
    if (&b == this) {
        return *this;
    }

    SpeciesThermo::operator=(b);
    m_sp.clear();
    // Copy SpeciesThermoInterpType objects from 'b'
    for (STIT_map::const_iterator iter = b.m_sp.begin();
         iter != b.m_sp.end();
         iter++) {
        for (size_t k = 0; k < iter->second.size(); k++) {
            size_t i = iter->second[k].first;
            shared_ptr<SpeciesThermoInterpType> spec(
                iter->second[k].second->duplMyselfAsSpeciesThermoInterpType());
            m_sp[iter->first].push_back(std::make_pair(i, spec));
        }
    }

    m_tpoly = b.m_tpoly;
    m_speciesLoc = b.m_speciesLoc;
    m_tlow_max = b.m_tlow_max;
    m_thigh_min = b.m_thigh_min;
    m_p0 = b.m_p0;

    return *this;
}

SpeciesThermo*
GeneralSpeciesThermo::duplMyselfAsSpeciesThermo() const
{
    return new GeneralSpeciesThermo(*this);
}

void GeneralSpeciesThermo::install(const std::string& name,
                                   size_t index,
                                   int type,
                                   const doublereal* c,
                                   doublereal minTemp_,
                                   doublereal maxTemp_,
                                   doublereal refPressure_)
{
    warn_deprecated("GeneralSpeciesThermo::install",
            "Use newSpeciesThermoInterpType and "
            "GeneralSpeciesThermo::install_STIT instead");
    if (minTemp_ <= 0.0) {
        throw CanteraError("GeneralSpeciesThermo::install",
                           "T_min must be positive");
    }

    /*
     * Create the necessary object
     */
    shared_ptr<SpeciesThermoInterpType> sp(newSpeciesThermoInterpType(type,
        minTemp_, maxTemp_, refPressure_, c));
    sp->validate(name);
    install_STIT(index, sp);
}

void GeneralSpeciesThermo::install_STIT(size_t index,
                                        shared_ptr<SpeciesThermoInterpType> stit_ptr)
{
    if (!stit_ptr) {
        throw CanteraError("GeneralSpeciesThermo::install_STIT",
                           "zero pointer");
    }
    AssertThrow(m_speciesLoc.find(index) == m_speciesLoc.end(),
                "Index position isn't null, duplication of assignment: " + int2str(index));
    int type = stit_ptr->reportType();
    m_speciesLoc[index] = std::make_pair(type, m_sp[type].size());
    m_sp[type].push_back(std::make_pair(index, stit_ptr));
    if (m_sp[type].size() == 1) {
        m_tpoly[type].resize(stit_ptr->temperaturePolySize());
    }

    // Calculate max and min T
    m_tlow_max = std::max(stit_ptr->minTemp(), m_tlow_max);
    m_thigh_min = std::min(stit_ptr->maxTemp(), m_thigh_min);
    markInstalled(index);
}

void GeneralSpeciesThermo::installPDSShandler(size_t k, PDSS* PDSS_ptr,
        VPSSMgr* vpssmgr_ptr)
{
    shared_ptr<SpeciesThermoInterpType> stit_ptr(new STITbyPDSS(k, vpssmgr_ptr, PDSS_ptr));
    install_STIT(k, stit_ptr);
}

void GeneralSpeciesThermo::update_one(size_t k, doublereal t, doublereal* cp_R,
                                      doublereal* h_RT, doublereal* s_R) const
{
    const SpeciesThermoInterpType* sp_ptr = provideSTIT(k);
    if (sp_ptr) {
        sp_ptr->updatePropertiesTemp(t, cp_R+k, h_RT+k, s_R+k);
    }
}

void GeneralSpeciesThermo::update(doublereal t, doublereal* cp_R,
                                  doublereal* h_RT, doublereal* s_R) const
{
    STIT_map::const_iterator iter = m_sp.begin();
    tpoly_map::iterator jter = m_tpoly.begin();
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

int GeneralSpeciesThermo::reportType(size_t index) const
{
    const SpeciesThermoInterpType* sp = provideSTIT(index);
    if (sp) {
        return sp->reportType();
    }
    return -1;
}

void GeneralSpeciesThermo::reportParams(size_t index, int& type,
        doublereal* const c, doublereal& minTemp_, doublereal& maxTemp_,
        doublereal& refPressure_) const
{
    const SpeciesThermoInterpType* sp = provideSTIT(index);
    size_t n;
    if (sp) {
        sp->reportParameters(n, type, minTemp_, maxTemp_,
                             refPressure_, c);
        if (n != index) {
            throw CanteraError("GeneralSpeciesThermo::reportParams",
                               "Internal error encountered");
        }
    } else {
        type = -1;
    }
}

doublereal GeneralSpeciesThermo::minTemp(size_t k) const
{
    if (k != npos) {
        const SpeciesThermoInterpType* sp = provideSTIT(k);
        if (sp) {
            return sp->minTemp();
        }
    }
    return m_tlow_max;
}

doublereal GeneralSpeciesThermo::maxTemp(size_t k) const
{
    if (k != npos) {
        const SpeciesThermoInterpType* sp = provideSTIT(k);
        if (sp) {
            return sp->maxTemp();
        }
    }
    return m_thigh_min;
}

doublereal GeneralSpeciesThermo::refPressure(size_t k) const
{
    if (k != npos) {
        const SpeciesThermoInterpType* sp = provideSTIT(k);
        if (sp) {
            return sp->refPressure();
        }
    }
    return m_p0;
}

SpeciesThermoInterpType* GeneralSpeciesThermo::provideSTIT(size_t k)
{
    try {
        const std::pair<int, size_t>& loc = getValue(m_speciesLoc, k);
        return getValue(m_sp, loc.first)[loc.second].second.get();
    } catch (std::out_of_range&) {
        return 0;
    }
}

const SpeciesThermoInterpType* GeneralSpeciesThermo::provideSTIT(size_t k) const
{
    try {
        const std::pair<int, size_t>& loc = getValue(m_speciesLoc, k);
        return getValue(m_sp, loc.first)[loc.second].second.get();
    } catch (std::out_of_range&) {
        return 0;
    }
}


doublereal GeneralSpeciesThermo::reportOneHf298(const size_t k) const
{
    const SpeciesThermoInterpType* sp_ptr = provideSTIT(k);
    doublereal h = -1.0;
    if (sp_ptr) {
        h = sp_ptr->reportHf298(0);
    }
    return h;
}

void GeneralSpeciesThermo::modifyOneHf298(const size_t k, const doublereal Hf298New)
{
    SpeciesThermoInterpType* sp_ptr = provideSTIT(k);
    if (sp_ptr) {
        sp_ptr->modifyOneHf298(k, Hf298New);
    }
}

}
