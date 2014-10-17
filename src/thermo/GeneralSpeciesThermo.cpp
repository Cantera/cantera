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
    SpeciesThermo(),
    m_tlow_max(0.0),
    m_thigh_min(1.0E30),
    m_p0(OneAtm)
{
}

GeneralSpeciesThermo::GeneralSpeciesThermo(const GeneralSpeciesThermo& b) :
    m_tpoly(b.m_tpoly),
    m_tlow_max(b.m_tlow_max),
    m_thigh_min(b.m_thigh_min),
    m_p0(b.m_p0)
{
    clear();
    // Copy SpeciesThermoInterpTypes from 'b'
    for (STIT_map::const_iterator iter = b.m_sp.begin();
         iter != b.m_sp.end();
         iter++) {
        for (size_t k = 0; k < iter->second.size(); k++) {
            const SpeciesThermoInterpType* spec = iter->second[k];
            m_sp[iter->first].push_back(spec->duplMyselfAsSpeciesThermoInterpType());
        }
    }
}

GeneralSpeciesThermo&
GeneralSpeciesThermo::operator=(const GeneralSpeciesThermo& b)
{
    if (&b == this) {
        return *this;
    }
    clear();
    // Copy SpeciesThermoInterpType objects from 'b'
    for (STIT_map::const_iterator iter = b.m_sp.begin();
         iter != b.m_sp.end();
         iter++) {
        for (size_t k = 0; k < iter->second.size(); k++) {
            const SpeciesThermoInterpType* spec = iter->second[k];
            m_sp[iter->first].push_back(spec->duplMyselfAsSpeciesThermoInterpType());
        }
    }

    m_tpoly = b.m_tpoly;
    m_tlow_max = b.m_tlow_max;
    m_thigh_min = b.m_thigh_min;
    m_p0 = b.m_p0;

    return *this;
}

GeneralSpeciesThermo::~GeneralSpeciesThermo()
{
    clear();
}

SpeciesThermo*
GeneralSpeciesThermo::duplMyselfAsSpeciesThermo() const
{
    return new GeneralSpeciesThermo(*this);
}

void GeneralSpeciesThermo::clear()
{
    for (STIT_map::const_iterator iter = m_sp.begin();
         iter != m_sp.end();
         iter++) {
        for (size_t k = 0; k < iter->second.size(); k++) {
            delete iter->second[k];
        }
    }
    m_sp.clear();
}

void GeneralSpeciesThermo::install(const std::string& name,
                                   size_t index,
                                   int type,
                                   const doublereal* c,
                                   doublereal minTemp_,
                                   doublereal maxTemp_,
                                   doublereal refPressure_)
{
    if (minTemp_ <= 0.0) {
        throw CanteraError("GeneralSpeciesThermo::install",
                           "T_min must be positive");
    }

    /*
     * Create the necessary object
     */
    SpeciesThermoInterpType* sp = newSpeciesThermoInterpType(type,
        minTemp_, maxTemp_, refPressure_, c);
    sp->setIndex(index);
    sp->validate(name);
    install_STIT(sp);
}

void GeneralSpeciesThermo::install_STIT(SpeciesThermoInterpType* stit_ptr)
{
    if (!stit_ptr) {
        throw CanteraError("GeneralSpeciesThermo::install_STIT",
                           "zero pointer");
    }
    size_t index = stit_ptr->speciesIndex();
    AssertThrow(m_speciesLoc.find(index) == m_speciesLoc.end(),
                "Index position isn't null, duplication of assignment: " + int2str(index));

    int type = stit_ptr->reportType();
    m_speciesLoc[index] = std::make_pair(type, m_sp[type].size());
    m_sp[type].push_back(stit_ptr);
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
    STITbyPDSS* stit_ptr = new STITbyPDSS(k, vpssmgr_ptr, PDSS_ptr);
    install_STIT(stit_ptr);
}

void GeneralSpeciesThermo::update_one(size_t k, doublereal t, doublereal* cp_R,
                                      doublereal* h_RT, doublereal* s_R) const
{
    const SpeciesThermoInterpType* sp_ptr = provideSTIT(k);
    if (sp_ptr) {
        sp_ptr->updatePropertiesTemp(t, cp_R, h_RT, s_R);
    }
}

void GeneralSpeciesThermo::update(doublereal t, doublereal* cp_R,
                                  doublereal* h_RT, doublereal* s_R) const
{
    STIT_map::const_iterator iter = m_sp.begin();
    tpoly_map::iterator jter = m_tpoly.begin();
    for (; iter != m_sp.end(); iter++, jter++) {
        const std::vector<SpeciesThermoInterpType*>& species = iter->second;
        double* tpoly = &jter->second[0];
        species[0]->updateTemperaturePoly(t, tpoly);
        for (size_t k = 0; k < species.size(); k++) {
            species[k]->updateProperties(tpoly, cp_R, h_RT, s_R);
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
        return getValue(m_sp, loc.first)[loc.second];
    } catch (std::out_of_range&) {
        return 0;
    }
}

const SpeciesThermoInterpType* GeneralSpeciesThermo::provideSTIT(size_t k) const
{
    try {
        const std::pair<int, size_t>& loc = getValue(m_speciesLoc, k);
        return getValue(m_sp, loc.first)[loc.second];
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
