/**
 *  @file GeneralSpeciesThermo.cpp
 *  Declarations for a completely general species thermodynamic property
 *  manager for a phase (see \ref spthermo and
 * \link Cantera::GeneralSpeciesThermo GeneralSpeciesThermo\endlink).
 */
// Copyright 2001-2004  California Institute of Technology

#include "cantera/thermo/GeneralSpeciesThermo.h"
#include "cantera/thermo/NasaPoly1.h"
#include "NasaPoly2.h"
#include "ShomatePoly.h"
#include "ConstCpPoly.h"
#include "cantera/thermo/Mu0Poly.h"
#include "cantera/thermo/AdsorbateThermo.h"
#include "cantera/thermo/SpeciesThermoFactory.h"
#include <iostream>

using namespace std;

namespace Cantera
{
GeneralSpeciesThermo::GeneralSpeciesThermo() :
    SpeciesThermo(),
    m_tlow_max(0.0),
    m_thigh_min(1.0E30),
    m_p0(OneAtm),
    m_kk(0)
{
    m_tlow_max = 0.0;
    m_thigh_min = 1.0E30;
}

GeneralSpeciesThermo::
GeneralSpeciesThermo(const GeneralSpeciesThermo& b) :
    m_tlow_max(b.m_tlow_max),
    m_thigh_min(b.m_thigh_min),
    m_kk(b.m_kk)
{
    m_sp.resize(m_kk, 0);
    for (size_t k = 0; k < m_kk; k++) {
        SpeciesThermoInterpType* bk = b.m_sp[k];
        if (bk) {
            m_sp[k] = bk->duplMyselfAsSpeciesThermoInterpType();
        }
    }
}

GeneralSpeciesThermo&
GeneralSpeciesThermo::operator=(const GeneralSpeciesThermo& b)
{
    if (&b != this) {
        m_tlow_max = b.m_tlow_max;
        m_thigh_min = b.m_thigh_min;

        for (size_t k = 0; k < m_kk; k++) {
            SpeciesThermoInterpType* sp = m_sp[k];
            if (sp) {
                delete sp;
                m_sp[k] = 0;
            }
        }
        m_kk = b.m_kk;
        m_sp.resize(m_kk, 0);
        for (size_t k = 0; k < m_kk; k++) {
            SpeciesThermoInterpType* bk = b.m_sp[k];
            if (bk) {
                m_sp[k] = bk->duplMyselfAsSpeciesThermoInterpType();
            }
        }
    }
    return *this;
}

GeneralSpeciesThermo::~GeneralSpeciesThermo()
{
    for (size_t k = 0; k < m_kk; k++) {
        SpeciesThermoInterpType* sp = m_sp[k];
        if (sp) {
            delete sp;
            m_sp[k] = 0;
        }
    }
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
    /*
     * Resize the arrays if necessary, filling the empty
     * slots with the zero pointer.
     */

    if (minTemp_ <= 0.0) {
        throw CanteraError("Error in GeneralSpeciesThermo.cpp",
                           " Cannot take 0 tmin as input. \n\n");
    }

    if (index >= m_kk) {
        m_sp.resize(index+1, 0);
        m_kk = index+1;
    }

    //int nfreq = 3;
    /*
     * Create the necessary object
     */

    switch (type) {
    case NASA1:
        m_sp[index] = new NasaPoly1(index, minTemp_, maxTemp_,
                                    refPressure_, c);
        break;
    case SHOMATE1:
        m_sp[index] = new ShomatePoly(index, minTemp_, maxTemp_,
                                      refPressure_, c);
        break;
    case CONSTANT_CP:
    case SIMPLE:
        m_sp[index] = new ConstCpPoly(index, minTemp_, maxTemp_,
                                      refPressure_, c);
        break;
    case MU0_INTERP:
        m_sp[index] = new Mu0Poly(index, minTemp_, maxTemp_,
                                  refPressure_, c);
        break;
    case SHOMATE2:
        m_sp[index] = new ShomatePoly2(index, minTemp_, maxTemp_,
                                       refPressure_, c);
        break;
    case NASA2:
        m_sp[index] = new NasaPoly2(index, minTemp_, maxTemp_,
                                    refPressure_, c);
        break;

    case STAT:
        m_sp[index] = new StatMech(index, minTemp_, maxTemp_,
                                   refPressure_, c, name);
        break;

    case ADSORBATE:
        m_sp[index] = new Adsorbate(index, minTemp_, maxTemp_,
                                    refPressure_, c);
        break;
    default:
        throw UnknownSpeciesThermoModel(
            "GeneralSpeciesThermo::install",
            "unknown species type", int2str(type));
        break;
    }
    if (!m_sp[index]) {
        cout << "Null m_sp... index = " << index << endl;
        cout << "type = " << type << endl;
    }
    m_tlow_max = max(minTemp_, m_tlow_max);
    m_thigh_min = min(maxTemp_, m_thigh_min);
}

void GeneralSpeciesThermo::install_STIT(SpeciesThermoInterpType* stit_ptr)
{
    /*
     * Resize the arrays if necessary, filling the empty
     * slots with the zero pointer.
     */
    if (!stit_ptr) {
        throw CanteraError("GeneralSpeciesThermo::install_STIT",
                           "zero pointer");
    }
    size_t index = stit_ptr->speciesIndex();
    if (index >= m_kk) {
        m_sp.resize(index+1, 0);
        m_kk = index+1;
    }
    AssertThrow(m_sp[index] == 0,
                "Index position isn't null, duplication of assignment: " + int2str(index));
    /*
     *  Now, simply assign the position
     */
    m_sp[index] = stit_ptr;

    /*
     * Calculate max and min
     */
    m_tlow_max = max(stit_ptr->minTemp(), m_tlow_max);
    m_thigh_min = min(stit_ptr->maxTemp(), m_thigh_min);
}

void GeneralSpeciesThermo::installPDSShandler(size_t k, PDSS* PDSS_ptr,
        VPSSMgr* vpssmgr_ptr)
{
    STITbyPDSS* stit_ptr = new STITbyPDSS(k, vpssmgr_ptr, PDSS_ptr);
    install_STIT(stit_ptr);
}

void GeneralSpeciesThermo::
update_one(size_t k, doublereal t, doublereal* cp_R,
           doublereal* h_RT, doublereal* s_R) const
{
    SpeciesThermoInterpType* sp_ptr = m_sp[k];
    if (sp_ptr) {
        sp_ptr->updatePropertiesTemp(t, cp_R, h_RT, s_R);
    }
}

void GeneralSpeciesThermo::
update(doublereal t, doublereal* cp_R,
       doublereal* h_RT, doublereal* s_R) const
{
    vector<SpeciesThermoInterpType*>::const_iterator _begin, _end;
    _begin  = m_sp.begin();
    _end    = m_sp.end();
    SpeciesThermoInterpType* sp_ptr = 0;
    for (; _begin != _end; ++_begin) {
        sp_ptr = *(_begin);
        if (sp_ptr) {
            sp_ptr->updatePropertiesTemp(t, cp_R, h_RT, s_R);
        }
        // else {
        //   writelog("General::update: sp_ptr is NULL!\n");
        //}
    }
}

int GeneralSpeciesThermo::reportType(size_t index) const
{
    SpeciesThermoInterpType* sp = m_sp[index];
    if (sp) {
        return sp->reportType();
    }
    return -1;
}

void GeneralSpeciesThermo::
reportParams(size_t index, int& type, doublereal* const c,
             doublereal& minTemp_, doublereal& maxTemp_, doublereal& refPressure_) const
{
    warn_deprecated("GeneralSpeciesThermo::reportParams");
    SpeciesThermoInterpType* sp = m_sp[index];
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
    if (k == npos) {
        return m_tlow_max;
    } else {
        SpeciesThermoInterpType* sp = m_sp[k];
        if (sp) {
            return sp->minTemp();
        }
    }
    return m_tlow_max;
}

doublereal GeneralSpeciesThermo::maxTemp(size_t k) const
{
    if (k == npos) {
        return m_thigh_min;
    } else {
        SpeciesThermoInterpType* sp = m_sp[k];
        if (sp) {
            return sp->maxTemp();
        }
    }
    return m_thigh_min;
}

doublereal GeneralSpeciesThermo::refPressure(size_t k) const
{
    if (k == npos) {
        return m_p0;
    } else {
        SpeciesThermoInterpType* sp = m_sp[k];
        if (sp) {
            return sp->refPressure();
        }
    }
    return m_p0;
}

SpeciesThermoInterpType* GeneralSpeciesThermo::provideSTIT(size_t k)
{
    return m_sp[k];
}

#ifdef H298MODIFY_CAPABILITY

doublereal GeneralSpeciesThermo::reportOneHf298(int k) const
{
    SpeciesThermoInterpType* sp_ptr = m_sp[k];
    doublereal h = -1.0;
    if (sp_ptr) {
        h = sp_ptr->reportHf298(0);
    }
    return h;
}

void GeneralSpeciesThermo::modifyOneHf298(const int k, const doublereal Hf298New)
{
    SpeciesThermoInterpType* sp_ptr = m_sp[k];
    if (sp_ptr) {
        sp_ptr->modifyOneHf298(k, Hf298New);
    }
}

#endif

}
