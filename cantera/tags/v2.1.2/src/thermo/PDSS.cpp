/**
 * @file PDSS.cpp
 * Implementation of a pressure dependent standard state
 * virtual function
 * (see class \link Cantera::PDSS PDSS\endlink).
 */
/*
 * Copyright (2006) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
#include "cantera/base/ct_defs.h"
#include "cantera/base/xml.h"
#include "cantera/base/ctml.h"
#include "cantera/thermo/PDSS.h"

#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/SpeciesThermo.h"

#include "cantera/thermo/VPStandardStateTP.h"

namespace Cantera
{
PDSS::PDSS() :
    m_pdssType(cPDSS_UNDEF),
    m_temp(-1.0),
    m_pres(-1.0),
    m_p0(-1.0),
    m_minTemp(-1.0),
    m_maxTemp(10000.0),
    m_tp(0),
    m_vpssmgr_ptr(0),
    m_mw(0.0),
    m_spindex(npos),
    m_spthermo(0),
    m_h0_RT_ptr(0),
    m_cp0_R_ptr(0),
    m_s0_R_ptr(0),
    m_g0_RT_ptr(0),
    m_V0_ptr(0),
    m_hss_RT_ptr(0),
    m_cpss_R_ptr(0),
    m_sss_R_ptr(0),
    m_gss_RT_ptr(0),
    m_Vss_ptr(0)
{
}

PDSS::PDSS(VPStandardStateTP* tp, size_t spindex) :
    m_pdssType(cPDSS_UNDEF),
    m_temp(-1.0),
    m_pres(-1.0),
    m_p0(-1.0),
    m_minTemp(-1.0),
    m_maxTemp(10000.0),
    m_tp(tp),
    m_vpssmgr_ptr(0),
    m_mw(0.0),
    m_spindex(spindex),
    m_spthermo(0),
    m_h0_RT_ptr(0),
    m_cp0_R_ptr(0),
    m_s0_R_ptr(0),
    m_g0_RT_ptr(0),
    m_V0_ptr(0),
    m_hss_RT_ptr(0),
    m_cpss_R_ptr(0),
    m_sss_R_ptr(0),
    m_gss_RT_ptr(0),
    m_Vss_ptr(0)
{
    if (tp) {
        m_spthermo = &(tp->speciesThermo());
    }
    if (tp) {
        m_vpssmgr_ptr = tp->provideVPSSMgr();
    }
}

PDSS::PDSS(const PDSS& b) :
    m_pdssType(cPDSS_UNDEF),
    m_temp(-1.0),
    m_pres(-1.0),
    m_p0(-1.0),
    m_minTemp(-1.0),
    m_maxTemp(10000.0),
    m_tp(0),
    m_vpssmgr_ptr(0),
    m_mw(b.m_mw),
    m_spindex(b.m_spindex),
    m_spthermo(b.m_spthermo),
    m_h0_RT_ptr(b.m_h0_RT_ptr),
    m_cp0_R_ptr(b.m_cp0_R_ptr),
    m_s0_R_ptr(b.m_s0_R_ptr),
    m_g0_RT_ptr(b.m_g0_RT_ptr),
    m_V0_ptr(b.m_V0_ptr),
    m_hss_RT_ptr(b.m_hss_RT_ptr),
    m_cpss_R_ptr(b.m_cpss_R_ptr),
    m_sss_R_ptr(b.m_sss_R_ptr),
    m_gss_RT_ptr(b.m_gss_RT_ptr),
    m_Vss_ptr(b.m_Vss_ptr)
{
    /*
     * Use the assignment operator to do the brunt
     * of the work for the copy constructor.
     */
    *this = b;
}

PDSS& PDSS::operator=(const PDSS& b)
{
    if (&b == this) {
        return *this;
    }

    m_pdssType     = b.m_pdssType;
    m_temp         = b.m_temp;
    m_pres         = b.m_pres;
    m_p0           = b.m_p0;
    m_minTemp      = b.m_minTemp;
    m_maxTemp      = b.m_maxTemp;

    // Pointers which are zero, are properly assigned in the
    // function, initAllPtrs(). which must be called after the
    // assignment operation.

    m_tp           = 0;
    m_vpssmgr_ptr  = 0;
    m_mw           = b.m_mw;
    m_spindex      = b.m_spindex;
    m_spthermo     = 0;
    m_cp0_R_ptr    = 0;
    m_h0_RT_ptr    = 0;
    m_s0_R_ptr     = 0;
    m_g0_RT_ptr    = 0;
    m_V0_ptr       = 0;
    m_cpss_R_ptr   = 0;
    m_hss_RT_ptr   = 0;
    m_sss_R_ptr    = 0;
    m_gss_RT_ptr   = 0;
    m_Vss_ptr      = 0;

    // Here we just fill these in so that local copies within the VPSS object work.
    m_tp           = b.m_tp;
    m_vpssmgr_ptr  = b.m_vpssmgr_ptr;
    m_spthermo     = b.m_spthermo;
    m_cp0_R_ptr    = b.m_cp0_R_ptr;
    m_h0_RT_ptr    = b.m_h0_RT_ptr;
    m_s0_R_ptr     = b.m_s0_R_ptr;
    m_g0_RT_ptr    = b.m_g0_RT_ptr;
    m_V0_ptr       = b.m_V0_ptr;
    m_cpss_R_ptr   = b.m_cpss_R_ptr;
    m_hss_RT_ptr   = b.m_hss_RT_ptr;
    m_sss_R_ptr    = b.m_sss_R_ptr;
    m_gss_RT_ptr   = b.m_gss_RT_ptr;
    m_Vss_ptr      = b.m_Vss_ptr;

    return *this;
}

PDSS::~PDSS()
{
}

PDSS* PDSS::duplMyselfAsPDSS() const
{
    return new PDSS(*this);
}

PDSS_enumType PDSS::reportPDSSType() const
{
    return m_pdssType;
}

void PDSS::initThermoXML(const XML_Node& phaseNode, const std::string& id)
{
    AssertThrow(m_tp != 0, "PDSS::initThermoXML()");
    m_p0 =  m_vpssmgr_ptr->refPressure(m_spindex);
    m_minTemp = m_vpssmgr_ptr->minTemp(m_spindex);
    m_maxTemp = m_vpssmgr_ptr->maxTemp(m_spindex);
}

void PDSS::initThermo()
{
    AssertThrow(m_tp != 0, "PDSS::initThermo()");
    m_vpssmgr_ptr = m_tp->provideVPSSMgr();
    m_vpssmgr_ptr->initThermo();
    initPtrs();
    m_mw = m_tp->molecularWeight(m_spindex);
}

void PDSS::initAllPtrs(VPStandardStateTP* tp, VPSSMgr* vpssmgr_ptr,
                       SpeciesThermo* spthermo)
{
    m_tp = tp;
    m_vpssmgr_ptr = vpssmgr_ptr;
    m_spthermo = spthermo;
    initPtrs();
}

void PDSS::initPtrs()
{
    AssertThrow(m_vpssmgr_ptr->mPDSS_h0_RT.size() != 0, "PDSS::initPtrs()");
    m_h0_RT_ptr  = &(m_vpssmgr_ptr->mPDSS_h0_RT[0]);
    m_cp0_R_ptr  = &(m_vpssmgr_ptr->mPDSS_cp0_R[0]);
    m_s0_R_ptr   = &(m_vpssmgr_ptr->mPDSS_s0_R[0]);
    m_g0_RT_ptr  = &(m_vpssmgr_ptr->mPDSS_g0_RT[0]);
    m_V0_ptr     = &(m_vpssmgr_ptr->mPDSS_V0[0]);

    m_hss_RT_ptr  = &(m_vpssmgr_ptr->mPDSS_hss_RT[0]);
    m_cpss_R_ptr  = &(m_vpssmgr_ptr->mPDSS_cpss_R[0]);
    m_sss_R_ptr   = &(m_vpssmgr_ptr->mPDSS_sss_R[0]);
    m_gss_RT_ptr  = &(m_vpssmgr_ptr->mPDSS_gss_RT[0]);
    m_Vss_ptr     = &(m_vpssmgr_ptr->mPDSS_Vss[0]);
}

doublereal PDSS::enthalpy_mole() const
{
    err("enthalpy_mole()");
    return 0.0;
}

doublereal PDSS::enthalpy_RT() const
{
    double RT = GasConstant * m_temp;
    return enthalpy_mole()/RT;
}

doublereal PDSS::intEnergy_mole() const
{
    err("intEnergy_mole()");
    return 0.0;
}

doublereal PDSS::entropy_mole() const
{
    err("entropy_mole()");
    return 0.0;
}

doublereal PDSS::entropy_R() const
{
    return entropy_mole()/GasConstant;
}

doublereal PDSS::gibbs_mole() const
{
    err("gibbs_mole()");
    return 0.0;
}

doublereal PDSS::gibbs_RT() const
{
    double RT = GasConstant * m_temp;
    return gibbs_mole()/RT;
}

doublereal PDSS::cp_mole() const
{
    err("cp_mole()");
    return 0.0;
}

doublereal PDSS::cp_R() const
{
    return cp_mole()/GasConstant;
}

doublereal PDSS::molarVolume() const
{
    err("molarVolume()");
    return 0.0;
}

doublereal PDSS::density() const
{
    err("density()");
    return 0.0;
}

doublereal PDSS::cv_mole() const
{
    err("cv_mole()");
    return 0.0;
}

doublereal PDSS::gibbs_RT_ref() const
{
    err("gibbs_RT_ref()");
    return 0.0;
}

doublereal PDSS::enthalpy_RT_ref() const
{
    err("enthalpy_RT_ref()");
    return 0.0;
}

doublereal PDSS::entropy_R_ref() const
{
    err("entropy_RT_ref()");
    return 0.0;
}

doublereal PDSS::cp_R_ref() const
{
    err("entropy_RT_ref()");
    return 0.0;
}

doublereal PDSS::molarVolume_ref() const
{
    err("molarVolume_ref()");
    return 0.0;
}

doublereal PDSS::
enthalpyDelp_mole() const
{
    doublereal RT = m_temp * GasConstant;
    doublereal tmp = enthalpy_RT_ref();
    return enthalpy_mole() - RT * tmp;
}

doublereal PDSS::entropyDelp_mole() const
{
    doublereal tmp = entropy_R_ref();
    return entropy_mole() - GasConstant * tmp;

}

doublereal PDSS::gibbsDelp_mole() const
{
    doublereal RT = m_temp * GasConstant;
    doublereal tmp = gibbs_RT_ref();
    return gibbs_mole() - RT * tmp;
}

doublereal PDSS::cpDelp_mole() const
{
    doublereal tmp = cp_R_ref();
    return cp_mole() - GasConstant * tmp;
}

doublereal PDSS::pressure() const
{
    return m_pres;
}

doublereal PDSS::thermalExpansionCoeff() const
{
    throw CanteraError("PDSS::thermalExpansionCoeff()", "unimplemented");
    return 0.0;
}

doublereal PDSS::critTemperature() const
{
    err("critTemperature()");
    return 0.0;
}

doublereal PDSS::critPressure() const
{
    err("critPressure()");
    return 0.0;
}

doublereal PDSS::critDensity() const
{
    err("critDensity()");
    return 0.0;
}

void PDSS::setPressure(doublereal pres)
{
    m_pres = pres;
}

doublereal PDSS::temperature() const
{
    return m_temp;
}

void PDSS::setTemperature(doublereal temp)
{
    m_temp = temp;
}

doublereal PDSS::molecularWeight() const
{
    return m_mw;
}
void PDSS::setMolecularWeight(doublereal mw)
{
    m_mw = mw;
}

void PDSS::setState_TP(doublereal temp, doublereal pres)
{
    err("setState_TP()");
}

void PDSS::setState_TR(doublereal temp, doublereal rho)
{
    err("setState_TR()");
}

doublereal PDSS::satPressure(doublereal t)
{
    err("satPressure()");
    return 0.0;
}

void PDSS::err(const std::string& msg) const
{
    throw CanteraError("PDSS::" + msg, "unimplemented");
}

void PDSS::reportParams(size_t& kindex, int& type,
                        doublereal* const c,
                        doublereal& minTemp_,
                        doublereal& maxTemp_,
                        doublereal& refPressure_) const
{
    kindex = m_spindex;
    type = m_pdssType;
    minTemp_ = m_minTemp;
    maxTemp_ = m_maxTemp;
    refPressure_ = m_p0;
}
}
