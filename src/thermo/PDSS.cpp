/**
 * @file PDSS.cpp
 * Implementation of a pressure dependent standard state
 * virtual function
 * (see class \link Cantera::PDSS PDSS\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/base/ctml.h"
#include "cantera/thermo/PDSS.h"
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
        m_spthermo = &tp->speciesThermo();
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
    warn_deprecated("PDSS copy constructor", "To be removed after"
        " Cantera 2.3 for all classes derived from PDSS.");
    // Use the assignment operator to do the brunt of the work for the copy
    // constructor.
    *this = b;
}

PDSS& PDSS::operator=(const PDSS& b)
{
    warn_deprecated("PDSS assignment operator", "To be removed after"
        " Cantera 2.3 for all classes derived from PDSS.");
    if (&b == this) {
        return *this;
    }

    m_pdssType = b.m_pdssType;
    m_temp = b.m_temp;
    m_pres = b.m_pres;
    m_p0 = b.m_p0;
    m_minTemp = b.m_minTemp;
    m_maxTemp = b.m_maxTemp;

    // Pointers which are zero, are properly assigned in the function,
    // initAllPtrs(). which must be called after the assignment operation.
    m_tp = 0;
    m_vpssmgr_ptr = 0;
    m_mw = b.m_mw;
    m_spindex = b.m_spindex;
    m_spthermo = 0;
    m_cp0_R_ptr = 0;
    m_h0_RT_ptr = 0;
    m_s0_R_ptr = 0;
    m_g0_RT_ptr = 0;
    m_V0_ptr = 0;
    m_cpss_R_ptr = 0;
    m_hss_RT_ptr = 0;
    m_sss_R_ptr = 0;
    m_gss_RT_ptr = 0;
    m_Vss_ptr = 0;

    // Here we just fill these in so that local copies within the VPSS object work.
    m_tp = b.m_tp;
    m_vpssmgr_ptr = b.m_vpssmgr_ptr;
    m_spthermo = b.m_spthermo;
    m_cp0_R_ptr = b.m_cp0_R_ptr;
    m_h0_RT_ptr = b.m_h0_RT_ptr;
    m_s0_R_ptr = b.m_s0_R_ptr;
    m_g0_RT_ptr = b.m_g0_RT_ptr;
    m_V0_ptr = b.m_V0_ptr;
    m_cpss_R_ptr = b.m_cpss_R_ptr;
    m_hss_RT_ptr = b.m_hss_RT_ptr;
    m_sss_R_ptr = b.m_sss_R_ptr;
    m_gss_RT_ptr = b.m_gss_RT_ptr;
    m_Vss_ptr = b.m_Vss_ptr;

    return *this;
}

PDSS* PDSS::duplMyselfAsPDSS() const
{
    warn_deprecated("PDSS::duplMyselfAsPDSS",
        "To be removed after Cantera 2.3.");
    return new PDSS(*this);
}

PDSS_enumType PDSS::reportPDSSType() const
{
    warn_deprecated("PDSS::reportPDSSType", "To be removed after Cantera 2.3.");
    return m_pdssType;
}

void PDSS::initThermoXML(const XML_Node& phaseNode, const std::string& id)
{
    AssertThrow(m_tp != 0, "PDSS::initThermoXML()");
    m_p0 = m_vpssmgr_ptr->refPressure(m_spindex);
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
                       MultiSpeciesThermo* spthermo)
{
    warn_deprecated("PDSS::initAllPtrs", "To be removed after Cantera 2.3 "
        "for all classes derived from PDSS.");
    m_tp = tp;
    m_vpssmgr_ptr = vpssmgr_ptr;
    m_spthermo = spthermo;
    initPtrs();
}

void PDSS::initPtrs()
{
    AssertThrow(m_vpssmgr_ptr->mPDSS_h0_RT.size() != 0, "PDSS::initPtrs()");
    m_h0_RT_ptr = &m_vpssmgr_ptr->mPDSS_h0_RT[0];
    m_cp0_R_ptr = &m_vpssmgr_ptr->mPDSS_cp0_R[0];
    m_s0_R_ptr = &m_vpssmgr_ptr->mPDSS_s0_R[0];
    m_g0_RT_ptr = &m_vpssmgr_ptr->mPDSS_g0_RT[0];
    m_V0_ptr = &m_vpssmgr_ptr->mPDSS_V0[0];

    m_hss_RT_ptr = &m_vpssmgr_ptr->mPDSS_hss_RT[0];
    m_cpss_R_ptr = &m_vpssmgr_ptr->mPDSS_cpss_R[0];
    m_sss_R_ptr = &m_vpssmgr_ptr->mPDSS_sss_R[0];
    m_gss_RT_ptr = &m_vpssmgr_ptr->mPDSS_gss_RT[0];
    m_Vss_ptr = &m_vpssmgr_ptr->mPDSS_Vss[0];
}

doublereal PDSS::enthalpy_mole() const
{
    throw NotImplementedError("PDSS::enthalpy_mole()");
}

doublereal PDSS::enthalpy_RT() const
{
    throw NotImplementedError("PDSS::enthalpy_RT()");
}

doublereal PDSS::intEnergy_mole() const
{
    throw NotImplementedError("PDSS::intEnergy_mole()");
}

doublereal PDSS::entropy_mole() const
{
    throw NotImplementedError("PDSS::entropy_mole()");
}

doublereal PDSS::entropy_R() const
{
    throw NotImplementedError("PDSS::entropy_R()");
}

doublereal PDSS::gibbs_mole() const
{
    throw NotImplementedError("PDSS::gibbs_mole()");
}

doublereal PDSS::gibbs_RT() const
{
    throw NotImplementedError("PDSS::gibbs_RT()");
}

doublereal PDSS::cp_mole() const
{
    throw NotImplementedError("PDSS::cp_mole()");
}

doublereal PDSS::cp_R() const
{
    throw NotImplementedError("PDSS::cp_R()");
}

doublereal PDSS::molarVolume() const
{
    throw NotImplementedError("PDSS::molarVolume()");
}

doublereal PDSS::density() const
{
    throw NotImplementedError("PDSS::density()");
}

doublereal PDSS::cv_mole() const
{
    throw NotImplementedError("PDSS::cv_mole()");
}

doublereal PDSS::gibbs_RT_ref() const
{
    throw NotImplementedError("PDSS::gibbs_RT_ref()");
}

doublereal PDSS::enthalpy_RT_ref() const
{
    throw NotImplementedError("PDSS::enthalpy_RT_ref()");
}

doublereal PDSS::entropy_R_ref() const
{
    throw NotImplementedError("PDSS::entropy_RT_ref()");
}

doublereal PDSS::cp_R_ref() const
{
    throw NotImplementedError("PDSS::entropy_RT_ref()");
}

doublereal PDSS::molarVolume_ref() const
{
    throw NotImplementedError("PDSS::molarVolume_ref()");
}

doublereal PDSS::enthalpyDelp_mole() const
{
    return enthalpy_mole() - m_temp * GasConstant * enthalpy_RT_ref();
}

doublereal PDSS::entropyDelp_mole() const
{
    return entropy_mole() - GasConstant * entropy_R_ref();
}

doublereal PDSS::gibbsDelp_mole() const
{
    return gibbs_mole() - m_temp * GasConstant * gibbs_RT_ref();
}

doublereal PDSS::cpDelp_mole() const
{
    return cp_mole() - GasConstant * cp_R_ref();
}

doublereal PDSS::pressure() const
{
    return m_pres;
}

doublereal PDSS::thermalExpansionCoeff() const
{
    throw NotImplementedError("PDSS::thermalExpansionCoeff()");
}

doublereal PDSS::critTemperature() const
{
    throw NotImplementedError("PDSS::critTemperature()");
}

doublereal PDSS::critPressure() const
{
    throw NotImplementedError("PDSS::critPressure()");
}

doublereal PDSS::critDensity() const
{
    throw NotImplementedError("PDSS::critDensity()");
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
    throw NotImplementedError("PDSS::setState_TP()");
}

void PDSS::setState_TR(doublereal temp, doublereal rho)
{
    throw NotImplementedError("PDSS::setState_TR()");
}

doublereal PDSS::satPressure(doublereal t)
{
    throw NotImplementedError("PDSS::satPressure()");
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

// PDSS_Molar methods

doublereal PDSS_Molar::enthalpy_RT() const
{
    return enthalpy_mole() / (GasConstant * temperature());
}

doublereal PDSS_Molar::entropy_R() const
{
    return entropy_mole() / GasConstant;
}

doublereal PDSS_Molar::gibbs_RT() const
{
    return gibbs_mole() / (GasConstant * temperature());
}

doublereal PDSS_Molar::cp_R() const
{
    return cp_mole() / GasConstant;
}

// PDSS_Nondimensional methods

doublereal PDSS_Nondimensional::enthalpy_mole() const
{
    return enthalpy_RT() * GasConstant * temperature();
}

doublereal PDSS_Nondimensional::entropy_mole() const
{
    return entropy_R() * GasConstant;
}

doublereal PDSS_Nondimensional::gibbs_mole() const
{
    return gibbs_RT() * GasConstant * temperature();
}

doublereal PDSS_Nondimensional::cp_mole() const
{
    return cp_R() * GasConstant;
}

}
