/**
 *  @file VPSSMgr.cpp
 * Definition file for a virtual base class that manages
 * the calculation of standard state properties for all of the
 * species in a single phase, assuming a variable P and T standard state
 * (see \ref mgrpdssthermocalc and
 * class \link Cantera::VPSSMgr VPSSMgr\endlink).
 */
/*
 * Copyright (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#include "cantera/thermo/VPSSMgr.h"
#include "cantera/thermo/VPStandardStateTP.h"
#include "cantera/thermo/SpeciesThermoFactory.h"
#include "cantera/thermo/PDSS.h"
#include "cantera/thermo/GeneralSpeciesThermo.h"

using namespace std;

namespace Cantera
{
class SpeciesThermo;

VPSSMgr::VPSSMgr(VPStandardStateTP* vptp_ptr, SpeciesThermo* spthermo) :
    m_kk(0),
    m_vptp_ptr(vptp_ptr),
    m_spthermo(spthermo),
    m_tlast(-1.0),
    m_plast(-1.0),
    m_p0(-1.0),
    m_minTemp(-1.0),
    m_maxTemp(1.0E8),
    m_useTmpRefStateStorage(false),
    m_useTmpStandardStateStorage(false)
{
    if (!m_vptp_ptr) {
        throw CanteraError("VPSSMgr",
                           "null pointer for VPStandardStateTP is not permissible");
    }
}

VPSSMgr::~VPSSMgr()
{
}

VPSSMgr::VPSSMgr(const VPSSMgr& right) :
    m_kk(0),
    m_vptp_ptr(0),
    m_spthermo(0),
    //  m_Tnow(300.),
    //   m_Pnow(OneAtm),
    m_tlast(-1.0),
    m_plast(-1.0),
    m_p0(-1.0),
    m_minTemp(-1.0),
    m_maxTemp(1.0E8),
    m_useTmpRefStateStorage(false),
    m_useTmpStandardStateStorage(false)
{
    *this = right;
}

VPSSMgr&
VPSSMgr::operator=(const VPSSMgr& right)
{
    if (&right == this) {
        return *this;
    }
    m_kk                          = right.m_kk;
    /*
     * What we are doing here is to make a shallow copy of the VPStandardStateTP
     * pointer in the "new" VPSSMgr object using the value from the "old"
     * VPSSMgr object. This is not appropriate if we are making a copy of a ThermoPhase
     * object and the VPSSMgr objects are owned by the ThermoPhase object.
     *
     * The new object will want to have a different value of m_vptp_ptr than the
     * value this is being copied here. It will want to refer to the copy of the
     * VPStandardStateTP object being made that will own the new VPSSMgr object.
     * However, the assignment object is not the place to carry out this fixup.
     *
     * We will have to "fix" up the shallow copies later.
     */
    m_vptp_ptr                    = right.m_vptp_ptr;
    m_spthermo                    = right.m_spthermo;
    m_tlast                       = -1.0;
    m_plast                       = -1.0;
    m_p0                          = right.m_p0;
    m_minTemp                     = right.m_minTemp;
    m_maxTemp                     = right.m_maxTemp;
    m_useTmpRefStateStorage       = right.m_useTmpRefStateStorage;
    m_h0_RT                       = right.m_h0_RT;
    m_cp0_R                       = right.m_cp0_R;
    m_g0_RT                       = right.m_g0_RT;
    m_s0_R                        = right.m_s0_R;
    m_V0                          = right.m_V0;
    m_useTmpStandardStateStorage  = right.m_useTmpStandardStateStorage;
    m_hss_RT                      = right.m_hss_RT;
    m_cpss_R                      = right.m_cpss_R;
    m_gss_RT                      = right.m_gss_RT;
    m_sss_R                       = right.m_sss_R;
    m_Vss                         = right.m_Vss;

    mPDSS_h0_RT                   = right.mPDSS_h0_RT;
    mPDSS_cp0_R                   = right.mPDSS_cp0_R;
    mPDSS_g0_RT                   = right.mPDSS_g0_RT;
    mPDSS_s0_R                    = right.mPDSS_s0_R;
    mPDSS_V0                      = right.mPDSS_V0;
    mPDSS_hss_RT                  = right.mPDSS_hss_RT;
    mPDSS_cpss_R                  = right.mPDSS_cpss_R;
    mPDSS_gss_RT                  = right.mPDSS_gss_RT;
    mPDSS_sss_R                   = right.mPDSS_sss_R;
    mPDSS_Vss                     = right.mPDSS_Vss;

    return *this;
}

VPSSMgr* VPSSMgr::duplMyselfAsVPSSMgr() const
{
    return new VPSSMgr(*this);
}

void VPSSMgr::initAllPtrs(VPStandardStateTP* vp_ptr,
                          SpeciesThermo* sp_ptr)
{
    m_vptp_ptr = vp_ptr;
    m_spthermo = sp_ptr;

    // Take care of STITTbyPDSS objects

    // Go see if the SpeciesThermo type is a GeneralSpeciesThermo
    GeneralSpeciesThermo* gst = dynamic_cast<GeneralSpeciesThermo*>(sp_ptr);
    if (gst) {
        for (size_t k = 0; k < m_kk; k++) {
            SpeciesThermoInterpType* st = gst->provideSTIT(k);
            STITbyPDSS* stpd = dynamic_cast<STITbyPDSS*>(st);
            if (stpd) {
                PDSS* PDSS_ptr = vp_ptr->providePDSS(k);
                stpd->initAllPtrs(k, this, PDSS_ptr);
            }
        }
    }

}

// Standard States

void
VPSSMgr::getStandardChemPotentials(doublereal* mu) const
{
    if (m_useTmpStandardStateStorage) {
        std::copy(m_gss_RT.begin(), m_gss_RT.end(), mu);
        doublereal _rt = GasConstant * m_tlast;
        scale(mu, mu+m_kk, mu, _rt);
    } else {
        err("getStandardChemPotentials");
    }
}

void
VPSSMgr::getGibbs_RT(doublereal* grt) const
{
    if (m_useTmpStandardStateStorage) {
        std::copy(m_gss_RT.begin(), m_gss_RT.end(), grt);
    } else {
        err("getGibbs_RT");
    }
}

void
VPSSMgr::getEnthalpy_RT(doublereal* hrt) const
{
    if (m_useTmpStandardStateStorage) {
        std::copy(m_hss_RT.begin(), m_hss_RT.end(), hrt);
    } else {
        err("getEnthalpy_RT");
    }
}

void
VPSSMgr::getEntropy_R(doublereal* sr) const
{
    if (m_useTmpStandardStateStorage) {
        std::copy(m_sss_R.begin(), m_sss_R.end(), sr);
    } else {
        err("getEntropy_RT");
    }
}

void
VPSSMgr::getIntEnergy_RT(doublereal* urt) const
{
    if (m_useTmpStandardStateStorage) {
        std::copy(m_hss_RT.begin(), m_hss_RT.end(), urt);
        doublereal pRT = m_plast / (GasConstant * m_tlast);
        for (size_t k = 0; k < m_kk; k++) {
            urt[k] -= pRT * m_Vss[k];
        }
    } else {
        err("getEntropy_RT");
    }
}

void
VPSSMgr::getCp_R(doublereal* cpr) const
{
    if (m_useTmpStandardStateStorage) {
        std::copy(m_cpss_R.begin(), m_cpss_R.end(), cpr);
    } else {
        err("getCp_R");
    }
}

void
VPSSMgr::getStandardVolumes(doublereal* vol) const
{
    if (m_useTmpStandardStateStorage) {
        std::copy(m_Vss.begin(), m_Vss.end(), vol);
    } else {
        err("getStandardVolumes");
    }
}
const vector_fp&
VPSSMgr::getStandardVolumes() const
{
    if (m_useTmpStandardStateStorage) {
        return m_Vss;
    } else {
        err("getStandardVolumes");
    }
}

/*****************************************************************/
void
VPSSMgr::getEnthalpy_RT_ref(doublereal* hrt) const
{
    if (m_useTmpRefStateStorage) {
        std::copy(m_h0_RT.begin(), m_h0_RT.end(), hrt);
    } else {
        err("getEnthalpy_RT_ref");
    }
}

void
VPSSMgr::getGibbs_RT_ref(doublereal* grt) const
{
    if (m_useTmpRefStateStorage) {
        std::copy(m_g0_RT.begin(), m_g0_RT.end(), grt);
    } else {
        err("getGibbs_RT_ref");
    }
}

void
VPSSMgr::getGibbs_ref(doublereal* g) const
{
    if (m_useTmpRefStateStorage) {
        std::copy(m_g0_RT.begin(), m_g0_RT.end(), g);
        doublereal _rt = GasConstant * m_tlast;
        scale(g, g+m_kk, g, _rt);
    } else {
        err("getGibbs_ref");
    }
}

void
VPSSMgr::getEntropy_R_ref(doublereal* sr) const
{
    if (m_useTmpRefStateStorage) {
        std::copy(m_s0_R.begin(), m_s0_R.end(), sr);
    } else {
        err("getEntropy_R_ref");
    }
}

void
VPSSMgr::getCp_R_ref(doublereal* cpr) const
{
    if (m_useTmpRefStateStorage) {
        std::copy(m_cp0_R.begin(), m_cp0_R.end(), cpr);
    } else {
        err("getCp_R_ref");
    }
}

void
VPSSMgr::getStandardVolumes_ref(doublereal* vol) const
{
    getStandardVolumes(vol);
    //err("getStandardVolumes_ref");
}

/*****************************************************************/

void VPSSMgr::setState_P(doublereal pres)
{
    if (m_plast != pres) {
        m_plast = pres;
        updateStandardStateThermo();
    }
}

void VPSSMgr::setState_T(doublereal temp)
{
    if (m_tlast != temp) {
        m_tlast = temp;
        updateRefStateThermo();
        updateStandardStateThermo();
    }
}

void VPSSMgr::setState_TP(doublereal temp, doublereal pres)
{
    if (m_tlast != temp) {
        m_tlast = temp;
        m_plast = pres;
        updateRefStateThermo();
        updateStandardStateThermo();
    } else if (m_plast != pres) {
        m_plast = pres;
        updateStandardStateThermo();
    }
}

void VPSSMgr::updateStandardStateThermo()
{
    _updateStandardStateThermo();
}

void VPSSMgr::updateRefStateThermo() const
{
    _updateRefStateThermo();
}

void VPSSMgr::_updateStandardStateThermo()
{
    for (size_t k = 0; k < m_kk; k++) {
        PDSS* kPDSS = m_vptp_ptr->providePDSS(k);
        kPDSS->setState_TP(m_tlast, m_plast);
    }
    err("_updateStandardStateThermo()");
}

void VPSSMgr::_updateRefStateThermo() const
{
    if (m_spthermo) {
        m_spthermo->update(m_tlast, &m_cp0_R[0], &m_h0_RT[0], &m_s0_R[0]);
        for (size_t k = 0; k < m_kk; k++) {
            m_g0_RT[k] = m_h0_RT[k] - m_s0_R[k];
        }
    }
}


/*****************************************************************/

void
VPSSMgr::initThermo()
{
    initLengths();
}

void
VPSSMgr::initLengths()
{
    m_kk = m_vptp_ptr->nSpecies();
    m_h0_RT.resize(m_kk, 0.0);
    m_cp0_R.resize(m_kk, 0.0);
    m_g0_RT.resize(m_kk, 0.0);
    m_s0_R.resize(m_kk, 0.0);
    m_V0.resize(m_kk, 0.0);
    m_hss_RT.resize(m_kk, 0.0);
    m_cpss_R.resize(m_kk, 0.0);
    m_gss_RT.resize(m_kk, 0.0);
    m_sss_R.resize(m_kk, 0.0);
    m_Vss.resize(m_kk, 0.0);

    // Storage used by the PDSS objects to store their answers.
    mPDSS_h0_RT.resize(m_kk, 0.0);
    mPDSS_cp0_R.resize(m_kk, 0.0);
    mPDSS_g0_RT.resize(m_kk, 0.0);
    mPDSS_s0_R.resize(m_kk, 0.0);
    mPDSS_V0.resize(m_kk, 0.0);
    mPDSS_hss_RT.resize(m_kk, 0.0);
    mPDSS_cpss_R.resize(m_kk, 0.0);
    mPDSS_gss_RT.resize(m_kk, 0.0);
    mPDSS_sss_R.resize(m_kk, 0.0);
    mPDSS_Vss.resize(m_kk, 0.0);
}

void VPSSMgr::initThermoXML(XML_Node& phaseNode, const std::string& id)
{
    const PDSS* kPDSS = m_vptp_ptr->providePDSS(0);
    m_p0 = kPDSS->refPressure();
    for (size_t i = 0; i < m_kk; i++) {
        const PDSS* kPDSS = m_vptp_ptr->providePDSS(i);
        doublereal mint = kPDSS->minTemp();
        if (mint > m_minTemp) {
            m_minTemp = mint;
        }
        mint = kPDSS->maxTemp();
        if (mint < m_maxTemp) {
            m_maxTemp = mint;
        }
    }
#ifdef DEBUG_MODE
    // Add a check to see that all references pressures are the same
    double m_p0_k;
    if (m_spthermo) {
        for (size_t k = 0; k < m_kk; k++) {
            m_p0_k = m_spthermo->refPressure(k);
            if (m_p0 != m_p0_k) {
                //throw CanteraError("VPSSMgr::initThermoXML",
                //    "inconsistent ref pressures" + fp2str(m_p0) + " "
                //    + fp2str(m_p0_k));
                // writelog("VPSSMgr::initThermoXML:"
                //    "inconsistent ref pressures: " + fp2str(m_p0) + " "
                //    + fp2str(m_p0_k) + " for SpeciesThermo k = " + int2str(k) + "\n");
            }
        }
    }

    for (size_t k = 0; k < m_kk; k++) {
        const PDSS* kPDSS = m_vptp_ptr->providePDSS(k);
        m_p0_k = kPDSS->refPressure();
        if (m_p0 != m_p0_k) {
            //throw CanteraError("VPSSMgr::initThermoXML",
            //    "inconsistent ref pressures" + fp2str(m_p0) + " "
            //    + fp2str(m_p0_k));
            //writelog("VPSSMgr::initThermoXML"
            //    "inconsistent ref pressures: " + fp2str(m_p0) + " "
            //    + fp2str(m_p0_k) + " for PDSS k = " + int2str(k) + "\n");
        }
    }
#endif
}

void VPSSMgr::installSTSpecies(size_t k,  const XML_Node& s,
                               const XML_Node* phaseNode_ptr)
{

    SpeciesThermoFactory*  f = SpeciesThermoFactory::factory();
    f->installThermoForSpecies(k, s, m_vptp_ptr, *m_spthermo, phaseNode_ptr);
    if (m_p0 < 0.0) {
        m_p0 = m_spthermo->refPressure(k);
    }
}

PDSS* VPSSMgr::createInstallPDSS(size_t k, const XML_Node& s,
                                 const XML_Node* phaseNode_ptr)
{
    err("VPSSMgr::createInstallPDSS");
    return (PDSS*) 0;
}

/*****************************************************************/
doublereal VPSSMgr::minTemp(size_t k) const
{
    if (k != npos) {
        const PDSS* kPDSS = m_vptp_ptr->providePDSS(k);
        return kPDSS->minTemp();
    }
    return m_minTemp;
}

doublereal VPSSMgr::maxTemp(size_t k) const
{
    if (k != npos) {
        const PDSS* kPDSS = m_vptp_ptr->providePDSS(k);
        return kPDSS->maxTemp();
    }
    return m_maxTemp;
}

doublereal VPSSMgr::refPressure(size_t k) const
{
    if (k != npos) {
        const PDSS* kPDSS = m_vptp_ptr->providePDSS(k);
        return kPDSS->refPressure();
    }
    return m_p0;
}

PDSS_enumType VPSSMgr::reportPDSSType(int index) const
{
    err("reportPDSSType()");
    return cPDSS_UNDEF;
}


VPSSMgr_enumType VPSSMgr::reportVPSSMgrType() const
{
    err("reportVPSSType()");
    return cVPSSMGR_UNDEF;
}

/*****************************************************************/

void VPSSMgr::err(const std::string& msg) const
{
    throw CanteraError("VPSSMgr::" + msg, "unimplemented");
}
}
