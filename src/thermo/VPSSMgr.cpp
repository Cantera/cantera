/**
 *  @file VPSSMgr.cpp
 * Definition file for a virtual base class that manages
 * the calculation of standard state properties for all of the
 * species in a single phase, assuming a variable P and T standard state
 * (see \ref mgrpdssthermocalc and
 * class \link Cantera::VPSSMgr VPSSMgr\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/VPSSMgr.h"
#include "cantera/thermo/VPStandardStateTP.h"
#include "cantera/thermo/SpeciesThermoFactory.h"
#include "cantera/thermo/PDSS.h"
#include "cantera/thermo/MultiSpeciesThermo.h"
#include "cantera/base/utilities.h"
#include "cantera/base/xml.h"

using namespace std;

namespace Cantera
{
VPSSMgr::VPSSMgr(VPStandardStateTP* vptp_ptr, MultiSpeciesThermo* spthermo) :
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

VPSSMgr::VPSSMgr(const VPSSMgr& right) :
    m_kk(0),
    m_vptp_ptr(0),
    m_spthermo(0),
    m_tlast(-1.0),
    m_plast(-1.0),
    m_p0(-1.0),
    m_minTemp(-1.0),
    m_maxTemp(1.0E8),
    m_useTmpRefStateStorage(false),
    m_useTmpStandardStateStorage(false)
{
    warn_deprecated("VPSSMgr copy constructor", "To be removed after"
        " Cantera 2.3 for all classes derived from VPSSMgr.");
    *this = right;
}

VPSSMgr& VPSSMgr::operator=(const VPSSMgr& right)
{
    warn_deprecated("VPSSMgr assignment operator", "To be removed after"
        " Cantera 2.3 for all classes derived from VPSSMgr.");
    if (&right == this) {
        return *this;
    }
    m_kk = right.m_kk;

    // What we are doing here is to make a shallow copy of the VPStandardStateTP
    // pointer in the "new" VPSSMgr object using the value from the "old"
    // VPSSMgr object. This is not appropriate if we are making a copy of a
    // ThermoPhase object and the VPSSMgr objects are owned by the ThermoPhase
    // object.
    //
    // The new object will want to have a different value of m_vptp_ptr than the
    // value this is being copied here. It will want to refer to the copy of the
    // VPStandardStateTP object being made that will own the new VPSSMgr object.
    // However, the assignment object is not the place to carry out this fixup.
    //
    // We will have to "fix" up the shallow copies later.
    m_vptp_ptr = right.m_vptp_ptr;
    m_spthermo = right.m_spthermo;
    m_tlast = -1.0;
    m_plast = -1.0;
    m_p0 = right.m_p0;
    m_minTemp = right.m_minTemp;
    m_maxTemp = right.m_maxTemp;
    m_useTmpRefStateStorage = right.m_useTmpRefStateStorage;
    m_h0_RT = right.m_h0_RT;
    m_cp0_R = right.m_cp0_R;
    m_g0_RT = right.m_g0_RT;
    m_s0_R = right.m_s0_R;
    m_V0 = right.m_V0;
    m_useTmpStandardStateStorage = right.m_useTmpStandardStateStorage;
    m_hss_RT = right.m_hss_RT;
    m_cpss_R = right.m_cpss_R;
    m_gss_RT = right.m_gss_RT;
    m_sss_R = right.m_sss_R;
    m_Vss = right.m_Vss;
    mPDSS_h0_RT = right.mPDSS_h0_RT;
    mPDSS_cp0_R = right.mPDSS_cp0_R;
    mPDSS_g0_RT = right.mPDSS_g0_RT;
    mPDSS_s0_R = right.mPDSS_s0_R;
    mPDSS_V0 = right.mPDSS_V0;
    mPDSS_hss_RT = right.mPDSS_hss_RT;
    mPDSS_cpss_R = right.mPDSS_cpss_R;
    mPDSS_gss_RT = right.mPDSS_gss_RT;
    mPDSS_sss_R = right.mPDSS_sss_R;
    mPDSS_Vss = right.mPDSS_Vss;

    return *this;
}

VPSSMgr* VPSSMgr::duplMyselfAsVPSSMgr() const
{
    warn_deprecated("VPSSMgr::duplMyselfAsVPSSMgr", "To be removed after"
        " Cantera 2.3 for all classes derived from VPSSMgr.");
    return new VPSSMgr(*this);
}

void VPSSMgr::initAllPtrs(VPStandardStateTP* vp_ptr,
                          MultiSpeciesThermo* sp_ptr)
{
    warn_deprecated("VPSSMgr::initAllPtrs", "To be removed after Cantera 2.3 "
        "for all classes derived from VPSSMgr.");
    m_vptp_ptr = vp_ptr;
    m_spthermo = sp_ptr;

    // Take care of STITTbyPDSS objects
    for (size_t k = 0; k < m_kk; k++) {
        SpeciesThermoInterpType* st = m_spthermo->provideSTIT(k);
        STITbyPDSS* stpd = dynamic_cast<STITbyPDSS*>(st);
        if (stpd) {
            PDSS* PDSS_ptr = vp_ptr->providePDSS(k);
            stpd->initAllPtrs(k, this, PDSS_ptr);
        }
    }
}

// Standard States

void VPSSMgr::getStandardChemPotentials(doublereal* mu) const
{
    if (m_useTmpStandardStateStorage) {
        std::copy(m_gss_RT.begin(), m_gss_RT.end(), mu);
        scale(mu, mu+m_kk, mu, GasConstant * m_tlast);
    } else {
        throw NotImplementedError("VPSSMgr::getStandardChemPotentials");
    }
}

void VPSSMgr::getGibbs_RT(doublereal* grt) const
{
    if (m_useTmpStandardStateStorage) {
        std::copy(m_gss_RT.begin(), m_gss_RT.end(), grt);
    } else {
        throw NotImplementedError("VPSSMgr::getGibbs_RT");
    }
}

void VPSSMgr::getEnthalpy_RT(doublereal* hrt) const
{
    if (m_useTmpStandardStateStorage) {
        std::copy(m_hss_RT.begin(), m_hss_RT.end(), hrt);
    } else {
        throw NotImplementedError("VPSSMgr::getEnthalpy_RT");
    }
}

void VPSSMgr::getEntropy_R(doublereal* sr) const
{
    if (m_useTmpStandardStateStorage) {
        std::copy(m_sss_R.begin(), m_sss_R.end(), sr);
    } else {
        throw NotImplementedError("VPSSMgr::getEntropy_RT");
    }
}

void VPSSMgr::getIntEnergy_RT(doublereal* urt) const
{
    if (m_useTmpStandardStateStorage) {
        std::copy(m_hss_RT.begin(), m_hss_RT.end(), urt);
        for (size_t k = 0; k < m_kk; k++) {
            urt[k] -= m_plast / (GasConstant * m_tlast) * m_Vss[k];
        }
    } else {
        throw NotImplementedError("VPSSMgr::getEntropy_RT");
    }
}

void VPSSMgr::getCp_R(doublereal* cpr) const
{
    if (m_useTmpStandardStateStorage) {
        std::copy(m_cpss_R.begin(), m_cpss_R.end(), cpr);
    } else {
        throw NotImplementedError("VPSSMgr::getCp_R");
    }
}

void VPSSMgr::getStandardVolumes(doublereal* vol) const
{
    if (m_useTmpStandardStateStorage) {
        std::copy(m_Vss.begin(), m_Vss.end(), vol);
    } else {
        throw NotImplementedError("VPSSMgr::getStandardVolumes");
    }
}
const vector_fp& VPSSMgr::getStandardVolumes() const
{
    if (!m_useTmpStandardStateStorage) {
        throw NotImplementedError("VPSSMgr::getStandardVolumes");
    }
    return m_Vss;
}

/*****************************************************************/

void VPSSMgr::getEnthalpy_RT_ref(doublereal* hrt) const
{
    if (m_useTmpRefStateStorage) {
        std::copy(m_h0_RT.begin(), m_h0_RT.end(), hrt);
    } else {
        throw NotImplementedError("VPSSMgr::getEnthalpy_RT_ref");
    }
}

void VPSSMgr::getGibbs_RT_ref(doublereal* grt) const
{
    if (m_useTmpRefStateStorage) {
        std::copy(m_g0_RT.begin(), m_g0_RT.end(), grt);
    } else {
        throw NotImplementedError("VPSSMgr::getGibbs_RT_ref");
    }
}

void VPSSMgr::getGibbs_ref(doublereal* g) const
{
    if (m_useTmpRefStateStorage) {
        std::copy(m_g0_RT.begin(), m_g0_RT.end(), g);
        scale(g, g+m_kk, g, GasConstant * m_tlast);
    } else {
        throw NotImplementedError("VPSSMgr::getGibbs_ref");
    }
}

void VPSSMgr::getEntropy_R_ref(doublereal* sr) const
{
    if (m_useTmpRefStateStorage) {
        std::copy(m_s0_R.begin(), m_s0_R.end(), sr);
    } else {
        throw NotImplementedError("VPSSMgr::getEntropy_R_ref");
    }
}

void VPSSMgr::getCp_R_ref(doublereal* cpr) const
{
    if (m_useTmpRefStateStorage) {
        std::copy(m_cp0_R.begin(), m_cp0_R.end(), cpr);
    } else {
        throw NotImplementedError("VPSSMgr::getCp_R_ref");
    }
}

void VPSSMgr::getStandardVolumes_ref(doublereal* vol) const
{
    getStandardVolumes(vol);
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
    throw NotImplementedError("VPSSMgr::_updateStandardStateThermo()");
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

void VPSSMgr::initThermo()
{
    initLengths();
}

void VPSSMgr::initLengths()
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
        m_minTemp = std::max(m_minTemp, kPDSS->minTemp());
        m_maxTemp = std::min(m_maxTemp, kPDSS->maxTemp());
    }
}

void VPSSMgr::installSTSpecies(size_t k, const XML_Node& s,
                               const XML_Node* phaseNode_ptr)
{
    shared_ptr<SpeciesThermoInterpType> stit(newSpeciesThermoInterpType(s.child("thermo")));
    stit->validate(s["name"]);
    m_spthermo->install_STIT(k, stit);
    if (m_p0 < 0.0) {
        m_p0 = m_spthermo->refPressure(k);
    }
}

PDSS* VPSSMgr::createInstallPDSS(size_t k, const XML_Node& s,
                                 const XML_Node* phaseNode_ptr)
{
    throw NotImplementedError("VPSSMgr::VPSSMgr::createInstallPDSS");
}

/*****************************************************************/

doublereal VPSSMgr::minTemp(size_t k) const
{
    if (k != npos) {
        return m_vptp_ptr->providePDSS(k)->minTemp();
    }
    return m_minTemp;
}

doublereal VPSSMgr::maxTemp(size_t k) const
{
    if (k != npos) {
        return m_vptp_ptr->providePDSS(k)->maxTemp();
    }
    return m_maxTemp;
}

doublereal VPSSMgr::refPressure(size_t k) const
{
    if (k != npos) {
        return m_vptp_ptr->providePDSS(k)->refPressure();
    }
    return m_p0;
}

PDSS_enumType VPSSMgr::reportPDSSType(int index) const
{
    throw NotImplementedError("VPSSMgr::reportPDSSType()");
}

VPSSMgr_enumType VPSSMgr::reportVPSSMgrType() const
{
    throw NotImplementedError("VPSSMgr::reportVPSSType()");
}

}
