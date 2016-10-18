/**
 *  @file VPSSMgr_IdealGas.cpp
 * Definition file for a derived class that handles the calculation
 * of standard state thermo properties for
 *  a set of species which have an Ideal Gas dependence
 * (see \ref thermoprops and
 * class \link Cantera::VPSSMgr_IdealGas VPSSMgr_IdealGas\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/VPSSMgr_IdealGas.h"
#include "cantera/base/ctml.h"
#include "cantera/thermo/SpeciesThermoFactory.h"
#include "cantera/thermo/MultiSpeciesThermo.h"
#include "cantera/thermo/PDSS_IdealGas.h"

using namespace std;

namespace Cantera
{

VPSSMgr_IdealGas::VPSSMgr_IdealGas(VPStandardStateTP* vp_ptr, MultiSpeciesThermo* spth) :
    VPSSMgr(vp_ptr, spth)
{
    m_useTmpRefStateStorage = true;
    m_useTmpStandardStateStorage = true;
}

VPSSMgr_IdealGas::VPSSMgr_IdealGas(const VPSSMgr_IdealGas& right) :
    VPSSMgr(right.m_vptp_ptr, right.m_spthermo)
{
    m_useTmpRefStateStorage = true;
    m_useTmpStandardStateStorage = true;
    *this = right;
}

VPSSMgr_IdealGas& VPSSMgr_IdealGas::operator=(const VPSSMgr_IdealGas& b)
{
    if (&b == this) {
        return *this;
    }
    VPSSMgr::operator=(b);
    return *this;
}

VPSSMgr* VPSSMgr_IdealGas::duplMyselfAsVPSSMgr() const
{
    return new VPSSMgr_IdealGas(*this);
}

void VPSSMgr_IdealGas::getIntEnergy_RT(doublereal* urt) const
{
    getEnthalpy_RT(urt);
    for (size_t k = 0; k < m_kk; k++) {
        urt[k] -= 1.0;
    }
}

void VPSSMgr_IdealGas::getStandardVolumes(doublereal* vol) const
{
    copy(m_Vss.begin(), m_Vss.end(), vol);
}

void VPSSMgr_IdealGas::_updateStandardStateThermo()
{
    doublereal pp = log(m_plast / m_p0);
    doublereal v = temperature() *GasConstant /m_plast;

    for (size_t k = 0; k < m_kk; k++) {
        m_hss_RT[k] = m_h0_RT[k];
        m_cpss_R[k] = m_cp0_R[k];
        m_sss_R[k] = m_s0_R[k] - pp;
        m_gss_RT[k] = m_hss_RT[k] - m_sss_R[k];
        m_Vss[k] = v;
    }
}

PDSS* VPSSMgr_IdealGas::createInstallPDSS(size_t k, const XML_Node& speciesNode,
                                          const XML_Node* const phaseNode_ptr)
{
    const XML_Node* ss = speciesNode.findByName("standardState");
    if (ss && ss->attrib("model") != "ideal_gas") {
        throw CanteraError("VPSSMgr_IdealGas::createInstallPDSS",
                           "standardState model for species isn't "
                           "ideal_gas: " + speciesNode["name"]);
    }
    if (m_Vss.size() < k+1) {
        m_Vss.resize(k+1, 0.0);
    }

    shared_ptr<SpeciesThermoInterpType> stit(
        newSpeciesThermoInterpType(speciesNode.child("thermo")));
    stit->validate(speciesNode["name"]);
    m_spthermo->install_STIT(k, stit);

    PDSS* kPDSS = new PDSS_IdealGas(m_vptp_ptr, k, speciesNode,
                                    *phaseNode_ptr, true);

    m_p0 = m_spthermo->refPressure(k);
    return kPDSS;
}

PDSS_enumType VPSSMgr_IdealGas::reportPDSSType(int k) const
{
    warn_deprecated("VPSSMgr_IdealGas::reportPDSSType",
        "To be removed after Cantera 2.3.");
    return cPDSS_IDEALGAS;
}

VPSSMgr_enumType VPSSMgr_IdealGas::reportVPSSMgrType() const
{
    warn_deprecated("VPSSMgr_IdealGas::reportVPSSMgrType",
        "To be removed after Cantera 2.3.");
    return cVPSSMGR_IDEALGAS;
}

}
