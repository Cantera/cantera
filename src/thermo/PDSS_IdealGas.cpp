/**
 * @file PDSS_IdealGas.cpp
 * Implementation of a pressure dependent standard state
 * virtual function.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/base/ctml.h"
#include "cantera/thermo/PDSS_IdealGas.h"
#include "cantera/thermo/VPStandardStateTP.h"

using namespace std;

namespace Cantera
{
PDSS_IdealGas::PDSS_IdealGas(VPStandardStateTP* tp, int spindex) :
    PDSS(tp, spindex)
{
}

PDSS_IdealGas::PDSS_IdealGas(VPStandardStateTP* tp, size_t spindex, const XML_Node& speciesNode,
                             const XML_Node& phaseRoot, bool spInstalled) :
    PDSS(tp, spindex)
{
    if (!spInstalled) {
        throw CanteraError("PDSS_IdealGas", "sp installing not done yet");
    }
    constructPDSSXML(tp, spindex, phaseRoot, "");
}

void PDSS_IdealGas::constructPDSSXML(VPStandardStateTP* tp, size_t spindex,
                                     const XML_Node& phaseNode, const std::string& id)
{
}

void PDSS_IdealGas::initThermo()
{
    PDSS::initThermo();
    m_p0 = m_tp->speciesThermo().refPressure(m_spindex);
    m_minTemp = m_spthermo->minTemp(m_spindex);
    m_maxTemp = m_spthermo->maxTemp(m_spindex);
}

doublereal PDSS_IdealGas::enthalpy_RT() const
{
    return m_h0_RT;
}

doublereal PDSS_IdealGas::intEnergy_mole() const
{
    return (m_h0_RT - 1.0) * GasConstant * m_temp;
}

doublereal PDSS_IdealGas::entropy_R() const
{
    return m_s0_R - log(m_pres/m_p0);
}

doublereal PDSS_IdealGas::gibbs_RT() const
{
    return m_g0_RT + log(m_pres/m_p0);
}

doublereal PDSS_IdealGas::cp_R() const
{
    return m_cp0_R;
}

doublereal PDSS_IdealGas::molarVolume() const
{
    return GasConstant * m_temp / m_pres;
}

doublereal PDSS_IdealGas::density() const
{
    return m_pres * m_mw / (GasConstant * m_temp);
}

doublereal PDSS_IdealGas::cv_mole() const
{
    return cp_mole() - GasConstant;
}

doublereal PDSS_IdealGas::gibbs_RT_ref() const
{
    return m_g0_RT;
}

doublereal PDSS_IdealGas::enthalpy_RT_ref() const
{
    return m_h0_RT;
}

doublereal PDSS_IdealGas::entropy_R_ref() const
{
    return m_s0_R;
}

doublereal PDSS_IdealGas::cp_R_ref() const
{
    return cp_R();
}

doublereal PDSS_IdealGas::molarVolume_ref() const
{
    return GasConstant * m_temp / m_p0;
}

doublereal PDSS_IdealGas::pressure() const
{
    throw CanteraError("PDSS_IdealGas::pressure()", "unimplemented");
}

void PDSS_IdealGas::setPressure(doublereal p)
{
    m_sss_R = m_s0_R + log(m_pres/m_p0);
    m_gss_RT = m_hss_RT - m_sss_R;
    m_Vss = GasConstant * m_temp / m_pres;
}

doublereal PDSS_IdealGas::temperature() const
{
    return m_temp;
}

void PDSS_IdealGas::setTemperature(doublereal temp)
{
    m_temp = temp;
    m_spthermo->update_single(m_spindex, temp, &m_cp0_R, &m_h0_RT, &m_s0_R);
    m_g0_RT = m_h0_RT - m_s0_R;
    m_V0 = GasConstant * m_temp / m_p0;
    m_hss_RT = m_h0_RT;
    m_cpss_R = m_cp0_R;
    m_sss_R = m_s0_R + log(m_pres/m_p0);
    m_gss_RT = m_hss_RT - m_sss_R;
    m_Vss = GasConstant * m_temp / m_pres;
}

void PDSS_IdealGas::setState_TP(doublereal temp, doublereal pres)
{
    m_pres = pres;
    setTemperature(temp);
}

void PDSS_IdealGas::setState_TR(doublereal temp, doublereal rho)
{
    m_pres = GasConstant * temp * rho / m_mw;
    setTemperature(temp);
}

}
