/**
 * @file PDSS_IdealGas.cpp
 * Implementation of a pressure dependent standard state
 * virtual function.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/ctml.h"
#include "cantera/thermo/PDSS_IdealGas.h"
#include "cantera/thermo/VPStandardStateTP.h"

using namespace std;

namespace Cantera
{

PDSS_IdealGas::PDSS_IdealGas()
{
}

void PDSS_IdealGas::initThermo()
{
    PDSS::initThermo();
    m_p0 = m_spthermo->refPressure();
    m_minTemp = m_spthermo->minTemp();
    m_maxTemp = m_spthermo->maxTemp();
}

doublereal PDSS_IdealGas::intEnergy_mole() const
{
    return (m_h0_RT - 1.0) * GasConstant * m_temp;
}

doublereal PDSS_IdealGas::cv_mole() const
{
    return cp_mole() - GasConstant;
}

doublereal PDSS_IdealGas::pressure() const
{
    throw NotImplementedError("PDSS_IdealGas::pressure");
}

void PDSS_IdealGas::setPressure(doublereal p)
{
    m_sss_R = m_s0_R - log(m_pres/m_p0);
    m_gss_RT = m_hss_RT - m_sss_R;
    m_Vss = GasConstant * m_temp / m_pres;
}

void PDSS_IdealGas::setTemperature(doublereal temp)
{
    m_temp = temp;
    m_spthermo->updatePropertiesTemp(temp, &m_cp0_R, &m_h0_RT, &m_s0_R);
    m_g0_RT = m_h0_RT - m_s0_R;
    m_V0 = GasConstant * m_temp / m_p0;
    m_hss_RT = m_h0_RT;
    m_cpss_R = m_cp0_R;
    m_sss_R = m_s0_R - log(m_pres/m_p0);
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
