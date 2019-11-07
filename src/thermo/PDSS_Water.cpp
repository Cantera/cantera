/**
 * @file PDSS_Water.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/ctml.h"
#include "cantera/thermo/PDSS_Water.h"
#include "cantera/thermo/WaterPropsIAPWS.h"
#include "cantera/base/stringUtils.h"
#include "cantera/thermo/Elements.h"

namespace Cantera
{
PDSS_Water::PDSS_Water() :
    m_waterProps(&m_sub),
    m_dens(1000.0),
    m_iState(WATER_LIQUID),
    EW_Offset(0.0),
    SW_Offset(0.0),
    m_allowGasPhase(false)
{
    m_minTemp = 200.;
    m_maxTemp = 10000.;
    m_mw = 2*getElementWeight("H") + getElementWeight("O");

    // Set the baseline
    doublereal T = 298.15;
    m_p0 = OneAtm;
    doublereal presLow = 1.0E-2;
    doublereal oneBar = 1.0E5;
    doublereal dens = 1.0E-9;
    m_dens = m_sub.density(T, presLow, WATER_GAS, dens);
    m_pres = presLow;
    SW_Offset = 0.0;
    doublereal s = entropy_mole();
    s -= GasConstant * log(oneBar/presLow);
    if (s != 188.835E3) {
        SW_Offset = 188.835E3 - s;
    }
    s = entropy_mole();
    s -= GasConstant * log(oneBar/presLow);

    doublereal h = enthalpy_mole();
    if (h != -241.826E6) {
        EW_Offset = -241.826E6 - h;
    }
    h = enthalpy_mole();

    // Set the initial state of the system to 298.15 K and 1 bar.
    setTemperature(298.15);
    m_dens = m_sub.density(298.15, OneAtm, WATER_LIQUID);
    m_pres = OneAtm;
}

doublereal PDSS_Water::enthalpy_mole() const
{
    return m_sub.enthalpy() + EW_Offset;
}

doublereal PDSS_Water::intEnergy_mole() const
{
    return m_sub.intEnergy() + EW_Offset;
}

doublereal PDSS_Water::entropy_mole() const
{
    return m_sub.entropy() + SW_Offset;
}

doublereal PDSS_Water::gibbs_mole() const
{
    return m_sub.Gibbs() + EW_Offset - SW_Offset*m_temp;
}

doublereal PDSS_Water::cp_mole() const
{
    return m_sub.cp();
}

doublereal PDSS_Water::cv_mole() const
{
    return m_sub.cv();
}

doublereal PDSS_Water::molarVolume() const
{
    return m_sub.molarVolume();
}

doublereal PDSS_Water::gibbs_RT_ref() const
{
    doublereal T = m_temp;
    m_sub.density(T, m_p0, m_iState);
    doublereal h = m_sub.enthalpy();
    m_sub.setState_TR(m_temp, m_dens);
    return (h + EW_Offset - SW_Offset*T)/(T * GasConstant);
}

doublereal PDSS_Water::enthalpy_RT_ref() const
{
    doublereal T = m_temp;
    m_sub.density(T, m_p0, m_iState);
    doublereal h = m_sub.enthalpy();
    m_sub.setState_TR(m_temp, m_dens);
    return (h + EW_Offset)/(T * GasConstant);
}

doublereal PDSS_Water::entropy_R_ref() const
{
    doublereal T = m_temp;
    m_sub.density(T, m_p0, m_iState);
    doublereal s = m_sub.entropy();
    m_sub.setState_TR(m_temp, m_dens);
    return (s + SW_Offset)/GasConstant;
}

doublereal PDSS_Water::cp_R_ref() const
{
    doublereal T = m_temp;
    m_sub.density(T, m_p0, m_iState);
    doublereal cp = m_sub.cp();
    m_sub.setState_TR(m_temp, m_dens);
    return cp/GasConstant;
}

doublereal PDSS_Water::molarVolume_ref() const
{
    doublereal T = m_temp;
    m_sub.density(T, m_p0, m_iState);
    doublereal mv = m_sub.molarVolume();
    m_sub.setState_TR(m_temp, m_dens);
    return mv;
}

doublereal PDSS_Water::pressure() const
{
    m_pres = m_sub.pressure();
    return m_pres;
}

void PDSS_Water::setPressure(doublereal p)
{
    // In this routine we must be sure to only find the water branch of the
    // curve and not the gas branch
    doublereal T = m_temp;
    doublereal dens = m_dens;
    int waterState = WATER_LIQUID;
    if (T > m_sub.Tcrit()) {
        waterState = WATER_SUPERCRIT;
    }

    doublereal dd = m_sub.density(T, p, waterState, dens);
    if (dd <= 0.0) {
        throw CanteraError("PDSS_Water:setPressure",
            "Failed to set water SS state: T = {} K and p = {} Pa", T, p);
    }
    m_dens = dd;
    m_pres = p;

    // We are only putting the phase check here because of speed considerations.
    m_iState = m_sub.phaseState(true);
    if (!m_allowGasPhase && m_iState != WATER_SUPERCRIT && m_iState != WATER_LIQUID && m_iState != WATER_UNSTABLELIQUID) {
        throw CanteraError("PDSS_Water::setPressure",
                           "Water State isn't liquid or crit");
    }
}

doublereal PDSS_Water::thermalExpansionCoeff() const
{
    return m_sub.coeffThermExp();
}

doublereal PDSS_Water::dthermalExpansionCoeffdT() const
{
    doublereal pres = pressure();
    doublereal dens_save = m_dens;
    doublereal tt = m_temp - 0.04;
    doublereal dd = m_sub.density(tt, pres, m_iState, m_dens);
    if (dd < 0.0) {
        throw CanteraError("PDSS_Water::dthermalExpansionCoeffdT",
            "unable to solve for the density at T = {}, P = {}", tt, pres);
    }
    doublereal vald = m_sub.coeffThermExp();
    m_sub.setState_TR(m_temp, dens_save);
    doublereal val2 = m_sub.coeffThermExp();
    return (val2 - vald) / 0.04;
}

doublereal PDSS_Water::isothermalCompressibility() const
{
    return m_sub.isothermalCompressibility();
}

doublereal PDSS_Water::critTemperature() const
{
    return m_sub.Tcrit();
}

doublereal PDSS_Water::critPressure() const
{
    return m_sub.Pcrit();
}

doublereal PDSS_Water::critDensity() const
{
    return m_sub.Rhocrit();
}

void PDSS_Water::setDensity(doublereal dens)
{
    m_dens = dens;
    m_sub.setState_TR(m_temp, m_dens);
}

doublereal PDSS_Water::density() const
{
    return m_dens;
}

void PDSS_Water::setTemperature(doublereal temp)
{
    m_temp = temp;
    m_sub.setState_TR(temp, m_dens);
}

void PDSS_Water::setState_TP(doublereal temp, doublereal pres)
{
    m_temp = temp;
    setPressure(pres);
}

void PDSS_Water::setState_TR(doublereal temp, doublereal dens)
{
    m_temp = temp;
    m_dens = dens;
    m_sub.setState_TR(m_temp, m_dens);
}

doublereal PDSS_Water::pref_safe(doublereal temp) const
{
    if (temp < m_sub.Tcrit()) {
        doublereal pp = m_sub.psat_est(temp);
        if (pp > OneAtm) {
            return pp;
        }
    } else  {
        return m_sub.Pcrit();
    }
    return OneAtm;
}

doublereal PDSS_Water::satPressure(doublereal t)
{
    doublereal pp = m_sub.psat(t, WATER_LIQUID);
    m_dens = m_sub.density();
    m_temp = t;
    return pp;
}

}
