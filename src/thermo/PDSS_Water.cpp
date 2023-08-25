/**
 * @file PDSS_Water.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/PDSS_Water.h"
#include "cantera/thermo/WaterPropsIAPWS.h"
#include "cantera/thermo/Elements.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/global.h"

namespace Cantera
{
PDSS_Water::PDSS_Water() :
    m_waterProps(&m_sub)
{
    m_minTemp = 200.;
    m_maxTemp = 10000.;
    m_mw = 2*getElementWeight("H") + getElementWeight("O");

    // Set the baseline
    double T = 298.15;
    m_p0 = OneAtm;
    double presLow = 1.0E-2;
    double oneBar = 1.0E5;
    double dens = 1.0E-9;
    m_dens = m_sub.density(T, presLow, WATER_GAS, dens);
    m_pres = presLow;
    SW_Offset = 0.0;
    double s = entropy_mole();
    s -= GasConstant * log(oneBar/presLow);
    if (s != 188.835E3) {
        SW_Offset = 188.835E3 - s;
    }
    s = entropy_mole();
    s -= GasConstant * log(oneBar/presLow);

    double h = enthalpy_mole();
    if (h != -241.826E6) {
        EW_Offset = -241.826E6 - h;
    }
    h = enthalpy_mole();

    // Set the initial state of the system to 298.15 K and 1 bar.
    setTemperature(298.15);
    m_dens = m_sub.density(298.15, OneAtm, WATER_LIQUID);
    m_pres = OneAtm;
}

double PDSS_Water::enthalpy_mole() const
{
    return m_sub.enthalpy_mass() * m_mw + EW_Offset;
}

double PDSS_Water::intEnergy_mole() const
{
    return m_sub.intEnergy_mass() * m_mw + EW_Offset;
}

double PDSS_Water::entropy_mole() const
{
    return m_sub.entropy_mass() * m_mw + SW_Offset;
}

double PDSS_Water::gibbs_mole() const
{
    return m_sub.gibbs_mass() * m_mw+ EW_Offset - SW_Offset * m_temp;
}

double PDSS_Water::cp_mole() const
{
    return m_sub.cp_mass() * m_mw;
}

double PDSS_Water::cv_mole() const
{
    return m_sub.cv_mass() * m_mw;
}

double PDSS_Water::molarVolume() const
{
    return m_mw / m_sub.density();
}

double PDSS_Water::gibbs_RT_ref() const
{
    m_sub.density(m_temp, m_p0, m_iState);
    double h = m_sub.enthalpy_mass() * m_mw;
    m_sub.setState_TD(m_temp, m_dens);
    return (h + EW_Offset - SW_Offset * m_temp) / (m_temp * GasConstant);
}

double PDSS_Water::enthalpy_RT_ref() const
{
    m_sub.density(m_temp, m_p0, m_iState);
    double h = m_sub.enthalpy_mass() * m_mw;
    m_sub.setState_TD(m_temp, m_dens);
    return (h + EW_Offset) / (m_temp * GasConstant);
}

double PDSS_Water::entropy_R_ref() const
{
    m_sub.density(m_temp, m_p0, m_iState);
    double s = m_sub.entropy_mass() * m_mw;
    m_sub.setState_TD(m_temp, m_dens);
    return (s + SW_Offset) / GasConstant;
}

double PDSS_Water::cp_R_ref() const
{
    m_sub.density(m_temp, m_p0, m_iState);
    double cp = m_sub.cp_mass() * m_mw;
    m_sub.setState_TD(m_temp, m_dens);
    return cp / GasConstant;
}

double PDSS_Water::molarVolume_ref() const
{
    m_sub.density(m_temp, m_p0, m_iState);
    double mv = m_mw / m_sub.density();
    m_sub.setState_TD(m_temp, m_dens);
    return mv;
}

double PDSS_Water::pressure() const
{
    m_pres = m_sub.pressure();
    return m_pres;
}

void PDSS_Water::setPressure(double p)
{
    // In this routine we must be sure to only find the water branch of the
    // curve and not the gas branch
    double T = m_temp;
    double dens = m_dens;
    int waterState = WATER_LIQUID;
    if (T > m_sub.Tcrit()) {
        waterState = WATER_SUPERCRIT;
    }

    double dd = m_sub.density(T, p, waterState, dens);
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

double PDSS_Water::thermalExpansionCoeff() const
{
    return m_sub.coeffThermExp();
}

double PDSS_Water::dthermalExpansionCoeffdT() const
{
    double pres = pressure();
    double dens_save = m_dens;
    double tt = m_temp - 0.04;
    double dd = m_sub.density(tt, pres, m_iState, m_dens);
    if (dd < 0.0) {
        throw CanteraError("PDSS_Water::dthermalExpansionCoeffdT",
            "unable to solve for the density at T = {}, P = {}", tt, pres);
    }
    double vald = m_sub.coeffThermExp();
    m_sub.setState_TD(m_temp, dens_save);
    double val2 = m_sub.coeffThermExp();
    return (val2 - vald) / 0.04;
}

double PDSS_Water::isothermalCompressibility() const
{
    return m_sub.isothermalCompressibility();
}

double PDSS_Water::critTemperature() const
{
    return m_sub.Tcrit();
}

double PDSS_Water::critPressure() const
{
    return m_sub.Pcrit();
}

double PDSS_Water::critDensity() const
{
    return m_sub.Rhocrit();
}

void PDSS_Water::setDensity(double dens)
{
    m_dens = dens;
    m_sub.setState_TD(m_temp, m_dens);
}

double PDSS_Water::density() const
{
    return m_dens;
}

void PDSS_Water::setTemperature(double temp)
{
    m_temp = temp;
    m_sub.setState_TD(temp, m_dens);
}

void PDSS_Water::setState_TP(double temp, double pres)
{
    m_temp = temp;
    setPressure(pres);
}

double PDSS_Water::pref_safe(double temp) const
{
    if (temp < m_sub.Tcrit()) {
        double pp = m_sub.psat_est(temp);
        if (pp > OneAtm) {
            return pp;
        }
    } else  {
        return m_sub.Pcrit();
    }
    return OneAtm;
}

double PDSS_Water::satPressure(double t)
{
    double pp = m_sub.psat(t, WATER_LIQUID);
    m_dens = m_sub.density();
    m_temp = t;
    return pp;
}

void PDSS_Water::getParameters(AnyMap& eosNode) const
{
    eosNode["model"] = "liquid-water-IAPWS95";
}

}
