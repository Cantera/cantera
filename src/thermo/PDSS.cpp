/**
 * @file PDSS.cpp
 * Implementation of a pressure dependent standard state
 * virtual function
 * (see class @link Cantera::PDSS PDSS@endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/global.h"
#include "cantera/thermo/PDSS.h"
#include "cantera/thermo/VPStandardStateTP.h"

namespace Cantera
{

double PDSS::enthalpy_mole() const
{
    throw NotImplementedError("PDSS::enthalpy_mole");
}

double PDSS::enthalpy_RT() const
{
    throw NotImplementedError("PDSS::enthalpy_RT");
}

double PDSS::intEnergy_mole() const
{
    throw NotImplementedError("PDSS::intEnergy_mole");
}

double PDSS::entropy_mole() const
{
    throw NotImplementedError("PDSS::entropy_mole");
}

double PDSS::entropy_R() const
{
    throw NotImplementedError("PDSS::entropy_R");
}

double PDSS::gibbs_mole() const
{
    throw NotImplementedError("PDSS::gibbs_mole");
}

double PDSS::gibbs_RT() const
{
    throw NotImplementedError("PDSS::gibbs_RT");
}

double PDSS::cp_mole() const
{
    throw NotImplementedError("PDSS::cp_mole");
}

double PDSS::cp_R() const
{
    throw NotImplementedError("PDSS::cp_R");
}

double PDSS::molarVolume() const
{
    throw NotImplementedError("PDSS::molarVolume");
}

double PDSS::density() const
{
    throw NotImplementedError("PDSS::density");
}

double PDSS::cv_mole() const
{
    throw NotImplementedError("PDSS::cv_mole");
}

double PDSS::gibbs_RT_ref() const
{
    throw NotImplementedError("PDSS::gibbs_RT_ref");
}

double PDSS::enthalpy_RT_ref() const
{
    throw NotImplementedError("PDSS::enthalpy_RT_ref");
}

double PDSS::entropy_R_ref() const
{
    throw NotImplementedError("PDSS::entropy_RT_ref");
}

double PDSS::cp_R_ref() const
{
    throw NotImplementedError("PDSS::entropy_RT_ref");
}

double PDSS::molarVolume_ref() const
{
    throw NotImplementedError("PDSS::molarVolume_ref");
}

double PDSS::pressure() const
{
    return m_pres;
}

double PDSS::thermalExpansionCoeff() const
{
    throw NotImplementedError("PDSS::thermalExpansionCoeff");
}

double PDSS::critTemperature() const
{
    throw NotImplementedError("PDSS::critTemperature");
}

double PDSS::critPressure() const
{
    throw NotImplementedError("PDSS::critPressure");
}

double PDSS::critDensity() const
{
    throw NotImplementedError("PDSS::critDensity");
}

void PDSS::setPressure(double pres)
{
    m_pres = pres;
}

double PDSS::temperature() const
{
    return m_temp;
}

void PDSS::setTemperature(double temp)
{
    m_temp = temp;
}

double PDSS::molecularWeight() const
{
    return m_mw;
}
void PDSS::setMolecularWeight(double mw)
{
    m_mw = mw;
}

void PDSS::setState_TP(double temp, double pres)
{
    throw NotImplementedError("PDSS::setState_TP");
}

double PDSS::satPressure(double t)
{
    throw NotImplementedError("PDSS::satPressure");
}

// PDSS_Molar methods

double PDSS_Molar::enthalpy_RT() const
{
    return enthalpy_mole() / (GasConstant * temperature());
}

double PDSS_Molar::entropy_R() const
{
    return entropy_mole() / GasConstant;
}

double PDSS_Molar::gibbs_RT() const
{
    return gibbs_mole() / (GasConstant * temperature());
}

double PDSS_Molar::cp_R() const
{
    return cp_mole() / GasConstant;
}

// PDSS_Nondimensional methods

PDSS_Nondimensional::PDSS_Nondimensional()
    : m_h0_RT(0.0)
    , m_cp0_R(0.0)
    , m_s0_R(0.0)
    , m_g0_RT(0.0)
    , m_V0(0.0)
    , m_hss_RT(0.0)
    , m_cpss_R(0.0)
    , m_sss_R(0.0)
    , m_gss_RT(0.0)
    , m_Vss(0.0)
{
}

double PDSS_Nondimensional::enthalpy_mole() const
{
    return enthalpy_RT() * GasConstant * temperature();
}

double PDSS_Nondimensional::entropy_mole() const
{
    return entropy_R() * GasConstant;
}

double PDSS_Nondimensional::gibbs_mole() const
{
    return gibbs_RT() * GasConstant * temperature();
}

double PDSS_Nondimensional::cp_mole() const
{
    return cp_R() * GasConstant;
}

double PDSS_Nondimensional::gibbs_RT_ref() const
{
    return m_g0_RT;
}

double PDSS_Nondimensional::enthalpy_RT_ref() const
{
    return m_h0_RT;
}

double PDSS_Nondimensional::entropy_R_ref() const
{
    return m_s0_R;
}

double PDSS_Nondimensional::cp_R_ref() const
{
    return m_cp0_R;
}

double PDSS_Nondimensional::molarVolume_ref() const
{
    return m_V0;
}

double PDSS_Nondimensional::enthalpy_RT() const
{
    return m_hss_RT;
}

double PDSS_Nondimensional::entropy_R() const
{
    return m_sss_R;
}

double PDSS_Nondimensional::gibbs_RT() const
{
    return m_gss_RT;
}

double PDSS_Nondimensional::cp_R() const
{
    return m_cpss_R;
}

double PDSS_Nondimensional::molarVolume() const
{
    return m_Vss;
}

double PDSS_Nondimensional::density() const
{
    return m_mw / m_Vss;
}

}
