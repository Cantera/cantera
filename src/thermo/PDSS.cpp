/**
 * @file PDSS.cpp
 * Implementation of a pressure dependent standard state
 * virtual function
 * (see class \link Cantera::PDSS PDSS\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/ctml.h"
#include "cantera/thermo/PDSS.h"
#include "cantera/thermo/VPStandardStateTP.h"

namespace Cantera
{
PDSS::PDSS() :
    m_temp(-1.0),
    m_pres(-1.0),
    m_p0(-1.0),
    m_minTemp(-1.0),
    m_maxTemp(10000.0),
    m_mw(0.0)
{
}

doublereal PDSS::enthalpy_mole() const
{
    throw NotImplementedError("PDSS::enthalpy_mole");
}

doublereal PDSS::enthalpy_RT() const
{
    throw NotImplementedError("PDSS::enthalpy_RT");
}

doublereal PDSS::intEnergy_mole() const
{
    throw NotImplementedError("PDSS::intEnergy_mole");
}

doublereal PDSS::entropy_mole() const
{
    throw NotImplementedError("PDSS::entropy_mole");
}

doublereal PDSS::entropy_R() const
{
    throw NotImplementedError("PDSS::entropy_R");
}

doublereal PDSS::gibbs_mole() const
{
    throw NotImplementedError("PDSS::gibbs_mole");
}

doublereal PDSS::gibbs_RT() const
{
    throw NotImplementedError("PDSS::gibbs_RT");
}

doublereal PDSS::cp_mole() const
{
    throw NotImplementedError("PDSS::cp_mole");
}

doublereal PDSS::cp_R() const
{
    throw NotImplementedError("PDSS::cp_R");
}

doublereal PDSS::molarVolume() const
{
    throw NotImplementedError("PDSS::molarVolume");
}

doublereal PDSS::density() const
{
    throw NotImplementedError("PDSS::density");
}

doublereal PDSS::cv_mole() const
{
    throw NotImplementedError("PDSS::cv_mole");
}

doublereal PDSS::gibbs_RT_ref() const
{
    throw NotImplementedError("PDSS::gibbs_RT_ref");
}

doublereal PDSS::enthalpy_RT_ref() const
{
    throw NotImplementedError("PDSS::enthalpy_RT_ref");
}

doublereal PDSS::entropy_R_ref() const
{
    throw NotImplementedError("PDSS::entropy_RT_ref");
}

doublereal PDSS::cp_R_ref() const
{
    throw NotImplementedError("PDSS::entropy_RT_ref");
}

doublereal PDSS::molarVolume_ref() const
{
    throw NotImplementedError("PDSS::molarVolume_ref");
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
    throw NotImplementedError("PDSS::thermalExpansionCoeff");
}

doublereal PDSS::critTemperature() const
{
    throw NotImplementedError("PDSS::critTemperature");
}

doublereal PDSS::critPressure() const
{
    throw NotImplementedError("PDSS::critPressure");
}

doublereal PDSS::critDensity() const
{
    throw NotImplementedError("PDSS::critDensity");
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
    throw NotImplementedError("PDSS::setState_TP");
}

void PDSS::setState_TR(doublereal temp, doublereal rho)
{
    throw NotImplementedError("PDSS::setState_TR");
}

doublereal PDSS::satPressure(doublereal t)
{
    throw NotImplementedError("PDSS::satPressure");
}

void PDSS::reportParams(size_t& kindex, int& type,
                        doublereal* const c,
                        doublereal& minTemp_,
                        doublereal& maxTemp_,
                        doublereal& refPressure_) const
{
    kindex = npos;
    type = 0;
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
