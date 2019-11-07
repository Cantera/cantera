/**
 * @file PDSS_ConstVol.cpp
 * Implementation of a pressure dependent standard state
 * virtual function.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/ctml.h"
#include "cantera/thermo/PDSS_ConstVol.h"
#include "cantera/thermo/VPStandardStateTP.h"

using namespace std;

namespace Cantera
{

PDSS_ConstVol::PDSS_ConstVol()
{
}

void PDSS_ConstVol::setParametersFromXML(const XML_Node& speciesNode)
{
    PDSS::setParametersFromXML(speciesNode);

    const XML_Node* ss = speciesNode.findByName("standardState");
    if (!ss) {
        throw CanteraError("PDSS_ConstVol::setParametersFromXML",
                           "no standardState Node for species '{}'",
                           speciesNode.name());
    }
    if (ss->attrib("model") != "constant_incompressible") {
        throw CanteraError("PDSS_ConstVol::setParametersFromXML",
                           "standardState model for species '{}' isn't "
                           "'constant_incompressible'", speciesNode.name());
    }

    setMolarVolume(getFloat(*ss, "molarVolume", "toSI"));
}

void PDSS_ConstVol::initThermo()
{
    PDSS::initThermo();
    if (m_input.hasKey("molar-volume")) {
        setMolarVolume(m_input.convert("molar-volume", "m^3/kmol"));
    }
    m_minTemp = m_spthermo->minTemp();
    m_maxTemp = m_spthermo->maxTemp();
    m_p0 = m_spthermo->refPressure();
    m_V0 = m_constMolarVolume;
    m_Vss = m_constMolarVolume;
}

doublereal PDSS_ConstVol::intEnergy_mole() const
{
    doublereal pV = (m_pres * m_Vss);
    return m_h0_RT * GasConstant * m_temp - pV;
}

doublereal PDSS_ConstVol::cv_mole() const
{
    return (cp_mole() - m_V0);
}

void PDSS_ConstVol::setPressure(doublereal p)
{
    m_pres = p;
    doublereal del_pRT = (m_pres - m_p0) / (GasConstant * m_temp);
    m_hss_RT = m_h0_RT + del_pRT * m_Vss;
    m_gss_RT = m_hss_RT - m_sss_R;
}

void PDSS_ConstVol::setTemperature(doublereal temp)
{
    m_temp = temp;
    m_spthermo->updatePropertiesTemp(temp, &m_cp0_R, &m_h0_RT, &m_s0_R);
    m_g0_RT = m_h0_RT - m_s0_R;

    doublereal del_pRT = (m_pres - m_p0) / (GasConstant * m_temp);

    m_hss_RT = m_h0_RT + del_pRT * m_Vss;
    m_cpss_R = m_cp0_R;
    m_sss_R = m_s0_R;
    m_gss_RT = m_hss_RT - m_sss_R;
}

void PDSS_ConstVol::setState_TP(doublereal temp, doublereal pres)
{
    setTemperature(temp);
    setPressure(pres);
}

void PDSS_ConstVol::setState_TR(doublereal temp, doublereal rho)
{
    doublereal rhoStored = m_mw / m_constMolarVolume;
    if (fabs(rhoStored - rho) / (rhoStored + rho) > 1.0E-4) {
        throw CanteraError("PDSS_ConstVol::setState_TR",
                           "Inconsistent supplied rho");
    }
    setTemperature(temp);
}

doublereal PDSS_ConstVol::satPressure(doublereal t)
{
    return 1.0E-200;
}

}
