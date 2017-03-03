/**
 * @file PDSS_SSVol.cpp
 * Implementation of a pressure dependent standard state
 * virtual function.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/base/ctml.h"
#include "cantera/thermo/PDSS_SSVol.h"
#include "cantera/thermo/VPStandardStateTP.h"

using namespace std;

namespace Cantera
{

PDSS_SSVol::PDSS_SSVol()
    : volumeModel_(SSVolume_Model::constant)
    , m_constMolarVolume(-1.0)
{
    TCoeff_[0] = 0.0;
    TCoeff_[1] = 0.0;
    TCoeff_[2] = 0.0;
}

void PDSS_SSVol::setParametersFromXML(const XML_Node& speciesNode)
{
    PDSS::setParametersFromXML(speciesNode);

    const XML_Node* ss = speciesNode.findByName("standardState");
    if (!ss) {
        throw CanteraError("PDSS_SSVol::constructPDSSXML",
                           "no standardState Node for species " + speciesNode.name());
    }
    std::string model = ss->attrib("model");
    if (model == "constant_incompressible" || model == "constant") {
        volumeModel_ = SSVolume_Model::constant;
        m_constMolarVolume = getFloat(*ss, "molarVolume", "toSI");
    } else if (model == "temperature_polynomial") {
        volumeModel_ = SSVolume_Model::tpoly;
        size_t num = getFloatArray(*ss, TCoeff_, true, "toSI", "volumeTemperaturePolynomial");
        if (num != 4) {
            throw CanteraError("PDSS_SSVol::constructPDSSXML",
                               " Didn't get 4 density polynomial numbers for species " + speciesNode.name());
        }
    } else if (model == "density_temperature_polynomial") {
        volumeModel_ = SSVolume_Model::density_tpoly;
        size_t num = getFloatArray(*ss, TCoeff_, true, "toSI", "densityTemperaturePolynomial");
        if (num != 4) {
            throw CanteraError("PDSS_SSVol::constructPDSSXML",
                               " Didn't get 4 density polynomial numbers for species " + speciesNode.name());
        }
    } else {
        throw CanteraError("PDSS_SSVol::constructPDSSXML",
                           "standardState model for species isn't constant_incompressible: " + speciesNode.name());
    }
}

void PDSS_SSVol::initThermo()
{
    PDSS::initThermo();
    m_minTemp = m_spthermo->minTemp();
    m_maxTemp = m_spthermo->maxTemp();
    m_p0 = m_spthermo->refPressure();
    m_V0 = m_constMolarVolume;
    m_Vss = m_constMolarVolume;
}

doublereal PDSS_SSVol::intEnergy_mole() const
{
    doublereal pV = m_pres * m_Vss;
    return m_h0_RT * GasConstant * m_temp - pV;
}

doublereal PDSS_SSVol::cv_mole() const
{
    return (cp_mole() - m_V0);
}

void PDSS_SSVol::calcMolarVolume()
{
    if (volumeModel_ == SSVolume_Model::constant) {
        m_Vss = m_constMolarVolume;
    } else if (volumeModel_ == SSVolume_Model::tpoly) {
        m_Vss = TCoeff_[0] + m_temp * (TCoeff_[1] + m_temp * (TCoeff_[2] + m_temp * TCoeff_[3]));
        dVdT_ = TCoeff_[1] + 2.0 * m_temp * TCoeff_[2] + 3.0 * m_temp * m_temp * TCoeff_[3];
        d2VdT2_ = 2.0 * TCoeff_[2] + 6.0 * m_temp * TCoeff_[3];
    } else if (volumeModel_ == SSVolume_Model::density_tpoly) {
        doublereal dens = TCoeff_[0] + m_temp * (TCoeff_[1] + m_temp * (TCoeff_[2] + m_temp * TCoeff_[3]));
        m_Vss = m_mw / dens;
        doublereal dens2 = dens * dens;
        doublereal ddensdT = TCoeff_[1] + 2.0 * m_temp * TCoeff_[2] + 3.0 * m_temp * m_temp * TCoeff_[3];
        doublereal d2densdT2 = 2.0 * TCoeff_[2] + 6.0 * m_temp * TCoeff_[3];
        dVdT_ = - m_mw / dens2 * ddensdT;
        d2VdT2_ = 2.0 * m_mw / (dens2 * dens) * ddensdT * ddensdT - m_mw / dens2 * d2densdT2;
    } else {
        throw CanteraError("PDSS_SSVol::calcMolarVolume", "unimplemented");
    }
}

void PDSS_SSVol::setPressure(doublereal p)
{
    m_pres = p;
    doublereal deltaP = m_pres - m_p0;
    if (fabs(deltaP) < 1.0E-10) {
        m_hss_RT = m_h0_RT;
        m_sss_R = m_s0_R;
        m_gss_RT = m_hss_RT - m_sss_R;
        m_cpss_R = m_cp0_R;
    } else {
        doublereal del_pRT = deltaP / (GasConstant * m_temp);
        doublereal sV_term = - deltaP / GasConstant * dVdT_;
        m_hss_RT = m_h0_RT + sV_term + del_pRT * m_Vss;
        m_sss_R = m_s0_R + sV_term;
        m_gss_RT = m_hss_RT - m_sss_R;
        m_cpss_R = m_cp0_R - m_temp * deltaP * d2VdT2_;
    }
}

void PDSS_SSVol::setTemperature(doublereal temp)
{
    m_temp = temp;
    m_spthermo->updatePropertiesTemp(temp, &m_cp0_R, &m_h0_RT, &m_s0_R);
    calcMolarVolume();
    m_g0_RT = m_h0_RT - m_s0_R;
    doublereal deltaP = m_pres - m_p0;
    if (fabs(deltaP) < 1.0E-10) {
        m_hss_RT = m_h0_RT;
        m_sss_R = m_s0_R;
        m_gss_RT = m_hss_RT - m_sss_R;
        m_cpss_R = m_cp0_R;
    } else {
        doublereal del_pRT = deltaP / (GasConstant * m_temp);
        doublereal sV_term = - deltaP / GasConstant * dVdT_;
        m_hss_RT = m_h0_RT + sV_term + del_pRT * m_Vss;
        m_sss_R = m_s0_R + sV_term;
        m_gss_RT = m_hss_RT - m_sss_R;
        m_cpss_R = m_cp0_R - m_temp * deltaP * d2VdT2_;
    }
}

void PDSS_SSVol::setState_TP(doublereal temp, doublereal pres)
{
    m_pres = pres;
    setTemperature(temp);
}

void PDSS_SSVol::setState_TR(doublereal temp, doublereal rho)
{
    doublereal rhoStored = m_mw / m_constMolarVolume;
    if (fabs(rhoStored - rho) / (rhoStored + rho) > 1.0E-4) {
        throw CanteraError("PDSS_SSVol::setState_TR",
                           "Inconsistent supplied rho");
    }
    setTemperature(temp);
}

doublereal PDSS_SSVol::satPressure(doublereal t)
{
    return 1.0E-200;
}

}
