/**
 * @file PDSS_SSVol.cpp
 * Implementation of a pressure dependent standard state
 * virtual function.
 */
/*
 * Copyright (2006) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
#include "cantera/base/ct_defs.h"
#include "cantera/base/xml.h"
#include "cantera/base/ctml.h"
#include "cantera/thermo/PDSS_SSVol.h"
#include "cantera/thermo/ThermoFactory.h"

#include "cantera/thermo/VPStandardStateTP.h"

#include <fstream>

using namespace std;

namespace Cantera
{
PDSS_SSVol::PDSS_SSVol(VPStandardStateTP* tp, size_t spindex) :
    PDSS(tp, spindex),
    volumeModel_(cSSVOLUME_CONSTANT),
    m_constMolarVolume(-1.0)
{
    m_pdssType = cPDSS_SSVOL;
    TCoeff_[0] = 0.0;
    TCoeff_[1] = 0.0;
    TCoeff_[2] = 0.0;
}

PDSS_SSVol::PDSS_SSVol(VPStandardStateTP* tp,
                       size_t spindex, const std::string& inputFile, const std::string& id) :
    PDSS(tp, spindex),
    volumeModel_(cSSVOLUME_CONSTANT),
    m_constMolarVolume(-1.0)
{

    m_pdssType = cPDSS_SSVOL;
    constructPDSSFile(tp, spindex, inputFile, id);
}

PDSS_SSVol::PDSS_SSVol(VPStandardStateTP* tp, size_t spindex,
                       const XML_Node& speciesNode,
                       const XML_Node& phaseRoot,
                       bool spInstalled) :
    PDSS(tp, spindex),
    volumeModel_(cSSVOLUME_CONSTANT),
    m_constMolarVolume(-1.0)
{
    m_pdssType = cPDSS_SSVOL;
    constructPDSSXML(tp, spindex, speciesNode,  phaseRoot, spInstalled) ;
}

PDSS_SSVol::PDSS_SSVol(const PDSS_SSVol& b) :
    PDSS(b),
    volumeModel_(cSSVOLUME_CONSTANT),
    m_constMolarVolume(-1.0)
{
    /*
     * Use the assignment operator to do the brunt
     * of the work for the copy constructor.
     */
    *this = b;
}

PDSS_SSVol& PDSS_SSVol::operator=(const PDSS_SSVol& b)
{
    if (&b == this) {
        return *this;
    }
    PDSS::operator=(b);
    volumeModel_ = b.volumeModel_;
    m_constMolarVolume = b.m_constMolarVolume;
    TCoeff_ = b.TCoeff_;
    return *this;
}

PDSS* PDSS_SSVol::duplMyselfAsPDSS() const
{
    return new PDSS_SSVol(*this);
}

void PDSS_SSVol::constructPDSSXML(VPStandardStateTP* tp, size_t spindex,
                                  const XML_Node& speciesNode,
                                  const XML_Node& phaseNode, bool spInstalled)
{
    PDSS::initThermo();
    SpeciesThermo& sp = m_tp->speciesThermo();
    m_p0 = sp.refPressure(m_spindex);

    if (!spInstalled) {
        throw CanteraError("PDSS_SSVol::constructPDSSXML", "spInstalled false not handled");
    }

    const XML_Node* ss = speciesNode.findByName("standardState");
    if (!ss) {
        throw CanteraError("PDSS_SSVol::constructPDSSXML",
                           "no standardState Node for species " + speciesNode.name());
    }
    std::string model = (*ss)["model"];
    if (model == "constant_incompressible" || model == "constant") {
        volumeModel_ = cSSVOLUME_CONSTANT;
        m_constMolarVolume = ctml::getFloat(*ss, "molarVolume", "toSI");
    } else if (model == "temperature_polynomial") {
        volumeModel_ = cSSVOLUME_TPOLY;
        size_t num = ctml::getFloatArray(*ss, TCoeff_, true, "toSI", "volumeTemperaturePolynomial");
        if (num != 4) {
            throw CanteraError("PDSS_SSVol::constructPDSSXML",
                               " Didn't get 4 density polynomial numbers for species " + speciesNode.name());
        }
    } else if (model == "density_temperature_polynomial") {
        volumeModel_ = cSSVOLUME_DENSITY_TPOLY;
        size_t num = ctml::getFloatArray(*ss, TCoeff_, true, "toSI", "densityTemperaturePolynomial");
        if (num != 4) {
            throw CanteraError("PDSS_SSVol::constructPDSSXML",
                               " Didn't get 4 density polynomial numbers for species " + speciesNode.name());
        }
    } else {
        throw CanteraError("PDSS_SSVol::constructPDSSXML",
                           "standardState model for species isn't constant_incompressible: " + speciesNode.name());
    }
    std::string id = "";
}

void PDSS_SSVol::constructPDSSFile(VPStandardStateTP* tp, size_t spindex,
                                   const std::string& inputFile, const std::string& id)
{
    if (inputFile.size() == 0) {
        throw CanteraError("PDSS_SSVol::initThermo",
                           "input file is null");
    }
    std::string path = findInputFile(inputFile);
    ifstream fin(path.c_str());
    if (!fin) {
        throw CanteraError("PDSS_SSVol::initThermo","could not open "
                           +path+" for reading.");
    }
    /*
     * The phase object automatically constructs an XML object.
     * Use this object to store information.
     */

    XML_Node* fxml = new XML_Node();
    fxml->build(fin);
    XML_Node* fxml_phase = findXMLPhase(fxml, id);
    if (!fxml_phase) {
        throw CanteraError("PDSS_SSVol::initThermo",
                           "ERROR: Can not find phase named " +
                           id + " in file named " + inputFile);
    }

    XML_Node& speciesList = fxml_phase->child("speciesArray");
    XML_Node* speciesDB = get_XML_NameID("speciesData", speciesList["datasrc"],
                                         &(fxml_phase->root()));
    const vector<string>&sss = tp->speciesNames();
    const XML_Node* s =  speciesDB->findByAttr("name", sss[spindex]);

    constructPDSSXML(tp, spindex, *s, *fxml_phase, true);
    delete fxml;
}

void PDSS_SSVol::initThermoXML(const XML_Node& phaseNode, const std::string& id)
{
    PDSS::initThermoXML(phaseNode, id);
    m_minTemp = m_spthermo->minTemp(m_spindex);
    m_maxTemp = m_spthermo->maxTemp(m_spindex);
    m_p0 = m_spthermo->refPressure(m_spindex);
    m_mw = m_tp->molecularWeight(m_spindex);
}

void PDSS_SSVol::initThermo()
{
    PDSS::initThermo();
    SpeciesThermo& sp = m_tp->speciesThermo();
    m_p0 = sp.refPressure(m_spindex);
    m_V0_ptr[m_spindex] = m_constMolarVolume;
    m_Vss_ptr[m_spindex] = m_constMolarVolume;
}

doublereal
PDSS_SSVol::enthalpy_mole() const
{
    doublereal val = enthalpy_RT();
    doublereal RT = GasConstant * m_temp;
    return val * RT;
}

doublereal
PDSS_SSVol::enthalpy_RT() const
{
    return m_hss_RT_ptr[m_spindex];
}

doublereal
PDSS_SSVol::intEnergy_mole() const
{
    doublereal pVRT = (m_pres * m_Vss_ptr[m_spindex]) / (GasConstant * m_temp);
    doublereal val = m_h0_RT_ptr[m_spindex] - pVRT;
    doublereal RT = GasConstant * m_temp;
    return val * RT;
}

doublereal
PDSS_SSVol::entropy_mole() const
{
    doublereal val = entropy_R();
    return val * GasConstant;
}

doublereal
PDSS_SSVol::entropy_R() const
{
    return m_sss_R_ptr[m_spindex];
}

doublereal
PDSS_SSVol::gibbs_mole() const
{
    doublereal val = gibbs_RT();
    doublereal RT = GasConstant * m_temp;
    return val * RT;
}

doublereal
PDSS_SSVol::gibbs_RT() const
{
    return m_gss_RT_ptr[m_spindex];
}

doublereal
PDSS_SSVol::cp_mole() const
{
    doublereal val = m_cpss_R_ptr[m_spindex];
    return val * GasConstant;
}

doublereal
PDSS_SSVol::cp_R() const
{
    return m_cpss_R_ptr[m_spindex];
}

doublereal
PDSS_SSVol::cv_mole() const
{
    return (cp_mole() -  m_V0_ptr[m_spindex]);
}

doublereal
PDSS_SSVol::molarVolume() const
{
    return m_Vss_ptr[m_spindex];
}

doublereal
PDSS_SSVol::density() const
{
    doublereal val = m_Vss_ptr[m_spindex];
    return m_mw/val;
}

doublereal
PDSS_SSVol::gibbs_RT_ref() const
{
    return m_g0_RT_ptr[m_spindex];
}

doublereal PDSS_SSVol::enthalpy_RT_ref() const
{
    return m_h0_RT_ptr[m_spindex];
}

doublereal PDSS_SSVol::entropy_R_ref() const
{
    return m_s0_R_ptr[m_spindex];
}

doublereal PDSS_SSVol::cp_R_ref() const
{
    return m_cp0_R_ptr[m_spindex];
}

doublereal PDSS_SSVol::molarVolume_ref() const
{
    return m_V0_ptr[m_spindex];
}

void PDSS_SSVol::calcMolarVolume() const
{
    if (volumeModel_ == cSSVOLUME_CONSTANT) {
        m_Vss_ptr[m_spindex] = m_constMolarVolume;
    } else if (volumeModel_ == cSSVOLUME_TPOLY) {
        m_Vss_ptr[m_spindex] = TCoeff_[0] + m_temp * (TCoeff_[1] + m_temp * (TCoeff_[2] + m_temp * TCoeff_[3]));
        dVdT_   = TCoeff_[1] + 2.0 * m_temp * TCoeff_[2] + 3.0 * m_temp * m_temp * TCoeff_[3];
        d2VdT2_ =  2.0 * TCoeff_[2] + 6.0 * m_temp * TCoeff_[3];
    } else  if (volumeModel_ == cSSVOLUME_DENSITY_TPOLY) {
        doublereal dens =  TCoeff_[0] + m_temp * (TCoeff_[1] + m_temp * (TCoeff_[2] + m_temp * TCoeff_[3]));
        m_Vss_ptr[m_spindex] = m_mw / dens;
        doublereal dens2 = dens * dens;
        doublereal ddensdT =  TCoeff_[1] + 2.0 * m_temp * TCoeff_[2] + 3.0 * m_temp * m_temp * TCoeff_[3];
        doublereal d2densdT2 = 2.0 * TCoeff_[2] + 6.0 * m_temp * TCoeff_[3];
        dVdT_   = - m_mw / (dens2) * (ddensdT);
        d2VdT2_ = 2.0 * m_mw / (dens2 * dens) * ddensdT * ddensdT - m_mw / dens2 * d2densdT2;
    } else {
        throw CanteraError("PDSS_SSVol::calcMolarVolume", "unimplemented");
    }
}

doublereal PDSS_SSVol::critTemperature() const
{
    throw CanteraError("PDSS_SSVol::critTemperature()", "unimplemented");
    return 0.0;
}

doublereal PDSS_SSVol::critPressure() const
{
    throw CanteraError("PDSS_SSVol::critPressure()", "unimplemented");
    return 0.0;
}

doublereal PDSS_SSVol::critDensity() const
{
    throw CanteraError("PDSS_SSVol::critDensity()", "unimplemented");
    return 0.0;
}

void PDSS_SSVol::setPressure(doublereal p)
{
    m_pres = p;
    doublereal deltaP = m_pres - m_p0;
    if (fabs(deltaP) < 1.0E-10) {
        m_hss_RT_ptr[m_spindex] = m_h0_RT_ptr[m_spindex];
        m_sss_R_ptr[m_spindex] = m_s0_R_ptr[m_spindex];
        m_gss_RT_ptr[m_spindex] = m_hss_RT_ptr[m_spindex] - m_sss_R_ptr[m_spindex];
        m_cpss_R_ptr[m_spindex] = m_cp0_R_ptr[m_spindex];
    } else {
        doublereal del_pRT = deltaP / (GasConstant * m_temp);
        doublereal sV_term =  - deltaP / (GasConstant) * dVdT_;
        m_hss_RT_ptr[m_spindex] = m_h0_RT_ptr[m_spindex] + sV_term + del_pRT * (m_Vss_ptr[m_spindex]);
        m_sss_R_ptr[m_spindex] = m_s0_R_ptr[m_spindex] + sV_term;
        m_gss_RT_ptr[m_spindex] = m_hss_RT_ptr[m_spindex] - m_sss_R_ptr[m_spindex];
        m_cpss_R_ptr[m_spindex] = m_cp0_R_ptr[m_spindex] - m_temp * deltaP * d2VdT2_;
    }
}

void PDSS_SSVol::setTemperature(doublereal temp)
{
    m_temp = temp;
    m_spthermo->update_one(m_spindex, temp, m_cp0_R_ptr, m_h0_RT_ptr, m_s0_R_ptr);
    calcMolarVolume();
    m_g0_RT_ptr[m_spindex] =  m_h0_RT_ptr[m_spindex] -  m_s0_R_ptr[m_spindex];
    doublereal deltaP = m_pres - m_p0;
    if (fabs(deltaP) < 1.0E-10) {
        m_hss_RT_ptr[m_spindex] = m_h0_RT_ptr[m_spindex];
        m_sss_R_ptr[m_spindex] = m_s0_R_ptr[m_spindex];
        m_gss_RT_ptr[m_spindex] = m_hss_RT_ptr[m_spindex] - m_sss_R_ptr[m_spindex];
        m_cpss_R_ptr[m_spindex] = m_cp0_R_ptr[m_spindex];
    } else {
        doublereal del_pRT = deltaP / (GasConstant * m_temp);
        doublereal sV_term =  - deltaP / (GasConstant) * dVdT_;
        m_hss_RT_ptr[m_spindex] = m_h0_RT_ptr[m_spindex] + sV_term + del_pRT * (m_Vss_ptr[m_spindex]);
        m_sss_R_ptr[m_spindex] = m_s0_R_ptr[m_spindex] + sV_term;
        m_gss_RT_ptr[m_spindex] = m_hss_RT_ptr[m_spindex] - m_sss_R_ptr[m_spindex];
        m_cpss_R_ptr[m_spindex] = m_cp0_R_ptr[m_spindex] - m_temp * deltaP * d2VdT2_;
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
