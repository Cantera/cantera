/**
 * @file PDSS_Water.cpp
 *
 */
/*
 * Copyright (2006) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
#include "cantera/base/ctml.h"
#include "cantera/thermo/PDSS_Water.h"
#include "cantera/thermo/WaterPropsIAPWS.h"
#include "cantera/base/stringUtils.h"

#include <fstream>

namespace Cantera
{
PDSS_Water::PDSS_Water() :
    m_waterProps(&m_sub),
    m_dens(1000.0),
    m_iState(WATER_LIQUID),
    EW_Offset(0.0),
    SW_Offset(0.0),
    m_verbose(0),
    m_allowGasPhase(false)
{
    m_pdssType = cPDSS_WATER;
    m_spthermo = 0;
    constructSet();
    m_minTemp = 200.;
    m_maxTemp = 10000.;
}

PDSS_Water::PDSS_Water(VPStandardStateTP* tp, int spindex) :
    PDSS(tp, spindex),
    m_waterProps(&m_sub),
    m_dens(1000.0),
    m_iState(WATER_LIQUID),
    EW_Offset(0.0),
    SW_Offset(0.0),
    m_verbose(0),
    m_allowGasPhase(false)
{
    m_pdssType = cPDSS_WATER;
    m_spthermo = 0;
    constructSet();
    m_minTemp = 200.;
    m_maxTemp = 10000.;
}

PDSS_Water::PDSS_Water(VPStandardStateTP* tp, int spindex,
                       const std::string& inputFile, const std::string& id) :
    PDSS(tp, spindex),
    m_waterProps(&m_sub),
    m_dens(1000.0),
    m_iState(WATER_LIQUID),
    EW_Offset(0.0),
    SW_Offset(0.0),
    m_verbose(0),
    m_allowGasPhase(false)
{
    m_pdssType = cPDSS_WATER;
    constructPDSSFile(tp, spindex, inputFile, id);
    m_spthermo = 0;
    m_minTemp = 200.;
    m_maxTemp = 10000.;
}

PDSS_Water::PDSS_Water(VPStandardStateTP* tp, int spindex,
                       const XML_Node& speciesNode,
                       const XML_Node& phaseRoot, bool spInstalled) :
    PDSS(tp, spindex),
    m_waterProps(&m_sub),
    m_dens(1000.0),
    m_iState(WATER_LIQUID),
    EW_Offset(0.0),
    SW_Offset(0.0),
    m_verbose(0),
    m_allowGasPhase(false)
{
    m_pdssType = cPDSS_WATER;
    std::string id= "";
    constructPDSSXML(tp, spindex, phaseRoot, id) ;
    initThermo();
    m_spthermo = 0;
    m_minTemp = 200.;
    m_maxTemp = 10000.;
}

PDSS_Water::PDSS_Water(const PDSS_Water& b) :
    PDSS(),
    m_waterProps(&m_sub),
    m_dens(1000.0),
    m_iState(WATER_LIQUID),
    EW_Offset(b.EW_Offset),
    SW_Offset(b.SW_Offset),
    m_verbose(b.m_verbose),
    m_allowGasPhase(b.m_allowGasPhase)
{
    /*
     * Use the assignment operator to do the brunt
     * of the work for the copy constructor.
     */
    *this = b;
}

PDSS_Water& PDSS_Water::operator=(const PDSS_Water& b)
{
    if (&b == this) {
        return *this;
    }
    /*
     * Call the base class operator
     */
    PDSS::operator=(b);

    m_sub = b.m_sub;
    m_waterProps = b.m_waterProps;
    m_dens          = b.m_dens;
    m_iState        = b.m_iState;
    EW_Offset       = b.EW_Offset;
    SW_Offset       = b.SW_Offset;
    m_verbose       = b.m_verbose;
    m_allowGasPhase = b.m_allowGasPhase;

    return *this;
}

PDSS* PDSS_Water::duplMyselfAsPDSS() const
{
    return new PDSS_Water(*this);
}

void PDSS_Water::constructPDSSXML(VPStandardStateTP* tp, int spindex,
                                  const XML_Node& phaseNode, const std::string& id)
{
    constructSet();
}

void PDSS_Water::constructPDSSFile(VPStandardStateTP* tp, int spindex,
                                   const std::string& inputFile, const std::string& id)
{
    if (inputFile.size() == 0) {
        throw CanteraError("PDSS_Water::constructPDSSFile",
                           "input file is null");
    }
    std::string path = findInputFile(inputFile);
    std::ifstream fin(path.c_str());
    if (!fin) {
        throw CanteraError("PDSS_Water::initThermo","could not open "
                           +path+" for reading.");
    }
    /*
     * The phase object automatically constructs an XML object.
     * Use this object to store information.
     */

    XML_Node fxml;
    fxml.build(fin);
    XML_Node* fxml_phase = findXMLPhase(&fxml, id);
    if (!fxml_phase) {
        throw CanteraError("PDSS_Water::initThermo",
                           "ERROR: Can not find phase named " +
                           id + " in file named " + inputFile);
    }
    constructPDSSXML(tp, spindex, *fxml_phase, id);
}

void PDSS_Water::constructSet()
{
    /*
     * Calculate the molecular weight.
     *  hard coded to Cantera's elements and Water.
     */
    m_mw = 2 * 1.00794 + 15.9994;

    /*
     * Set the baseline
     */
    doublereal T = 298.15;

    m_p0 = OneAtm;

    doublereal presLow = 1.0E-2;
    doublereal oneBar = 1.0E5;
    doublereal dens = 1.0E-9;
    m_dens = m_sub.density(T, presLow, WATER_GAS, dens);
    m_pres = presLow;
    SW_Offset = 0.0;
    doublereal s = entropy_mole();
    s -=  GasConstant * log(oneBar/presLow);
    if (s != 188.835E3) {
        SW_Offset = 188.835E3 - s;
    }
    s = entropy_mole();
    s -=  GasConstant * log(oneBar/presLow);

    doublereal h = enthalpy_mole();
    if (h != -241.826E6) {
        EW_Offset = -241.826E6 - h;
    }
    h = enthalpy_mole();

    /*
     * Set the initial state of the system to 298.15 K and
     * 1 bar.
     */
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

doublereal  PDSS_Water::molarVolume() const
{
    return m_sub.molarVolume();
}

doublereal PDSS_Water::gibbs_RT_ref() const
{
    doublereal T = m_temp;
    m_sub.density(T, m_p0);
    doublereal h = m_sub.enthalpy();
    m_sub.setState_TR(m_temp, m_dens);
    return (h + EW_Offset - SW_Offset*T)/(T * GasConstant);
}

doublereal PDSS_Water::enthalpy_RT_ref() const
{
    doublereal T = m_temp;
    m_sub.density(T, m_p0);
    doublereal h = m_sub.enthalpy();
    m_sub.setState_TR(m_temp, m_dens);
    return (h + EW_Offset)/(T * GasConstant);
}

doublereal PDSS_Water::entropy_R_ref() const
{
    doublereal T = m_temp;
    m_sub.density(T, m_p0);
    doublereal s = m_sub.entropy();
    m_sub.setState_TR(m_temp, m_dens);
    return (s + SW_Offset)/GasConstant;
}

doublereal PDSS_Water::cp_R_ref() const
{
    doublereal T = m_temp;
    m_sub.density(T, m_p0);
    doublereal cp = m_sub.cp();
    m_sub.setState_TR(m_temp, m_dens);
    return cp/GasConstant;
}

doublereal PDSS_Water::molarVolume_ref() const
{
    doublereal T = m_temp;
    m_sub.density(T, m_p0);
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
        std::string stateString = "T = " +
                                  fp2str(T) + " K and p = " + fp2str(p) + " Pa";
        throw CanteraError("PDSS_Water:setPressure()",
                           "Failed to set water SS state: " + stateString);
    }
    m_dens = dd;
    m_pres = p;

    // We are only putting the phase check here because of speed considerations.
    m_iState = m_sub.phaseState(true);
    if (! m_allowGasPhase) {
        if (m_iState != WATER_SUPERCRIT && m_iState != WATER_LIQUID && m_iState != WATER_UNSTABLELIQUID) {
            throw CanteraError("PDSS_Water::setPressure",
                               "Water State isn't liquid or crit");
        }
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
                           "unable to solve for the density at T = " + fp2str(tt) + ", P = " + fp2str(pres));
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
