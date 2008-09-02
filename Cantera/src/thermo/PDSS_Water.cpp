/**
 * @file PDSS_Water.cpp
 *
 */
/*
 * Copywrite (2006) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
/*
 * $Id$
 */
#include "ct_defs.h"
#include "xml.h"
#include "ctml.h"
#include "PDSS_Water.h"
#include "WaterPropsIAPWS.h"
#include "ThermoFactory.h"
#include <math.h>


#include "VPStandardStateTP.h"

namespace Cantera {
  /**
   * Basic list of constructors and duplicators
   */
  PDSS_Water::PDSS_Water() :
    PDSS(),
    m_sub(0),
    m_temp(0.0),
    m_dens(1000.0),
    m_iState(-3000),
    EW_Offset(0.0),
    SW_Offset(0.0),
    m_verbose(0),
    m_allowGasPhase(false)
  {
    m_pdssType = cPDSS_WATER;
    m_sub = new WaterPropsIAPWS();  
    m_spthermo = 0;
    constructSet();
  }

  PDSS_Water::PDSS_Water(VPStandardStateTP *tp, int spindex) :
    PDSS(tp, spindex),
    m_sub(0),
    m_temp(0.0),
    m_dens(1000.0),
    m_iState(-3000),
    EW_Offset(0.0),
    SW_Offset(0.0),
    m_verbose(0),
    m_allowGasPhase(false)
  {
    m_pdssType = cPDSS_WATER;
    m_sub = new WaterPropsIAPWS();
    m_spthermo = 0;
    constructSet();
  }


  PDSS_Water::PDSS_Water(VPStandardStateTP *tp, int spindex, 
		       std::string inputFile, std::string id) :
    PDSS(tp, spindex),
    m_sub(0),
    m_temp(0.0),
    m_dens(1000.0),
    m_iState(-3000),
    EW_Offset(0.0),
    SW_Offset(0.0),
    m_verbose(0),
    m_allowGasPhase(false)
  {
    m_pdssType = cPDSS_WATER;
    m_sub = new WaterPropsIAPWS();  
    constructPDSSFile(tp, spindex, inputFile, id);
    m_spthermo = 0;
  }

  PDSS_Water::PDSS_Water(VPStandardStateTP *tp, int spindex,
			 const XML_Node& speciesNode, 
			 const XML_Node& phaseRoot, bool spInstalled) :
    PDSS(tp, spindex),
    m_sub(0),
    m_temp(0.0),
    m_dens(1000.0),
    m_iState(-3000),
    EW_Offset(0.0),
    SW_Offset(0.0),
    m_verbose(0),
    m_allowGasPhase(false)
  {
    m_pdssType = cPDSS_WATER;
    m_sub = new WaterPropsIAPWS();
    std::string id= "";
    constructPDSSXML(tp, spindex, phaseRoot, id) ;
    initThermo();
    m_spthermo = 0;
  }



  PDSS_Water::PDSS_Water(const PDSS_Water &b) :
    PDSS(),
    m_sub(0),
    m_temp(0.0),
    m_dens(1000.0),
    m_iState(-3000),
    EW_Offset(b.EW_Offset),
    SW_Offset(b.SW_Offset),
    m_verbose(b.m_verbose),
    m_allowGasPhase(b.m_allowGasPhase)
  {
    m_sub = new WaterPropsIAPWS();  
    /*
     * Use the assignment operator to do the brunt
     * of the work for the copy construtor.
     */
    *this = b;
  }

  /**
   * Assignment operator
   */
  PDSS_Water& PDSS_Water::operator=(const PDSS_Water&b) {
    if (&b == this) return *this;
    /*
     * Call the base class operator
     */
    PDSS::operator=(b);

    m_sub->operator=(*(b.m_sub));
    m_temp          = b.m_temp;
    m_dens          = b.m_dens;
    m_iState        = b.m_iState;
    EW_Offset       = b.EW_Offset;
    SW_Offset       = b.SW_Offset;
    m_verbose       = b.m_verbose;
    m_allowGasPhase = b.m_allowGasPhase;
    return *this;
  }

  PDSS_Water::~PDSS_Water() { 
    delete m_sub; 
  }

  PDSS *PDSS_Water::duplMyselfAsPDSS() const {
    PDSS_Water *kPDSS = new PDSS_Water(*this);
    return (PDSS *) kPDSS;
  }

  /**
   * constructPDSSXML:
   *
   * Initialization of a Debye-Huckel phase using an
   * xml file.
   *
   * This routine is a precursor to initThermo(XML_Node*)
   * routine, which does most of the work.
   *
   * @param infile XML file containing the description of the
   *        phase
   *
   * @param id  Optional parameter identifying the name of the
   *            phase. If none is given, the first XML
   *            phase element will be used.
   */
  void PDSS_Water::constructPDSSXML(VPStandardStateTP *tp, int spindex,
				   const XML_Node& phaseNode, std::string id) {
    constructSet();
  }
   
  /**
   * constructPDSSFile():
   *
   * Initialization of a Debye-Huckel phase using an
   * xml file.
   *
   * This routine is a precursor to initThermo(XML_Node*)
   * routine, which does most of the work.
   *
   * @param infile XML file containing the description of the
   *        phase
   *
   * @param id  Optional parameter identifying the name of the
   *            phase. If none is given, the first XML
   *            phase element will be used.
   */
  void PDSS_Water::constructPDSSFile(VPStandardStateTP *tp, int spindex,
				    std::string inputFile, std::string id) {

    if (inputFile.size() == 0) {
      throw CanteraError("aterTp::initThermo",
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

    XML_Node *fxml = new XML_Node();
    fxml->build(fin);
    XML_Node *fxml_phase = findXMLPhase(fxml, id);
    if (!fxml_phase) {
      throw CanteraError("PDSS_Water::initThermo",
			 "ERROR: Can not find phase named " +
			 id + " in file named " + inputFile);
    }	
    constructPDSSXML(tp, spindex, *fxml_phase, id);
    delete fxml;
  }



  void PDSS_Water::constructSet() {
    if (m_sub) delete m_sub;
    m_sub = new WaterPropsIAPWS();
    if (m_sub == 0) {
      throw CanteraError("PDSS_Water::initThermo",
			 "could not create new substance object.");
    }
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
    doublereal dd = m_sub->density(T, presLow, WATER_GAS, dens);
    setTemperature(T);
    m_dens = dd;
    SW_Offset = 0.0;
    doublereal s = entropy_mole();
    s -=  GasConstant * log(oneBar/presLow);
    if (s != 188.835E3) {
      SW_Offset = 188.835E3 - s;
    }
    s = entropy_mole();
    s -=  GasConstant * log(oneBar/presLow);
    //printf("s = %g\n", s);

    doublereal h = enthalpy_mole();
    if (h != -241.826E6) {
      EW_Offset = -241.826E6 - h;
    }
    h = enthalpy_mole();

    //printf("h = %g\n", h);


    /*
     * Set the initial state of the system to 298.15 K and 
     * 1 bar.
     */
    setTemperature(298.15);
    doublereal rho0 = m_sub->density(298.15, OneAtm, WATER_LIQUID);
    m_dens = rho0;
  }

  void PDSS_Water::initThermo() {
    PDSS::initThermo();
  }

  void PDSS_Water::
  initThermoXML(const XML_Node& phaseNode, std::string id) {
    PDSS::initThermoXML(phaseNode, id);
  }

  doublereal PDSS_Water::
  enthalpy_mole() const {
    doublereal T = m_temp;
    doublereal dens = m_dens;
    doublereal h = m_sub->enthalpy(T, dens);
    return (h + EW_Offset);
  }

  doublereal PDSS_Water::
  intEnergy_mole() const {
    doublereal T = m_dens;
    doublereal dens = m_temp;
    doublereal u = m_sub->intEnergy(T, dens);
    return (u + EW_Offset);            
  }

  doublereal PDSS_Water::
  entropy_mole() const {
    doublereal T = m_temp;
    doublereal dens = m_dens;
    doublereal s = m_sub->entropy(T, dens);
    return (s + SW_Offset); 
  }


  doublereal PDSS_Water::
  gibbs_mole() const {
    doublereal T = m_temp;
    doublereal dens = m_dens;
    doublereal g = m_sub->Gibbs(T, dens);
    return (g + EW_Offset - SW_Offset*T);
  }


  doublereal PDSS_Water::
  cp_mole() const {
    doublereal T = m_temp;
    doublereal dens = m_dens;
    doublereal cp = m_sub->cp(T, dens);
    return cp;            
  }

 
  doublereal PDSS_Water::
  cv_mole() const {
    doublereal T = m_temp;
    doublereal dens = m_dens;
    doublereal cv = m_sub->cv(T, dens);
    return cv;
  }

  doublereal 
  PDSS_Water::molarVolume() const {
    doublereal T = m_temp;
    doublereal dens = m_dens;
    doublereal mv = m_sub->molarVolume(T, dens);
    return (mv);
  }

  doublereal
  PDSS_Water::gibbs_RT_ref() const {
    doublereal T = m_temp;
    doublereal dens0 = m_sub->density(T, m_p0);
    doublereal h = m_sub->enthalpy(T, dens0);
    m_sub->setState(m_temp, m_dens);
    return ((h + EW_Offset - SW_Offset*T)/(T * GasConstant));
  }


  doublereal
  PDSS_Water::enthalpy_RT_ref() const {
    doublereal T = m_temp;
    doublereal dens0 = m_sub->density(T, m_p0);
    doublereal h = m_sub->enthalpy(T, dens0);
    m_sub->setState(m_temp, m_dens);
    return ((h + EW_Offset)/(T * GasConstant));
  }

  doublereal PDSS_Water::
  entropy_R_ref() const {
    doublereal T = m_temp;
    doublereal dens0 = m_sub->density(T, m_p0);
    doublereal s = m_sub->entropy(T, dens0);
    m_sub->setState(m_temp, m_dens);
    return ((s + SW_Offset)/GasConstant); 
  }

  doublereal PDSS_Water::
  cp_R_ref() const {
    doublereal T = m_temp;
    doublereal dens0 = m_sub->density(T, m_p0);
    doublereal cp = m_sub->cp(T, dens0);
    m_sub->setState(m_temp, m_dens);
    return (cp/GasConstant); 
  }

  doublereal PDSS_Water::
  molarVolume_ref() const {
    doublereal T = m_temp;
    doublereal dens0 = m_sub->density(T, m_p0);
    doublereal mv = m_sub->molarVolume(T, dens0);
    m_sub->setState(m_temp, m_dens);
    return (mv); 
  }


  /**
   * Calculate the pressure (Pascals), given the temperature and density
   *  Temperature: kelvin
   *  rho: density in kg m-3
   */
  doublereal PDSS_Water::
  pressure() const {
    doublereal T = m_temp;
    doublereal dens = m_dens;
    doublereal p = m_sub->pressure(T, dens);
    m_pres = p;
    return p;
  }


  void PDSS_Water::
  setPressure(doublereal p) {
    doublereal T = m_temp;
    doublereal dens = m_dens;
    int waterState = WATER_GAS;
    doublereal rc = m_sub->Rhocrit();
    if (dens > rc) {
      waterState = WATER_LIQUID;
    }
#ifdef DEBUG_HKM
    //printf("waterPDSS: set pres = %g t = %g, waterState = %d\n",
    //      p, T, waterState);
#endif
    doublereal dd = m_sub->density(T, p, waterState, dens);
    if (dd <= 0.0) {
      std::string stateString = "T = " +
	fp2str(T) + " K and p = " + fp2str(p) + " Pa";
      throw CanteraError("PDSS_Water:setPressure()", 
			 "Failed to set water SS state: " + stateString);
    }
    m_dens = dd;
    m_pres = p;
  }
 
  // Return the volumetric thermal expansion coefficient. Units: 1/K.
  /*
   * The thermal expansion coefficient is defined as
   * \f[
   * \beta = \frac{1}{v}\left(\frac{\partial v}{\partial T}\right)_P
   * \f]
   */
  doublereal PDSS_Water::thermalExpansionCoeff() const {
    doublereal pres = pressure();
    doublereal val = m_sub->coeffThermExp(m_temp, pres);
    return val; 
  }

  doublereal PDSS_Water::dthermalExpansionCoeffdT() const {
    doublereal pres = pressure();
    double tt = m_temp - 0.04;
    doublereal vald = m_sub->coeffThermExp(tt, pres);
    doublereal val2 = m_sub->coeffThermExp(m_temp, pres);
    doublereal val = (val2 - vald) / 0.04;
    return val; 
  }
   
  doublereal PDSS_Water::isothermalCompressibility() const {
    doublereal val = m_sub->isothermalCompressibility();
    return val; 
  }

  /// critical temperature 
  doublereal PDSS_Water::critTemperature() const { return m_sub->Tcrit(); }
        
  /// critical pressure
  doublereal PDSS_Water::critPressure() const { return m_sub->Pcrit(); }
        
  /// critical density
  doublereal PDSS_Water::critDensity() const { return m_sub->Rhocrit(); }
        
  void PDSS_Water::setDensity(doublereal dens) {
    m_dens = dens;
    m_sub->setState(m_temp, m_dens);
  }

  doublereal PDSS_Water::density() const {
    return m_dens;
  }

  void PDSS_Water::setTemperature(doublereal temp) {
    m_temp = temp;
    doublereal dd = m_dens;
    m_sub->setState(temp, dd);
  }

  void PDSS_Water::setState_TP(doublereal temp, doublereal pres) {
    m_temp = temp;
    setPressure(pres);
  }

  // saturation pressure
  doublereal PDSS_Water::satPressure(doublereal t){
    doublereal pp = m_sub->psat(t);
    doublereal dens = m_dens;
    m_temp = t;
    m_dens = dens;
    return pp;
  }
        


}
