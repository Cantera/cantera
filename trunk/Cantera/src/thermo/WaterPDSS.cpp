/**
 * @file WaterPDSS.cpp
 *
 */
/*
 * Copywrite (2006) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
/*
 * $Id: WaterPDSS.cpp,v 1.1 2006/07/04 00:01:53 hkmoffa Exp $
 */
#include "ct_defs.h"
#include "xml.h"
#include "ctml.h"
#include "WaterPDSS.h"
#include "WaterPropsIAPWS.h"
#include "importCTML.h"


#include "ThermoPhase.h"

namespace Cantera {
  /**
   * Basic list of constructors and duplicators
   */

  WaterPDSS::WaterPDSS(ThermoPhase *tp, int spindex) :
    PDSS(tp, spindex),
    m_sub(0),
    m_iState(-1),
    EW_Offset(0.0),
    SW_Offset(0.0),
    m_verbose(0),
    m_allowGasPhase(false)
  {
    constructPDSS(tp, spindex);
  }


  WaterPDSS::WaterPDSS(ThermoPhase *tp, int spindex, 
		       string inputFile, string id) :
    PDSS(tp, spindex),
    m_sub(0),
    m_iState(-1),
    m_mw(0.0),
    EW_Offset(0.0),
    SW_Offset(0.0),
    m_verbose(0),
    m_allowGasPhase(false)
  {
    constructPDSSFile(tp, spindex, inputFile, id);
  }

  WaterPDSS::WaterPDSS(ThermoPhase *tp, int spindex,
		       XML_Node& phaseRoot, string id) :
    PDSS(tp, spindex),
    m_sub(0),
    m_iState(-1),
    m_mw(0.0),
    EW_Offset(0.0),
    SW_Offset(0.0),
    m_verbose(0),
    m_allowGasPhase(false)
  {
    constructPDSSXML(tp, spindex, phaseRoot, id) ;
  }



  WaterPDSS::WaterPDSS(const WaterPDSS &b) :
    PDSS(b),
    m_sub(0),
    m_iState(-1),
    m_mw(b.m_mw),
    EW_Offset(b.EW_Offset),
    SW_Offset(b.SW_Offset),
    m_verbose(b.m_verbose),
    m_allowGasPhase(b.m_allowGasPhase)
  {
    m_sub = new WaterPropsIAPWS(*(b.m_sub));  
    /*
     * Use the assignment operator to do the brunt
     * of the work for the copy construtor.
     */
    *this = b;
  }

  /**
   * Assignment operator
   */
  WaterPDSS& WaterPDSS::operator=(const WaterPDSS&b) {
    if (&b == this) return *this;
    m_sub->operator=(*(b.m_sub));
    PDSS::operator=(b);
    m_verbose = b.m_verbose;
    m_allowGasPhase = b.m_allowGasPhase;
    return *this;
  }

  WaterPDSS::~WaterPDSS() { 
    delete m_sub; 
  }
  
  void WaterPDSS::constructPDSS(ThermoPhase *tp, int spindex) {
    initThermo();
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
  void WaterPDSS::constructPDSSXML(ThermoPhase *tp, int spindex,
				   XML_Node& phaseNode, string id) {
    initThermo();
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
  void WaterPDSS::constructPDSSFile(ThermoPhase *tp, int spindex,
				    string inputFile, string id) {

    if (inputFile.size() == 0) {
      throw CanteraError("WaterTp::initThermo",
			 "input file is null");
    }
    string path = findInputFile(inputFile);
    ifstream fin(path.c_str());
    if (!fin) {
      throw CanteraError("WaterPDSS::initThermo","could not open "
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
      throw CanteraError("WaterPDSS::initThermo",
			 "ERROR: Can not find phase named " +
			 id + " in file named " + inputFile);
    }	
    constructPDSSXML(tp, spindex, *fxml_phase, id);
    delete fxml;
  }

  void WaterPDSS::
  initThermoXML(XML_Node& phaseNode, string id) {
    initThermo();
  }

  void WaterPDSS::initThermo() {
    if (m_sub) delete m_sub;
    m_sub = new WaterPropsIAPWS();
    if (m_sub == 0) {
      throw CanteraError("WaterPDSS::initThermo",
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
    double rho0 = m_sub->density(298.15, OneAtm, WATER_LIQUID);
    m_dens = rho0;


  }

  void WaterPDSS::
  setParametersFromXML(const XML_Node& eosdata) {
 
  }

  /**
   * Return the molar enthalpy in units of J kmol-1
   */
  doublereal WaterPDSS::
  enthalpy_mole() const {
    double T = m_temp;
    double dens = m_dens;
    doublereal h = m_sub->enthalpy(T, dens);
    return (h + EW_Offset);
  }

  /**
   * Calculate the internal energy in mks units of
   * J kmol-1 
   */
  doublereal WaterPDSS::
  intEnergy_mole() const {
    double T = m_dens;
    double dens = m_temp;
    doublereal u = m_sub->intEnergy(T, dens);
    return (u + EW_Offset);            
  }

  /**
   * Calculate the entropy in mks units of 
   * J kmol-1 K-1
   */
  doublereal WaterPDSS::
  entropy_mole() const {
    double T = m_temp;
    double dens = m_dens;
    doublereal s = m_sub->entropy(T, dens);
    return (s + SW_Offset); 
  }

  /**
   * Calculate the Gibbs free energy in mks units of
   * J kmol-1 K-1.
   */
  doublereal WaterPDSS::
  gibbs_mole() const {
    double T = m_temp;
    double dens = m_dens;
    doublereal g = m_sub->Gibbs(T, dens);
    return (g + EW_Offset - SW_Offset*T);
  }

  /**
   * Calculate the constant pressure heat capacity
   * in mks units of J kmol-1 K-1
   */
  doublereal WaterPDSS::
  cp_mole() const {
    double T = m_temp;
    double dens = m_dens;
    doublereal cp = m_sub->cp(T, dens);
    return cp;            
  }

  /**
   * Calculate the constant volume heat capacity
   * in mks units of J kmol-1 K-1
   */
  doublereal WaterPDSS::
  cv_mole() const {
    double T = m_temp;
    double dens = m_dens;
    doublereal cv = m_sub->cv(T, dens);
    return cv;
  }

  /**
   * Calculate the pressure (Pascals), given the temperature and density
   *  Temperature: kelvin
   *  rho: density in kg m-3
   */
  doublereal WaterPDSS::
  pressure() const {
    double T = m_temp;
    double dens = m_dens;
    doublereal p = m_sub->pressure(T, dens);
    return p;
  }
        
  void WaterPDSS::
  setTempPressure(doublereal t, doublereal p) {
    m_temp = t;
    setPressure(p);
  }

  void WaterPDSS::
  setPressure(doublereal p) {
    double T = m_temp;
    double dens = m_dens;
    int waterState = WATER_GAS;
    double rc = m_sub->Rhocrit();
    if (dens > rc) {
      waterState = WATER_LIQUID;
    }
#ifdef DEBUG_HKM
    //printf("waterPDSS: set pres = %g t = %g, waterState = %d\n",
    //      p, T, waterState);
#endif
    doublereal dd = m_sub->density(T, p, waterState, dens);
    if (dd <= 0.0) {
      printf("throw an error\n");
      throw CanteraError("WaterPDSS:pressure", "Failed to set water state");
    }
    m_dens = dd;
  }
 

  /// critical temperature 
  doublereal WaterPDSS::critTemperature() const { return m_sub->Tcrit(); }
        
  /// critical pressure
  doublereal WaterPDSS::critPressure() const { return m_sub->Pcrit(); }
        
  /// critical density
  doublereal WaterPDSS::critDensity() const { return m_sub->Rhocrit(); }
        
  void WaterPDSS::setDensity(double dens) {
    m_dens = dens;
    m_sub->setState(m_temp, m_dens);
  }

  double WaterPDSS::density() const {
    return m_dens;
  }


  double WaterPDSS::temperature() const {
    return m_temp;
  }
 

  void WaterPDSS::setTemperature(double temp) {
    m_temp = temp;
    doublereal dd = m_dens;
    m_sub->setState(temp, dd);
  }

  doublereal WaterPDSS::molecularWeight() const {
    return m_mw;
  }
  void WaterPDSS::setMolecularWeight(double mw) {
    m_mw = mw;
  }

  void WaterPDSS::setState_TP(double temp, double pres) {
    m_temp = temp;
    setPressure(pres);
  }

  /// saturation pressure
  doublereal WaterPDSS::satPressure(doublereal t){
    doublereal pp = m_sub->psat(t);
    double dens = m_dens;
    m_temp = t;
    m_dens = dens;
    return pp;
  }
        


}
