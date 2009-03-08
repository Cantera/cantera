/**
 * @file PDSS.cpp
 *
 * Implementation of a pressure dependent standard state 
 * virtual function.
 */
/*
 * Copywrite (2006) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
/*
 * $Id: PDSS.cpp,v 1.1 2006/07/04 00:01:53 hkmoffa Exp $
 */

#include "ct_defs.h"
#include "xml.h"
#include "ctml.h"
#include "PDSS.h"
#include "importCTML.h"


#include "ThermoPhase.h"

namespace Cantera {
  /**
   * Basic list of constructors and duplicators
   */

  PDSS::PDSS(ThermoPhase *tp, int spindex) :
    m_temp(-1.0),
    m_dens(-1.0),
    m_tp(tp),
    m_spindex(spindex),
    m_mw(0.0)
  {
    constructPDSS(tp, spindex);
  }


  PDSS::PDSS(ThermoPhase *tp, int spindex, string inputFile, string id) :
    m_temp(-1.0),
    m_dens(-1.0),
    m_tp(tp),
    m_spindex(spindex),
    m_mw(0.0)
  {
    constructPDSSFile(tp, spindex, inputFile, id);
  }


  PDSS::PDSS(ThermoPhase *tp, int spindex, XML_Node& phaseRoot, string id) :
    m_temp(-1.0),
    m_dens(-1.0),
    m_tp(0),
    m_spindex(0),
    m_mw(0.0)
  {
    constructPDSSXML(tp, spindex, phaseRoot, id) ;
  }


  PDSS::PDSS(const PDSS &b) :
    m_temp(-1.0),
    m_dens(-1.0),
    m_tp(0),
    m_spindex(0),
    m_mw(b.m_mw)
  {
    /*
     * Use the assignment operator to do the brunt
     * of the work for the copy construtor.
     */
    *this = b;
  }

  /**
   * Assignment operator
   */
  PDSS& PDSS::operator=(const PDSS&b) {
    if (&b == this) return *this;
    m_tp = b.m_tp;
    m_spindex = b.m_spindex;
    m_temp = b.m_temp;
    m_dens = b.m_dens;
    m_mw = b.m_mw;
    return *this;
  }

  PDSS::~PDSS() { 
  }
  
  void PDSS::constructPDSS(ThermoPhase *tp, int spindex) {
    initThermo();
  }

   
  /**
   * constructPDSSXML:
   *
   * Initialization of a PDSS object using an
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
  void PDSS::constructPDSSXML(ThermoPhase *tp, int spindex, 
			      XML_Node& phaseNode, string id) {
    initThermo();
  }

   
  /**
   * constructPDSSFile():
   *
   * Initialization of a PDSS object using an
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
  void PDSS::constructPDSSFile(ThermoPhase *tp, int spindex,
			       string inputFile, string id) {

    if (inputFile.size() == 0) {
      throw CanteraError("PDSS::initThermo",
			 "input file is null");
    }
    string path = findInputFile(inputFile);
    ifstream fin(path.c_str());
    if (!fin) {
      throw CanteraError("PDSS::initThermo","could not open "
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
      throw CanteraError("PDSS::initThermo",
			 "ERROR: Can not find phase named " +
			 id + " in file named " + inputFile);
    }	
    constructPDSSXML(tp, spindex, *fxml_phase, id);
    delete fxml;
  }

  void PDSS::
  initThermoXML(XML_Node& phaseNode, string id) {
    initThermo();
  }

  void PDSS::initThermo() {
  }

  void PDSS::
  setParametersFromXML(const XML_Node& eosdata) {
  }

  /**
   * Return the molar enthalpy in units of J kmol-1
   */
  doublereal PDSS::
  enthalpy_mole() const {
    throw CanteraError("PDSS::enthalpy_mole()", "unimplemented");
    return (0.0);
  }

  /**
   * Calculate the internal energy in mks units of
   * J kmol-1 
   */
  doublereal PDSS::
  intEnergy_mole() const {
    throw CanteraError("PDSS::enthalpy_mole()", "unimplemented");
    return (0.0);
  }

  /**
   * Calculate the entropy in mks units of 
   * J kmol-1 K-1
   */
  doublereal PDSS::
  entropy_mole() const {

    throw CanteraError("PDSS::entropy_mole()", "unimplemented");
    return (0.0);
  }

  /**
   * Calculate the Gibbs free energy in mks units of
   * J kmol-1 K-1.
   */
  doublereal PDSS::
  gibbs_mole() const {
    throw CanteraError("PDSS::gibbs_mole()", "unimplemented");
    return (0.0);
  }

  /**
   * Calculate the constant pressure heat capacity
   * in mks units of J kmol-1 K-1
   */
  doublereal PDSS::
  cp_mole() const {
    throw CanteraError("PDSS::cp_mole()", "unimplemented");
    return (0.0);
  }

  /**
   * Calculate the constant volume heat capacity
   * in mks units of J kmol-1 K-1
   */
  doublereal PDSS::
  cv_mole() const {
    throw CanteraError("PDSS::cv_mole()", "unimplemented");
    return (0.0);
  }

   
  /**
   * Return the difference in enthalpy between current p
   * and ref p0, in mks units of
   * in units of J kmol-1
   */
  doublereal PDSS::
  enthalpyDelp_mole() const {
    throw CanteraError("PDSS::enthalpy_mole()", "unimplemented");
    return (0.0);
  }

  /**
   * Calculate difference in the internal energy between current p
   * and ref p0, in mks units of
   * J kmol-1 
   */
  doublereal PDSS::
  intEnergyDelp_mole() const {
    throw CanteraError("PDSS::enthalpyDelp_mole()", "unimplemented");
    return (0.0);
  }

  /**
   *  Return the difference in entropy between current p
   * and ref p0, in mks units of
   * J kmol-1 K-1
   */
  doublereal PDSS::
  entropyDelp_mole() const {

    throw CanteraError("PDSS::entropyDelp_mole()", "unimplemented");
    return (0.0);
  }

  /**
   * Calculate the difference in Gibbs free energy between current p and
   * the ref p0, in mks units of
   * J kmol-1 K-1.
   */
  doublereal PDSS::
  gibbsDelp_mole() const {
    throw CanteraError("PDSS::gibbsDelp_mole()", "unimplemented");
    return (0.0);
  }

  /**
   * Calculate the difference in the constant pressure heat capacity
   * between the current p and the ref p0,
   * in mks units of J kmol-1 K-1
   */
  doublereal PDSS::
  cpDelp_mole() const {
    throw CanteraError("PDSS::cpDelp_mole()", "unimplemented");
    return (0.0);
  }

  /**
   * Calculate the difference in constant volume heat capacity
   * between the current p and the ref p0
   * in mks units of J kmol-1 K-1
   */
  doublereal PDSS::
  cvDelp_mole() const {
    throw CanteraError("PDSS::cvDelp_mole()", "unimplemented");
    return (0.0);
  }

  /**
   * Calculate the pressure (Pascals), given the temperature and density
   *  Temperature: kelvin
   *  rho: density in kg m-3
   */
  doublereal PDSS::
  pressure() const {
    throw CanteraError("PDSS::pressure()", "unimplemented");
    return (0.0);
  }
        
  void PDSS::
  setPressure(doublereal p) {
    throw CanteraError("PDSS::pressure()", "unimplemented");
  }
 

  /// critical temperature 
  doublereal PDSS::critTemperature() const { 
    throw CanteraError("PDSS::critTemperature()", "unimplemented");
    return (0.0);
  }
        
  /// critical pressure
  doublereal PDSS::critPressure() const {
    throw CanteraError("PDSS::critPressure()", "unimplemented");
    return (0.0);
  }
        
  /// critical density
  doublereal PDSS::critDensity() const {
    throw CanteraError("PDSS::critDensity()", "unimplemented");
    return (0.0);
  }
        
  void PDSS::setDensity(double dens) {
    m_dens = dens;
  }

  /**
   * Return the density of the standard state
   *
   *  We assume that the storred density is current.
   * Note, this is the density of the standard state,
   * not of the mixture.
   */
  double PDSS::density() const {
    return m_dens;
  }

  /**
   * Return the temperature 
   *
   * Obtain the temperature from the owning ThermoPhase object
   * if you can. 
   */
  double PDSS::temperature() const {
    if (m_tp) {
      m_temp = m_tp->temperature();
    }
    return m_temp;
  }
 
  void PDSS::setTemperature(double temp) {
    m_temp = temp;
  }

  doublereal PDSS::molecularWeight() const {
    return m_mw;
  }
  void PDSS::setMolecularWeight(double mw) {
    m_mw = mw;
  }

  void PDSS::setState_TP(double temp, double pres) {
    throw CanteraError("PDSS::setState_TP()", "unimplemented");
  }

  /// saturation pressure
  doublereal PDSS::satPressure(doublereal t){
    throw CanteraError("PDSS::satPressure()", "unimplemented");
    return (0.0);
  }
    

}
