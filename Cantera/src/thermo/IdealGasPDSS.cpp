/**
 * @file IdealGasPDSS.cpp
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
 * $Id$
 */

#include "ct_defs.h"
#include "xml.h"
#include "ctml.h"
#include "IdealGasPDSS.h"
#include "importCTML.h"


#include "ThermoPhase.h"

using namespace std;

namespace Cantera {
  /**
   * Basic list of constructors and duplicators
   */

  IdealGasPDSS::IdealGasPDSS(ThermoPhase *tp, int spindex) :
    PDSS(tp, spindex)
  {
    constructPDSS(tp, spindex);
  }


  IdealGasPDSS::IdealGasPDSS(ThermoPhase *tp, int spindex, std::string inputFile, std::string id) :
    PDSS(tp, spindex)
  {
    constructPDSSFile(tp, spindex, inputFile, id);
  }


  IdealGasPDSS::IdealGasPDSS(ThermoPhase *tp, int spindex, XML_Node& phaseRoot, std::string id) :
    PDSS(tp, spindex)
  {
    constructPDSSXML(tp, spindex, phaseRoot, id) ;
  }


  IdealGasPDSS::IdealGasPDSS(const IdealGasPDSS &b) :
    PDSS(b)
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
  IdealGasPDSS& IdealGasPDSS::operator=(const IdealGasPDSS&b) {
    if (&b == this) return *this;
    PDSS::operator=(b);
    return *this;
  }

  IdealGasPDSS::~IdealGasPDSS() { 
  }
  
  void IdealGasPDSS::constructPDSS(ThermoPhase *tp, int spindex) {
    initThermo();
  }

   
  /**
   * constructPDSSXML:
   *
   * Initialization of a IdealGasPDSS object using an
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
  void IdealGasPDSS::constructPDSSXML(ThermoPhase *tp, int spindex, 
			      XML_Node& phaseNode, std::string id) {
    initThermo();
  }

   
  /**
   * constructPDSSFile():
   *
   * Initialization of a IdealGasPDSS object using an
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
  void IdealGasPDSS::constructPDSSFile(ThermoPhase *tp, int spindex,
			       std::string inputFile, std::string id) {

    if (inputFile.size() == 0) {
      throw CanteraError("IdealGasPDSS::initThermo",
			 "input file is null");
    }
    std::string path = findInputFile(inputFile);
    ifstream fin(path.c_str());
    if (!fin) {
      throw CanteraError("IdealGasPDSS::initThermo","could not open "
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
      throw CanteraError("IdealGasPDSS::initThermo",
			 "ERROR: Can not find phase named " +
			 id + " in file named " + inputFile);
    }	
    constructPDSSXML(tp, spindex, *fxml_phase, id);
    delete fxml;
  }

  void IdealGasPDSS::
  initThermoXML(XML_Node& phaseNode, std::string id) {
    initThermo();
  }

  void IdealGasPDSS::initThermo() {
  }

  void IdealGasPDSS::
  setParametersFromXML(const XML_Node& eosdata) {
  }

  /**
   * Return the molar enthalpy in units of J kmol-1
   */
  doublereal IdealGasPDSS::
  enthalpy_mole() const {
    throw CanteraError("IdealGasPDSS::enthalpy_mole()", "unimplemented");
    return (0.0);
  }

  /**
   * Calculate the internal energy in mks units of
   * J kmol-1 
   */
  doublereal IdealGasPDSS::
  intEnergy_mole() const {
    throw CanteraError("IdealGasPDSS::enthalpy_mole()", "unimplemented");
    return (0.0);
  }

  /**
   * Calculate the entropy in mks units of 
   * J kmol-1 K-1
   */
  doublereal IdealGasPDSS::
  entropy_mole() const {

    throw CanteraError("IdealGasPDSS::entropy_mole()", "unimplemented");
    return (0.0);
  }

  /**
   * Calculate the Gibbs free energy in mks units of
   * J kmol-1 K-1.
   */
  doublereal IdealGasPDSS::
  gibbs_mole() const {
    throw CanteraError("IdealGasPDSS::gibbs_mole()", "unimplemented");
    return (0.0);
  }

  /**
   * Calculate the constant pressure heat capacity
   * in mks units of J kmol-1 K-1
   */
  doublereal IdealGasPDSS::
  cp_mole() const {
    throw CanteraError("IdealGasPDSS::cp_mole()", "unimplemented");
    return (0.0);
  }

  /**
   * Calculate the constant volume heat capacity
   * in mks units of J kmol-1 K-1
   */
  doublereal IdealGasPDSS::
  cv_mole() const {
    throw CanteraError("IdealGasPDSS::cv_mole()", "unimplemented");
    return (0.0);
  }

   
  /**
   * Return the difference in enthalpy between current p
   * and ref p0, in mks units of
   * in units of J kmol-1
   */
  doublereal IdealGasPDSS::
  enthalpyDelp_mole() const {
    throw CanteraError("IdealGasPDSS::enthalpy_mole()", "unimplemented");
    return (0.0);
  }

  /**
   * Calculate difference in the internal energy between current p
   * and ref p0, in mks units of
   * J kmol-1 
   */
  doublereal IdealGasPDSS::
  intEnergyDelp_mole() const {
    throw CanteraError("IdealGasPDSS::enthalpyDelp_mole()", "unimplemented");
    return (0.0);
  }

  /**
   *  Return the difference in entropy between current p
   * and ref p0, in mks units of
   * J kmol-1 K-1
   */
  doublereal IdealGasPDSS::
  entropyDelp_mole() const {

    throw CanteraError("IdealGasPDSS::entropyDelp_mole()", "unimplemented");
    return (0.0);
  }

  /**
   * Calculate the difference in Gibbs free energy between current p and
   * the ref p0, in mks units of
   * J kmol-1 K-1.
   */
  doublereal IdealGasPDSS::
  gibbsDelp_mole() const {
    throw CanteraError("IdealGasPDSS::gibbsDelp_mole()", "unimplemented");
    return (0.0);
  }

  /**
   * Calculate the difference in the constant pressure heat capacity
   * between the current p and the ref p0,
   * in mks units of J kmol-1 K-1
   */
  doublereal IdealGasPDSS::
  cpDelp_mole() const {
    throw CanteraError("IdealGasPDSS::cpDelp_mole()", "unimplemented");
    return (0.0);
  }

  /**
   * Calculate the difference in constant volume heat capacity
   * between the current p and the ref p0
   * in mks units of J kmol-1 K-1
   */
  doublereal IdealGasPDSS::
  cvDelp_mole() const {
    throw CanteraError("IdealGasPDSS::cvDelp_mole()", "unimplemented");
    return (0.0);
  }

  /**
   * Calculate the pressure (Pascals), given the temperature and density
   *  Temperature: kelvin
   *  rho: density in kg m-3
   */
  doublereal IdealGasPDSS::
  pressure() const {
    throw CanteraError("IdealGasPDSS::pressure()", "unimplemented");
    return (0.0);
  }
        
  void IdealGasPDSS::
  setPressure(doublereal p) {
    throw CanteraError("IdealGasPDSS::pressure()", "unimplemented");
  }
 

  /// critical temperature 
  doublereal IdealGasPDSS::critTemperature() const { 
    throw CanteraError("IdealGasPDSS::critTemperature()", "unimplemented");
    return (0.0);
  }
        
  /// critical pressure
  doublereal IdealGasPDSS::critPressure() const {
    throw CanteraError("IdealGasPDSS::critPressure()", "unimplemented");
    return (0.0);
  }
        
  /// critical density
  doublereal IdealGasPDSS::critDensity() const {
    throw CanteraError("IdealGasPDSS::critDensity()", "unimplemented");
    return (0.0);
  }
        
  void IdealGasPDSS::setDensity(double dens) {
    m_dens = dens;
  }

  /**
   * Return the density of the standard state
   *
   *  We assume that the storred density is current.
   * Note, this is the density of the standard state,
   * not of the mixture.
   */
  double IdealGasPDSS::density() const {
    return m_dens;
  }

  /**
   * Return the temperature 
   *
   * Obtain the temperature from the owning ThermoPhase object
   * if you can. 
   */
  double IdealGasPDSS::temperature() const {
    if (m_tp) {
      m_temp = m_tp->temperature();
    }
    return m_temp;
  }
 
  void IdealGasPDSS::setTemperature(double temp) {
    m_temp = temp;
  }

  doublereal IdealGasPDSS::molecularWeight() const {
    return m_mw;
  }
  void IdealGasPDSS::setMolecularWeight(double mw) {
    m_mw = mw;
  }

  void IdealGasPDSS::setState_TP(double temp, double pres) {
    throw CanteraError("IdealGasPDSS::setState_TP()", "unimplemented");
  }

  /// saturation pressure
  doublereal IdealGasPDSS::satPressure(doublereal t){
    throw CanteraError("IdealGasPDSS::satPressure()", "unimplemented");
    return (0.0);
  }
    

}
