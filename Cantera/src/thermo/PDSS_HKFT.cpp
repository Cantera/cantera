/*
 * $Id$
 */
#include "ct_defs.h"
#include "xml.h"
#include "ctml.h"
#include "PDSS_HKFT.h"

#include "VPStandardStateTP.h"

using namespace std;

namespace Cantera {
  /**
   * Basic list of constructors and duplicators
   */



  PDSS_HKFT::PDSS_HKFT(VPStandardStateTP *tp, int spindex) :
    PDSS(tp, spindex)
  {
  }


  PDSS_HKFT::PDSS_HKFT(VPStandardStateTP *tp, int spindex, std::string inputFile, std::string id) :
    PDSS(tp, spindex)
  {
  }

  PDSS_HKFT::PDSS_HKFT(VPStandardStateTP *tp, int spindex, const XML_Node& speciesNode,
                       const XML_Node& phaseRoot, bool spInstalled) :
    PDSS(tp, spindex)
  {
  }

  PDSS_HKFT::PDSS_HKFT(const PDSS_HKFT &b) :
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
  PDSS_HKFT& PDSS_HKFT::operator=(const PDSS_HKFT&b) {
    if (&b == this) return *this;
    m_tp = b.m_tp;
    m_spindex = b.m_spindex;
    m_temp = b.m_temp;
    m_pres = b.m_pres;
    m_mw = b.m_mw;
    return *this;
  }
  
  /**
   * Destructor for the PDSS_HKFT class
   */
  PDSS_HKFT::~PDSS_HKFT() { 
  }
  




  /**
   * Return the molar enthalpy in units of J kmol-1
   */
  doublereal 
  PDSS_HKFT::enthalpy_mole() const {
    throw CanteraError("PDSS_HKFT::enthalpy_mole()", "unimplemented");
    return (0.0);
  }

  doublereal 
  PDSS_HKFT::enthalpy_RT() const {
    throw CanteraError("PDSS_HKFT::enthalpy_RT()", "unimplemented");
    return (0.0);
  }

  /**
   * Calculate the internal energy in mks units of
   * J kmol-1 
   */
  doublereal 
  PDSS_HKFT::intEnergy_mole() const {
    throw CanteraError("PDSS_HKFT::enthalpy_mole()", "unimplemented");
    return (0.0);
  }

  /**
   * Calculate the entropy in mks units of 
   * J kmol-1 K-1
   */
  doublereal
  PDSS_HKFT::entropy_mole() const {

    throw CanteraError("PDSS_HKFT::entropy_mole()", "unimplemented");
    return (0.0);
  }

  /**
   * Calculate the Gibbs free energy in mks units of
   * J kmol-1 K-1.
   */
  doublereal 
  PDSS_HKFT::gibbs_mole() const {
    throw CanteraError("PDSS_HKFT::gibbs_mole()", "unimplemented");
    return (0.0);
  }

  /**
   * Calculate the constant pressure heat capacity
   * in mks units of J kmol-1 K-1
   */
  doublereal 
  PDSS_HKFT::cp_mole() const {
    throw CanteraError("PDSS_HKFT::cp_mole()", "unimplemented");
    return (0.0);
  }

  /**
   * Calculate the constant volume heat capacity
   * in mks units of J kmol-1 K-1
   */
  doublereal 
  PDSS_HKFT::cv_mole() const {
    throw CanteraError("PDSS_HKFT::cv_mole()", "unimplemented");
    return (0.0);
  }

  doublereal 
  PDSS_HKFT::molarVolume() const {
    throw CanteraError("PDSS_HKFT::molarVolume()", "unimplemented");
    return (0.0);
  }

  doublereal 
  PDSS_HKFT::gibbs_RT_ref() const {
    throw CanteraError("PDSS_HKFT::gibbs_RT_ref()", "unimplemented");
    return (0.0);
  }

  doublereal 
  PDSS_HKFT::enthalpy_RT_ref() const {
    throw CanteraError("PDSS_HKFT::enthalpy_RT_ref()", "unimplemented");
    return (0.0);
  }

  doublereal 
  PDSS_HKFT::entropy_R_ref() const {
    throw CanteraError("PDSS_HKFT::entropy_RT_ref()", "unimplemented");
    return (0.0);
  }

  doublereal 
  PDSS_HKFT::cp_R_ref() const {
    throw CanteraError("PDSS_HKFT::cp_RT_ref()", "unimplemented");
    return (0.0);
  }
  
  doublereal
  PDSS_HKFT::molarVolume_ref() const {
    throw CanteraError("PDSS_HKFT::molarVolume_ref()", "unimplemented");
    return (0.0);
  }
  
  /**
   * Calculate the pressure (Pascals), given the temperature and density
   *  Temperature: kelvin
   *  rho: density in kg m-3
   */
  doublereal
  PDSS_HKFT::pressure() const {
    throw CanteraError("PDSS_HKFT::pressure()", "unimplemented");
    return (0.0);
  }
        
  void 
  PDSS_HKFT::setPressure(doublereal p) {
    throw CanteraError("PDSS_HKFT::pressure()", "unimplemented");
  }
 
  void PDSS_HKFT::setTemperature(doublereal temp) {
    m_temp = temp;
  }

  doublereal PDSS_HKFT::temperature() const {
    return m_temp;
  }


  void PDSS_HKFT::setState_TP(doublereal temp, doublereal pres) {
    setTemperature(temp);
    setPressure(pres);
  }

  /// critical temperature 
  doublereal 
  PDSS_HKFT::critTemperature() const { 
    throw CanteraError("PDSS_HKFT::critTemperature()", "unimplemented");
    return (0.0);
  }
        
  /// critical pressure
  doublereal PDSS_HKFT::critPressure() const {
    throw CanteraError("PDSS_HKFT::critPressure()", "unimplemented");
    return (0.0);
  }
        
  /// critical density
  doublereal PDSS_HKFT::critDensity() const {
    throw CanteraError("PDSS_HKFT::critDensity()", "unimplemented");
    return (0.0);
  }
        
 


  void PDSS_HKFT::initThermo() {
    PDSS::initThermo();
    SpeciesThermo &sp = m_tp->speciesThermo();
    m_p0 = sp.refPressure(m_spindex);
  
  }


  void PDSS_HKFT::initThermoXML(const XML_Node& phaseNode, std::string id) {
    PDSS::initThermoXML(phaseNode, id);
  }

  void PDSS_HKFT::constructPDSSXML(VPStandardStateTP *tp, int spindex,
				   const XML_Node& speciesNode, 
				   const XML_Node& phaseNode, bool spInstalled) {
    PDSS::initThermo();
    SpeciesThermo &sp = m_tp->speciesThermo();
    m_p0 = sp.refPressure(m_spindex);

    if (!spInstalled) {
      throw CanteraError("PDSS_HKFT::constructPDSSXML", "spInstalled false not handled");
    }

    const XML_Node *ss = speciesNode.findByName("standardState");
    if (!ss) {
      throw CanteraError("PDSS_HKFT::constructPDSSXML",
			 "no standardState Node for species " + speciesNode.name());
    }
    std::string model = (*ss)["model"];
    if (model != "constant_incompressible") {
      throw CanteraError("PDSS_HKFT::initThermoXML",
			 "standardState model for species isn't constant_incompressible: "
			 + speciesNode.name());
    }
  
  
    std::string id = "";
    initThermoXML(phaseNode, id);
  }

  void PDSS_HKFT::constructPDSSFile(VPStandardStateTP *tp, int spindex,
				    std::string inputFile, std::string id) {

    if (inputFile.size() == 0) {
      throw CanteraError("PDSS_HKFT::initThermo",
			 "input file is null");
    }
    std::string path = findInputFile(inputFile);
    ifstream fin(path.c_str());
    if (!fin) {
      throw CanteraError("PDSS_HKFT::initThermo","could not open "
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
      throw CanteraError("PDSS_HKFT::initThermo",
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

}
