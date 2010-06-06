/**
 *
 *  @file LatticePhase.cpp
 *  Definitions for a simple thermodynamics model of a bulk phase
 *  derived from ThermoPhase,
 *  assuming a lattice of solid atoms
 *  (see \ref thermoprops and class \link Cantera::LatticePhase LatticePhase\endlink).
 *
 */
/*
 * $Id$
 */

#include "config.h"
#ifdef WITH_LATTICE_SOLID

#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "ct_defs.h"
#include "mix_defs.h"
#include "LatticePhase.h"
#include "SpeciesThermo.h"
#include "ThermoFactory.h"

#include <cmath>
#include <string>

namespace Cantera {

  // Base Empty constructor
  LatticePhase::LatticePhase() :
    m_tlast(0.0) 
  {
  }
  
  // Copy Constructor
  /*
   * @param right Object to be copied
   */
  LatticePhase::LatticePhase(const LatticePhase &right) :
    m_tlast(0.0)
  {
    *this = operator=(right);
  }
  
  // Assignment operator
  /*
   * @param right Object to be copied
   */
  LatticePhase& LatticePhase::operator=(const LatticePhase& right) {
    if (&right != this) {
      ThermoPhase::operator=(right);
      m_mm         = right.m_mm;
      m_tmin       = right.m_tmin;
      m_tmax       = right.m_tmax;
      m_p0         = right.m_p0;
      m_tlast      = right.m_tlast;
      m_h0_RT      = right.m_h0_RT;
      m_cp0_R      = right.m_cp0_R;
      m_g0_RT      = right.m_g0_RT;
      m_s0_R       = right.m_s0_R;
      m_press      = right.m_press;
      m_vacancy    = right.m_vacancy;
      m_molar_density = right.m_molar_density;
    }
    return *this;
  }
  
  // Destructor
  LatticePhase::~LatticePhase() {
  }
 

  // Full constructor for a lattice phase
  /*
   * @param inputFile String name of the input file
   * @param id        string id of the phase name
   */
  LatticePhase::LatticePhase(std::string inputFile, std::string id) {
    constructPhaseFile(inputFile, id);
  }
  
  // Full constructor for a water phase
  /*
   * @param phaseRef  XML node referencing the lattice phase.
   * @param id        string id of the phase name
   */
  LatticePhase::LatticePhase(XML_Node& phaseRef, std::string id) {
    constructPhaseXML(phaseRef, id);
  }


  // Duplication function
  /*
   * This virtual function is used to create a duplicate of the
   * current phase. It's used to duplicate the phase when given
   * a ThermoPhase pointer to the phase.
   *
   * @return It returns a ThermoPhase pointer.
   */
  ThermoPhase *LatticePhase::duplMyselfAsThermoPhase() const {
    LatticePhase *igp = new LatticePhase(*this);
    return (ThermoPhase *) igp;
  }
  
  /*
   * @param infile XML file containing the description of the
   *        phase
   *
   * @param id  Optional parameter identifying the name of the
   *            phase. If none is given, the first XML
   *            phase element will be used.
   */
  void LatticePhase::constructPhaseXML(XML_Node& phaseNode, std::string idTarget) {
    std::string idattrib = phaseNode.id();
    if (idTarget != idattrib) {
      throw CanteraError("LatticePhase::constructPhaseXML","ids don't match");
    }

    /*
     * Call the Cantera importPhase() function. This will import
     * all of the species into the phase. This will also handle
     * all of the solvent and solute standard states.
     */
    bool m_ok = importPhase(phaseNode, this);
    if (!m_ok) {
      throw CanteraError("LatticePhase::constructPhaseXML","importPhase failed ");
    }
  }

  /*
   * constructPhaseFile
   *
   *
   * This routine is a precursor to constructPhaseXML(XML_Node*)
   * routine, which does most of the work.
   *
   * @param inputFile XML file containing the description of the
   *        phase
   *
   * @param id  Optional parameter identifying the name of the
   *            phase. If none is given, the first XML
   *            phase element will be used.
   */
  void LatticePhase::constructPhaseFile(std::string inputFile, std::string id) {

    if (inputFile.size() == 0) {
      throw CanteraError("LatticePhase::constructPhaseFile",
                         "input file is null");
    }
    std::string path = findInputFile(inputFile);
    std::ifstream fin(path.c_str());
    if (!fin) {
      throw CanteraError("LatticePhase::constructPhaseFile","could not open "
                         +path+" for reading.");
    }
    /*
     * The phase object automatically constructs an XML object.
     * Use this object to store information.
     */
    XML_Node &phaseNode_XML = xml();
    XML_Node *fxml = new XML_Node();
    fxml->build(fin);
    XML_Node *fxml_phase = findXMLPhase(fxml, id);
    if (!fxml_phase) {
      throw CanteraError("LatticePhase::constructPhaseFile",
                         "ERROR: Can not find phase named " +
                         id + " in file named " + inputFile);
    }
    fxml_phase->copy(&phaseNode_XML);
    constructPhaseXML(*fxml_phase, id);
    delete fxml;
  }

  
  doublereal LatticePhase::
  enthalpy_mole() const {
    doublereal p0 = m_spthermo->refPressure();
    return GasConstant * temperature() * 
      mean_X(&enthalpy_RT_ref()[0]) 
      + (pressure() - p0)/molarDensity();
  }

  doublereal LatticePhase::intEnergy_mole() const {
    doublereal p0 = m_spthermo->refPressure();
    return GasConstant * temperature() * 
      mean_X(&enthalpy_RT_ref()[0]) 
      - p0/molarDensity();
  }

  doublereal LatticePhase::entropy_mole() const {
    return GasConstant * (mean_X(&entropy_R_ref()[0]) -
			  sum_xlogx());
  }

  doublereal LatticePhase::gibbs_mole() const {
    return enthalpy_mole() - temperature() * entropy_mole();
  }

  doublereal LatticePhase::cp_mole() const {
    return GasConstant * mean_X(&cp_R_ref()[0]);
  }

  doublereal LatticePhase::cv_mole() const {
    return cp_mole();
  }
  
  

  void LatticePhase::setPressure(doublereal p) {
    m_press = p;
    setMolarDensity(m_molar_density);
  }

  void LatticePhase::getActivityConcentrations(doublereal* c) const {
    getMoleFractions(c);
  }

  void LatticePhase::getActivityCoefficients(doublereal* ac) const {
    for (int k = 0; k < m_kk; k++) {
      ac[k] = 1.0;
    }
  }

  doublereal LatticePhase::standardConcentration(int k) const {
    return 1.0;
  }

  doublereal LatticePhase::logStandardConc(int k) const {
    return 0.0;
  }

  void LatticePhase::getChemPotentials(doublereal* mu) const {
    doublereal vdp = ((pressure() - m_spthermo->refPressure())/
		      molarDensity());
    doublereal xx;
    doublereal rt = temperature() * GasConstant;
    const array_fp& g_RT = gibbs_RT_ref();
    for (int k = 0; k < m_kk; k++) {
      xx = fmaxx(SmallNumber, moleFraction(k));
      mu[k] = rt*(g_RT[k] + log(xx)) + vdp;
    }
  }

  void LatticePhase::getPartialMolarVolumes(doublereal* vbar) const {
    getStandardVolumes(vbar); 
  }

  void LatticePhase::getStandardChemPotentials(doublereal* mu0) const {
    const array_fp& gibbsrt = gibbs_RT_ref();
    scale(gibbsrt.begin(), gibbsrt.end(), mu0, _RT());
  }
  
  void LatticePhase::getPureGibbs(doublereal* gpure) const {
    const array_fp& gibbsrt = gibbs_RT_ref();
    scale(gibbsrt.begin(), gibbsrt.end(), gpure, _RT());
  }

  void LatticePhase::getEnthalpy_RT(doublereal* hrt) const {
    const array_fp& _h = enthalpy_RT_ref();
    std::copy(_h.begin(), _h.end(), hrt);
    doublereal tmp = (pressure() - m_p0) / (molarDensity() * GasConstant * temperature());
    for (int k = 0; k < m_kk; k++) {
      hrt[k] += tmp;
    }
  }

  void LatticePhase::getEntropy_R(doublereal* sr) const {
    const array_fp& _s = entropy_R_ref();
    std::copy(_s.begin(), _s.end(), sr);
  }

  void LatticePhase::getGibbs_RT(doublereal* grt) const {
    const array_fp& gibbsrt = gibbs_RT_ref();
    std::copy(gibbsrt.begin(), gibbsrt.end(), grt);
  }

  void LatticePhase::getCp_R(doublereal* cpr) const {
    const array_fp& _cpr = cp_R_ref();
    std::copy(_cpr.begin(), _cpr.end(), cpr);
  }

  void LatticePhase::getStandardVolumes(doublereal* vbar) const {
    doublereal vv = 1.0/m_molar_density;
    for (int k = 0; k < m_kk; k++) {
      vbar[k] = vv;
    }
  }

  void LatticePhase::initThermo() {
    m_kk = nSpecies();
    m_mm = nElements();
    doublereal tmin = m_spthermo->minTemp();
    doublereal tmax = m_spthermo->maxTemp();
    if (tmin > 0.0) m_tmin = tmin;
    if (tmax > 0.0) m_tmax = tmax;
    m_p0 = refPressure();

    int leng = m_kk;
    m_h0_RT.resize(leng);
    m_g0_RT.resize(leng);
    m_cp0_R.resize(leng);
    m_s0_R.resize(leng);
    setMolarDensity(m_molar_density);
  }


  void LatticePhase::_updateThermo() const {
    doublereal tnow = temperature();
    if (fabs(molarDensity() - m_molar_density)/m_molar_density > 0.0001) {
      throw CanteraError("_updateThermo","molar density changed from "
			 +fp2str(m_molar_density)+" to "+fp2str(molarDensity()));
    }
    if (m_tlast != tnow) {
      m_spthermo->update(tnow, &m_cp0_R[0], &m_h0_RT[0], 
			 &m_s0_R[0]);
      m_tlast = tnow;
      int k;
      for (k = 0; k < m_kk; k++) {
	m_g0_RT[k] = m_h0_RT[k] - m_s0_R[k];
      }
      m_tlast = tnow;
    }
  }

  void LatticePhase::setParameters(int n, doublereal* const c) {
    m_molar_density = c[0];
    setMolarDensity(m_molar_density);
  }

  void LatticePhase::getParameters(int &n, doublereal * const c) const {
    double d = molarDensity();
    c[0] = d;
    n = 1;
  }

  void LatticePhase::setParametersFromXML(const XML_Node& eosdata) {
    eosdata._require("model", "Lattice");
    m_molar_density = getFloat(eosdata, "site_density", "toSI");
    m_vacancy = getChildValue(eosdata, "vacancy_species");
  }
}

#endif
