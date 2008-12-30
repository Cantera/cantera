/**
 *  @file Phase.cpp
 *   Definition file for class, Phase, which contains functions for setting the
 *   state of a phase, and for referencing species by name
 *   (see \ref phases and class \link Cantera::Phase Phase\endlink).
 */

// Copyright 2001  California Institute of Technology

#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "ct_defs.h"
#include "Phase.h"
#include "vec_functions.h"
#include "ctexceptions.h"

using namespace std;

namespace Cantera {    

  Phase::Phase() : 
    Constituents(),
    State(),
    m_kk(-1), 
    m_ndim(3),
    m_index(-1), 
    m_xml(new XML_Node("phase")), 
    m_id("<phase>"),
    m_name("") 
  {
  }

  /*
   * Copy Constructor
   *
   * This function just does the default initialization, and
   * then calls the assignment operator.
   */
  Phase::Phase(const Phase &right) :
    Constituents(),
    State(),
    m_kk(-1),
    m_ndim(3),
    m_index(-1), 
    m_xml(new XML_Node("phase")), 
    m_id("<phase>"),
    m_name("") 
  {
    /*
     * Call the assignment operator.
     */
    *this = operator=(right);
  }
    
  /*
   * Assignment operator
   *
   * This operation is sort of complicated. We have to
   * call the assignment operator for the Constituents and
   * State operators that Phase inherits from. Then,
   * we have to copy our own data, making sure to do a 
   * deep copy on the XML_Node data owned by this object.
   */
  Phase &Phase::operator=(const Phase &right) {
    /*
     * Check for self assignment.
     */
    if (this == &right) return *this;
    /*
     * Now call the inherited-classes assignment operators.
     */
    (void) Constituents::operator=(right);
    (void) State::operator=(right);
    /*
     * Handle its own data
     */
    m_kk    = right.m_kk;
    m_ndim  = right.m_ndim;
    m_index = right.m_index;
    m_data  = right.m_data;
    /*
     * This is a little complicated. -> Because we delete m_xml
     * in the destructor, we own m_xml completely, and we need
     * to have our own individual copies of the XML data tree
     * in each object
     */
    if (m_xml) {
      delete m_xml;
      m_xml = 0;
    }
    if (right.m_xml) {
      m_xml   = new XML_Node();
      (right.m_xml)->copy(m_xml);
    }
    m_id    = right.m_id;
    m_name  = right.m_name;

    return *this;
  }

  // Destructor.
  Phase::~Phase() { 
    if (m_xml) {
      delete m_xml;
      m_xml = 0;
    }
  }

  XML_Node& Phase::xml() { 
    return *m_xml;
  }

  std::string Phase::id() const { 
    return m_id; 
  }

  void Phase::setID(std::string id) {
    m_id = id;
  } 

  std::string Phase::name() const {
    return m_name; 
  }

  void Phase::setName(std::string nm) { 
    m_name = nm; 
  }

  int Phase::index() const { 
    return m_index;
  }

  void Phase::setIndex(int m) { 
    m_index = m;
  }

  void Phase::saveState(vector_fp& state) const {
    state.resize(nSpecies() + 2);
    saveState(state.size(),&(state[0]));
  }
  void Phase::saveState(int lenstate, doublereal* state) const {
    state[0] = temperature();
    state[1] = density();
    getMassFractions(state + 2);
  }

  void Phase::restoreState(const vector_fp& state) {
    restoreState(state.size(),&state[0]);
  }

  void Phase::restoreState(int lenstate, const doublereal* state) {
    if (int(lenstate) >= nSpecies() + 2) {
      setMassFractions_NoNorm(state + 2);
      setTemperature(state[0]);
      setDensity(state[1]);
    }
    else {
      throw ArraySizeError("Phase::restoreState",
			   lenstate,nSpecies()+2);
    }
  }

  void Phase::setMoleFractionsByName(compositionMap& xMap) {
    int kk = nSpecies();
    doublereal x;
    vector_fp mf(kk, 0.0);
    for (int k = 0; k < kk; k++) {
      x = xMap[speciesName(k)];
      if (x > 0.0) mf[k] = x;
    }
    setMoleFractions(&mf[0]);
  }

  void Phase::setMoleFractionsByName(const std::string& x) {
    compositionMap xx;
    int kk = nSpecies();
    for (int k = 0; k < kk; k++) { 
      xx[speciesName(k)] = -1.0;
    }
    parseCompString(x, xx);
    setMoleFractionsByName(xx);
    //int kk = nSpecies();
    //vector_fp mf(kk);
    //for (int k = 0; k < kk; k++) { 
    //    mf[k] = xx[speciesName(k)];
    //}
    //setMoleFractions(mf.begin());
  }

  void Phase::setMassFractionsByName(compositionMap& yMap) {
    int kk = nSpecies();
    doublereal y;
    vector_fp mf(kk, 0.0);
    for (int k = 0; k < kk; k++) { 
      y = yMap[speciesName(k)];
      if (y > 0.0) mf[k] = y;
    }
    setMassFractions(&mf[0]);
  }

  void Phase::setMassFractionsByName(const std::string& y) {
    compositionMap yy;
    int kk = nSpecies();
    for (int k = 0; k < kk; k++) { 
      yy[speciesName(k)] = -1.0;
    }
    parseCompString(y, yy);
    setMassFractionsByName(yy);
  }

  /** Set the temperature (K), density (kg/m^3), and mole fractions. */
  void Phase::setState_TRX(doublereal t, doublereal dens, 
			   const doublereal* x) {
    setMoleFractions(x); setTemperature(t); setDensity(dens);
  }

  void Phase::setState_TNX(doublereal t, doublereal n, 
			   const doublereal* x) {
    setMoleFractions(x); setTemperature(t); setMolarDensity(n);
  }

  /** Set the temperature (K), density (kg/m^3), and mole fractions. */
  void Phase::setState_TRX(doublereal t, doublereal dens, 
			   compositionMap& x) {
    setMoleFractionsByName(x); setTemperature(t); setDensity(dens);
  }

  /** Set the temperature (K), density (kg/m^3), and mass fractions. */
  void Phase::setState_TRY(doublereal t, doublereal dens, 
			   const doublereal* y) {
    setMassFractions(y); setTemperature(t); setDensity(dens);
  }        

  /** Set the temperature (K), density (kg/m^3), and mass fractions. */
  void Phase::setState_TRY(doublereal t, doublereal dens, 
			   compositionMap& y) {
    setMassFractionsByName(y); setTemperature(t); setDensity(dens);
  }
    
  /** Set the temperature (K) and density (kg/m^3) */
  void Phase::setState_TR(doublereal t, doublereal rho) {
    setTemperature(t); setDensity(rho);
  }
    
  /** Set the temperature (K) and mole fractions.  */
  void Phase::setState_TX(doublereal t, doublereal* x) {
    setTemperature(t); setMoleFractions(x);
  }

  /** Set the temperature (K) and mass fractions.  */
  void Phase::setState_TY(doublereal t, doublereal* y) {
    setTemperature(t); setMassFractions(y);
  }

  /** Set the density (kg/m^3) and mole fractions.  */
  void Phase::setState_RX(doublereal rho, doublereal* x) {
    setMoleFractions(x); setDensity(rho);
  }

  /** Set the density (kg/m^3) and mass fractions.  */
  void Phase::setState_RY(doublereal rho, doublereal* y) {
    setMassFractions(y); setDensity(rho);
  }

  /*
   * Copy the vector of molecular weights into vector weights.
   */
  void Phase::getMolecularWeights(vector_fp& weights) const {
    const array_fp& mw = Constituents::molecularWeights();
    if (weights.size() < mw.size()) weights.resize(mw.size());
    copy(mw.begin(), mw.end(), weights.begin());
  }

  /*
   * Copy the vector of molecular weights into array weights.
   * @deprecated
   */
  void Phase::getMolecularWeights(int iwt, doublereal* weights) const {
    const array_fp& mw = Constituents::molecularWeights();
    copy(mw.begin(), mw.end(), weights);
  }

  /*
   * Copy the vector of molecular weights into array weights.
   */
  void Phase::getMolecularWeights(doublereal* weights) const {
    const array_fp& mw = Constituents::molecularWeights();
    copy(mw.begin(), mw.end(), weights);
  }

  /**
   * Return a const reference to the internal vector of
   * molecular weights.
   */
  const array_fp& Phase::molecularWeights() const {
    return Constituents::molecularWeights(); 
  }


  /**
   * Get the mole fractions by name. 
   */
  void Phase::getMoleFractionsByName(compositionMap& x) const {
    x.clear();
    int kk = nSpecies();
    for (int k = 0; k < kk; k++) {
      x[speciesName(k)] = State::moleFraction(k);
    }
  }

  doublereal Phase::moleFraction(int k) const {
    return State::moleFraction(k);
  }

  doublereal Phase::moleFraction(std::string name) const {
    int iloc = speciesIndex(name);
    if (iloc >= 0) return State::moleFraction(iloc);
    else return 0.0;
  }

  doublereal Phase::massFraction(int k) const {
    return State::massFraction(k);
  }

  doublereal Phase::massFraction(std::string name) const {
    int iloc = speciesIndex(name);
    if (iloc >= 0) return massFractions()[iloc];
    else return 0.0;
  }

  doublereal Phase::chargeDensity() const {
    int k;
    int nsp = nSpecies();
    doublereal cdens = 0.0;
    for (k = 0; k < nsp; k++) 
      cdens += charge(k)*State::moleFraction(k);
    cdens *= Faraday;
    return cdens;
  }

  /** 
   *  Finished adding species, prepare to use them for calculation
   *  of mixture properties.
   */
  void Phase::freezeSpecies() {
    Constituents::freezeSpecies();
    init(Constituents::molecularWeights());
    int kk = nSpecies();
    int nv = kk + 2;
    m_data.resize(nv,0.0);
    m_data[0] = 300.0;
    m_data[1] = 0.001;
    m_data[2] = 1.0;

    //setState_TRY(300.0, density(), &m_data[2]);

    m_kk = nSpecies();
  } 

  bool Phase::ready() const {
    return (m_kk > 0 && Constituents::ready() && State::ready());
  }

  //         int Phase::installUpdater_T(Updater* u) {
  //             return m_T_updater.install(u);
  //         }

  //         int Phase::installUpdater_C(Updater* u) {
  //             return m_C_updater.install(u);
  //         }
}
