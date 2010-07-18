/**
 *  @file LatticeSolidPhase.h
 *  Definitions for a simple thermodynamics model of a bulk solid phase
 *  derived from %ThermoPhase,
 *  assuming an ideal solution model based on a lattice of solid atoms
 *  (see \ref thermoprops and class \link Cantera::LatticeSolidPhase LatticeSolidPhase\endlink).
 */
/*
 * $Id$
 */

#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif


#include "ct_defs.h"
#ifdef WITH_LATTICE_SOLID

#include "mix_defs.h"
#include "LatticeSolidPhase.h"
#include "LatticePhase.h"
#include "SpeciesThermo.h"
#include "ThermoFactory.h"

#include <string>

using namespace std;
//======================================================================================================================
namespace Cantera {

  //====================================================================================================================
  // Base empty constructor
  LatticeSolidPhase::LatticeSolidPhase() : 
    m_mm(0),
    m_tlast(0.0),
    m_press(-1.0),
    m_molar_density(0.0),
    m_nlattice(0),
    m_lattice(0),
    m_x(0)
  {
  }
  //====================================================================================================================
  // Copy Constructor
  /*
   * @param right Object to be copied
   */
  LatticeSolidPhase::LatticeSolidPhase(const LatticeSolidPhase &right) :
    m_mm(0),
    m_tlast(0.0),
    m_press(-1.0),
    m_molar_density(0.0),
    m_nlattice(0),
    m_lattice(0),
    m_x(0)
  {
    *this = operator=(right);
  }
  //====================================================================================================================
  // Assignment operator
  /*
   * @param right Object to be copied
   */
  LatticeSolidPhase& 
  LatticeSolidPhase::operator=(const LatticeSolidPhase& right) {
    if (&right != this) {
      ThermoPhase::operator=(right);
      m_mm = right.m_mm;
      m_tlast = right.m_tlast;
      m_press = right.m_press;
      m_molar_density = right.m_molar_density;
      m_nlattice = right.m_nlattice;
      deepStdVectorPointerCopy<LatticePhase>(right.m_lattice, m_lattice);
      m_x = right.m_x;
    }
    return *this;
  }
  //====================================================================================================================
  // Destructor
  LatticeSolidPhase::~LatticeSolidPhase() {
    for (int n = 0; n < m_nlattice; n++) {
      delete m_lattice[n];
      m_lattice[n] = 0;
    }
  }
  //====================================================================================================================
  // Duplication function
  /*
   * This virtual function is used to create a duplicate of the
   * current phase. It's used to duplicate the phase when given
   * a ThermoPhase pointer to the phase.
   *
   * @return It returns a %ThermoPhase pointer.
   */
  ThermoPhase *LatticeSolidPhase::duplMyselfAsThermoPhase() const {
    LatticeSolidPhase *igp = new LatticeSolidPhase(*this);
    return (ThermoPhase *) igp;
  }
  //====================================================================================================================
  doublereal LatticeSolidPhase::
  enthalpy_mole() const {
    _updateThermo();
    doublereal ndens, sum = 0.0;
    int n;
    for (n = 0; n < m_nlattice; n++) {
      ndens = m_lattice[n]->molarDensity();
      sum += ndens * m_lattice[n]->enthalpy_mole();
    }
    return sum/molarDensity();
  }
  //====================================================================================================================
  doublereal LatticeSolidPhase::intEnergy_mole() const {
    _updateThermo();
    doublereal ndens, sum = 0.0;
    int n;
    for (n = 0; n < m_nlattice; n++) {
      ndens = m_lattice[n]->molarDensity();
      sum += ndens * m_lattice[n]->intEnergy_mole();
    }
    return sum/molarDensity();
  }
  //====================================================================================================================
  doublereal LatticeSolidPhase::entropy_mole() const {
    _updateThermo();
    doublereal ndens, sum = 0.0;
    int n;
    for (n = 0; n < m_nlattice; n++) {
      ndens = m_lattice[n]->molarDensity();
      sum += ndens * m_lattice[n]->entropy_mole();
    }
    return sum/molarDensity();
  }
  //====================================================================================================================
  doublereal LatticeSolidPhase::gibbs_mole() const {
    _updateThermo();
    doublereal ndens, sum = 0.0;
    for (int n = 0; n < m_nlattice; n++) {
      ndens = m_lattice[n]->molarDensity(); 
      sum += ndens * m_lattice[n]->gibbs_mole();
    }
    return sum/molarDensity();
  }
  //====================================================================================================================
  doublereal LatticeSolidPhase::cp_mole() const {
    _updateThermo();
    doublereal sum = 0.0;
    for (int n = 0; n < m_nlattice; n++) {
      doublereal ndens = m_lattice[n]->molarDensity(); 
      sum += ndens * m_lattice[n]->cp_mole();
    }
    return sum/molarDensity();
  }
  //====================================================================================================================
  void LatticeSolidPhase::getActivityConcentrations(doublereal* c) const {
    _updateThermo();
    int strt = 0;
    for (int n = 0; n < m_nlattice; n++) {
      m_lattice[n]->getMoleFractions(c+strt);
      strt += m_lattice[n]->nSpecies();
    }
  }
  //====================================================================================================================
  void LatticeSolidPhase::getActivityCoefficients(doublereal* ac) const {
    for (int k = 0; k < m_kk; k++) {
      ac[k] = 1.0;
    }
  }
  //====================================================================================================================
  doublereal LatticeSolidPhase::standardConcentration(int k) const {
    return 1.0;
  }
  //====================================================================================================================
  doublereal LatticeSolidPhase::logStandardConc(int k) const {
    return 0.0;
  }

  //====================================================================================================================
  // Set the mole fractions to the specified values, and then 
  // normalize them so that they sum to 1.0 for each of the subphases
  /*
   *
   * @param x  Input vector of mole fractions. There is no restriction
   *           on the sum of the mole fraction vector. Internally,
   *           this object will pass portions of this vector to the sublattices which assume that the portions
   *           individually sum to one.
   *           Length is m_kk.
   */
  void LatticeSolidPhase::setMoleFractions(const doublereal* const x) {
    int nsp, strt = 0;
    doublereal sum = 0.0;
    for (int n = 0; n < m_nlattice; n++) {
      nsp =  m_lattice[n]->nSpecies();
      m_lattice[n]->setMoleFractions(x+strt);
      for (int k = 0; k < nsp; k++) {
	sum += x[strt + k]; 
      }
      strt += nsp;
    }
    for (int k = 0; k < strt; k++) {
      m_x[k] = x[k] / sum;
    }
    State::setMoleFractions(DATA_PTR(m_x));
  }
  //====================================================================================================================
  // Get the species mole fraction vector.
  /*
   * On output the mole fraction vector will sum to one for each of the subphases which make up this phase.
   *
   * @param x On return, x contains the mole fractions. Must have a
   *          length greater than or equal to the number of species.
   */
  void LatticeSolidPhase::getMoleFractions(doublereal* const x) const {
    int nsp, strt = 0;
    State::getMoleFractions(x);
    doublereal sum;
    for (int n = 0; n < m_nlattice; n++) {
      nsp =  m_lattice[n]->nSpecies();
      sum = 0.0;
      for (int k = 0; k < nsp; k++) {
	sum += (x + strt)[k];
      }
      for (int k = 0; k < nsp; k++) {
	(x + strt)[k] /= sum;
      }
      /*
       * At this point we can check against the mole fraction vector of the underlying LatticePhase objects and
       * get the same answer.
       */
#ifdef DEBUG_MODE
      m_lattice[n]->getMoleFractions(&(m_x[strt]));
      for (int k = 0; k < nsp; k++) {
	if (fabs((x + strt)[k] - m_x[strt+k]) > 1.0E-14) {
	  throw CanteraError("LatticeSolidPhase::getMoleFractions()",
			     "internal error");
	}
      }
#endif
      strt += nsp;
    }
  }
  //====================================================================================================================
  void LatticeSolidPhase::getChemPotentials(doublereal* mu) const {
    _updateThermo();
    int strt = 0;
    for (int n = 0; n < m_nlattice; n++) {
      doublereal dratio = m_lattice[n]->molarDensity()/molarDensity();
      m_lattice[n]->getChemPotentials(mu+strt);
      scale(mu + strt, mu + strt + m_lattice[n]->nSpecies(), mu + strt, dratio);
      strt += m_lattice[n]->nSpecies();
    }
  }
  //====================================================================================================================
  void LatticeSolidPhase::getStandardChemPotentials(doublereal* mu0) const {
    _updateThermo();
    int strt = 0;
    for (int n = 0; n < m_nlattice; n++) {
      doublereal dratio = m_lattice[n]->molarDensity()/molarDensity();
      m_lattice[n]->getStandardChemPotentials(mu0+strt);
      scale(mu0 + strt, mu0 + strt + m_lattice[n]->nSpecies(), mu0 + strt, dratio);
      strt += m_lattice[n]->nSpecies();
    }
  }
  //====================================================================================================================
  // Add in species from Slave phases
  /*
   *  This hook is used for  cSS_CONVENTION_SLAVE phases
   *
   *  @param  phaseNode    XML_Node for the current phase
   */
  void LatticeSolidPhase::installSlavePhases(Cantera::XML_Node* phaseNode) 
  {
    int m, k;
    for (int n = 0; n < m_nlattice; n++) {
      LatticePhase *lp = m_lattice[n];
      int nsp =  lp->nSpecies();
      vector<doublereal> constArr(lp->nElements());
      for (k = 0; k < nsp; k++) {
	std::string sname = lp->speciesName(k);
	std::map<std::string, double> comp;
	lp->getAtoms(k, DATA_PTR(constArr));
	for (m = 0; m < lp->nElements(); m++) {
	  if (constArr[m] != 0.0) {
	    std::string ename = lp->elementName(m);
	    comp[ename] = constArr[m];
	  }
	}
	int nel = nElements();
	vector_fp ecomp(nel, 0.0);            
	for (m = 0; m < nel; m++) {
	  double anum = comp[elementName(m)];
	  if (anum  != 0.0) {
	    ecomp[m] = anum;
	  }
	}
	double chrg = lp->charge(k);
	double sz = lp->size(k);
        addUniqueSpecies(sname, &ecomp[0], chrg, sz);
      }
    }
  }
  //====================================================================================================================
  void LatticeSolidPhase::initThermo() {
    m_kk = nSpecies();
    m_mm = nElements();
    m_x.resize(m_kk);
    int nsp, k, loc = 0;
    doublereal ndens;
    m_molar_density = 0.0;
    for (int n = 0; n < m_nlattice; n++) {
      nsp = m_lattice[n]->nSpecies();
      ndens = m_lattice[n]->molarDensity();
      for (k = 0; k < nsp; k++) {
	m_x[loc] = ndens * m_lattice[n]->moleFraction(k);
	loc++;
      }
      m_molar_density += ndens;
    }
    setMoleFractions(DATA_PTR(m_x));

    //          const vector<string>& spnames = speciesNames();
    //          int n, k, kl, namesize;
    //          int nl = m_sitedens.size();
    //          string s;
    //          m_lattice.resize(m_kk,-1);
    //          vector_fp conc(m_kk, 0.0);

    //          compositionMap xx;
    //          for (n = 0; n < nl; n++) {
    //              for (k = 0; k < m_kk; k++) { 
    //                  xx[speciesName(k)] = -1.0;
    //              }
    //              parseCompString(m_sp[n], xx);
    //              for (k = 0; k < m_kk; k++) { 
    //                  if (xx[speciesName(k)] != -1.0) {
    //                      conc[k] = m_sitedens[n]*xx[speciesName(k)];
    //                      m_lattice[k] = n;
    //                  }
    //              }

    //          }
    //          for (k = 0; k < m_kk; k++) {
    //              if (m_lattice[k] == -1) {
    //                  throw CanteraError("LatticeSolidPhase::"
    //                      "setParametersFromXML","Species "+speciesName(k)
    //                      +" not a member of any lattice.");
    //              }                    
    //          }
    //          setMoleFractions(DATA_PTR(conc));
  }

  //====================================================================================================================
  void LatticeSolidPhase::_updateThermo() const {
    doublereal tnow = temperature();
    //        if (fabs(molarDensity() - m_molar_density)/m_molar_density > 0.0001) {
    //   throw CanteraError("_updateThermo","molar density changed from "
    //        +fp2str(m_molar_density)+" to "+fp2str(molarDensity()));
    //}
    if (m_tlast != tnow) {
      int n;
      getMoleFractions(DATA_PTR(m_x));
      int strt = 0;
      for (n = 0; n < m_nlattice; n++) {
	m_lattice[n]->setTemperature(tnow);
	m_lattice[n]->setMoleFractions(DATA_PTR(m_x) + strt);
	m_lattice[n]->setPressure(m_press);
	strt += m_lattice[n]->nSpecies();
      }
      m_tlast = tnow;
    }
  }
  //====================================================================================================================
  void LatticeSolidPhase::setLatticeMoleFractions(int nn, std::string x) {
    m_lattice[nn]->setMoleFractionsByName(x);
    int n, k, loc=0, nsp;
    doublereal ndens;
    for (n = 0; n < m_nlattice; n++) {
      nsp = m_lattice[n]->nSpecies();
      ndens = m_lattice[n]->molarDensity();
      for (k = 0; k < nsp; k++) {
	m_x[loc] = ndens * m_lattice[n]->moleFraction(k);
	loc++;
      }
    }
    setMoleFractions(DATA_PTR(m_x));
  }
  //====================================================================================================================
  void LatticeSolidPhase::setParametersFromXML(const XML_Node& eosdata) {
    eosdata._require("model","LatticeSolid");
    XML_Node& la = eosdata.child("LatticeArray");
    std::vector<XML_Node*> lattices;
    la.getChildren("phase",lattices);
    int n;
    int nl = lattices.size();
    m_nlattice = nl;
    for (n = 0; n < nl; n++) {
      XML_Node& i = *lattices[n];
      m_lattice.push_back((LatticePhase*)newPhase(i));
    }
  }
  //====================================================================================================================


 doublereal LatticeSolidPhase::err(std::string msg) const {
   throw CanteraError("LatticeSolidPhase","Unimplemented " + msg);
   return 0.0;
  }

} // End namespace Cantera
//======================================================================================================================
#endif  // End #define WITH_LATTICE_SOLID
//======================================================================================================================
