/**
 *  @file Kinetics.cpp
 *      Declarations for the base class for kinetics 
 *    managers (see \ref  kineticsmgr and class 
 *  \link Cantera::Kinetics Kinetics\endlink).
 *
 *      Kinetics managers calculate rates of progress of species due to homogeneous or heterogeneous kinetics.
 */
/*
 *  $Date: 2008/12/16 20:32:18 $
 *  $Revision: 1.3 $
 */

// Copyright 2001-2004  California Institute of Technology            


             
#include "InterfaceKinetics.h"
#include "SurfPhase.h"            
#include "StoichManager.h"
#include "RateCoeffMgr.h"
                                                              
#include "ImplicitSurfChem.h"
                    
#include <iostream>
using namespace std;
                                                            
                                  
namespace Cantera {

    
  Kinetics::Kinetics() : m_ii(0), m_thermo(0),
			 m_index(-1), m_surfphase(-1), m_rxnphase(-1), 
			 m_mindim(4) {}

  Kinetics::~Kinetics(){}


  //  Copy Constructor for the %Kinetics object.
  /* 
   * Currently, this is not fully implemented. If called it will
   * throw an exception.
   */
  Kinetics::Kinetics(const Kinetics &right) :
    m_ii(0), 
    m_thermo(0),
    m_index(-1), 
    m_surfphase(-1),
    m_rxnphase(-1), 
    m_mindim(4)
  {
    /*
     * Call the assignment operator
     */
    *this = operator=(right);
  }
  
  // Assignment operator
  /*
   *  This is NOT a virtual function.
   *
   * @param right    Reference to %Kinetics object to be copied into the
   *                 current one.
   */
  Kinetics& Kinetics::
  operator=(const Kinetics &right) {
    /*
     * Check for self assignment.
     */
    if (this == &right) return *this;
    
    m_ii                = right.m_ii;
    m_perturb           = right.m_perturb;
    m_reactants         = right.m_reactants;
    m_products          = right.m_products;
   
    m_thermo            = right.m_thermo; //  DANGER -> shallow pointer copy
    
    m_start             = right.m_start;
    m_phaseindex        = right.m_phaseindex;
    m_index             = right.m_index;
    m_surfphase         = right.m_surfphase;
    m_rxnphase          = right.m_rxnphase;
    m_mindim            = right.m_mindim;
    m_dummygroups       = right.m_dummygroups;

    return *this;
  }


  // Duplication routine for objects which inherit from
  // Kinetics
  /*
   *  This virtual routine can be used to duplicate %Kinetics objects
   *  inherited from %Kinetics even if the application only has
   *  a pointer to %Kinetics to work with.
   *
   *  These routines are basically wrappers around the derived copy
   *  constructor.
   */
  Kinetics *Kinetics::duplMyselfAsKinetics() const {
    Kinetics* tp = new Kinetics(*this);
    return tp;
  }



  int Kinetics::ID() const {
    return 0;
  }

  int Kinetics::type() const {
    return 0;
  }


  /**
   * Takes as input an array of properties for all species in the
   * mechanism and copies those values beloning to a particular
   * phase to the output array.
   * @param data Input data array.
   * @param phase Pointer to one of the phase objects participating
   * in this reaction mechanism
   * @param phase_data Output array where the values for the the
   * specified phase are to be written. 
   */
  void Kinetics::selectPhase(const doublereal* data, const thermo_t* phase,
			     doublereal* phase_data) {
    int n, nsp, np = nPhases();
    for (n = 0; n < np; n++) {
      if (phase == m_thermo[n]) {
	nsp = phase->nSpecies();
	copy(data + m_start[n], 
	     data + m_start[n] + nsp, phase_data);
	return;
      }
    }
    throw CanteraError("Kinetics::selectPhase", "Phase not found.");
  }


  /**
   * kineticsSpeciesName():
   *
   * Return the string name of the kth species in the kinetics
   * manager. k is an integer from 0 to ktot - 1, where ktot is
   * the number of species in the kinetics manager, which is the
   * sum of the number of species in all phases participating in
   * the kinetics manager.  If k is out of bounds, the string
   * "<unknown>" is returned. 
   */
  string Kinetics::kineticsSpeciesName(int k) const {
    int np = m_start.size();
    for (int n = np-1; n >= 0; n--) {
      if (k >= m_start[n]) {
	return thermo(n).speciesName(k - m_start[n]);
      }
    }
    return "<unknown>";
  }

  /**
   * kineticsSpeciesIndex():
   *
   * This routine will look up a species number based on
   * the input string nm. The lookup of species will
   * occur for all phases listed in the kinetics object,
   * unless the string ph refers to a specific phase of
   * the object. 
   *
   *  return
   *   - If a match is found, the position in the species list
   *   is returned. 
   *   - If a specific phase is specified and no match is found,
   *   the value -1 is returned.
   *   - If no match is found in any phase, the value -2 is returned.
   */
  int Kinetics::kineticsSpeciesIndex(std::string nm, std::string ph) const {
    int np = static_cast<int>(m_thermo.size());
    int k;
    string id;
    for (int n = 0; n < np; n++) {
      id = thermo(n).id();
      if (ph == id) {
	k = thermo(n).speciesIndex(nm);
	if (k < 0) return -1;
	return k + m_start[n];
      }
      else if (ph == "<any>") {
	/*
	 * Call the speciesIndex() member function of the
	 * ThermoPhase object to find a match.
	 */
	k = thermo(n).speciesIndex(nm);
	if (k >= 0) return k + m_start[n];
      }                    
    }
    return -2;
  }

  /**
   * This function looks up the string name of a species and
   * returns a reference to the ThermoPhase object of the
   * phase where the species resides.
   * Will throw an error if the species string doesn't match.
   */
  thermo_t& Kinetics::speciesPhase(std::string nm) {
    int np = static_cast<int>(m_thermo.size());
    int k;
    string id;
    for (int n = 0; n < np; n++) {
      k = thermo(n).speciesIndex(nm);
      if (k >= 0) return thermo(n);
    }
    throw CanteraError("speciesPhase", "unknown species "+nm);
  }

  /**
   * This function takes as an argument the kineticsSpecies index
   * (i.e., the list index in the list of species in the kinetics
   * manager) and returns the index of the phase owning the 
   * species.
   */
  int Kinetics::speciesPhaseIndex(int k) {
    int np = m_start.size();
    for (int n = np-1; n >= 0; n--) {
      if (k >= m_start[n]) {
	return n;
      }
    }
    throw CanteraError("speciesPhaseIndex", 
		       "illegal species index: "+int2str(k));
  }

  /*
   * Add a phase to the kinetics manager object. This must
   * be done before the function init() is called or 
   * before any reactions are input.
   * The following fields are updated:
   *  m_start -> vector of integers, containing the
   *             starting position of the species for
   *             each phase in the kinetics mechanism.
   *  m_surfphase -> index of the surface phase.
   *  m_thermo -> vector of pointers to ThermoPhase phases
   *              that participate in the kinetics 
   *              mechanism.
   *  m_phaseindex -> map containing the string id of each
   *              ThermoPhase phase as a key and the
   *              index of the phase within the kinetics
   *              manager object as the value.
   */
  void Kinetics::addPhase(thermo_t& thermo) {

    // if not the first thermo object, set the start position
    // to that of the last object added + the number of its species 
    if (m_thermo.size() > 0) {
      m_start.push_back(m_start.back() 
			+ m_thermo.back()->nSpecies());
    }
    // otherwise start at 0
    else {
      m_start.push_back(0);
    }

    // the phase with lowest dimensionality is assumed to be the
    // phase/interface at which reactions take place
    if (thermo.nDim() <= m_mindim) {
      m_mindim = thermo.nDim();
      m_rxnphase = nPhases();
    }

    // there should only be one surface phase
    int ptype = -100;
    if (type() == cEdgeKinetics) ptype = cEdge;
    else if (type() == cInterfaceKinetics) ptype = cSurf;
    if (thermo.eosType() == ptype) {
      //   if (m_surfphase >= 0) {
      //    throw CanteraError("Kinetics::addPhase",
      //        "cannot add more than one surface phase");
      // }
      m_surfphase = nPhases();
      m_rxnphase = nPhases();
    }
    m_thermo.push_back(&thermo);
    m_phaseindex[m_thermo.back()->id()] = nPhases();
  }

  
  //! Private function of the class Kinetics, indicating that a function
  //!  inherited from the base class hasn't had a definition assigned to it
  /*!
   * @param m String message
   */
  void Kinetics::err(std::string m) const {
    throw CanteraError("Kinetics::" + m, 
		       "The default Base class method was called, when "
		       "the inherited class's method should "
		       "have been called");
  }

}
