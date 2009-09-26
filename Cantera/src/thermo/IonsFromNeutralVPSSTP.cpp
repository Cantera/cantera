/**
 *  @file PseudoBinaryVPSSTP.cpp
 *   Definitions for intermediate ThermoPhase object for phases which
 *   employ excess gibbs free energy formulations
 *  (see \ref thermoprops 
 * and class \link Cantera::PseudoBinaryVPSSTP PseudoBinaryVPSSTP\endlink).
 *
 * Header file for a derived class of ThermoPhase that handles
 * variable pressure standard state methods for calculating
 * thermodynamic properties that are further based upon expressions
 * for the excess gibbs free energy expressed as a function of
 * the mole fractions.
 */
/*
 * Copywrite (2009) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
/*
 *  $Date: 2009/03/27 01:08:55 $
 *  $Revision: 1.2 $
 */


#include "IonsFromNeutralVPSSTP.h"
#include "ThermoFactory.h"

#include "PDSS_IonsFromNeutral.h"
#include "mix_defs.h"

#include <cmath>

using namespace std;

#ifndef MIN
# define MIN(x,y) (( (x) < (y) ) ? (x) : (y))
#endif

namespace Cantera {

  static  const double xxSmall = 1.0E-150;

  /*
   * Default constructor.
   *
   */
  IonsFromNeutralVPSSTP::IonsFromNeutralVPSSTP() :
    GibbsExcessVPSSTP(),
    ionSolnType_(cIonSolnType_SINGLEANION),
    numNeutralMoleculeSpecies_(0),
    indexSpecialSpecies_(-1),
    indexSecondSpecialSpecies_(-1),
    numCationSpecies_(0),
    numAnionSpecies_(0),
    numPassThroughSpecies_(0),
    neutralMoleculePhase_(0),
    IOwnNThermoPhase_(true)
  {
  }

 
  // Construct and initialize an IonsFromNeutralVPSSTP object
  // directly from an asci input file
  /*
   * Working constructors
   *
   *  The two constructors below are the normal way
   *  the phase initializes itself. They are shells that call
   *  the routine initThermo(), with a reference to the
   *  XML database to get the info for the phase.
   *
   * @param inputFile Name of the input file containing the phase XML data
   *                  to set up the object
   * @param id        ID of the phase in the input file. Defaults to the
   *                  empty string.
   * @param neutralPhase   The object takes a neutralPhase ThermoPhase
   *                       object as input. It can either take a pointer
   *                       to an existing object in the parameter list,
   *                       in which case it does not own the object, or
   *                       it can construct a neutral Phase as a slave
   *                       object, in which case, it does own the slave
   *                       object, for purposes of who gets to destroy
   *                       the object.
   *                       If this parameter is zero, then a slave
   *                       neutral phase object is created and used.
   */
  IonsFromNeutralVPSSTP::IonsFromNeutralVPSSTP(std::string inputFile, std::string id,
					       ThermoPhase *neutralPhase) :
    GibbsExcessVPSSTP(),
    ionSolnType_(cIonSolnType_SINGLEANION),
    numNeutralMoleculeSpecies_(0),
    indexSpecialSpecies_(-1),
    indexSecondSpecialSpecies_(-1),
    numCationSpecies_(0),
    numAnionSpecies_(0),
    numPassThroughSpecies_(0),
    neutralMoleculePhase_(neutralPhase),
    IOwnNThermoPhase_(true)
  {
    if (neutralPhase) {
      IOwnNThermoPhase_ = false;
    }
    constructPhaseFile(inputFile, id);
  }

  IonsFromNeutralVPSSTP::IonsFromNeutralVPSSTP(XML_Node& phaseRoot, std::string id,
					       ThermoPhase *neutralPhase) :
    GibbsExcessVPSSTP(),
    ionSolnType_(cIonSolnType_SINGLEANION),
    numNeutralMoleculeSpecies_(0),
    indexSpecialSpecies_(-1),
    indexSecondSpecialSpecies_(-1),
    numCationSpecies_(0),
    numAnionSpecies_(0),
    numPassThroughSpecies_(0),
    neutralMoleculePhase_(neutralPhase),
    IOwnNThermoPhase_(true)
  {
    if (neutralPhase) {
      IOwnNThermoPhase_ = false;
    }
    constructPhaseXML(phaseRoot, id);
  }



  /*
   * Copy Constructor:
   *
   *  Note this stuff will not work until the underlying phase
   *  has a working copy constructor
   */
  IonsFromNeutralVPSSTP::IonsFromNeutralVPSSTP(const IonsFromNeutralVPSSTP &b) :
    GibbsExcessVPSSTP(),
    ionSolnType_(cIonSolnType_SINGLEANION),
    numNeutralMoleculeSpecies_(0),
    indexSpecialSpecies_(-1),
    indexSecondSpecialSpecies_(-1),
    numCationSpecies_(0),
    numAnionSpecies_(0),
    numPassThroughSpecies_(0),
    neutralMoleculePhase_(0),
    IOwnNThermoPhase_(true)
  {
    *this = operator=(b);
  }

  /*
   * operator=()
   *
   *  Note this stuff will not work until the underlying phase
   *  has a working assignment operator
   */
  IonsFromNeutralVPSSTP& IonsFromNeutralVPSSTP::
  operator=(const IonsFromNeutralVPSSTP &b) {
    if (&b != this) {
      GibbsExcessVPSSTP::operator=(b);
    }  

    ionSolnType_                = b.ionSolnType_;
    numNeutralMoleculeSpecies_  = b.numNeutralMoleculeSpecies_;
    indexSpecialSpecies_        = b.indexSpecialSpecies_;
    indexSecondSpecialSpecies_  = b.indexSecondSpecialSpecies_;
    NeutralMolecMoleFractions_  = b.NeutralMolecMoleFractions_;
    cationList_                 = b.cationList_;
    numCationSpecies_           = b.numCationSpecies_;
    anionList_                  = b.anionList_;
    numAnionSpecies_            = b.numAnionSpecies_;
    passThroughList_            = b.passThroughList_;
    numPassThroughSpecies_      = b.numPassThroughSpecies_;
    /*
     *  This is a shallow copy. We need to figure this out
     */
    neutralMoleculePhase_       = b.neutralMoleculePhase_;
    if (neutralMoleculePhase_) {
      exit(-1);
    }
    IOwnNThermoPhase_           = b.IOwnNThermoPhase_;

    moleFractionsTmp_           = b.moleFractionsTmp_;

    return *this;
  }

  /**
   *
   * ~IonsFromNeutralVPSSTP():   (virtual)
   *
   * Destructor: does nothing:
   *
   */
  IonsFromNeutralVPSSTP::~IonsFromNeutralVPSSTP() {
    if (IOwnNThermoPhase_) {
      delete neutralMoleculePhase_;
      neutralMoleculePhase_=0;
    }
  }

  /*
   * This routine duplicates the current object and returns
   * a pointer to ThermoPhase.
   */
  ThermoPhase* 
  IonsFromNeutralVPSSTP::duplMyselfAsThermoPhase() const {
    IonsFromNeutralVPSSTP* mtp = new IonsFromNeutralVPSSTP(*this);
    return (ThermoPhase *) mtp;
  }

  /*
   *  -------------- Utilities -------------------------------
   */

 
  // Equation of state type flag.
  /*
   * The ThermoPhase base class returns
   * zero. Subclasses should define this to return a unique
   * non-zero value. Known constants defined for this purpose are
   * listed in mix_defs.h. The IonsFromNeutralVPSSTP class also returns
   * zero, as it is a non-complete class.
   */
  int IonsFromNeutralVPSSTP::eosType() const { 
    return cIonsFromNeutral;
  }

 

  /*
   * ------------ Molar Thermodynamic Properties ----------------------
   */
  /*
   * Molar enthalpy of the solution. Units: J/kmol.
   */
  doublereal IonsFromNeutralVPSSTP::enthalpy_mole() const {
    getPartialMolarEnthalpies(DATA_PTR(m_pp));
    return mean_X(DATA_PTR(m_pp));
  }

  /**
   * Molar internal energy of the solution. Units: J/kmol.
   *
   * This is calculated from the soln enthalpy and then
   * subtracting pV.
   */
  doublereal IonsFromNeutralVPSSTP::intEnergy_mole() const {
    double hh = enthalpy_mole();
    double pres = pressure();
    double molarV = 1.0/molarDensity();
    double uu = hh - pres * molarV;
    return uu;
  }

  /**
   *  Molar soln entropy at constant pressure. Units: J/kmol/K.
   *
   *  This is calculated from the partial molar entropies.
   */
  doublereal IonsFromNeutralVPSSTP::entropy_mole() const {
    getPartialMolarEntropies(DATA_PTR(m_pp));
    return mean_X(DATA_PTR(m_pp));
  }

  /// Molar Gibbs function. Units: J/kmol.
  doublereal IonsFromNeutralVPSSTP::gibbs_mole() const {
    getChemPotentials(DATA_PTR(m_pp));
    return mean_X(DATA_PTR(m_pp));
  }
 /** Molar heat capacity at constant pressure. Units: J/kmol/K.
   *
   * Returns the solution heat capacition at constant pressure.
   * This is calculated from the partial molar heat capacities.
   */
  doublereal IonsFromNeutralVPSSTP::cp_mole() const {
    getPartialMolarCp(DATA_PTR(m_pp));
    double val = mean_X(DATA_PTR(m_pp));
    return val;
  }

  /// Molar heat capacity at constant volume. Units: J/kmol/K.
  doublereal IonsFromNeutralVPSSTP::cv_mole() const {
    // Need to revisit this, as it is wrong
    getPartialMolarCp(DATA_PTR(m_pp));
    return mean_X(DATA_PTR(m_pp));
    //err("not implemented");
    //return 0.0;
  }

  /*
   * - Activities, Standard States, Activity Concentrations -----------
   */

  // This method returns an array of generalized concentrations
  /*
   * \f$ C^a_k\f$ are defined such that \f$ a_k = C^a_k /
   * C^0_k, \f$ where \f$ C^0_k \f$ is a standard concentration
   * defined below and \f$ a_k \f$ are activities used in the
   * thermodynamic functions.  These activity (or generalized)
   * concentrations are used
   * by kinetics manager classes to compute the forward and
   * reverse rates of elementary reactions. Note that they may
   * or may not have units of concentration --- they might be
   * partial pressures, mole fractions, or surface coverages,
   * for example.
   *
   *  Here we define the activity concentrations as equal
   *  to the activities, because the standard concentration is 1.
   *
   * @param c Output array of generalized concentrations. The
   *           units depend upon the implementation of the
   *           reaction rate expressions within the phase.
   */
  void IonsFromNeutralVPSSTP::getActivityConcentrations(doublereal* c) const {
    getActivities(c);
  }
 
  // Return the standard concentration for the kth species
  /*
   * The standard concentration \f$ C^0_k \f$ used to normalize
   * the activity (i.e., generalized) concentration. In many cases, this quantity
   * will be the same for all species in a phase - for example,
   * for an ideal gas \f$ C^0_k = P/\hat R T \f$. For this
   * reason, this method returns a single value, instead of an
   * array.  However, for phases in which the standard
   * concentration is species-specific (e.g. surface species of
   * different sizes), this method may be called with an
   * optional parameter indicating the species.
   *
   * @param k Optional parameter indicating the species. The default
   *          is to assume this refers to species 0.
   * @return 
   *   Returns the standard concentration. The units are by definition
   *   dependent on the ThermoPhase and kinetics manager representation.
   */
  doublereal IonsFromNeutralVPSSTP::standardConcentration(int k) const {
    return 1.0;
  }

  // Natural logarithm of the standard concentration of the kth species.
  /*
   * @param k    index of the species (defaults to zero)
   */
  doublereal IonsFromNeutralVPSSTP::logStandardConc(int k) const {
    return 0.0;
  }

  // Returns the units of the standard and generalized concentrations.
  /*
   * Note they have the same units, as their
   * ratio is defined to be equal to the activity of the kth
   * species in the solution, which is unitless.
   *
   * This routine is used in print out applications where the
   * units are needed. Usually, MKS units are assumed throughout
   * the program and in the XML input files.
   *
   * The base %ThermoPhase class assigns the default quantities
   * of (kmol/m3) for all species.
   * Inherited classes are responsible for overriding the default 
   * values if necessary.
   *
   * @param uA Output vector containing the units
   *  uA[0] = kmol units - default  = 1
   *  uA[1] = m    units - default  = -nDim(), the number of spatial
   *                                dimensions in the Phase class.
   *  uA[2] = kg   units - default  = 0;
   *  uA[3] = Pa(pressure) units - default = 0;
   *  uA[4] = Temperature units - default = 0;
   *  uA[5] = time units - default = 0
   * @param k species index. Defaults to 0.
   * @param sizeUA output int containing the size of the vector.
   *        Currently, this is equal to 6.
   */
  void IonsFromNeutralVPSSTP::getUnitsStandardConc(double *uA, int k,
						   int sizeUA) const {
    uA[0] = 0;
    uA[1] = 0;
    uA[2] = 0;
    uA[3] = 0;
    uA[4] = 0;
    uA[5] = 0;
  }

  // Get the array of non-dimensional molar-based activity coefficients at
  // the current solution temperature, pressure, and solution concentration.
  /*
   * @param ac Output vector of activity coefficients. Length: m_kk.
   */
  void IonsFromNeutralVPSSTP::getActivityCoefficients(doublereal* ac) const {

    // This stuff has moved to the setState routines
    //   calcNeutralMoleculeMoleFractions();
    //   neutralMoleculePhase_->setState_TPX(temperature(), pressure(), DATA_PTR(NeutralMolecMoleFractions_)); 
    //   neutralMoleculePhase_->getStandardChemPotentials(DATA_PTR(muNeutralMolecule_));

    /*
     * Update the activity coefficients
     */
    s_update_lnActCoeff();

    /*
     * take the exp of the internally storred coefficients.
     */
    for (int k = 0; k < m_kk; k++) {
      ac[k] = exp(lnActCoeff_Scaled_[k]);      
    }
  }

  /*
   * ---------  Partial Molar Properties of the Solution -------------------------------
   */

  // Get the species chemical potentials. Units: J/kmol.
  /*
   * This function returns a vector of chemical potentials of the
   * species in solution at the current temperature, pressure
   * and mole fraction of the solution.
   *
   * @param mu  Output vector of species chemical
   *            potentials. Length: m_kk. Units: J/kmol
   */
  void  
  IonsFromNeutralVPSSTP::getChemPotentials(doublereal* mu) const {
    int k, icat, jNeut;
    doublereal xx, fact2;
    /*
     *  Transfer the mole fractions to the slave neutral molecule
     *  phase 
     *   Note we may move this in the future.
     */
    //calcNeutralMoleculeMoleFractions();
    //neutralMoleculePhase_->setState_TPX(temperature(), pressure(), DATA_PTR(NeutralMolecMoleFractions_)); 

    /*
     * Get the standard chemical potentials of netural molecules
     */
    neutralMoleculePhase_->getStandardChemPotentials(DATA_PTR(muNeutralMolecule_));

    doublereal RT_ = GasConstant * temperature();

    switch (ionSolnType_) {
    case cIonSolnType_PASSTHROUGH:
      neutralMoleculePhase_->getChemPotentials(mu);
      break;
    case cIonSolnType_SINGLEANION:
      neutralMoleculePhase_->getActivityCoefficients(DATA_PTR(gammaNeutralMolecule_));

      fact2 = 2.0 * RT_ * log(2.0);

      // Do the cation list
      for (k = 0; k < (int) cationList_.size(); k++) {
	//! Get the id for the next cation
        icat = cationList_[k];
	jNeut = fm_invert_ionForNeutral[icat];
	xx = fmaxx(SmallNumber, moleFractions_[icat]);
	mu[icat] = muNeutralMolecule_[jNeut] + fact2 + RT_ * log(gammaNeutralMolecule_[jNeut] * xx);
      }

      // Do the anion list
      icat = anionList_[0];
      jNeut = fm_invert_ionForNeutral[icat];
      xx = fmaxx(SmallNumber, moleFractions_[icat]);
      mu[icat] = RT_ * log(xx);

      // Do the list of neutral molecules
      for (k = 0; k <  numPassThroughSpecies_; k++) {
	icat = passThroughList_[k];
	jNeut = fm_invert_ionForNeutral[icat];
	xx = fmaxx(SmallNumber, moleFractions_[icat]);
	mu[icat] = muNeutralMolecule_[jNeut] + RT_ * log( gammaNeutralMolecule_[jNeut] * xx);
      }
      break;
 
    case cIonSolnType_SINGLECATION:
      throw CanteraError("eosType", "Unknown type");
      break;
    case cIonSolnType_MULTICATIONANION:
      throw CanteraError("eosType", "Unknown type");
      break;
    default:
      throw CanteraError("eosType", "Unknown type");
      break;
    }
  }


  // Returns an array of partial molar enthalpies for the species
  // in the mixture.
  /*
   * Units (J/kmol)
   *
   * For this phase, the partial molar enthalpies are equal to the
   * standard state enthalpies modified by the derivative of the
   * molality-based activity coefficent wrt temperature
   *
   *  \f[
   * \bar h_k(T,P) = h^o_k(T,P) - R T^2 \frac{d \ln(\gamma_k)}{dT}
   * \f]
   *
   */
  void IonsFromNeutralVPSSTP::getPartialMolarEnthalpies(doublereal* hbar) const {
   /*
     * Get the nondimensional standard state enthalpies
     */
    getEnthalpy_RT(hbar);
    /*
     * dimensionalize it.
     */
    double T = temperature();
    double RT = GasConstant * T;
    for (int k = 0; k < m_kk; k++) {
      hbar[k] *= RT;
    }
    /*
     * Update the activity coefficients, This also update the
     * internally storred molalities.
     */
    s_update_lnActCoeff();
    s_update_dlnActCoeffdT();
    double RTT = RT * T;
    for (int k = 0; k < m_kk; k++) {
      hbar[k] -= RTT * dlnActCoeffdT_Scaled_[k];
    }
  }

  // Returns an array of partial molar entropies for the species
  // in the mixture.
  /*
   * Units (J/kmol)
   *
   * For this phase, the partial molar enthalpies are equal to the
   * standard state enthalpies modified by the derivative of the
   * activity coefficent wrt temperature
   *
   *  \f[
   * \bar s_k(T,P) = s^o_k(T,P) - R T^2 \frac{d \ln(\gamma_k)}{dT}
   * \f]
   *
   */
  void IonsFromNeutralVPSSTP::getPartialMolarEntropies(doublereal* sbar) const {
    double xx;
    /*
     * Get the nondimensional standard state entropies
     */
    getEntropy_R(sbar);
    double T = temperature();
    /*
     * Update the activity coefficients, This also update the
     * internally storred molalities.
     */
    s_update_lnActCoeff();
    s_update_dlnActCoeffdT();

    for (int k = 0; k < m_kk; k++) {
      xx = fmaxx(moleFractions_[k], xxSmall);
      sbar[k] += - lnActCoeff_Scaled_[k] -log(xx) - T * dlnActCoeffdT_Scaled_[k];
    }  
    /*
     * dimensionalize it.
     */
   for (int k = 0; k < m_kk; k++) {
      sbar[k] *= GasConstant;
    }
  }


  // This is temporary. We will get rid of this
  void IonsFromNeutralVPSSTP::setTemperature(doublereal t) {
    double p = pressure();
    IonsFromNeutralVPSSTP::setState_TP(t, p);
  }

 // This is temporary. We will get rid of this
  void IonsFromNeutralVPSSTP::setPressure(doublereal p) {
    double t = temperature();
    IonsFromNeutralVPSSTP::setState_TP(t, p);
  }

  // Set the temperature (K) and pressure (Pa)
  /*
   * Setting the pressure may involve the solution of a nonlinear equation.
   *
   * @param t    Temperature (K)
   * @param p    Pressure (Pa)
   */
  void IonsFromNeutralVPSSTP::setState_TP(doublereal t, doublereal p) {
    /*
     *  This is a two phase process. First, we calculate the standard states 
     *  within the neutral molecule phase.
     */
    neutralMoleculePhase_->setState_TP(t, p);
    VPStandardStateTP::setState_TP(t,p);

    /*
     * Calculate the partial molar volumes, and then the density of the fluid
     */
    calcDensity();
  }
  
  // Calculate ion mole fractions from neutral molecule 
  // mole fractions.
  /*
   *  @param mf Dump the mole fractions into this vector.
   */
  void IonsFromNeutralVPSSTP::calcIonMoleFractions(doublereal * const mf) const {
    int k;
    doublereal fmij;
    /*
     * Download the neutral mole fraction vector into the
     * vector, NeutralMolecMoleFractions_[]
     */
    neutralMoleculePhase_->getMoleFractions(DATA_PTR(NeutralMolecMoleFractions_)); 
  
    // Zero the mole fractions
    fbo_zero_dbl_1(mf, m_kk);
  
    /*
     *  Use the formula matrix to calculate the relative mole numbers.
     */
    for (int jNeut = 0; jNeut <  numNeutralMoleculeSpecies_; jNeut++) {
      for (k = 0; k < m_kk; k++) {
	fmij =  fm_neutralMolec_ions_[k + jNeut * m_kk];
	mf[k] += fmij * NeutralMolecMoleFractions_[jNeut];
      }
    }

    /*
     * Normalize the new mole fractions
     */
    doublereal sum = 0.0;
     for (k = 0; k < m_kk; k++) {
       sum += mf[k];
     }
     for (k = 0; k < m_kk; k++) {
       mf[k] /= sum; 
     }

  }

  // Calculate neutral molecule mole fractions
  /*
   *  This routine calculates the neutral molecule mole
   *  fraction given the vector of ion mole fractions,
   *  i.e., the mole fractions from this ThermoPhase.
   *  Note, this routine basically assumes that there
   *  is charge neutrality. If there isn't, then it wouldn't
   *  make much sense. 
   *
   *  for the case of  cIonSolnType_SINGLEANION, some slough
   *  in the charge neutrality is allowed. The cation number
   *  is followed, while the difference in charge neutrality
   *  is dumped into the anion mole number to fix the imbalance.
   */
  void IonsFromNeutralVPSSTP::calcNeutralMoleculeMoleFractions() const {
    int k, icat, jNeut;
    doublereal sumCat; 
    doublereal sumAnion;
    doublereal fmij;
    doublereal sum = 0.0;

    //! Zero the vector we are trying to find.
    for (k = 0; k < numNeutralMoleculeSpecies_; k++) {
      NeutralMolecMoleFractions_[k] = 0.0;
    }
#ifdef DEBUG_MODE
    sum = -1.0;
    for (k = 0; k < m_kk; k++) {
     sum += moleFractions_[k];
    }
    if (fabs(sum) > 1.0E-11) {
      throw CanteraError("IonsFromNeutralVPSSTP::calcNeutralMoleculeMoleFractions", 
			 "molefracts don't sum to one: " + fp2str(sum));
    }
#endif

    // bool fmSimple = true;

    switch (ionSolnType_) {

    case cIonSolnType_PASSTHROUGH:

      for (k = 0; k < m_kk; k++) {
	NeutralMolecMoleFractions_[k] = moleFractions_[k];
      }
      break;

    case cIonSolnType_SINGLEANION:

      sumCat = 0.0;
      sumAnion = 0.0;
      for (k = 0; k < numNeutralMoleculeSpecies_; k++) {
	NeutralMolecMoleFractions_[k] = 0.0;
      }

      for (k = 0; k < (int) cationList_.size(); k++) {
	//! Get the id for the next cation
        icat = cationList_[k];
	jNeut = fm_invert_ionForNeutral[icat];
	if (jNeut >= 0) {
	  fmij =  fm_neutralMolec_ions_[icat + jNeut * m_kk];
	  AssertTrace(fmij != 0.0);
	  NeutralMolecMoleFractions_[jNeut] += moleFractions_[icat] / fmij;
	}
      }

      for (k = 0; k <  numPassThroughSpecies_; k++) {
	icat = passThroughList_[k];
	jNeut = fm_invert_ionForNeutral[icat];
	fmij = fm_neutralMolec_ions_[ icat + jNeut * m_kk];
	NeutralMolecMoleFractions_[jNeut] += moleFractions_[icat] / fmij;
      }

#ifdef DEBUG_MODE
      for (k = 0; k < m_kk; k++) {
	moleFractionsTmp_[k] = moleFractions_[k];
      }
      for (jNeut = 0; jNeut <  numNeutralMoleculeSpecies_; jNeut++) {
	for (k = 0; k < m_kk; k++) {
	  fmij =  fm_neutralMolec_ions_[k + jNeut * m_kk];
	  moleFractionsTmp_[k] -= fmij * NeutralMolecMoleFractions_[jNeut];
	}
      }
      for (k = 0; k < m_kk; k++) {
	if (fabs(moleFractionsTmp_[k]) > 1.0E-13) {
	  //! Check to see if we have in fact found the inverse.
	  if (anionList_[0] != k) {
	    throw CanteraError("", "neutral molecule calc error");
	  } else {
	    //! For the single anion case, we will allow some slippage
	    if (fabs(moleFractionsTmp_[k]) > 1.0E-5) {
	      throw CanteraError("", "neutral molecule calc error - anion");
	    }
	  }
	}
      }
#endif      

      // Normalize the Neutral Molecule mole fractions
      sum = 0.0;
      for (k = 0; k < numNeutralMoleculeSpecies_; k++) {
        sum += NeutralMolecMoleFractions_[k];
      }
      for (k = 0; k < numNeutralMoleculeSpecies_; k++) {
        NeutralMolecMoleFractions_[k] /= sum;
      }

      break;

    case  cIonSolnType_SINGLECATION:

      throw CanteraError("eosType", "Unknown type");
     
      break;
     
    case  cIonSolnType_MULTICATIONANION:

      throw CanteraError("eosType", "Unknown type");
      break;

    default:

      throw CanteraError("eosType", "Unknown type");
      break;

    } 
  }

  void IonsFromNeutralVPSSTP::setMassFractions(const doublereal* const y) {
    GibbsExcessVPSSTP::setMassFractions(y);
    calcNeutralMoleculeMoleFractions();
    neutralMoleculePhase_->setMoleFractions(DATA_PTR(NeutralMolecMoleFractions_));
  }

  void IonsFromNeutralVPSSTP::setMassFractions_NoNorm(const doublereal* const y) {
    GibbsExcessVPSSTP::setMassFractions_NoNorm(y);
    calcNeutralMoleculeMoleFractions();
   neutralMoleculePhase_->setMoleFractions(DATA_PTR(NeutralMolecMoleFractions_));
  }

  void IonsFromNeutralVPSSTP::setMoleFractions(const doublereal* const y) {
    GibbsExcessVPSSTP::setMoleFractions(y);
    calcNeutralMoleculeMoleFractions();
    neutralMoleculePhase_->setMoleFractions(DATA_PTR(NeutralMolecMoleFractions_));
  }

  void IonsFromNeutralVPSSTP::setMoleFractions_NoNorm(const doublereal* const y) {
    GibbsExcessVPSSTP::setMoleFractions_NoNorm(y);
    calcNeutralMoleculeMoleFractions();
    neutralMoleculePhase_->setMoleFractions(DATA_PTR(NeutralMolecMoleFractions_));
  }


  void IonsFromNeutralVPSSTP::setConcentrations(const doublereal* const y) {
    GibbsExcessVPSSTP::setConcentrations(y);
    calcNeutralMoleculeMoleFractions();
    neutralMoleculePhase_->setMoleFractions(DATA_PTR(NeutralMolecMoleFractions_));
  }

  /*
   * ------------ Partial Molar Properties of the Solution ------------
   */


  doublereal IonsFromNeutralVPSSTP::err(std::string msg) const {
    throw CanteraError("IonsFromNeutralVPSSTP","Base class method "
		       +msg+" called. Equation of state type: "+int2str(eosType()));
    return 0;
  }
  /*
   *   Import, construct, and initialize a phase
   *   specification from an XML tree into the current object.
   *
   * This routine is a precursor to constructPhaseXML(XML_Node*)
   * routine, which does most of the work.
   *
   * @param infile XML file containing the description of the
   *        phase
   *
   * @param id  Optional parameter identifying the name of the
   *            phase. If none is given, the first XML
   *            phase element will be used.
   */
 void IonsFromNeutralVPSSTP::constructPhaseFile(std::string inputFile, std::string id) {

    if (inputFile.size() == 0) {
      throw CanteraError("MargulesVPSSTP:constructPhaseFile",
                         "input file is null");
    }
    string path = findInputFile(inputFile);
    std::ifstream fin(path.c_str());
    if (!fin) {
      throw CanteraError("MargulesVPSSTP:constructPhaseFile","could not open "
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
      throw CanteraError("MargulesVPSSTP:constructPhaseFile",
                         "ERROR: Can not find phase named " +
                         id + " in file named " + inputFile);
    }
    fxml_phase->copy(&phaseNode_XML);
    constructPhaseXML(*fxml_phase, id);
    delete fxml;
  }

  /*
   *   Import, construct, and initialize a HMWSoln phase
   *   specification from an XML tree into the current object.
   *
   *   Most of the work is carried out by the cantera base
   *   routine, importPhase(). That routine imports all of the
   *   species and element data, including the standard states
   *   of the species.
   *
   *   Then, In this routine, we read the information
   *   particular to the specification of the activity
   *   coefficient model for the Pitzer parameterization.
   *
   *   We also read information about the molar volumes of the
   *   standard states if present in the XML file.
   *
   * @param phaseNode This object must be the phase node of a
   *             complete XML tree
   *             description of the phase, including all of the
   *             species data. In other words while "phase" must
   *             point to an XML phase object, it must have
   *             sibling nodes "speciesData" that describe
   *             the species in the phase.
   * @param id   ID of the phase. If nonnull, a check is done
   *             to see if phaseNode is pointing to the phase
   *             with the correct id.
   */
  void IonsFromNeutralVPSSTP::constructPhaseXML(XML_Node& phaseNode, std::string id) {
    string stemp;
    if (id.size() > 0) {
      string idp = phaseNode.id();
      if (idp != id) {
        throw CanteraError("IonsFromNeutralVPSSTP::constructPhaseXML",
                           "phasenode and Id are incompatible");
      }
    }

    /*
     * Find the Thermo XML node
     */
    if (!phaseNode.hasChild("thermo")) {
      throw CanteraError("IonsFromNeutralVPSSTP::constructPhaseXML",
                         "no thermo XML node");
    }
    XML_Node& thermoNode = phaseNode.child("thermo");


   
    /*
     * Make sure that the thermo model is IonsFromNeutralMolecule
     */
    stemp = thermoNode.attrib("model");
    string formString = lowercase(stemp);
    if (formString != "ionsfromneutralmolecule") {
      throw CanteraError("IonsFromNeutralVPSSTP::constructPhaseXML",
                         "model name isn't IonsFromNeutralMolecule: " + formString);
    }

    /*
     * Find the Neutral Molecule Phase
     */
    if (!thermoNode.hasChild("neutralMoleculePhase")) {
      throw CanteraError("IonsFromNeutralVPSSTP::constructPhaseXML",
                         "no neutralMoleculePhase XML node");
    }
    XML_Node& neutralMoleculeNode = thermoNode.child("neutralMoleculePhase");

    string nsource = neutralMoleculeNode["datasrc"];
    XML_Node *neut_ptr = get_XML_Node(nsource, 0);
    if (!neut_ptr) {
      throw CanteraError("IonsFromNeutralVPSSTP::constructPhaseXML",
                         "neut_ptr = 0");
    }

    /*
     *  Create the neutralMolecule ThermoPhase if we haven't already
     */
    if (!neutralMoleculePhase_) {
      neutralMoleculePhase_  = newPhase(*neut_ptr);
    }

    /*
     * Call the Cantera importPhase() function. This will import
     * all of the species into the phase. This will also handle
     * all of the solvent and solute standard states
     */
    bool m_ok = importPhase(phaseNode, this);
    if (!m_ok) {
      throw CanteraError("IonsFromNeutralVPSSTP::constructPhaseXML",
			 "importPhase failed ");
    }

  }


  /*
   * @internal Initialize. This method is provided to allow
   * subclasses to perform any initialization required after all
   * species have been added. For example, it might be used to
   * resize internal work arrays that must have an entry for
   * each species.  The base class implementation does nothing,
   * and subclasses that do not require initialization do not
   * need to overload this method.  When importing a CTML phase
   * description, this method is called just prior to returning
   * from function importPhase.
   *
   * @see importCTML.cpp
   */
  void IonsFromNeutralVPSSTP::initThermo() {
    initLengths();
    GibbsExcessVPSSTP::initThermo();
  }


  //   Initialize lengths of local variables after all species have
  //   been identified.
  void  IonsFromNeutralVPSSTP::initLengths() {
    m_kk = nSpecies();
    numNeutralMoleculeSpecies_ =  neutralMoleculePhase_->nSpecies();
    moleFractions_.resize(m_kk);
    fm_neutralMolec_ions_.resize(numNeutralMoleculeSpecies_ * m_kk);
    fm_invert_ionForNeutral.resize(m_kk);
    NeutralMolecMoleFractions_.resize(numNeutralMoleculeSpecies_);
    cationList_.resize(m_kk);
    anionList_.resize(m_kk);
    passThroughList_.resize(m_kk);
    moleFractionsTmp_.resize(m_kk);
    muNeutralMolecule_.resize(numNeutralMoleculeSpecies_);
    gammaNeutralMolecule_.resize(numNeutralMoleculeSpecies_);
    dlnActCoeffdT_NeutralMolecule_.resize(numNeutralMoleculeSpecies_);
  }

  static double factorOverlap(const std::vector<std::string>&  elnamesVN ,
			      const std::vector<double>& elemVectorN,
			      const int  nElementsN,
			      const std::vector<std::string>&  elnamesVI ,
			      const std::vector<double>& elemVectorI,
			      const int  nElementsI)
  {
    double fMax = 1.0E100;
    for (int mi = 0; mi < nElementsI; mi++) {
      if (elnamesVI[mi] != "E") {
	if (elemVectorI[mi] > 1.0E-13) {
	  double eiNum = elemVectorI[mi];
	  for (int mn = 0; mn < nElementsN; mn++) {
	    if (elnamesVI[mi] == elnamesVN[mn]) {
	      if (elemVectorN[mn] <= 1.0E-13) {
		return 0.0;
	      }
	      fMax = MIN(fMax, elemVectorN[mn]/eiNum);
	    }
	  }
	}
      }
    }
    return fMax;
  }

  /*
   * initThermoXML()                (virtual from ThermoPhase)
   *   Import and initialize a ThermoPhase object
   *
   * @param phaseNode This object must be the phase node of a
   *             complete XML tree
   *             description of the phase, including all of the
   *             species data. In other words while "phase" must
   *             point to an XML phase object, it must have
   *             sibling nodes "speciesData" that describe
   *             the species in the phase.
   * @param id   ID of the phase. If nonnull, a check is done
   *             to see if phaseNode is pointing to the phase
   *             with the correct id. 
   */
  void IonsFromNeutralVPSSTP::initThermoXML(XML_Node& phaseNode, std::string id) {
    int k;
    /*
     *   variables that need to be populated
     *
     *    cationList_
     *      numCationSpecies_;
     */
 
    numCationSpecies_ = 0;
    cationList_.clear();
    for (k = 0; k < m_kk; k++) {
      if (charge(k) > 0) {
	cationList_.push_back(k);
	numCationSpecies_++;
      }
    }

    numAnionSpecies_ = 0;
    anionList_.clear();
    for (k = 0; k < m_kk; k++) {
      if (charge(k) < 0) {
	anionList_.push_back(k);
	numAnionSpecies_++;
      }
    }

    numPassThroughSpecies_= 0;
    passThroughList_.clear();
    for (k = 0; k < m_kk; k++) {
      if (charge(k) == 0) {
       passThroughList_.push_back(k);
       numPassThroughSpecies_++;
      }
    }

    PDSS_IonsFromNeutral *speciesSS = 0;
    indexSpecialSpecies_ = -1;
    for (k = 0; k < m_kk; k++) {
      speciesSS = dynamic_cast<PDSS_IonsFromNeutral *>(providePDSS(k));
      if (!speciesSS) {
	throw CanteraError("initThermoXML", "Dynamic cast failed");
      }
      if (speciesSS->specialSpecies_ == 1) {
	indexSpecialSpecies_ = k;
      }
      if (speciesSS->specialSpecies_ == 2) {
	indexSecondSpecialSpecies_ = k;
      }
    }


    int nElementsN =  neutralMoleculePhase_->nElements();
    const std::vector<std::string>&  elnamesVN = neutralMoleculePhase_->elementNames();
    std::vector<double> elemVectorN(nElementsN);
    std::vector<double> elemVectorN_orig(nElementsN);

    int nElementsI =  nElements();
    const std::vector<std::string>&  elnamesVI = elementNames();
    std::vector<double> elemVectorI(nElementsI);

    vector<doublereal> fm_tmp(m_kk);
    for (int jNeut = 0; jNeut <  numNeutralMoleculeSpecies_; jNeut++) {
      fm_invert_ionForNeutral[jNeut] = -1;
    }
    for (int jNeut = 0; jNeut <  numNeutralMoleculeSpecies_; jNeut++) {
      for (int m = 0; m < nElementsN; m++) {
	 elemVectorN[m] = neutralMoleculePhase_->nAtoms(jNeut, m);
      }
      elemVectorN_orig = elemVectorN;
      fvo_zero_dbl_1(fm_tmp, m_kk);

      for (int m = 0; m < nElementsI; m++) {
	 elemVectorI[m] = nAtoms(indexSpecialSpecies_, m);
      }
      double fac = factorOverlap(elnamesVN, elemVectorN, nElementsN,
				 elnamesVI ,elemVectorI, nElementsI);
      if (fac > 0.0) {
	for (int m = 0; m < nElementsN; m++) {
	  std::string mName = elnamesVN[m];
	  for (int mi = 0; mi < nElementsI; mi++) {
	    std::string eName = elnamesVI[mi];
	    if (mName == eName) {
	      elemVectorN[m] -= fac * elemVectorI[mi];
	    }
	     
	  }
	}
      }
      fm_neutralMolec_ions_[indexSpecialSpecies_  + jNeut * m_kk ] += fac;
    
    
      for (k = 0; k < m_kk; k++) {
	for (int m = 0; m < nElementsI; m++) {
	  elemVectorI[m] = nAtoms(k, m);
	}
	double fac = factorOverlap(elnamesVN, elemVectorN, nElementsN,
				   elnamesVI ,elemVectorI, nElementsI);
	if (fac > 0.0) {
	  for (int m = 0; m < nElementsN; m++) {
	    std::string mName = elnamesVN[m];
	    for (int mi = 0; mi < nElementsI; mi++) {
	      std::string eName = elnamesVI[mi];
	      if (mName == eName) {
		elemVectorN[m] -= fac * elemVectorI[mi];
	      }
	     
	    }
	  }
	  bool notTaken = true;
	  for (int iNeut = 0; iNeut < jNeut; iNeut++) {
	    if (fm_invert_ionForNeutral[iNeut] == k) {
	      notTaken = false;
	    }
	  }
	  if (notTaken) {
	    fm_invert_ionForNeutral[jNeut] = k;
	  }
	}
	fm_neutralMolec_ions_[k  + jNeut * m_kk] += fac;
      }

      // Ok check the work
      for (int m = 0; m < nElementsN; m++) {
	if (fabs(elemVectorN[m]) > 1.0E-13) {
	  throw CanteraError("IonsFromNeutralVPSSTP::initThermoXML", 
			     "Simple formula matrix generation failed");
	}
      }

  
    }
    /*
     * This includes the setStateFromXML calls
     */
    GibbsExcessVPSSTP::initThermoXML(phaseNode, id);

    /*
     * There is one extra step here. We assure ourselves that we
     * have charge conservation.
     */
  }

  // Update the activity coefficients
  /*
   * This function will be called to update the internally storred
   * natural logarithm of the activity coefficients
   *
   *   he = X_A X_B(B + C(X_A - X_B))
   */
  void IonsFromNeutralVPSSTP::s_update_lnActCoeff() const {
    int k, icat, jNeut;
    doublereal fmij;
    /*
     * Get the activity coefficiens of the neutral molecules
     */
    neutralMoleculePhase_->getActivityCoefficients(DATA_PTR(gammaNeutralMolecule_));

    switch (ionSolnType_) {
    case cIonSolnType_PASSTHROUGH:
      break;
    case cIonSolnType_SINGLEANION:
   
      // Do the cation list
      for (k = 0; k < (int) cationList_.size(); k++) {
	//! Get the id for the next cation
        icat = cationList_[k];
	jNeut = fm_invert_ionForNeutral[icat];
	fmij =  fm_neutralMolec_ions_[icat + jNeut * m_kk];
        lnActCoeff_Scaled_[icat] = fmij * log(gammaNeutralMolecule_[jNeut]);
      }

      // Do the anion list
      icat = anionList_[0];
      jNeut = fm_invert_ionForNeutral[icat];
      lnActCoeff_Scaled_[icat]= 0.0;

      // Do the list of neutral molecules
      for (k = 0; k <  numPassThroughSpecies_; k++) {
	icat = passThroughList_[k];
	jNeut = fm_invert_ionForNeutral[icat];
	lnActCoeff_Scaled_[icat] = log(gammaNeutralMolecule_[jNeut]);
      }
      break;
 
    case cIonSolnType_SINGLECATION:
      throw CanteraError("IonsFromNeutralVPSSTP::s_update_lnActCoeff", "Unimplemented type");
      break;
    case cIonSolnType_MULTICATIONANION:
      throw CanteraError("IonsFromNeutralVPSSTP::s_update_lnActCoeff", "Unimplemented type");
      break;
    default:
      throw CanteraError("IonsFromNeutralVPSSTP::s_update_lnActCoeff", "Unimplemented type");
      break;
    }

  }


  // Update the temperatture derivative of the ln activity coefficients
  /*
   * This function will be called to update the internally storred
   * temperature derivative of the natural logarithm of the activity coefficients
   */
  void IonsFromNeutralVPSSTP::s_update_dlnActCoeffdT() const {
    int k, icat, jNeut;
    doublereal fmij;
    /*
     * Get the activity coefficients of the neutral molecules
     */
    GibbsExcessVPSSTP *geThermo = dynamic_cast<GibbsExcessVPSSTP *>(neutralMoleculePhase_);
    if (!geThermo) {
      fvo_zero_dbl_1(dlnActCoeffdT_Scaled_, m_kk);
      return;
    }

    geThermo->getdlnActCoeffdT(DATA_PTR(dlnActCoeffdT_NeutralMolecule_));

    switch (ionSolnType_) {
    case cIonSolnType_PASSTHROUGH:
      break;
    case cIonSolnType_SINGLEANION:
   
      // Do the cation list
      for (k = 0; k < (int) cationList_.size(); k++) {
	//! Get the id for the next cation
        icat = cationList_[k];
	jNeut = fm_invert_ionForNeutral[icat];
	fmij =  fm_neutralMolec_ions_[icat + jNeut * m_kk];
        dlnActCoeffdT_Scaled_[icat] = fmij * dlnActCoeffdT_NeutralMolecule_[jNeut];
      }

      // Do the anion list
      icat = anionList_[0];
      jNeut = fm_invert_ionForNeutral[icat];
      dlnActCoeffdT_Scaled_[icat]= 0.0;

      // Do the list of neutral molecules
      for (k = 0; k <  numPassThroughSpecies_; k++) {
	icat = passThroughList_[k];
	jNeut = fm_invert_ionForNeutral[icat];
	dlnActCoeffdT_Scaled_[icat] = dlnActCoeffdT_NeutralMolecule_[jNeut];
      }
      break;
 
    case cIonSolnType_SINGLECATION:
      throw CanteraError("IonsFromNeutralVPSSTP::s_update_lnActCoeff", "Unimplemented type");
      break;
    case cIonSolnType_MULTICATIONANION:
      throw CanteraError("IonsFromNeutralVPSSTP::s_update_lnActCoeff", "Unimplemented type");
      break;
    default:
      throw CanteraError("IonsFromNeutralVPSSTP::s_update_lnActCoeff", "Unimplemented type");
      break;
    }

  }

  /**
   * Format a summary of the mixture state for output.
   */           
  std::string IonsFromNeutralVPSSTP::report(bool show_thermo) const {


    char p[800];
    string s = "";
    try {
      if (name() != "") {
	sprintf(p, " \n  %s:\n", name().c_str());
	s += p;
      }
      sprintf(p, " \n       temperature    %12.6g  K\n", temperature());
      s += p;
      sprintf(p, "          pressure    %12.6g  Pa\n", pressure());
      s += p;
      sprintf(p, "           density    %12.6g  kg/m^3\n", density());
      s += p;
      sprintf(p, "  mean mol. weight    %12.6g  amu\n", meanMolecularWeight());
      s += p;

      doublereal phi = electricPotential();
      sprintf(p, "         potential    %12.6g  V\n", phi);
      s += p;

      int kk = nSpecies();
      array_fp x(kk);
      array_fp molal(kk);
      array_fp mu(kk);
      array_fp muss(kk);
      array_fp acMolal(kk);
      array_fp actMolal(kk);
      getMoleFractions(&x[0]);
   
      getChemPotentials(&mu[0]);
      getStandardChemPotentials(&muss[0]);
      getActivities(&actMolal[0]);
 

      if (show_thermo) {
        sprintf(p, " \n");
        s += p;
        sprintf(p, "                          1 kg            1 kmol\n");
        s += p;
        sprintf(p, "                       -----------      ------------\n");
        s += p;
        sprintf(p, "          enthalpy    %12.6g     %12.4g     J\n", 
		enthalpy_mass(), enthalpy_mole());
        s += p;
        sprintf(p, "   internal energy    %12.6g     %12.4g     J\n", 
		intEnergy_mass(), intEnergy_mole());
        s += p;
        sprintf(p, "           entropy    %12.6g     %12.4g     J/K\n", 
		entropy_mass(), entropy_mole());
        s += p;
        sprintf(p, "    Gibbs function    %12.6g     %12.4g     J\n", 
		gibbs_mass(), gibbs_mole());
        s += p;
        sprintf(p, " heat capacity c_p    %12.6g     %12.4g     J/K\n", 
		cp_mass(), cp_mole());
        s += p;
        try {
	  sprintf(p, " heat capacity c_v    %12.6g     %12.4g     J/K\n", 
		  cv_mass(), cv_mole());
	  s += p;
        }
        catch(CanteraError) {
	  sprintf(p, " heat capacity c_v    <not implemented>       \n");
	  s += p;
        }
      }
  
    } catch (CanteraError) {
      ;
    }
    return s;
  }

 
}

