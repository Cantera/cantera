/**
 *  @file GibbsExcessVPSSTP.cpp
 *   Definitions for intermediate ThermoPhase object for phases which
 *   employ excess gibbs free energy formulations
 *  (see \ref thermoprops 
 * and class \link Cantera::GibbsExcessVPSSTP GibbsExcessVPSSTP\endlink).
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
 *  $Date: 2009/03/03 21:08:31 $
 *  $Revision: 1.3 $
 */


#include "GibbsExcessVPSSTP.h"
using namespace std;

namespace Cantera {

  /*
   * Default constructor.
   *
   */
  GibbsExcessVPSSTP::GibbsExcessVPSSTP() :
    VPStandardStateTP()
  {
  }

  /*
   * Copy Constructor:
   *
   *  Note this stuff will not work until the underlying phase
   *  has a working copy constructor
   */
  GibbsExcessVPSSTP::GibbsExcessVPSSTP(const GibbsExcessVPSSTP &b) :
    VPStandardStateTP()
  {
    *this = operator=(b);
  }

  /*
   * operator=()
   *
   *  Note this stuff will not work until the underlying phase
   *  has a working assignment operator
   */
  GibbsExcessVPSSTP& GibbsExcessVPSSTP::
  operator=(const GibbsExcessVPSSTP &b) {
    if (&b != this) {
      VPStandardStateTP::operator=(b);
    }

    moleFractions_       = b.moleFractions_;
    lnActCoeff_Scaled_   = b.lnActCoeff_Scaled_;
    m_pp                 = b.m_pp;

    return *this;
  }

  /**
   *
   * ~GibbsExcessVPSSTP():   (virtual)
   *
   * Destructor: does nothing:
   *
   */
  GibbsExcessVPSSTP::~GibbsExcessVPSSTP() {
  }

  /*
   * This routine duplicates the current object and returns
   * a pointer to ThermoPhase.
   */
  ThermoPhase* 
  GibbsExcessVPSSTP::duplMyselfAsThermoPhase() const {
    GibbsExcessVPSSTP* mtp = new GibbsExcessVPSSTP(*this);
    return (ThermoPhase *) mtp;
  }

  /*
   *  -------------- Utilities -------------------------------
   */

  void GibbsExcessVPSSTP::setMassFractions(const doublereal* const y) {
#if DEBUG_MODE
    checkMFSum(y);
#endif
    State::setMassFractions(y);
    getMoleFractions(DATA_PTR(moleFractions_));
  }

  void GibbsExcessVPSSTP::setMassFractions_NoNorm(const doublereal* const y) {
#if DEBUG_MODE
    checkMFSum(y);
#endif
    State::setMassFractions_NoNorm(y);
    getMoleFractions(DATA_PTR(moleFractions_));
  }

 void GibbsExcessVPSSTP::setMoleFractions(const doublereal* const x) {
#if DEBUG_MODE
    checkMFSum(x);
#endif
    State::setMoleFractions(x);
    getMoleFractions(DATA_PTR(moleFractions_));
  }

  void GibbsExcessVPSSTP::setMoleFractions_NoNorm(const doublereal* const x) {
#if DEBUG_MODE
    checkMFSum(x);
#endif
    State::setMoleFractions_NoNorm(x);
    getMoleFractions(DATA_PTR(moleFractions_));
  }


  void GibbsExcessVPSSTP::setConcentrations(const doublereal* const c) {
    State::setConcentrations(c);
    getMoleFractions(DATA_PTR(moleFractions_));
  }


  // Equation of state type flag.
  /*
   * The ThermoPhase base class returns
   * zero. Subclasses should define this to return a unique
   * non-zero value. Known constants defined for this purpose are
   * listed in mix_defs.h. The GibbsExcessVPSSTP class also returns
   * zero, as it is a non-complete class.
   */
  int GibbsExcessVPSSTP::eosType() const { 
    return 0;
  }

 

  /*
   * ------------ Molar Thermodynamic Properties ----------------------
   */

  /*
   *
   * ------------ Mechanical Properties ------------------------------
   *
   */

  /*
   * Set the pressure at constant temperature. Units: Pa.
   * This method sets a constant within the object.
   * The mass density is not a function of pressure.
   */
  void GibbsExcessVPSSTP::setPressure(doublereal p) {
#ifdef DEBUG_MODE
    //printf("setPressure: %g\n", p);
#endif
    /*
     * Store the current pressure
     */
    m_Pcurrent = p;
    /*
     * update the standard state thermo
     * -> This involves calling the water function and setting the pressure
     */
    updateStandardStateThermo();
  
    /*
     * Calculate all of the other standard volumes
     * -> note these are constant for now
     */
    calcDensity();
  }

  void GibbsExcessVPSSTP::calcDensity() {
    double *vbar = &m_pp[0];
    getPartialMolarVolumes(vbar);
   
    doublereal vtotal = 0.0;
    for (int i = 0; i < m_kk; i++) {
      vtotal += vbar[i] * moleFractions_[i];
    }
    doublereal dd = meanMolecularWeight() / vtotal;
    State::setDensity(dd);
  }

 


  /*
   * - Activities, Standard States, Activity Concentrations -----------
   */


  doublereal GibbsExcessVPSSTP::standardConcentration(int k) const {
    err("standardConcentration");
    return -1.0;
  }

  doublereal GibbsExcessVPSSTP::logStandardConc(int k) const {
    err("logStandardConc");
    return -1.0;
  }

  void GibbsExcessVPSSTP::getActivities(doublereal* ac) const {
    getActivityCoefficients(ac);
    getMoleFractions(DATA_PTR(moleFractions_));
    for (int k = 0; k < m_kk; k++) {
      ac[k] *= moleFractions_[k];
    }
  }


  void GibbsExcessVPSSTP::getElectrochemPotentials(doublereal* mu) const {
    getChemPotentials(mu);
    double ve = Faraday * electricPotential();
    for (int k = 0; k < m_kk; k++) {
      mu[k] += ve*charge(k);
    }
  }
  
  /*
   * ------------ Partial Molar Properties of the Solution ------------
   */

  // Return an array of partial molar volumes for the
  // species in the mixture. Units: m^3/kmol.
  /*
   *  Frequently, for this class of thermodynamics representations,
   *  the excess Volume due to mixing is zero. Here, we set it as
   *  a default. It may be overriden in derived classes.
   *
   *  @param vbar   Output vector of speciar partial molar volumes.
   *                Length = m_kk. units are m^3/kmol.
   */
  void GibbsExcessVPSSTP::getPartialMolarVolumes(doublereal* vbar) const {
    /*
     * Get the standard state values in m^3 kmol-1
     */
    getStandardVolumes(vbar);
  }
  


  doublereal GibbsExcessVPSSTP::err(std::string msg) const {
    throw CanteraError("GibbsExcessVPSSTP","Base class method "
		       +msg+" called. Equation of state type: "+int2str(eosType()));
    return 0;
  }

    //@}
    /// @name  Properties of the Standard State of the Species in the Solution
    //@{

     

    //@}
    /// @name Thermodynamic Values for the Species Reference States
    //@{

  double GibbsExcessVPSSTP::checkMFSum(const doublereal * const x) const {
    doublereal norm = accumulate(x, x + m_kk, 0.0);
    if (fabs(norm - 1.0) > 1.0E-9) {
      throw CanteraError("GibbsExcessVPSSTP::checkMFSun",
			 "MF sum exceeded tolerance of 1.0E-9:" + fp2str(norm));
    }
    return norm;
  }

  /*
   * Returns the units of the standard and general concentrations
   * Note they have the same units, as their divisor is 
   * defined to be equal to the activity of the kth species
   * in the solution, which is unitless.
   *
   * This routine is used in print out applications where the
   * units are needed. Usually, MKS units are assumed throughout
   * the program and in the XML input files. 
   *
   * On return uA contains the powers of the units (MKS assumed)
   * of the standard concentrations and generalized concentrations
   * for the kth species.
   *
   *  uA[0] = kmol units - default  = 1
   *  uA[1] = m    units - default  = -nDim(), the number of spatial
   *                                dimensions in the Phase class.
   *  uA[2] = kg   units - default  = 0;
   *  uA[3] = Pa(pressure) units - default = 0;
   *  uA[4] = Temperature units - default = 0;
   *  uA[5] = time units - default = 0
   */
  void GibbsExcessVPSSTP::getUnitsStandardConc(double *uA, int k, int sizeUA) const {
    for (int i = 0; i < sizeUA; i++) {
      if (i == 0) uA[0] = 1.0;
      if (i == 1) uA[1] = -nDim();
      if (i == 2) uA[2] = 0.0;
      if (i == 3) uA[3] = 0.0;
      if (i == 4) uA[4] = 0.0;
      if (i == 5) uA[5] = 0.0;
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
  void GibbsExcessVPSSTP::initThermo() {
    initLengths();
    VPStandardStateTP::initThermo();

  }


  //   Initialize lengths of local variables after all species have
  //   been identified.
  void  GibbsExcessVPSSTP::initLengths() {
    m_kk = nSpecies();
    moleFractions_.resize(m_kk);
    lnActCoeff_Scaled_.resize(m_kk);
    m_pp.resize(m_kk);
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
  void GibbsExcessVPSSTP::initThermoXML(XML_Node& phaseNode, std::string id) {


    VPStandardStateTP::initThermoXML(phaseNode, id);
  }
  
 /**
   * Format a summary of the mixture state for output.
   */           
  std::string GibbsExcessVPSSTP::report(bool show_thermo) const {


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

