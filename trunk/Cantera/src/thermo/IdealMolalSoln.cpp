/**
 *
 *  @file IdealMolalSoln.cpp
 */
/*
 * Copywrite (2006) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
/*
 *  $Author: hkmoffa $
 *  $Date: 2006/08/14 19:18:34 $
 *  $Revision: 1.6 $
 */


#ifndef MAX
#define MAX(x,y)    (( (x) > (y) ) ? (x) : (y))
#endif

#include "IdealMolalSoln.h"
#include "importCTML.h"

namespace Cantera {

  /**
   * Default constructor
   */
  IdealMolalSoln::IdealMolalSoln() :
    MolalityVPSSTP(),
    m_Pcurrent(OneAtm),
    m_formGC(2)
  {
  }

  /**
   * Copy Constructor:
   *
   *  Note this stuff will not work until the underlying phase
   *  has a working copy constructor
   */
  IdealMolalSoln::IdealMolalSoln(const IdealMolalSoln &b) :
    MolalityVPSSTP(b)
  {
    /*
     * Use the assignment operator to do the brunt
     * of the work for the copy construtor.
     */
    *this = b;
  }

  /**
   * operator=()
   *
   *  Note this stuff will not work until the underlying phase
   *  has a working assignment operator
   */
  IdealMolalSoln& IdealMolalSoln::
  operator=(const IdealMolalSoln &b) {
    if (&b != this) {
      MolalityVPSSTP::operator=(b);
      m_speciesMolarVolume  = b.m_speciesMolarVolume;
      m_Pcurrent            = b.m_Pcurrent;
      m_formGC              = b.m_formGC;
      m_expg0_RT            = b.m_expg0_RT;
      m_pe                  = b.m_pe;
      m_pp                  = b.m_pp;
      m_tmpV                = b.m_tmpV;
    }
    return *this;
  }

  IdealMolalSoln::IdealMolalSoln(string inputFile, string id) :
    MolalityVPSSTP()
  {
    constructPhaseFile(inputFile, id);
  }

  IdealMolalSoln::IdealMolalSoln(XML_Node& root, string id) :
    MolalityVPSSTP()
  {
    constructPhaseXML(root, id);
  }

  /**
   *
   * ~IdealMolalSoln():   (virtual)
   *
   * Destructor: does nothing:
   *
   */
  IdealMolalSoln::~IdealMolalSoln() {
  }

  /**
   *
   */
  ThermoPhase* IdealMolalSoln::duplMyselfAsThermoPhase() {
    IdealMolalSoln* mtp = new IdealMolalSoln(*this);
    return (ThermoPhase *) mtp;
  }

  //
  // -------- Molar Thermodynamic Properties of the Solution --------------- 
  //
  /**
   * Molar enthalpy of the solution: Units: J/kmol.
   *
   * Returns the amount of enthalpy per mole of solution.
   * For an ideal molal solution,
   * \f[
   * \bar{h}(T, P, X_k) = \sum_k X_k \bar{h}_k(T)  
   * \f]
   * The formula is written in terms of the partial molar enthalpies.
   * \f$ \bar{h}_k(T, p, m_k) \f$.
   * See the partial molar enthalpy function, getPartialMolarEnthalpies(),
   * for details.
   *
   * Units: J/kmol
   */
  doublereal IdealMolalSoln::enthalpy_mole() const {
    getPartialMolarEnthalpies(DATA_PTR(m_tmpV));
    getMoleFractions(DATA_PTR(m_pp));
    double val = mean_X(DATA_PTR(m_tmpV));
    return val;
  }

  /**
   * Molar internal energy of the solution: Units: J/kmol.
   *
   * Returns the amount of internal energy per mole of solution.
   * For an ideal molal solution,
   * \f[
   * \bar{u}(T, P, X_k) = \sum_k X_k \bar{u}_k(T)  
   * \f]
   * The formula is written in terms of the partial molar internal energy.
   * \f$ \bar{u}_k(T, p, m_k) \f$.
   */
  doublereal IdealMolalSoln::intEnergy_mole() const {
    getPartialMolarEnthalpies(DATA_PTR(m_tmpV));
    return mean_X(DATA_PTR(m_tmpV));
  }

  /**
   * Molar entropy of the solution: Units J/kmol/K.
   *
   * Returns the amount of entropy per mole of solution.
   * For an ideal molal solution,
   * \f[
   * \bar{s}(T, P, X_k) = \sum_k X_k \bar{s}_k(T)  
   * \f]
   * The formula is written in terms of the partial molar entropies.
   * \f$ \bar{s}_k(T, p, m_k) \f$.
   * See the partial molar entropies function, getPartialMolarEntropies(),
   * for details.
   *
   * Units: J/kmol/K.
   */
  doublereal IdealMolalSoln::entropy_mole() const {
    getPartialMolarEntropies(DATA_PTR(m_tmpV));
    return mean_X(DATA_PTR(m_tmpV));
  }

  /**
   * Molar Gibbs function for the solution: Units J/kmol. 
   *
   * Returns the gibbs free energy of the solution per mole
   * of the solution.
   *
   * \f[
   * \bar{g}(T, P, X_k) = \sum_k X_k \mu_k(T)  
   * \f]
   *
   * Units: J/kmol
   */
  doublereal IdealMolalSoln::gibbs_mole() const {
    getChemPotentials(DATA_PTR(m_tmpV));
    return mean_X(DATA_PTR(m_tmpV));
  }

  /**
   * Molar heat capacity at constant pressure: Units: J/kmol/K. 
   *  * \f[
   * \bar{c}_p(T, P, X_k) = \sum_k X_k \bar{c}_{p,k}(T)  
   * \f]
   *
   * Units: J/kmol/K
   */
  doublereal IdealMolalSoln::cp_mole() const {
    getPartialMolarCp(DATA_PTR(m_tmpV));
    double val = mean_X(DATA_PTR(m_tmpV));
    return val;
  }

  /**
   * Molar heat capacity at constant volume: Units: J/kmol/K. 
   * NOT IMPLEMENTED.
   * Units: J/kmol/K
   */
  doublereal IdealMolalSoln::cv_mole() const {
    return err("not implemented");
  }

  //
  // ------- Mechanical Equation of State Properties ------------------------
  //

  /**
   * Pressure. Units: Pa.
   * For this incompressible system, we return the internally storred
   * independent value of the pressure.
   */
  doublereal IdealMolalSoln::pressure() const {
    return m_Pcurrent;
  }

  /**
   * The isothermal compressibility. Units: 1/Pa.
   * The isothermal compressibility is defined as
   * \f[
   * \kappa_T = -\frac{1}{v}\left(\frac{\partial v}{\partial P}\right)_T
   * \f]
   *
   *  It's equal to zero for this model, since the molar volume
   *  doesn't change with pressure or temperature.
   */
  doublereal IdealMolalSoln::isothermalCompressibility() const {
    return 0.0;
  }

  /**
   * The thermal expansion coefficient. Units: 1/K.
   * The thermal expansion coefficient is defined as
   *
   * \f[
   * \beta = \frac{1}{v}\left(\frac{\partial v}{\partial T}\right)_P
   * \f]
   *
   *  It's equal to zero for this model, since the molar volume
   *  doesn't change with pressure or temperature.
   */
  doublereal IdealMolalSoln::thermalExpansionCoeff() const {
    return 0.0;
  }
    
  /**
   * Overwritten setDensity() function is necessary because the
   * density is not an indendent variable.
   *
   * This function will now throw an error condition
   *
   * @internal May have to adjust the strategy here to make
   * the eos for these materials slightly compressible, in order
   * to create a condition where the density is a function of
   * the pressure.
   *
   * This function will now throw an error condition.
   *
   *  NOTE: This is an overwritten function from the State.h
   *        class
   */
  void IdealMolalSoln::setDensity(doublereal rho) {
    double dens = density();
    if (rho != dens) {
      throw CanteraError("Idea;MolalSoln::setDensity",
			 "Density is not an independent variable");
    }
  }

  /**
   * Overwritten setMolarDensity() function is necessary because the
   * density is not an indendent variable.
   *
   * This function will now throw an error condition.
   *
   *  NOTE: This is a virtual function, overwritten function from the State.h
   *        class
   */
  void IdealMolalSoln::setMolarDensity(doublereal conc) {
    double concI = State::molarDensity();
    if (conc != concI) {
      throw CanteraError("IdealMolalSoln::setMolarDensity",
			 "molarDensity/denisty is not an independent variable");
    }
  }

  //
  // ------- Activities and Activity Concentrations
  //

  /**
   * This method returns an array of activity concentrations \f$ C^a_k\f$.
   * \f$ C^a_k\f$ are defined such that 
   * \f$ a_k = C^a_k / C^s_k, \f$ where \f$ C^s_k \f$ 
   * is a standard concentration
   * defined below.  These activity concentrations are used
   * by kinetics manager classes to compute the forward and
   * reverse rates of elementary reactions. 
   *
   * @param c Array of activity concentrations. The 
   *           units depend upon the implementation of the
   *           reaction rate expressions within the phase.
   */
  void IdealMolalSoln::getActivityConcentrations(doublereal* c) const {
    if (m_formGC != 1) {
      double c_solvent = standardConcentration();
      getActivities(c);
      for (int k = 0; k < m_kk; k++) {
	c[k] *= c_solvent;
      }
    } else {
      getActivities(c);
      for (int k = 0; k < m_kk; k++) {
	double c0 = standardConcentration(k);
	c[k] *= c0;
      }
    }
  }

  /**
   * The standard concentration \f$ C^s_k \f$ used to normalize
   * the activity concentration. In many cases, this quantity
   * will be the same for all species in a phase - for example,
   * for an ideal gas \f$ C^s_k = P/\hat R T \f$. For this
   * reason, this method returns a single value, instead of an
   * array.  However, for phases in which the standard
   * concentration is species-specific (e.g. surface species of
   * different sizes), this method may be called with an
   * optional parameter indicating the species.
   *
   */
  doublereal IdealMolalSoln::standardConcentration(int k) const {
    double c0 = 1.0, mvSolvent;
    switch (m_formGC) {
    case 0:
      break;
    case 1:
      c0 = 1.0 /m_speciesMolarVolume[m_indexSolvent];
      break;
    case 2:
      mvSolvent = m_speciesMolarVolume[m_indexSolvent];
      c0 = 1.0 / mvSolvent;
      break;
    }
    return c0;
  }
    
  /**
   * Returns the natural logarithm of the standard 
   * concentration of the kth species
   */
  doublereal IdealMolalSoln::logStandardConc(int k) const {
    double c0 = standardConcentration(k);
    return log(c0);
  }
    
  /**
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
  void IdealMolalSoln::getUnitsStandardConc(double *uA, int k, int sizeUA) {
    int eos = eosType();
    if (eos == 0) {
      for (int i = 0; i < sizeUA; i++) {
        uA[i] = 0.0;
      }
    } else {
      for (int i = 0; i < sizeUA; i++) {
        if (i == 0) uA[0] = 1.0;
        if (i == 1) uA[1] = -nDim();
        if (i == 2) uA[2] = 0.0;
        if (i == 3) uA[3] = 0.0;
        if (i == 4) uA[4] = 0.0;
        if (i == 5) uA[5] = 0.0;
      }
    }
  }

  /**
   * Get the array of non-dimensional molality-based
   * activities at the current solution temperature, 
   * pressure, and solution concentration.
   *
   *  The max against 8.689E-3 is to limit the activity
   *  coefficient to be greater than 1.0E-50.
   */
  void IdealMolalSoln::getActivities(doublereal* ac) const {
    /*
     * Update the molality array, m_molalities()
     *   This requires an update due to mole fractions
     */
    calcMolalities();
    for (int k = 0; k < m_kk; k++) {
      ac[k] = m_molalities[k];
    }
    double xmolSolvent = moleFraction(m_indexSolvent);
    xmolSolvent = MAX(8.689E-3, xmolSolvent);
    ac[m_indexSolvent] = 
      exp((xmolSolvent - 1.0)/xmolSolvent);
  }

  /**
   * Get the array of non-dimensional Molality based
   * activity coefficients at
   * the current solution temperature, pressure, and
   * solution concentration.
   * See Denbigh
   * (note solvent activity coefficient is on the molar scale).
   *
   *  The max against 5.0E-3 (1/200) is to limit the activity
   *  coefficient to be greater than 1.0E-50.
   */
  void IdealMolalSoln::
  getMolalityActivityCoefficients(doublereal* acMolality) const {
    for (int k = 0; k < m_kk; k++) {
      acMolality[k] = 1.0;
    }
    double xmolSolvent = moleFraction(m_indexSolvent);
    xmolSolvent = MAX(8.689E-3, xmolSolvent);
    acMolality[m_indexSolvent] = 
      exp((xmolSolvent - 1.0)/xmolSolvent) / xmolSolvent;
  }

  //
  // ------ Partial Molar Properties of the Solution -----------------
  //

  /**
   * Get the species chemical potentials: Units: J/kmol.
   *
   * This function returns a vector of chemical potentials of the 
   * species in solution.
   *
   * \f[
   *    \mu_k = \mu^{o}_k(T,P) + R T \ln(\frac{m_k}{m^\Delta})
   * \f]
   * \f[
   *    \mu_w = \mu^{o}_w(T,P) +
   *            R T ((X_w - 1.0) / X_w)
   * \f]
   *
   * \f$ w \f$ refers to the solvent species. 
   * \f$ X_w \f$ is the mole fraction of the solvent.
   * \f$ m_k \f$ is the molality of the kth solute.
   * \f$ m^\Delta is 1 gmol solute per kg solvent. \f$
   *
   * Units: J/kmol.
   */
  void IdealMolalSoln::getChemPotentials(doublereal* mu) const{
    double xx;
    const double xxSmall = 1.0E-150; 
    /*
     * First get the standard chemical potentials
     *  -> this requires updates of standard state as a function
     *     of T and P
     * These are defined at unit molality.
     */
    getStandardChemPotentials(mu);
    /*
     * Update the molality array, m_molalities()
     *   This requires an update due to mole fractions
     */
    calcMolalities();
    /*
     *   
     */
    doublereal RT = GasConstant * temperature();
    for (int k = 0; k < m_kk; k++) {
      if (k != m_indexSolvent) {
	xx = MAX(m_molalities[k], xxSmall);
	mu[k] += RT * log(xx);
      }
    }
    /*
     * Do the solvent 
     *  -> see my notes
     */
    double xmolSolvent = moleFraction(m_indexSolvent);
    xx = MAX(xmolSolvent, xxSmall);
    mu[m_indexSolvent] += 
      (RT * (xmolSolvent - 1.0) / xx);
  }

  /**
   * Returns an array of partial molar enthalpies for the species
   * in the mixture: Units (J/kmol).
   *
   * For this phase, the partial molar enthalpies are equal to the
   * SS species enthalpies
   *  \f[
   * \bar{h}_k(T,P) = \hat h^{0}_k(T,P) 
   * \f]
   *     
   *   note hbar = ubar + T d(ubar/dT)
   *   see note about partial molar entropies.
   */
  void IdealMolalSoln::getPartialMolarEnthalpies(doublereal* hbar) const {
    getEnthalpy_RT(hbar);
    doublereal RT = _RT();
    for (int k = 0; k < m_kk; k++) {
      hbar[k] *= RT;
    }
  }

  /**
   * Returns an array of partial molar entropies of the species in the
   * solution: Units: J/kmol.
   *
   * Maxwell's equations provide an insight in how to calculate this
   * (p.215 Smith and Van Ness)
   * \f[
   *      \frac{d(\mu_k)}{dT} = -\bar{s}_i
   * \f]
   * For this phase, the partial molar entropies are equal to the
   * standard state species entropies plus the ideal molal solution contribution.
   * 
   * \f[
   *   \bar{s}_k(T,P) =  s^0_k(T) - R log( m_k )
   * \f]
   * \f[
   *   \bar{s}_w(T,P) =  s^0_w(T) - R ((X_w - 1.0) / X_w)
   * \f]
   *
   * The subscript, w, refers to the solvent species. \f$ X_w \f$ is
   * the mole fraction of solvent.
   * The reference-state pure-species entropies,\f$ s^0_k(T) \f$,
   * at the reference pressure, \f$ P_{ref} \f$, are computed by the
   * species thermodynamic
   * property manager. They are polynomial functions of temperature.
   * @see SpeciesThermo
   */
  void IdealMolalSoln::
  getPartialMolarEntropies(doublereal* sbar) const {
    getEntropy_R(sbar);
    doublereal R = GasConstant;
    doublereal mm;
    calcMolalities();
    for (int k = 0; k < m_kk; k++) {
      if (k != m_indexSolvent) {
	mm = fmaxx(SmallNumber, m_molalities[k]);
	sbar[k] -= R * log(mm);
      }
    }
    double xmolSolvent = moleFraction(m_indexSolvent);
    sbar[m_indexSolvent] -= (R * (xmolSolvent - 1.0) / xmolSolvent);
  }

  /**
   * Returns an array of partial molar volumes of the species
   * in the solution: Units: m^3 kmol-1.
   *
   * For this solution, the partial molar volumes are equal to the
   * constant species molar volumes.
   *
   * Units: m^3 kmol-1.
   */
  void IdealMolalSoln::getPartialMolarVolumes(doublereal* vbar) const {
    getStandardVolumes(vbar);
  }

  /**
   * Partial molar heat capacity of the solution: Units: J/kmol/K.
   *
   *   The kth partial molar heat capacity is equal to 
   *   the temperature derivative of the partial molar
   *   enthalpy of the kth species in the solution at constant
   *   P and composition (p. 220 Smith and Van Ness).
   *  \f[
   *  \bar{Cp}_k(T,P) =  {Cp}^0_k(T) 
   * \f]
   *
   *   For this solution, this is equal to the reference state
   *   heat capacities.
   *
   *  Units: J/kmol/K
   */
  void IdealMolalSoln::getPartialMolarCp(doublereal* cpbar) const {
    /*
     * Get the nondimensional gibbs standard state of the
     * species at the T and P of the solution.
     */
    getCp_R(cpbar);
	
    for (int k = 0; k < m_kk; k++) {
      cpbar[k] *= GasConstant;
    }
  }
	
  /*
   * -------- Properties of the Standard State of the Species
   *           in the Solution ------------------
   */

  /** 
   *  Get the standard state chemical potentials of the species.
   *  This is the array of chemical potentials at unit activity 
   *  (Mole fraction scale)
   *  \f$ \mu^0_k(T,P) \f$.
   *  We define these here as the chemical potentials of the pure
   *  species at the temperature and pressure of the solution.
   *  This function is used in the evaluation of the 
   *  equilibrium constant Kc. Therefore, Kc will also depend
   *  on T and P. This is the norm for liquid and solid systems.
   *
   *  units = J / kmol
   */
  void IdealMolalSoln::getStandardChemPotentials(doublereal* mu) const {
    getGibbs_ref(mu);
    doublereal pref;
    doublereal delta_p;
    for (int k = 0; k < m_kk; k++) {
      pref = m_spthermo->refPressure(k);
      delta_p = m_Pcurrent - pref;
      mu[k] += delta_p * m_speciesMolarVolume[k];
    }
  }
    
  /**
   * Get the nondimensional gibbs function for the species
   * standard states at the current T and P of the solution.
   *
   *  \f[
   *  \mu^0_k(T,P) = \mu^{ref}_k(T) + (P - P_{ref}) * V_k
   * \f]
   * where \f$V_k\f$ is the molar volume of pure species <I>k</I>.
   * \f$ \mu^{ref}_k(T)\f$ is the chemical potential of pure
   * species <I>k</I> at the reference pressure, \f$P_{ref}\f$.
   *
   * @param grt Vector of length m_kk, which on return sr[k]
   *           will contain the nondimensional 
   *           standard state gibbs function for species k. 
   */
  void IdealMolalSoln::getGibbs_RT(doublereal* grt) const {
    getPureGibbs(grt);
    doublereal invRT = 1.0 / _RT();
    for (int k = 0; k < m_kk; k++) {
      grt[k] *= invRT;
    }
  }
    
  /**
   * Get the Gibbs functions for the pure species
   * at the current <I>T</I> and <I>P</I> of the solution.
   * We assume an incompressible constant partial molar
   * volume here:
   * \f[
   *  \mu^0_k(T,p) = \mu^{ref}_k(T) + (P - P_{ref}) * V_k
   * \f]
   * where \f$V_k\f$ is the molar volume of pure species <I>k</I>.
   * \f$ u^{ref}_k(T)\f$ is the chemical potential of pure
   * species <I>k</I> at the reference pressure, \f$P_{ref}\f$.
   *
   * Units: J/kmol
   */
  void IdealMolalSoln::getPureGibbs(doublereal* gpure) const {
    getGibbs_ref(gpure);
    doublereal pref;
    doublereal delta_p;
    for (int k = 0; k < m_kk; k++) {
      pref = m_spthermo->refPressure(k);
      delta_p = m_Pcurrent - pref;
      gpure[k] += delta_p * m_speciesMolarVolume[k];
    }
  }

  /**
   * Get the array of nondimensional Enthalpy functions for the ss
   * species at the current <I>T</I> and <I>P</I> of the solution.
   * We assume an incompressible constant partial molar
   * volume here:
   * \f[
   *  h^0_k(T,P) = h^{ref}_k(T) + (P - P_{ref}) * V_k
   * \f]
   * where \f$V_k\f$ is the molar volume of SS species <I>k</I>.
   * \f$ h^{ref}_k(T)\f$ is the enthalpy of the SS
   * species <I>k</I> at the reference pressure, \f$P_{ref}\f$.
   *
   * Units: dimensionless.
   */
  void IdealMolalSoln::
  getEnthalpy_RT(doublereal* hrt) const {
    getEnthalpy_RT_ref(hrt);
    doublereal pref;
    doublereal delta_p;
    double RT = _RT();
    for (int k = 0; k < m_kk; k++) {
      pref = m_spthermo->refPressure(k);
      delta_p = m_Pcurrent - pref;
      hrt[k] += delta_p/ RT * m_speciesMolarVolume[k];
    }
  }
    
  /**
   * Get the nondimensional Entropies for the species
   * standard states: Units: J/kmol/K
   *
   * Note, this is equal to the reference state entropies
   * due to the zero volume expansivity.
   * i.e., (dS/dp)_T = (dV/dT)_P = 0.0
   *
   * \f[
   *  S^0_k(T,P) = S^{ref}_k(T)
   * \f]
   *
   * Units: dimensionless
   *
   * @param sr Vector of length m_kk, which on return sr[k]
   *           will contain the nondimensional
   *           standard state entropy of species k.
   */
  void IdealMolalSoln::
  getEntropy_R(doublereal* sr) const {
    getEntropy_R_ref(sr);
  }

  /**
   * Get the nondimensional heat capacity at constant pressure
   * function for the species
   * standard states: Units J/kmol/K
   * \f[</I>
   *  Cp^0_k(T,P) = Cp^{ref}_k(T)
   * \f]
   * where \f$V_k\f$ is the molar volume of pure species <I>k</I>.
   * \f$ Cp^{ref}_k(T)\f$ is the constant pressure heat capacity
   * of species <I>k</I> at the reference pressure, \f$p_{ref}\f$.
   *
   * @param cpr Vector of length m_kk, which on return cpr[k]
   *           will contain the nondimensional 
   *           constant pressure heat capacity for species k. 
   */
  void IdealMolalSoln::getCp_R(doublereal* cpr) const {
    getCp_R_ref(cpr); 
  }
    
  /**
   * Get the molar volumes of each species in their standard
   * states at the current
   * <I>T</I> and <I>P</I> of the solution.
   *
   * \f[
   *  V^0_k(T,P) = V^{ref}_k()
   * \f]

   * Units = m^3 / kmol
   */
  void IdealMolalSoln::getStandardVolumes(doublereal *vol) const {
    copy(m_speciesMolarVolume.begin(),
	 m_speciesMolarVolume.end(), vol);
  }
    

  /*
   * ------ Thermodynamic Values for the Species Reference States ---
   */

  // -> This is handled by VPStandardStatesTP

  /*
   *  -------------- Utilities -------------------------------
   */

  /**
   *  Initialization routine for an IdealMolalSoln phase.
   *
   * This is a virtual routine. This routine will call initThermo()
   * for the parent class as well.
   */
  void IdealMolalSoln::initThermo() {
    initLengths();
    MolalityVPSSTP::initThermo();
  }

  /**
   * Initialization of an IdealMolalSoln phase using an
   * xml file
   *
   * This routine is a precursor to constructPhaseFile(XML_Node*)
   * routine, which does most of the work.
   *
   * @param inputFile XML file containing the description of the
   *        phase
   *
   * @param id  Optional parameter identifying the name of the
   *            phase. If none is given, the first XML
   *            phase element will be used.
   */
  void IdealMolalSoln::constructPhaseFile(string inputFile, string id) {

    if (inputFile.size() == 0) {
      throw CanteraError("IdealMolalSoln::constructPhaseFile",
			 "input file is null");
    }
    string path = findInputFile(inputFile);
    ifstream fin(path.c_str());
    if (!fin) {
      throw CanteraError("IdealMolalSoln::constructPhaseFile",
			 "could not open "
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
      throw CanteraError("IdealMolalSoln::constructPhaseFile",
			 "ERROR: Can not find phase named " +
			 id + " in file named " + inputFile);
    }
    fxml_phase->copy(&phaseNode_XML);
    constructPhaseXML(*fxml_phase, id);
    delete fxml;
  }
   	
  /**
   *   Import and initialize an IdealMolalSoln phase 
   *   specification in an XML tree into the current object.
   *   Here we read an XML description of the phase.
   *   We import descriptions of the elements that make up the
   *   species in a phase.
   *   We import information about the species, including their
   *   reference state thermodynamic polynomials. We then freeze
   *   the state of the species.
   *
   *   Then, we read the species molar volumes from the xml 
   *   tree to finish the initialization.
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
  void IdealMolalSoln::constructPhaseXML(XML_Node& phaseNode, string id) {
    if (id.size() > 0) {
      string idp = phaseNode.id();
      if (idp != id) {
	throw CanteraError("IdealMolalSoln::constructPhaseXML", 
			   "phasenode and Id are incompatible");
      }
    }

    /*
     * Find the Thermo XML node 
     */
    if (!phaseNode.hasChild("thermo")) {
      throw CanteraError("IdealMolalSoln::constructPhaseXML",
			 "no thermo XML node");
    }
 
    /*
     * Call the Cantera importPhase() function. This will import
     * all of the species into the phase. Then, it will call
     * initThermoXML() below.
     */
    bool m_ok = importPhase(phaseNode, this);
    if (!m_ok) {
      throw CanteraError("IdealMolalSoln::constructPhaseXML","importPhase failed "); 
    }
  }

  /**
   *   Import and initialize an IdealMolalSoln phase 
   *   specification in an XML tree into the current object.
   *
   *   This routine is called from importPhase() to finish
   *   up the initialization of the thermo object. It reads in the
   *   species molar volumes.
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
  void IdealMolalSoln::initThermoXML(XML_Node& phaseNode, string id) {

    /*
     * Initialize the whole thermo object, using a virtual function.
     */
    initThermo();

    if (id.size() > 0) {
      string idp = phaseNode.id();
      if (idp != id) {
	throw CanteraError("IdealMolalSoln::initThermo", 
			   "phasenode and Id are incompatible");
      }
    }

    /*
     * Find the Thermo XML node 
     */
    if (!phaseNode.hasChild("thermo")) {
      throw CanteraError("IdealMolalSoln::initThermo",
			 "no thermo XML node");
    }
    XML_Node& thermoNode = phaseNode.child("thermo");

    /*
     * Possible change the form of the standard concentrations
     */
    if (thermoNode.hasChild("standardConc")) {
      XML_Node& scNode = thermoNode.child("standardConc");
      m_formGC = 2;
      string formString = scNode.attrib("model");
      if (formString != "") {
	if (formString == "unity") {
	  m_formGC = 0;
	} else if (formString == "molar_volume") {
	  m_formGC = 1;
	} else if (formString == "solvent_volume") {
	  m_formGC = 2;
	} else {
	  throw CanteraError("IdealMolalSoln::initThermo",
			     "Unknown standardConc model: " + formString);
	}
      }
    }

    /*
     * Get the Name of the Solvent:
     *      <solvent> solventName </solvent>
     */
    string solventName = "";
    if (thermoNode.hasChild("solvent")) {
      XML_Node& scNode = thermoNode.child("solvent");
      vector<string> nameSolventa;
      getStringArray(scNode, nameSolventa);
      int nsp = static_cast<int>(nameSolventa.size());
      if (nsp != 1) {
	throw CanteraError("IdealMolalSoln::initThermoXML",
			   "badly formed solvent XML node");
      }
      solventName = nameSolventa[0];
    }


    /*
     * Reconcile the solvent name and index.
     */
    for (int k = 0; k < m_kk; k++) {
      string sname = speciesName(k);
      if (solventName == sname) {
	m_indexSolvent = k;
	break;
      }
    }
    if (m_indexSolvent == -1) {
      cout << "IdealMolalSoln::initThermo: Solvent Name not found" 
	   << endl;
      throw CanteraError("IdealMolalSoln::initThermo",
			 "Solvent name not found");
    }
    if (m_indexSolvent != 0) {
      throw CanteraError("IdealMolalSoln::initThermo",
			 "Solvent " + solventName +
			 " should be first species");
    }

    /*
     * Now go get the molar volumes
     */
    XML_Node& speciesList = phaseNode.child("speciesArray");
    XML_Node* speciesDB =
      get_XML_NameID("speciesData", speciesList["datasrc"],
		     &phaseNode.root());
    const vector<string>&sss = speciesNames();

    for (int k = 0; k < m_kk; k++) {
      XML_Node* s =  speciesDB->findByAttr("name", sss[k]);
      XML_Node *ss = s->findByName("standardState");
      m_speciesMolarVolume[k] = getFloat(*ss, "molarVolume", "-");
    }

    /*
     * Set the state
     */
    if (phaseNode.hasChild("state")) {
      XML_Node& stateNode = phaseNode.child("state");
      setStateFromXML(stateNode);
    }
    
  }

  /**
   * @internal
   * Set equation of state parameters. The number and meaning of
   * these depends on the subclass. 
   * @param n number of parameters
   * @param c array of \i n coefficients
   * 
   */
  void IdealMolalSoln::setParameters(int n, doublereal* c) {
  }
  void IdealMolalSoln::getParameters(int &n, doublereal * const c) {
  }
  /**
   * Set equation of state parameter values from XML
   * entries. This method is called by function importPhase in
   * file importCTML.cpp when processing a phase definition in
   * an input file. It should be overloaded in subclasses to set
   * any parameters that are specific to that particular phase
   * model.
   *
   * @param eosdata An XML_Node object corresponding to
   * the "thermo" entry for this phase in the input file.
   *
   * HKM -> Right now, the parameters are set elsewhere (initThermo)
   *        It just didn't seem to fit.
   */
  void IdealMolalSoln::setParametersFromXML(const XML_Node& eosdata) {
  }

  /*
   * ----------- Critical State Properties --------------------------
   */

  /*
   * ------------ Private and Restricted Functions ------------------
   */

  /**
   * Bail out of functions with an error exit if they are not
   * implemented.
   */
  doublereal IdealMolalSoln::err(string msg) const {
    throw CanteraError("IdealMolalSoln",
		       "Unfinished func called: " + msg );
    return 0.0;
  }

  /**
   * This internal function adjusts the lengths of arrays.
   *
   * This function is not virtual nor is it inherited
   */
  void IdealMolalSoln::initLengths() {
    m_kk = nSpecies();
    /*
     * Obtain the limits of the temperature from the species
     * thermo handler's limits.
     */
    int leng = m_kk;
    m_expg0_RT.resize(leng);
    m_pe.resize(leng, 0.0);
    m_pp.resize(leng);
    m_speciesMolarVolume.resize(leng);
    m_tmpV.resize(leng);
  }
 
}

