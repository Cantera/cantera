/**
 *  @file DebyeHuckel.h
 *
 */
/*
 * Copywrite (2006) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
/*
 * $Id: DebyeHuckel.h,v 1.5 2006/07/13 20:05:11 hkmoffa Exp $
 */

#ifndef CT_DEBYEHUCKEL_H
#define CT_DEBYEHUCKEL_H

#include "MolalityVPSSTP.h"
#include "electrolytes.h"
#include "Array.h"

namespace Cantera {

  /**
   * @defgroup thermoprops Thermodynamic Properties
   *
   * These classes are used to compute thermodynamic properties.
   */

  /**
   * DebyeHuckel.h
   *
   * Major Parameters:
   *
   * m_formDH = Form of the Debye-Huckel expression
   *
   *  DHFORM_DILUTE_LIMIT = 0
   *
   *      This form assumes a dilute limit to DH, and is mainly
   *      for informational purposes:
   *
   *   ln(gamma_k)/RT = -z_k**2 * alpha * sqrt(I)
   *
   *              where  I = 1/2 sum_k( molality_k * z_k**2)
   * 
   *  DHFORM_BDOT_AK       = 1
   *
   *      This form assumes Bethke's format for the DH coefficient
   *
   *   ln(gamma_k)/RT = -z_k**2 * alpha * sqrt(I) / (1 + B * a_k * sqrt(I))
   *                        + bdot_k * I
   *      
   *         (note, this particular form where a_k can differ in 
   *          multielectrolyte
   *          solutions has problems wrt a gibbs-duhem analysis. However
   *          we include it here because there is a lot of data fit to it)
   *
   *  DHFORM_BDOT_AUNIFORM = 2
   *
   *      This form assumes Bethke's format for the DH coefficient
   *
   *   ln(gamma_k)/RT = -z_k**2 * alpha * sqrt(I) / (1 + B * a * sqrt(I))
   *                        + bdot_k * I
   *      
   *         The value of a is determined at the beginning of the 
   *         calculation, and not changed.
   *
   *  DHFORM_BETAIJ        = 3
   * 
   *      This form assumes a linear expansion in a virial coefficient form
   *      It is used extensively in Newmann's book, and is the beginning of
   *      more complex treatments for stronger electrolytes, like Pitzer
   *      and HMW treatments.
   *
   *  ln(gamma_k)/RT = -z_k**2 * alpha * sqrt(I) / (1 + B * a * sqrt(I))
   *                        + 2* sum_j (beta_jk m_j)
   *  
   *  DHFORM_PITZER_BETAIJ  = 4
   * 
   *      This form assumes an activity coefficient formulation consistent
   *      with a truncated form of Pitzer's formulation.
   *
   *  ln(gamma_k)/RT = -z_k**2 * alpha * sqrt(I) / (1 + B * a * sqrt(I))
   *       -2 * z_k**2 * alpha * ln(1 + B * a * sqrt(I)) / (B * a)
   *                        + 2 * sum_j (beta_jk m_j)
   *  
   */
#define DHFORM_DILUTE_LIMIT  0
#define DHFORM_BDOT_AK       1
#define DHFORM_BDOT_ACOMMON  2
#define DHFORM_BETAIJ        3
#define DHFORM_PITZER_BETAIJ 4
 
  /*
   *  Acceptable ways to calculate the value of A_Debye
   */
#define    A_DEBYE_CONST  0
#define    A_DEBYE_WATER  1

  class WaterProps;
  class WaterPDSS;

  /**
   * Definition of the DebyeHuckel object
   */
  class DebyeHuckel : public MolalityVPSSTP {

  public:
        
    /// Constructors 
    DebyeHuckel();
    DebyeHuckel(const DebyeHuckel &);
    DebyeHuckel& operator=(const	DebyeHuckel&);

    DebyeHuckel(string inputFile, string id = "");
    DebyeHuckel(XML_Node& phaseRef, string id = "");

    /// Destructor. 
    virtual ~DebyeHuckel();


    ThermoPhase *duplMyselfAsThermoPhase();

    /**
     *   
     * @name  Utilities  
     * @{
     */

    /** 
     * Equation of state type flag. The base class returns
     * zero. Subclasses should define this to return a unique
     * non-zero value. Constants defined for this purpose are
     * listed in mix_defs.h.
     */
    virtual int eosType() const;

    /**
     * @} 
     * @name  Molar Thermodynamic Properties of the Solution --------------
     * @{
     */

    /// Molar enthalpy. Units: J/kmol. 
    /**
     * Molar enthalpy of the solution. Units: J/kmol.
     *      (HKM -> Bump up to Parent object)
     */
    virtual doublereal enthalpy_mole() const;

    /// Molar internal energy. Units: J/kmol. 
    /**
     * Molar internal energy of the solution. Units: J/kmol.
     *      (HKM -> Bump up to Parent object)
     */
    virtual doublereal intEnergy_mole() const;

    /// Molar entropy. Units: J/kmol/K. 
    /**
     * Molar entropy of the solution. Units: J/kmol/K.
     * For an ideal, constant partial molar volume solution mixture with
     * pure species phases which exhibit zero volume expansivity:
     * \f[
     * \hat s(T, P, X_k) = \sum_k X_k \hat s^0_k(T)  
     *      - \hat R  \sum_k X_k log(X_k)
     * \f]
     * The reference-state pure-species entropies 
     * \f$ \hat s^0_k(T,p_{ref}) \f$ are computed by the
     *  species thermodynamic 
     * property manager. The pure species entropies are independent of 
     * temperature since the volume expansivities are equal to zero.
     * @see SpeciesThermo
     *
     *      (HKM -> Bump up to Parent object)
     */
    virtual doublereal entropy_mole() const;

    /// Molar Gibbs function. Units: J/kmol. 
    /*
     *      (HKM -> Bump up to Parent object)
     */
    virtual doublereal gibbs_mole() const;

    /// Molar heat capacity at constant pressure. Units: J/kmol/K. 
    /*
     *     
     */
    virtual doublereal cp_mole() const;

    /// Molar heat capacity at constant volume. Units: J/kmol/K.
    /*
     *      (HKM -> Bump up to Parent object)
     */
    virtual doublereal cv_mole() const;

    //@}
    /** @name Mechanical Equation of State Properties -------------------------
     //@{
     *
     *   In this equation of state implementation, the density is a 
     *   function only of the mole fractions. Therefore, it can't be 
     *   an independent variable. Instead, the pressure is used as the
     *   independent variable. Functions which try to set the thermodynamic
     *   state by calling setDensity() may cause an exception to be
     *   thrown.  
     */

    /**
     * Pressure. Units: Pa.
     * For this incompressible system, we return the internally storred
     * independent value of the pressure.
     */ 
    virtual doublereal pressure() const;

    /**
     * Set the pressure at constant temperature. Units: Pa.
     * This method sets a constant within the object.
     * The mass density is not a function of pressure.
     */
    virtual void setPressure(doublereal p);

    /**
     * Calculate the density of the mixture using the partial 
     * molar volumes and mole fractions as input
     *
     * The formula for this is
     *
     * \f[ 
     * \rho = \frac{\sum_k{X_k W_k}}{\sum_k{X_k V_k}} 
     * \f]
     *
     * where \f$X_k\f$ are the mole fractions, \f$W_k\f$ are
     * the molecular weights, and \f$V_k\f$ are the pure species
     * molar volumes.
     *
     * Note, the basis behind this formula is that in an ideal
     * solution the partial molar volumes are equal to the pure
     * species molar volumes. We have additionally specified
     * in this class that the pure species molar volumes are
     * independent of temperature and pressure.
     *
     * NOTE: This is a non-virtual function, which is not a 
     *       member of the ThermoPhase base class. 
     */
    void calcDensity();

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
    void setDensity(doublereal rho);

    /**
     * Overwritten setMolarDensity() function is necessary because the
     * density is not an indendent variable.
     *
     * This function will now throw an error condition.
     *
     *  NOTE: This is a virtual function overwritten from the State.h
     *        class
     */
    virtual void setMolarDensity(doublereal conc);

   /**
     * Overwritten setTemperature(double) from State.h. This
     * function sets the temperature, and makes sure that
     * the value propagates to underlying objects.
     */
    virtual void setTemperature(doublereal temp);

    /**
     * The isothermal compressibility. Units: 1/Pa.
     * The isothermal compressibility is defined as
     * \f[
     * \kappa_T = -\frac{1}{v}\left(\frac{\partial v}{\partial P}\right)_T
     * \f]
     */
    virtual doublereal isothermalCompressibility() const;

    /**
     * The thermal expansion coefficient. Units: 1/K.
     * The thermal expansion coefficient is defined as
     *
     * \f[
     * \beta = \frac{1}{v}\left(\frac{\partial v}{\partial T}\right)_P
     * \f]
     */
    virtual doublereal thermalExpansionCoeff() const;

    /**
     * @} 
     * @name Potential Energy
     * 
     * Species may have an additional potential energy due to the
     * presence of external gravitation or electric fields. These
     * methods allow specifying a potential energy for individual
     * species.
     * @{
     */

    /**
     * Set the potential energy of species k to pe.
     * Units: J/kmol.
     * This function must be reimplemented in inherited classes
     * of ThermoPhase.
     */
    virtual void setPotentialEnergy(int k, doublereal pe) {
      err("setPotentialEnergy");
    }

    /**
     * Get the potential energy of species k.
     * Units: J/kmol.
     * This function must be reimplemented in inherited classes
     * of ThermoPhase.
     */
    virtual doublereal potentialEnergy(int k) const {
      return err("potentialEnergy");
    }

    /**
     * Set the electric potential of this phase (V).
     * This is used by classes InterfaceKinetics and EdgeKinetics to
     * compute the rates of charge-transfer reactions, and in computing
     * the electrochemical potentials of the species.
     */
    void setElectricPotential(doublereal v) {
      m_phi = v;
    }

    /// The electric potential of this phase (V).
    doublereal electricPotential() const { return m_phi; }


    /**
     * @}
     * @name Activities, Standard States,  and Activity Concentrations
     *
     * The activity \f$a_k\f$ of a species in solution is
     * related to the chemical potential by \f[ \mu_k = \mu_k^0(T)
     * + \hat R T \log a_k. \f] The quantity \f$\mu_k^0(T,P)\f$ is
     * the chemical potential at unit activity, which depends only
     * on temperature and the pressure.
     * Activity is assumed to be molality-based here.
     * @{
     */

    /**
     * This method returns an array of generalized concentrations
     * \f$ C_k\f$ that are defined such that 
     * \f$ a_k = C_k / C^0_k, \f$ where \f$ C^0_k \f$ 
     * is a standard concentration
     * defined below.  These generalized concentrations are used
     * by kinetics manager classes to compute the forward and
     * reverse rates of elementary reactions. 
     *
     * @param c Array of generalized concentrations. The 
     *          units depend upon the implementation of the
     *          reaction rate expressions within the phase.
     */
    virtual void getActivityConcentrations(doublereal* c) const;

    /**
     * The standard concentration \f$ C^0_k \f$ used to normalize
     * the generalized concentration. In many cases, this quantity
     * will be the same for all species in a phase - for example,
     * for an ideal gas \f$ C^0_k = P/\hat R T \f$. For this
     * reason, this method returns a single value, instead of an
     * array.  However, for phases in which the standard
     * concentration is species-specific (e.g. surface species of
     * different sizes), this method may be called with an
     * optional parameter indicating the species.
     */
    virtual doublereal standardConcentration(int k=0) const;

    /**
     * Returns the natural logarithm of the standard 
     * concentration of the kth species
     */
    virtual doublereal logStandardConc(int k=0) const;

    /**
     * Returns the units of the standard and generalized
     * concentrations Note they have the same units, as their
     * ratio is defined to be equal to the activity of the kth
     * species in the solution, which is unitless.
     *
     * This routine is used in print out applications where the
     * units are needed. Usually, MKS units are assumed throughout
     * the program and in the XML input files.
     *
     *  uA[0] = kmol units - default  = 1
     *  uA[1] = m    units - default  = -nDim(), the number of spatial
     *                                dimensions in the Phase class.
     *  uA[2] = kg   units - default  = 0;
     *  uA[3] = Pa(pressure) units - default = 0;
     *  uA[4] = Temperature units - default = 0;
     *  uA[5] = time units - default = 0
     */
    virtual void getUnitsStandardConc(double *uA, int k = 0,
				      int sizeUA = 6);

    /**
     * Get the array of non-dimensional molality-based activities at
     * the current solution temperature, pressure, and
     * solution concentration.
     * (note solvent is on molar scale).
     */
    virtual void getActivities(doublereal* ac) const;

    /**
     * Get the array of non-dimensional molality-based 
     * activity coefficients at
     * the current solution temperature, pressure, and
     * solution concentration.
     * (note solvent is on molar scale. The solvent molar
     *  based activity coefficient is returned).
     */
    virtual void 
    getMolalityActivityCoefficients(doublereal* acMolality) const;

    //@}
    /// @name  Partial Molar Properties of the Solution -----------------
    //@{

    /**
     * Get the species chemical potentials. Units: J/kmol.
     *
     * This function returns a vector of chemical potentials of the 
     * species in solution.
     * \f[
     *    \mu_k = \mu^{ref}_k(T) + V_k * (p - p_o) + R T ln(X_k)
     * \f]
     *  or another way to phrase this is
     * \f[
     *    \mu_k = \mu^o_k(T,p) + R T ln(X_k) 
     * \f]
     *  where \f$ \mu^o_k(T,p) = \mu^{ref}_k(T) + V_k * (p - p_o)\f$
     */
    virtual void getChemPotentials(doublereal* mu) const;


    /**
     * Get the species electrochemical potentials. 
     * These are partial molar quantities.
     * This method adds a term \f$ Fz_k \phi_k \f$ to the 
     * to each chemical potential.
     *
     * Units: J/kmol
     */
    void getElectrochemPotentials(doublereal* mu) const {
      getChemPotentials(mu);
      double ve = Faraday * electricPotential();
      for (int k = 0; k < m_kk; k++) {
	mu[k] += ve*charge(k);
      }
    }

    /**
     * Returns an array of partial molar enthalpies for the species
     * in the mixture.
     * Units (J/kmol)
     * For this phase, the partial molar enthalpies are equal to the
     * pure species enthalpies
     *  \f[
     * \bar h_k(T,P) = \hat h^{ref}_k(T) + (P - P_{ref}) \hat V^0_k
     * \f]
     * The reference-state pure-species enthalpies, 
     * \f$ \hat h^{ref}_k(T) \f$,
     * at the reference pressure,\f$ P_{ref} \f$,
     * are computed by the species thermodynamic 
     * property manager. They are polynomial functions of temperature.
     * @see SpeciesThermo
     */
    virtual void getPartialMolarEnthalpies(doublereal* hbar) const;

    /**
     * getPartialMolarEntropies()        (virtual, const)
     *
     * Returns an array of partial molar entropies of the species in the
     * solution. Units: J/kmol.
     *
     * Maxwell's equations provide an insight in how to calculate this
     * (p.215 Smith and Van Ness)
     *
     *      d(chemPot_i)/dT = -sbar_i
     *      
     *
     * For this phase, the partial molar entropies are equal to the
     * SS species entropies plus the ideal solution contribution.following
     * contribution:
     *  \f[
     * \bar s_k(T,P) =  \hat s^0_k(T) - R log(M0 * molality[k])
     * \f]
     * \f[
     * \bar s_solvent(T,P) =  \hat s^0_solvent(T) 
     *             - R ((xmolSolvent - 1.0) / xmolSolvent)
     * \f]
     *
     * The reference-state pure-species entropies,\f$ \hat s^0_k(T) \f$,
     * at the reference pressure, \f$ P_{ref} \f$,  are computed by the
     * species thermodynamic
     * property manager. They are polynomial functions of temperature.
     * @see SpeciesThermo
     */
    virtual void getPartialMolarEntropies(doublereal* sbar) const;
     
    /**
     * returns an array of partial molar volumes of the species
     * in the solution. Units: m^3 kmol-1.
     *
     * For this solution, thepartial molar volumes are equal to the
     * constant species molar volumes.
     */
    virtual void getPartialMolarVolumes(doublereal* vbar) const;

    virtual void getPartialMolarCp(doublereal* cpbar) const;


    //@}

    /// @name  Properties of the Standard State of the Species
    //          in the Solution --
    //@{

     
    /**
     *  Get the standard state chemical potentials of the species.
     *  This is the array of chemical potentials at unit activity 
     *  \f$ \mu^0_k(T,P) \f$.
     *  Activity is molality based in this object.
     *  We define these here as the chemical potentials of the pure
     *  species at the temperature and pressure of the solution.
     *  This function is used in the evaluation of the 
     *  equilibrium constant Kc. Therefore, Kc will also depend
     *  on T and P. This is the norm for liquid and solid systems.
     *
     *  units = J / kmol
     */
    virtual void getStandardChemPotentials(doublereal* mu) const;

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
    virtual void getGibbs_RT(doublereal* grt) const;

    /**
     * Get the nondimensional Gibbs functions for the standard
     * state of the species at the current T and P.
     */
    virtual void getPureGibbs(doublereal* gpure) const;

    /**
     *
     * getEnthalpy_RT()        (virtual, const)
     *
     * Get the array of nondimensional Enthalpy functions for the 
     * standard states 
     * species at the current <I>T</I> and <I>P</I> of the solution.
     * We assume an incompressible constant partial molar
     * volume here:
     * \f[
     *  h^0_k(T,P) = h^{ref}_k(T) + (P - P_{ref}) * V_k
     * \f]
     * where \f$V_k\f$ is the molar volume of SS species <I>k<\I>.
     * \f$ h^{ref}_k(T)\f$ is the enthalpy of the SS
     * species <I>k<\I> at the reference pressure, \f$P_{ref}\f$.
     */
    virtual void getEnthalpy_RT(doublereal* hrt) const;

    /**
     * Get the nondimensional Entropies for the species
     * standard states at the current T and P of the solution.
     *
     * Note, this is equal to the reference state entropies
     * due to the zero volume expansivity:
     * i.e., (dS/dp)_T = (dV/dT)_P = 0.0
     *
     * @param sr Vector of length m_kk, which on return sr[k]
     *           will contain the nondimensional
     *           standard state entropy of species k.
     */
    virtual void getEntropy_R(doublereal* sr) const;

    /**
     * Get the nondimensional heat capacity at constant pressure
     * function for the species
     * standard states at the current T and P of the solution.
     * \f[
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
    virtual void getCp_R(doublereal* cpr) const;

    /**
     * Get the molar volumes of each species in their standard
     * states at the current
     * <I>T</I> and <I>P</I> of the solution.
     * units = m^3 / kmol
     */
    virtual void getStandardVolumes(doublereal *vol) const;

    //@}
    /// @name Thermodynamic Values for the Species Reference States ---
    //@{


    ///////////////////////////////////////////////////////
    //
    //  The methods below are not virtual, and should not
    //  be overloaded.
    //
    //////////////////////////////////////////////////////

    /**
     * @name Specific Properties
     * @{
     */


    /**
     * @name Setting the State
     *
     * These methods set all or part of the thermodynamic
     * state.
     * @{
     */

    //@}

    /**
     * @name Chemical Equilibrium
     * Chemical equilibrium.
     * @{
     */

    /**
     * This method is used by the ChemEquil equilibrium solver.
     * It sets the state such that the chemical potentials satisfy
     * \f[ \frac{\mu_k}{\hat R T} = \sum_m A_{k,m}
     * \left(\frac{\lambda_m} {\hat R T}\right) \f] where 
     * \f$ \lambda_m \f$ is the element potential of element m. The
     * temperature is unchanged.  Any phase (ideal or not) that
     * implements this method can be equilibrated by ChemEquil.
     */ 
    virtual void setToEquilState(const doublereal* lambda_RT) {
      err("setToEquilState");
    }

    // called by function 'equilibrate' in ChemEquil.h to transfer
    // the element potentials to this object
    void setElementPotentials(const vector_fp& lambda) {
      m_lambda = lambda;
    }

    void getElementPotentials(doublereal* lambda) {
      copy(m_lambda.begin(), m_lambda.end(), lambda);
    }

    //@}


    /**
     * @internal
     * Set equation of state parameters. The number and meaning of
     * these depends on the subclass. 
     * @param n number of parameters
     * @param c array of \i n coefficients
     * 
     */
    virtual void setParameters(int n, doublereal* c);
    virtual void getParameters(int &n, doublereal * const c);

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
     */
    virtual void setParametersFromXML(const XML_Node& eosdata);
 
    //---------------------------------------------------------
    /// @name Critical state properties.
    /// These methods are only implemented by some subclasses.
        
    //@{
        
    /// Critical temperature (K).
    virtual doublereal critTemperature() const {
      err("critTemperature"); return -1.0;
    }
        
    /// Critical pressure (Pa).
    virtual doublereal critPressure() const {
      err("critPressure"); return -1.0;
    }
        
    /// Critical density (kg/m3).
    virtual doublereal critDensity() const {
      err("critDensity"); return -1.0;
    }                
        
    //@}
        
    /// @name Saturation properties.
    /// These methods are only implemented by subclasses that 
    /// implement full liquid-vapor equations of state.
    ///
    virtual doublereal satTemperature(doublereal p) const {
      err("satTemperature"); return -1.0;
    }
        
    virtual doublereal satPressure(doublereal t) const {
      err("satPressure"); return -1.0;
    }
        
    virtual doublereal vaporFraction() const {
      err("vaprFraction"); return -1.0;
    }
        
    virtual void setState_Tsat(doublereal t, doublereal x) {
      err("setState_sat"); 
    }

    virtual void setState_Psat(doublereal p, doublereal x) {
      err("setState_sat"); 
    }

    //@}


    /*
     *  -------------- Utilities -------------------------------
     */

    /**
     * @internal Install a species thermodynamic property
     * manager. The species thermodynamic property manager
     * computes properties of the pure species for use in
     * constructing solution properties. It is meant for internal
     * use, and some classes derived from ThermoPhase may not use
     * any species thermodynamic property manager.
     */
    void setSpeciesThermo(SpeciesThermo* spthermo) 
    { m_spthermo = spthermo; }

    /**
     * Return a reference to the species thermodynamic property
     * manager.  @todo This method will fail if no species thermo
     * manager has been installed.
     */
    SpeciesThermo& speciesThermo() { return *m_spthermo; }


    /**
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
    virtual void initThermo();

    /*
     * Initialization of a DebyeHuckel phase using an
     * xml file
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
    virtual void constructPhaseFile(string infile, string id="");

    /*
     *   Import and initialize a DebyeHuckel phase 
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
    virtual void constructPhaseXML(XML_Node& phaseNode, string id="");


    virtual void  initThermoXML(XML_Node& phaseNode, string id);

    /**
     * Report the molar volume of species k
     *
     * units - \f$ m^3 kmol^-1 \f$
     */
    double speciesMolarVolume(int k) const;

    /**
     * Fill in a return vector containing the species molar volumes
     * units - \f$ m^3 kmol^-1 \f$
     */
    //void   getSpeciesMolarVolumes(double *smv) const;


    /**
     *  Value of the Debye Huckel constant as a function of temperature
     * and pressure.
     *
     *            A_Debye = (F e B_Debye) / (8 Pi epsilon R T)
     *
     *            Units = sqrt(kg/gmol)
     */
    virtual double A_Debye_TP(double temperature = -1.0, 
			      double pressure = -1.0) const;

    /**
     * Value of the derivative of the Debye Huckel constant with 
     * respect to temperature as a function of temperature
     * and pressure.
     *
     *            A_Debye = (F e B_Debye) / (8 Pi epsilon R T)
     *
     *            Units = sqrt(kg/gmol)
     */
    virtual double dA_DebyedT_TP(double temperature = -1.0, 
				 double pressure = -1.0) const;

    /**
     * Value of the 2nd derivative of the Debye Huckel constant with 
     * respect to temperature as a function of temperature
     * and pressure.
     *
     *            A_Debye = (F e B_Debye) / (8 Pi epsilon R T)
     *
     *            Units = sqrt(kg/gmol)
     */
    virtual double d2A_DebyedT2_TP(double temperature = -1.0, 
				   double pressure = -1.0) const;

    /**
     * Value of the derivative of the Debye Huckel constant with 
     * respect to pressure, as a function of temperature
     * and pressure.
     *
     *      A_Debye = (F e B_Debye) / (8 Pi epsilon R T)
     *
     *  Units = sqrt(kg/gmol)
     */
    virtual double dA_DebyedP_TP(double temperature = -1.0, 
				 double pressure = -1.0) const;

    /*
     * AionicRadius()
     *
     *      Reports the ionic radius of the kth species
     */
    double AionicRadius(int k = 0) const;

    /**
     *
     * formDH():
     *
     *  Returns the form of the Debye-Huckel parameterization used
     */
    int formDH() const { return m_formDH; }

    Array2D& get_Beta_ij() { return m_Beta_ij; }

  private:
    /*   Static function that implements the non-polar species
     *   salt-out modifications.
     *   Returns the calculated activity coefficients.
     */
    double _nonpolarActCoeff(double IionicMolality) const;

    /**
     *      Formula for the osmotic coefficient that occurs in
     *      the GWB. It is originally from Helgeson for a variable
     *      NaCl brine. It's to be used with extreme caution.
     */
    double _osmoticCoeffHelgesonFixedForm() const;
    double _lnactivityWaterHelgesonFixedForm() const;

    //@}

         
  protected:

    /**
     * This is the form of the Debye-Huckel parameterization
     * used in this model.
     * The options are described at the top of this document,
     * and in the general documentation.
     * The list is repeated here:
     *
     * DHFORM_DILUTE_LIMIT  = 0       (default)
     * DHFORM_BDOT_AK       = 1
     * DHFORM_BDOT_AUNIFORM = 2
     * DHFORM_BETAIJ        = 3
     * DHFORM_PITZER_BETAIJ = 4
     */
    int m_formDH;

    /**
     * Format for the generalized concentration:
     *
     *  0 = unity
     *  1 = molar_volume
     *  2 = solvent_volume    (default)
     *
     * The generalized concentrations can have three different forms
     * depending on the value of the member attribute m_formGC, which
     * is supplied in the constructor.
     *                          <TABLE>
     *  <TR><TD> m_formGC </TD><TD> GeneralizedConc </TD><TD> StandardConc </TD></TR>
     *  <TR><TD> 0        </TD><TD> X_k             </TD><TD> 1.0          </TD></TR>
     *  <TR><TD> 1        </TD><TD> X_k / V_k       </TD><TD> 1.0 / V_k    </TD></TR>
     *  <TR><TD> 2        </TD><TD> X_k / V_N       </TD><TD> 1.0 / V_N    </TD></TR>
     *                         </TABLE>
     *
     * The value and form of the generalized concentration will affect
     * reaction rate constants involving species in this phase.
     *
     * (HKM Note: Using option #1 may lead to spurious results and
     *  has been included only with warnings. The reason is that it
     *  molar volumes of electrolytes may often be negative. The
     *  molar volume of H+ is defined to be zero too. Either options
     *  0 or 2 are the appropriate choice. Option 0 leads to 
     *  bulk reaction rate constants which have units of s-1.
     *  Option 2 leads to bulk reaction rate constants for 
     *  bimolecular rxns which have units of m-3 kmol-1 s-1.)
     */
    int m_formGC;

    /**
     *  Current pressure in Pascal
     */
    double m_Pcurrent;

	
    vector_int  m_electrolyteSpeciesType;

    /**
     * Species molar volumes \f$ m^3 kmol^-1 \f$
     *  -> m_speciesSize in Constituents.h
     */
    //array_fp m_speciesMolarVolume;

    /**
     *  a_k = Size of the ionic species in the DH formulation
     *        units = meters
     */
    array_fp m_Aionic;

    /**
     * Current value of the ionic strength on the molality scale
     */
    mutable double m_IionicMolality;

    /**
     * Maximum value of the ionic strength allowed in the
     * calculation of the activity coefficients.
     */
    double m_maxIionicStrength;

    /**
     * If true, then the fixed for of Helgeson's activity
     * for water is used instead of the rigoruous form 
     * obtained from Gibbs-Duhem relation. This should be
     * used with caution, and is really only included as a 
     * validation exercise.
     */
  public:
    bool m_useHelgesonFixedForm;
  protected:
    /**
     * Stoichiometric ionic strength on the molality scale
     */
    mutable double m_IionicMolalityStoich;

  public:
    /**
     * Form of the constant outside the Debye-Huckel term
     * called A. It's normally a function of temperature 
     * and pressure. However, it can be set from the
     * input file in order to aid in numerical comparisons.
     * Acceptable forms:
     *
     *       A_DEBYE_CONST  0
     *       A_DEBYE_WATER  1
     *
     * The A_DEBYE_WATER form may be used for water solvents
     * with needs to cover varying temperatures and pressures.
     * Note, the dielectric constant of water is a relatively
     * strong function of T, and its variability must be
     * accounted for,
     */
    int m_form_A_Debye;

  protected:
    /**
     * A_Debye -> this expression appears on the top of the
     *            ln actCoeff term in the general Debye-Huckel
     *            expression
     *            It depends on temperature and pressure.
     *            
     *            A_Debye = (F e B_Debye) / (8 Pi epsilon R T)
     *
     *            Units = sqrt(kg/gmol)
     *
     *            Nominal value(298K, atm) = 1.172576 sqrt(kg/gmol)
     *                  based on:
     *                    epsilon/epsilon_0 = 78.54
     *                           (water at 25C)
     *                    epsilon_0 = 8.854187817E12 C2 N-1 m-2
     *                    e = 8.314472E3 kg m2 s-2 kmol-1 K-1
     *                    F = 9.6485309E7 C kmol-1
     *                    R = 8.314472E3 kg m2 s-2 kmol-1 K-1
     *                    T = 298.15 K
     *                    B_Debye = 3.28640E9 sqrt(kg/gmol)/m
     *
     *            note in Pitzer's nomenclature, A_phi = A_Debye/3.0
     */
    mutable double m_A_Debye;

    /**
     * B_Debye -> this expression appears on the bottom of the
     *            ln actCoeff term in the general Debye-Huckel
     *            expression
     *            It depends on temperature
     *            
     *            B_Bebye = F / sqrt( epsilon R T / 2 )
     *
     *            Units = sqrt(kg/gmol) / m
     *
     *            Nominal value = 3.28640E9 sqrt(kg/gmol) / m
     *                  based on:
     *                    epsilon/epsilon_0 = 78.54
     *                           (water at 25C)
     *                    epsilon_0 = 8.854187817E12 C2 N-1 m-2
     *                    e = 8.314472E3 kg m2 s-2 kmol-1 K-1
     *                    F = 9.6485309E7 C kmol-1
     *                    R = 8.314472E3 kg m2 s-2 kmol-1 K-1
     *                    T = 298.15 K
     */
    double m_B_Debye;

    /**
     *  B_Dot ->  This expression is an extension of the 
     *            Debye-Huckel expression used in some formulations
     *            to extend DH to higher molalities.
     *            B_dot is specific to the major ionic pair.
     */
    array_fp  m_B_Dot;

    /**
     * m_npActCoeff -> These are coefficients to describe
     *  the increase in activity coeff for non-polar molecules
     *  due to the electrolyte becoming stronger (the so-called
     *  salt-out effect)
     */
    array_fp m_npActCoeff;

    /**
     *  Water standard state -> derived from the
     *  equation of state for water.
     */
    WaterPDSS *m_waterSS;
    double m_densWaterSS;

    /**
     *  Pointer to the water property calculator
     */
    WaterProps *m_waterProps;

    /**
     * Vector containing the species reference exp(-G/RT) functions
     * at T = m_tlast
     */
    mutable vector_fp      m_expg0_RT;

    /**
     * Vector of potential energies for the species.
     */
    mutable vector_fp      m_pe;

    /**
     * Temporary array used in equilibrium calculations
     */
    mutable vector_fp      m_pp;

    /**
     * vector of size m_kk, used as a temporary holding area.
     */
    mutable vector_fp      m_tmpV;

    /**
     * Stoichiometric species charge -> This is for calculations
     * of the ionic strength which ignore ion-ion pairing into
     * neutral molecules. The Stoichiometric species charge is the
     * charge of one of the ion that would occur if the species broke
     * into two charged ion pairs.
     *  NaCl ->   m_speciesCharge_Stoich = -1;
     *  HSO4- -> H+ + SO42-              = -2
     *      -> The other charge is calculated.
     * For species that aren't ion pairs, its equal to the
     * m_speciesCharge[] value.
     */
    vector_fp  m_speciesCharge_Stoich;

    /**
     *  Array of 2D data used in the DHFORM_BETAIJ formulation
     *  Beta_ij.value(i,j) is the coefficient of the jth species
     *  for the specification of the chemical potential of the ith
     *  species.
     */
    Array2D m_Beta_ij;

    /**
     *  Logarithm of the activity coefficients on the molality
     *  scale.
     *       mutable because we change this if the composition
     *       or temperature or pressure changes.
     */
    mutable array_fp m_lnActCoeffMolal;
    mutable array_fp m_dlnActCoeffMolaldT;
    mutable array_fp m_d2lnActCoeffMolaldT2;
    mutable array_fp m_dlnActCoeffMolaldP;
  
  private:
    doublereal err(string msg) const;


    void initLengths();

    /*
     * This function will be called to update the internally storred
     * natural logarithm of the molality activity coefficients 
     */
    void s_update_lnMolalityActCoeff() const;

    void s_update_dlnMolalityActCoeff_dT() const;
    void s_update_d2lnMolalityActCoeff_dT2() const;
    void s_update_dlnMolalityActCoeff_dP() const;
  };

}
        
#endif





