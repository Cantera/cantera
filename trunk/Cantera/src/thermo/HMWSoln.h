/**
 *  @file HMWSoln.h
 *
 *    Header file for Pitzer activity coefficient implementation 
 */
/*
 * Copywrite (2006) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
/*
 *  $Id: HMWSoln.h,v 1.4 2006/07/12 15:41:38 hkmoffa Exp $
 */

#ifndef CT_HMWSOLN_H
#define CT_HMWSOLN_H

#include "MolalityVPSSTP.h"
#include "electrolytes.h"

namespace Cantera {

  /**
   * @defgroup thermoprops Thermodynamic Properties
   *
   * These classes are used to compute thermodynamic properties.
   */

  /**
   * HMWSoln.h
   *
   * Major Parameters:
   *   The form of the Pitzer expression refers to the 
   *   form of the Gibbs free energy expression. The temperature
   *   dependence of the Pitzer coefficients are handled by
   *   another parameter.
   *
   * m_formPitzer = Form of the Pitzer expression
   *
   *  PITZERFORM_BASE = 0
   *
   *   Only one form is supported atm. This parameter is included for
   *   future expansion.
   *  
   */
#define PITZERFORM_BASE 0


  /*
   * Formulations for the temperature dependence of the Pitzer 
   * coefficients. Note, the temperature dependence of the
   * Gibbs free energy also depends on the temperature dependence
   * of the standard state and the temperature dependence of the
   * Debye-Huckel constant, which includes the dielectric constant
   * and the density. Therefore, this expression defines only part
   * of the temperature dependence for the mixture thermodynamic
   * functions.
   *
   *  PITZER_TEMP_CONSTANT
   *     All coefficients are considered constant wrt temperature
   *  PITZER_TEMP_LINEAR
   *     All coefficients are assumed to have a linear dependence
   *     wrt to temperature.
   *  PITZER_TEMP_COMPLEX1
   *     All coefficnets are assumed to have a complex functional
   *     based dependence wrt temperature;  See:
   *    (Silvester, Pitzer, J. Phys. Chem. 81, 19 1822 (1977)).
   *
   *       beta0 = q0 + q3(1/T - 1/Tr) + q4(ln(T/Tr)) +
   *               q1(T - Tr) + q2(T**2 - Tr**2)
   */
#define PITZER_TEMP_CONSTANT   0
#define PITZER_TEMP_LINEAR     1
#define PITZER_TEMP_COMPLEX1   2

  /*
   *  Acceptable ways to calculate the value of A_Debye
   */
#define    A_DEBYE_CONST  0
#define    A_DEBYE_WATER  1

  class WaterProps;
  class WaterPDSS;

  /**
   * Definition of the HMWSoln object
   */
  class HMWSoln : public MolalityVPSSTP {

  public:
        
    /// Constructors 
    HMWSoln();

    HMWSoln(const HMWSoln &);
    HMWSoln& operator=(const	HMWSoln&);

    HMWSoln(string inputFile, string id = "");
    HMWSoln(XML_Node& phaseRef, string id = "");

    /**
     *  This is a special constructor, used to replicate test problems
     *  during the initial verification of the object
     */
    HMWSoln(int testProb);

    /// Destructor. 
    virtual ~HMWSoln();


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

    /**
     * Excess molar enthalpy of the solution from 
     * the mixing process. Units: J/ kmol.
     *
     * Note this is kmol of the total solution.
     */
    virtual doublereal relative_enthalpy() const;

    /**
     * Excess molar enthalpy of the solution from 
     * the mixing process on a molality basis.
     *  Units: J/ (kmol add salt).
     *
     * Note this is kmol of the guessed at salt composition
     */
    virtual doublereal relative_molal_enthalpy() const;


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
    /// @name Mechanical Equation of State Properties ---------------------
    //@{
    /**
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
     *  NOTE: This is an overwritten function from the State.h
     *        class
     */
    void setMolarDensity(doublereal rho);

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
        
    virtual doublereal satPressure(doublereal t) const;
        
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

    /*
     * constructPhaseFile()                 (virtual from HMWSoln)
     *
     *   Import, construct, and initialize a HMWSoln phase 
     *   specification from an XML tree into the current object.
     *
     * This routine is a precursor to constructPhaseXML(XML_Node*)
     * routine, which does most of the work.
     */
    virtual void constructPhaseFile(string inputFile, string id);

    /*
     * constructPhaseXML                    (virtual from HMWSoln)
     *
     * This is the main routine for constructing the phase.
     *
     *   Most of the work is carried out by the cantera base
     *   routine, importPhase(). That routine imports all of the
     *   species and element data, including the standard states
     *   of the species.
     *
     *   Then, In this routine, we read the information 
     *   particular to the specification of the activity 
     *   coefficient model for the Pitzer parameterization.
     */
    virtual void constructPhaseXML(XML_Node& phaseNode, string id);

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
     * initThermoXML()                 (virtual from ThermoPhase)
     *
     *  This gets called from importPhase(). It processes the XML file
     *  after the species are set up. This is the main routine for
     *  reading in activity coefficient parameters.
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
    virtual void initThermoXML(XML_Node& phaseNode, string id);

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
     * Value of the Debye Huckel constant as a function of temperature
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

    /**
     * Return Pitzer's definition of A_L. This is basically the
     * derivative of the A_phi multiplied by 4 R T**2
     *
     *            A_Debye = (F e B_Debye) / (8 Pi epsilon R T)
     *            dA_phidT = d(A_Debye)/dT / 3.0
     *            A_L = dA_phidT * (4 * R * T * T)
     *
     *            Units = sqrt(kg/gmol) (RT)
     * 
     */
    double ADebye_L(double temperature = -1.0,
		    double pressure = -1.0) const;


    /**
     * Return Pitzer's definition of A_J. This is basically the
     * temperature derivative of A_L, and the second derivative
     * of A_phi
     *
     *            A_Debye = (F e B_Debye) / (8 Pi epsilon R T)
     *            dA_phidT = d(A_Debye)/dT / 3.0
     *            A_J = 2 A_L/T + 4 * R * T * T * d2(A_phi)/dT2
     *
     *            Units = sqrt(kg/gmol) (R)
     */
    double ADebye_J(double temperature = -1.0,
		    double pressure = -1.0) const;
    /**
     * Return Pitzer's definition of A_V. This is the
     * derivative wrt pressure of A_phi multiplied by - 4 R T
     *
     *            A_Debye = (F e B_Debye) / (8 Pi epsilon R T)
     *            dA_phidT = d(A_Debye)/dP / 3.0
     *            A_V = - dA_phidP * (4 * R * T)
     *
     *            Units = sqrt(kg/gmol) (RT) / Pascal
     * 
     */
    double ADebye_V(double temperature = -1.0,
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

    /*
     * AionicRadius()
     *
     *      Reports the ionic radius of the kth species
     */
    double AionicRadius(int k = 0) const;

    /**
     *
     * formPitzer():
     *
     *  Returns the form of the Pitzer parameterization used
     */
    int formPitzer() const { return m_formPitzer; }

    /**
     * Print out all of the input coefficients.
     */
    void printCoeffs () const;


    //@}
         
  protected:

    /**
     * This is the form of the Pitzer parameterization
     * used in this model.
     * The options are described at the top of this document,
     * and in the general documentation.
     * The list is repeated here:
     *
     * PITZERFORM_BASE  = 0    (only one supported atm)
     *
     */
    int m_formPitzer;

    /**
     * This is the form of the temperature dependence of Pitzer
     * parameterization used in the model.
     *
     *       PITZER_TEMP_CONSTANT   0
     *       PITZER_TEMP_LINEAR     1
     *       PITZER_TEMP_COMPLEX1   2
     */
    int m_formPitzerTemp;

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
     *  Current pressure in Pascal. This is now the independent variable
     *  as it must be for multicomponent solutions.
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
     * Associated Salts, if present in the mechanism,
     * don't contribute to the value of the ionic strength
     * in this version of the Ionic strength.
     */
    mutable double m_IionicMolality;

    /**
     * Maximum value of the ionic strength allowed in the
     * calculation of the activity coefficients.
     */
    double m_maxIionicStrength;

    /**
     * Reference Temperature for the Pitzer formulations.
     */
    double m_TempPitzerRef;

  protected:
    /**
     * Stoichiometric ionic strength on the molality scale.
     * This differs from m_IionicMolality in the sense that
     * associated salts are treated as unassociated salts,
     * when calculating the Ionic strength by this method.
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
     *            It depends on temperature. And, therefore,
     *            most be recalculated whenever T or P changes.
     *            This variable is a local copy of the calculation.
     *            
     *            A_Debye = (F e B_Debye) / (8 Pi epsilon R T)
     *
     *                 where B_Debye = F / sqrt(epsilon R T/2) 
     *                                 (dw/1000)^(1/2)
     *
     *            A_Debye = (1/ (8 Pi)) (2 Pi * Na * dw/1000)^(1/2)
     *                       (e * e / (epsilon * kb * T))^(3/2) 
     *
     *            Units = sqrt(kg/gmol)
     *
     *            Nominal value = 1.172576 sqrt(kg/gmol)
     *                  based on:
     *                    epsilon/epsilon_0 = 78.54
     *                           (water at 25C)
     *                    epsilon_0 = 8.854187817E-12 C2 N-1 m-2
     *                    e = 1.60217653 E-19 C
     *                    F = 9.6485309E7 C kmol-1
     *                    R = 8.314472E3 kg m2 s-2 kmol-1 K-1
     *                    T = 298.15 K
     *                    B_Debye = 3.28640E9 sqrt(kg/gmol)/m
     *                    dw = C_0 * M_0 (density of water) (kg/m3)
     *                       = 1.0E3 at 25C
     */
    mutable double m_A_Debye;

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
     *  Array of 2D data used in the Pitzer/HMW formulation.
     *  Beta0_ij[i][j] is the value of the Beta0 coefficient
     *  for the ij salt. It will be nonzero iff i and j are
     *  both charged and have opposite sign. The array is also
     *  symmetric.
     *     counterIJ where counterIJ = m_counterIJ[i][j]
     *  is used to access this array.
     */
    mutable vector_fp  m_Beta0MX_ij;
    mutable vector_fp  m_Beta0MX_ij_L;
    mutable vector_fp  m_Beta0MX_ij_LL;
    mutable vector_fp  m_Beta0MX_ij_P;
    mutable Array2D    m_Beta0MX_ij_coeff;

    /**
     *  Array of 2D data used in the Pitzer/HMW formulation.
     *  Beta1_ij[i][j] is the value of the Beta1 coefficient
     *  for the ij salt. It will be nonzero iff i and j are
     *  both charged and have opposite sign. The array is also
     *  symmetric. 
     *     counterIJ where counterIJ = m_counterIJ[i][j]
     *  is used to access this array.
     */
    mutable vector_fp m_Beta1MX_ij;
    mutable vector_fp m_Beta1MX_ij_L;
    mutable vector_fp m_Beta1MX_ij_LL;
    mutable vector_fp m_Beta1MX_ij_P;
    mutable Array2D   m_Beta1MX_ij_coeff;

    /**
     *  Array of 2D data used in the Pitzer/HMW formulation.
     *  Beta2_ij[i][j] is the value of the Beta2 coefficient
     *  for the ij salt. It will be nonzero iff i and j are
     *  both charged and have opposite sign, and i and j
     *  both have charges of 2 or more. The array is also
     *  symmetric.
     *  counterIJ where counterIJ = m_counterIJ[i][j]
     *  is used to access this array.
     */
    vector_fp m_Beta2MX_ij;
    vector_fp m_Beta2MX_ij_L;
    vector_fp m_Beta2MX_ij_LL;
    vector_fp m_Beta2MX_ij_P;

    /**
     *  Array of 2D data used in the Pitzer/HMW formulation.
     *  Alpha1MX_ij[i][j] is the value of the alpha1 coefficient
     *  for the ij interaction. It will be nonzero iff i and j are
     *  both charged and have opposite sign, and i and j
     *  both have charges of 2 or more. The array is also
     *  symmetric.
     *  counterIJ where counterIJ = m_counterIJ[i][j]
     *  is used to access this array.
     */
    vector_fp m_Alpha1MX_ij;

    /**
     *  Array of 2D data used in the Pitzer/HMW formulation.
     *  CphiMX_ij[i][j] is the value of the Cphi coefficient
     *  for the ij interaction. It will be nonzero iff i and j are
     *  both charged and have opposite sign, and i and j
     *  both have charges of 2 or more. The array is also
     *  symmetric.
     *  counterIJ where counterIJ = m_counterIJ[i][j]
     *  is used to access this array.
     */
    mutable vector_fp m_CphiMX_ij;
    mutable vector_fp m_CphiMX_ij_L;
    mutable vector_fp m_CphiMX_ij_LL;
    mutable vector_fp m_CphiMX_ij_P;
    mutable Array2D   m_CphiMX_ij_coeff;

    /**
     *  Array of 2D data used in the Pitzer/HMW formulation.
     *  Theta_ij[i][j] is the value of the theta coefficient
     *  for the ij interaction. It will be nonzero for charged
     *  ions with the same sign. It is symmetric.
     *  counterIJ where counterIJ = m_counterIJ[i][j]
     *  is used to access this array.
     *
     *  HKM Recent Pitzer papers have used a functional form
     *      for Theta_ij, which depends on the ionic strength.
     */
    vector_fp m_Theta_ij;
    vector_fp m_Theta_ij_L;
    vector_fp m_Theta_ij_LL;
    vector_fp m_Theta_ij_P;

    /**
     * Array of 3D data sed in the Pitzer/HMW formulation.
     * Psi_ijk[n] is the value of the psi coefficient for the
     * ijk interaction where
     *
     *   n = k + j * m_kk + i * m_kk * m_kk;
     *
     * It is potentially nonzero everywhere. 
     * The first two coordinates are symmetric wrt cations,
     * and the last two coordinates are symmetric wrt anions.
     */
    vector_fp m_Psi_ijk;
    vector_fp m_Psi_ijk_L;
    vector_fp m_Psi_ijk_LL;
    vector_fp m_Psi_ijk_P;

    /*
     * Array of 2D data used in the Pitzer/HMW formulation.
     * Lambda_ij[i][j] represents the lambda coefficient for the
     * ij interaction. This is a general interaction representing
     * neutral species. The neutral species occupy the first
     * index, i.e., i. The charged species occupy the j coordinate.
     * neutral, neutral interactions are also included here.
     */
    Array2D   m_Lambda_ij;
    Array2D   m_Lambda_ij_L;
    Array2D   m_Lambda_ij_LL;
    Array2D   m_Lambda_ij_P;

    /**
     *  Logarithm of the activity coefficients on the molality
     *  scale.
     *       mutable because we change this if the composition
     *       or temperature or pressure changes.
     */
    mutable vector_fp m_lnActCoeffMolal;
    mutable vector_fp m_dlnActCoeffMolaldT;
    mutable vector_fp m_d2lnActCoeffMolaldT2;
    mutable vector_fp m_dlnActCoeffMolaldP;

    /*
     * -------- Temporary Variables Used in the Activity Coeff Calc
     */

    /*
     * Set up a counter variable for keeping track of symmetric binary
     * interactions amongst the solute species.
     *
     * n = m_kk*i + j
     * m_CounterIJ[n] = counterIJ
     */
    mutable array_int m_CounterIJ;

    /**
     *  This is elambda, MEC
     */
    mutable double elambda[17];

    /**
     *  This is elambda1, MEC
     */
    mutable double elambda1[17];

    /**
     *  Various temporary arrays used in the calculation of
     *  the Pitzer activity coefficents.
     *  The subscript, L, denotes the same quantity's derivative
     *  wrt temperature
     */
    mutable vector_fp m_gfunc_IJ;
    mutable vector_fp m_hfunc_IJ;
    mutable vector_fp m_BMX_IJ;
    mutable vector_fp m_BMX_IJ_L;
    mutable vector_fp m_BMX_IJ_LL;
    mutable vector_fp m_BMX_IJ_P;
    mutable vector_fp m_BprimeMX_IJ;
    mutable vector_fp m_BprimeMX_IJ_L;
    mutable vector_fp m_BprimeMX_IJ_LL;
    mutable vector_fp m_BprimeMX_IJ_P;
    mutable vector_fp m_BphiMX_IJ;
    mutable vector_fp m_BphiMX_IJ_L;
    mutable vector_fp m_BphiMX_IJ_LL;
    mutable vector_fp m_BphiMX_IJ_P;
    mutable vector_fp m_Phi_IJ;
    mutable vector_fp m_Phi_IJ_L;
    mutable vector_fp m_Phi_IJ_LL;
    mutable vector_fp m_Phi_IJ_P;
    mutable vector_fp m_Phiprime_IJ;
    mutable vector_fp m_PhiPhi_IJ;
    mutable vector_fp m_PhiPhi_IJ_L;
    mutable vector_fp m_PhiPhi_IJ_LL;
    mutable vector_fp m_PhiPhi_IJ_P;
    mutable vector_fp m_CMX_IJ;
    mutable vector_fp m_CMX_IJ_L;
    mutable vector_fp m_CMX_IJ_LL;
    mutable vector_fp m_CMX_IJ_P;

    mutable vector_fp m_gamma;

  private:
    doublereal err(string msg) const;


    void initLengths();

    /*
     * This function will be called to update the internally storred
     * natural logarithm of the molality activity coefficients 
     */
    void s_update_lnMolalityActCoeff() const;
  public:
    void s_Pitzer_dlnMolalityActCoeff_dT() const;
    void s_Pitzer_dlnMolalityActCoeff_dP() const;
  private:
    /**
     * This function calculates the temperature derivative of the
     * natural logarithm of the molality activity coefficients.
     */
    void s_update_dlnMolalityActCoeff_dT() const;

    /**
     * This function calcultes the temperature second derivative
     * of the natural logarithm of the molality activity 
     * coefficients.
     */
    void s_update_d2lnMolalityActCoeff_dT2() const;
    /**
     * This function calculates the pressure derivative of the
     * natural logarithm of the molality activity coefficients.
     */
    void s_update_dlnMolalityActCoeff_dP() const;

    /**
     * This function calculates the temperature derivatives
     * of the Pitzer coefficients
     */
    void s_updatePitzerCoeffWRTemp(int doDerivs = 2) const;

    /**
     * This function does the main pitzer coefficient 
     * calculation
     */
    void s_updatePitzerSublnMolalityActCoeff() const;
    /*
     * Calculate the lambda interactions. 
     * 
     *
     * Calculate E-lambda terms for charge combinations of like sign,
     *   using method of Pitzer (1975).
     */
    void calc_lambdas(double is) const;

    /**
     *  Calculate etheta and etheta_prime
     *
     *  This interaction will be nonzero for species with the 
     *  same charge. this routine is not to be called for 
     *  neutral species; it core dumps or error exits.
     *
     * MEC implementation routine.
     *
     *  @param z1 charge of the first molecule
     *  @param z2 charge of the second molecule
     *  @param etheta return pointer containing etheta
     *  @param etheta_prime Return pointer containing etheta_prime.
     *
     *  This routine uses the internal variables, 
     *   elambda[] and elambda1[].
     *
     *  There is no prohibition against calling 
     *
     */
    void calc_thetas(int z1, int z2,
		     double *etheta, double *etheta_prime) const;


    void counterIJ_setup(void) const;
    void readXMLBinarySalt(XML_Node &BinSalt);
    void readXMLThetaAnion(XML_Node &BinSalt);
    void readXMLThetaCation(XML_Node &BinSalt);
    void readXMLPsiCommonAnion(XML_Node &BinSalt);
    void readXMLPsiCommonCation(XML_Node &BinSalt);
    void readXMLLambdaNeutral(XML_Node &BinSalt);


  public:
    /*
     * Turn on copious debug printing
     */
    mutable int m_debugCalc;

  };

}
        
#endif

