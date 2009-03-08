/**
 *  @file IdealMolalSoln.h
 *
 * Header file for a derived class of ThermoPhase that handles
 * variable pressure standard state methods for calculating
 * thermodynamic properties that are further based upon 
 * activities on the molality scale. The Ideal molal
 * solution assumes that all molality-based activity
 * coefficients are equal to one.
 */
/*
 * Copywrite (2006) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
/*
 *  $Author: hkmoffa $
 *  $Date: 2006/07/06 15:39:34 $
 *  $Revision: 1.3 $
 */

#ifndef CT_IDEALMOLALSOLN_H
#define CT_IDEALMOLALSOLN_H

#include "MolalityVPSSTP.h"

namespace Cantera {

  /**  \addtogroup thermoprops */
  /* @{ */
 

  /**
   * This phase is based upon the mixing-rule assumption that 
   * all molality-based activity coefficients are equal
   * to one.
   *
   * This is a full instanteation of a ThermoPhase object.
   * The assumption is that the molality-based activity
   * coefficient is equal to one. This also implies that
   * the osmotic coefficient is equal to one.
   *
   * Note, this does not mean that the solution is an
   * ideal solution. In fact, there is a singularity in
   * the formulation as
   * the solvent concentration goes to zero.
   *
   * The mechanical equation of state is currently assumed to
   * be that of an incompressible solution. This may change
   * in the future. Each species has its own molar volume.
   * The molar volume is a constant. 
   *
   * Class IdealMolalSoln represents a condensed phase.
   * The phase and the pure species phases which
   * comprise the standard states of the species are assumed to have
   * zero volume expansivity and zero isothermal compressibility.
   * Each species does, however, have constant but distinct partial
   * molar volumes equal to their pure species molar volumes.
   * The class derives from class ThermoPhase,
   * and overloads the virtual methods defined there with ones that
   * use expressions appropriate for incompressible mixtures.
   *
   * The standard concentrations can have three different forms
   * depending on the value of the member attribute m_formGC, which
   * is supplied in the XML file.
   *
   *                          <TABLE>
   *  <TR><TD> m_formGC </TD><TD> ActivityConc </TD><TD> StandardConc </TD></TR>
   *  <TR><TD> 0        </TD><TD> \f$ {m_k}/ { m^{\Delta}}\f$     </TD><TD> \f$ 1.0        \f$ </TD></TR>
   *  <TR><TD> 1        </TD><TD> \f$  m_k / (m^{\Delta} V_k)\f$  </TD><TD> \f$ 1.0 / V_k  \f$ </TD></TR>
   *  <TR><TD> 2        </TD><TD> \f$  m_k / (m^{\Delta} V^0_0)\f$</TD><TD> \f$ 1.0 / V^0_0\f$ </TD></TR>
   *                         </TABLE>
   *
   * \f$ V^0_0 \f$ is the solvent standard molar volume. \f$ m^{\Delta} \f$ is a constant equal to a 
   * molality of \f$ 1.0 \quad\mbox{gm kmol}^{-1} \f$.
   *
   * The current default is to have mformGC = 2.
   * 
   * The value and form of the activity concentration will affect
   * reaction rate constants involving species in this phase.
   */
  class IdealMolalSoln : public MolalityVPSSTP {

  public:
        
    /// Constructors 
    IdealMolalSoln();
    IdealMolalSoln(const IdealMolalSoln &);
    IdealMolalSoln& operator=(const	IdealMolalSoln&);

    IdealMolalSoln(string inputFile, string id = "");
    IdealMolalSoln(XML_Node& phaseRef, string id = "");

    /// Destructor. 
    virtual ~IdealMolalSoln();


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
    virtual int eosType() const { return 0; }

    /**
     * @} 
     * @name  Molar Thermodynamic Properties of the Solution --------------- 
     * @{
     */

    //Molar enthalpy. Units: J/kmol. 
    virtual doublereal enthalpy_mole() const;

    // Molar internal energy. Units: J/kmol.
    virtual doublereal intEnergy_mole() const;

    // Molar entropy. Units: J/kmol/K. 
    virtual doublereal entropy_mole() const;

    // Molar Gibbs function. Units: J/kmol.
    virtual doublereal gibbs_mole() const;

    // Molar heat capacity at constant pressure. Units: J/kmol/K. 
    virtual doublereal cp_mole() const;

    // Molar heat capacity at constant volume. Units: J/kmol/K. 
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


    /*
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
    virtual void setPressure(doublereal p) {
      m_Pcurrent = p;
    }

    /*
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

    /*
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

    /*
     * Overwritten setMolarDensity() function is necessary because the
     * density is not an indendent variable.
     *
     * This function will now throw an error condition.
     *
     *  NOTE: This is an overwritten function from the State.h
     *        class
     */
    void setMolarDensity(doublereal rho);

    /*
     * The isothermal compressibility. Units: 1/Pa.
     * The isothermal compressibility is defined as
     * \f[
     * \kappa_T = -\frac{1}{v}\left(\frac{\partial v}{\partial P}\right)_T
     * \f]
     */
    virtual doublereal isothermalCompressibility() const;

    /*
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

    /*
     * Set the potential energy of species k to pe.
     * Units: J/kmol.
     * This function must be reimplemented in inherited classes
     * of ThermoPhase.
     */
    virtual void setPotentialEnergy(int k, doublereal pe) {
      err("setPotentialEnergy");
    }

    /*
     * Get the potential energy of species k.
     * Units: J/kmol.
     * This function must be reimplemented in inherited classes
     * of ThermoPhase.
     */
    virtual doublereal potentialEnergy(int k) const {
      return err("potentialEnergy");
    }

    /*
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
     * @name Activities and Activity Concentrations
     *
     * The activity \f$a_k\f$ of a species in solution is
     * related to the chemical potential by \f[ \mu_k = \mu_k^0(T)
     * + \hat R T \log a_k. \f] The quantity \f$\mu_k^0(T)\f$ is
     * the chemical potential at unit activity, which depends only
     * on temperature and the pressure.
     * @{
     */

    /*
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

    /*
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

    /*
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

    /*
     * Get the array of non-dimensional activities at
     * the current solution temperature, pressure, and
     * solution concentration.
     *
     * (note solvent is on molar scale)
     */
    virtual void getActivities(doublereal* ac) const;

    /*
     * Get the array of non-dimensional molality-based
     * activity coefficients at the current solution temperature, 
     * pressure, and solution concentration.
     * 
     * 
     * (note solvent is on molar scale. The solvent molar
     *  based activity coefficient is returned).
     */
    virtual void 
    getMolalityActivityCoefficients(doublereal* acMolality) const;

    //@}
    /// @name  Partial Molar Properties of the Solution -----------------
    //@{

    // Get the species chemical potentials: Units: J/kmol.
    // This is also the partial molar gibbs free energies.
    virtual void getChemPotentials(doublereal* mu) const;

    /**
     * Get the species electrochemical potentials: Units: J/kmol.
     *
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

    /*
     * Returns an array of partial molar enthalpies for the species
     * in the mixture.
     * Units (J/kmol)
     * For this phase, the partial molar enthalpies are equal to the
     * pure species enthalpies
     *  \f[
     * \bar h_k(T,P) = \hat h^{ref}_k(T) + (P - P_{ref}) \hat V^0_k
     * \f]
     * The reference-state pure-species enthalpies, \f$ \hat h^{ref}_k(T) \f$,
     * at the reference pressure,\f$ P_{ref} \f$,
     * are computed by the species thermodynamic 
     * property manager. They are polynomial functions of temperature.
     * @see SpeciesThermo
     */
    virtual void getPartialMolarEnthalpies(doublereal* hbar) const;

    /*
     * getPartialMolarEntropies()        (virtual, const)
     *
     * Returns an array of partial molar entropies of the species in the
     * solution. Units: J/kmol.
     *
     */
    virtual void getPartialMolarEntropies(doublereal* sbar) const;
     
    // partial molar volumes of the species Units: m^3 kmol-1.
    virtual void getPartialMolarVolumes(doublereal* vbar) const;

    /*
     * Partial molar heat capacity of the solution:
     *   The kth partial molar heat capacity is  equal to 
     *   the temperature derivative of the partial molar
     *   enthalpy of the kth species in the solution at constant
     *   P and composition (p. 220 Smith and Van Ness).
     *
     *     Cp = -T d2(chemPot_i)/dT2
     */
    virtual void getPartialMolarCp(doublereal* cpbar) const;

    //@}
    /// @name  Properties of the Standard State of the Species
    //          in the Solution --
    //@{
     
    /*
     *  Get the standard state chemical potentials of the species.
     *  This is the array of chemical potentials at unit activity 
     *  \f$ \mu^0_k(T,P) \f$.
     *  We define these here as the chemical potentials of the pure
     *  species at the temperature and pressure of the solution.
     *  This function is used in the evaluation of the 
     *  equilibrium constant Kc. Therefore, Kc will also depend
     *  on T and P. This is the norm for liquid and solid systems.
     *
     *  units = J / kmol
     */
    virtual void getStandardChemPotentials(doublereal* mu) const;

    /*
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

    /*
     * Get the nondimensional Gibbs functions for the standard
     * state of the species at the current T and P.
     */
    virtual void getPureGibbs(doublereal* gpure) const;

    /*
     *
     * getEnthalpy_RT()        (virtual, const)
     *
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
     */
    virtual void getEnthalpy_RT(doublereal* hrt) const;

    /*
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

    /*
     * Get the nondimensional heat capacity at constant pressure
     * function for the species
     * standard states at the current <I>T</I> and P of the solution.
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

    /*
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
     * 
     * Not implemented.
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
     * @param c array of <I>n</I> coefficients
     * 
     */
    virtual void setParameters(int n, doublereal* c);
    virtual void getParameters(int &n, doublereal * const c);

    /*
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
        
    /**
     *   Critical temperature (K).
     *   Not implemented for this phase type.
     */
    virtual doublereal critTemperature() const {
      err("critTemperature"); return -1.0;
    }
        
    /**
     *  Critical pressure (Pa).
     *
     * Not implemented for this phase type.
     */
    virtual doublereal critPressure() const {
      err("critPressure"); return -1.0;
    }
        
    /**
     *  Critical density (kg/m3).
     *  Not implemented for this phase type.
     */
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

    /*
     * initThermo()                  (virtual from ThermoPhase)
     *
     *  This internal routine is responsible for setting up
     *  the internal storage.
     */
    virtual void initThermo();

    /*
     * constructPhaseFile                   (virtual from here)
     *
     * Initialization of an IdealSolidSolnPhase phase using an
     * xml file identified by its path
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
    virtual void constructPhaseFile(string infile, string id="");

    /*
     * constructPhaseXML                    (virtual from here)
     *
     * This is the main routine for constructing the phase.
     * It processes the XML file, and then it calls importPhase().
     * Then, initThermoXML() is called after importPhase().
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
    virtual void constructPhaseXML(XML_Node& phaseNode, string id);

    /*
     *  initThermoXML                      (virtual from ThermoPhase) 
     * 
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
    virtual void initThermoXML(XML_Node& phaseNode, string id="");

    /*
     * Report the molar volume of species k
     *
     * units - \f$ m^3 kmol^-1 \f$
     */
    double speciesMolarVolume(int k) const;

    /*
     * Fill in a return vector containing the species molar volumes
     * units - \f$ m^3 kmol^-1 \f$
     */
    void   getSpeciesMolarVolumes(double *smv) const;
    //@}

  protected:
    /**
     * Species molar volume \f$ m^3 kmol^-1 \f$
     */
    array_fp   m_speciesMolarVolume;

    /**
     *  Current pressure in Pascal.
     *
     * This is an independent variable in the problem.
     */
    double m_Pcurrent;

    /**
     * The standard concentrations can have three different forms
     * depending on the value of the member attribute m_formGC, which
     * is supplied in the XML file.
     *
     *                          <TABLE>
     *  <TR><TD> m_formGC </TD><TD> ActivityConc </TD><TD> StandardConc </TD></TR>
     *  <TR><TD> 0        </TD><TD> \f$ {m_k}/ { m^{\Delta}}\f$     </TD><TD> \f$ 1.0        \f$ </TD></TR>
     *  <TR><TD> 1        </TD><TD> \f$  m_k / (m^{\Delta} V_k)\f$  </TD><TD> \f$ 1.0 / V_k  \f$ </TD></TR>
     *  <TR><TD> 2        </TD><TD> \f$  m_k / (m^{\Delta} V^0_0)\f$</TD><TD> \f$ 1.0 / V^0_0\f$ </TD></TR>
     *                         </TABLE>
     */
     int m_formGC;

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
  
  private:
    doublereal err(string msg) const;


    void initLengths();
  };

  /* @} */
}
        
#endif





