/**
 *
 *  @file IdealGasPhase.h
 *   `
 *   ThermoPhase object for the ideal gas equation of state.
 */

/*  $Author$
 *  $Date$
 *  $Revision$
 *
 *  Copyright 2001 California Institute of Technology
 *
 */


#ifndef CT_IDEALGASPHASE_H
#define CT_IDEALGASPHASE_H

//#include "ct_defs.h"
#include "mix_defs.h"
#include "ThermoPhase.h"
#include "SpeciesThermo.h"
#include "utilities.h"

namespace Cantera {

    
  //!Class IdealGasPhase represents low-density gases that obey the
  //!  ideal gas equation of state. 
  /*!
   *
   * %IdealGasPhase derives from class ThermoPhase,
   * and overloads the virtual methods defined there with ones that
   * use expressions appropriate for ideal gas mixtures.
   *
   * This class is optimized for speed of execution.
   *
   * @ingroup thermoprops
   */
  class IdealGasPhase : public ThermoPhase  {

  public:

    //! Empty Constructor
    IdealGasPhase();

    //! Destructor
    virtual ~IdealGasPhase() {}
    
    //! Equation of state flag.
    /*!
     *  Returns the value cIdealGas, defined  in mix_defs.h.
     */
    virtual int eosType() const { return cIdealGas; }

    /**
     * @name Molar Thermodynamic Properties of the Solution ------------------------------
     * @{
     */

    
    //! Return the Molar enthalpy. Units: J/kmol.
    /*!
     * For an ideal gas mixture,
     * \f[
     * \hat h(T) = \sum_k X_k \hat h^0_k(T),
     * \f]
     * and is a function only of temperature.
     * The standard-state pure-species enthalpies 
     * \f$ \hat h^0_k(T) \f$ are computed by the species thermodynamic 
     * property manager.
     *
     * \see SpeciesThermo
     */
    virtual doublereal enthalpy_mole() const {
      return GasConstant * temperature() * 
	mean_X(&enthalpy_RT_ref()[0]);
    }

    /**
     * Molar internal energy. J/kmol. For an ideal gas mixture,
     * \f[
     * \hat u(T) = \sum_k X_k \hat h^0_k(T) - \hat R T,
     * \f]
     * and is a function only of temperature.
     * The reference-state pure-species enthalpies 
     * \f$ \hat h^0_k(T) \f$ are computed by the species thermodynamic 
     * property manager.
     * @see SpeciesThermo
     */
    virtual doublereal intEnergy_mole() const {
      return GasConstant * temperature()
	* ( mean_X(&enthalpy_RT_ref()[0]) - 1.0);
    }

    /**
     * Molar entropy. Units: J/kmol/K.
     * For an ideal gas mixture,
     * \f[
     * \hat s(T, P) = \sum_k X_k \hat s^0_k(T) - \hat R \log (P/P^0).
     * \f]
     * The reference-state pure-species entropies 
     * \f$ \hat s^0_k(T) \f$ are computed by the species thermodynamic 
     * property manager.
     * @see SpeciesThermo
     */
    virtual doublereal entropy_mole() const {
      return GasConstant * (mean_X(&entropy_R_ref()[0]) -
			    sum_xlogx() - std::log(pressure()/m_spthermo->refPressure()));
    }

    /**
     * Molar Gibbs free Energy for an ideal gas.
     * Units =  J/kmol.
     */
    virtual doublereal gibbs_mole() const {
      return enthalpy_mole() - temperature() * entropy_mole();
    }


    /**
     * Molar heat capacity at constant pressure. Units: J/kmol/K.
     * For an ideal gas mixture, 
     * \f[
     * \hat c_p(t) = \sum_k \hat c^0_{p,k}(T).
     * \f]
     * The reference-state pure-species heat capacities  
     * \f$ \hat c^0_{p,k}(T) \f$ are computed by the species thermodynamic 
     * property manager.
     * @see SpeciesThermo
     */
    virtual doublereal cp_mole() const {
      return GasConstant * mean_X(&cp_R_ref()[0]);
    }

    /**
     * Molar heat capacity at constant volume. Units: J/kmol/K.
     * For an ideal gas mixture,
     * \f[ \hat c_v = \hat c_p - \hat R. \f]
     */
    virtual doublereal cv_mole() const {
      return cp_mole() - GasConstant;
    }

    //@}

    /**
     * @name Mechanical Equation of State ------------------------------------------------
     * @{
     */

    /**
     * Pressure. Units: Pa.
     * For an ideal gas mixture, 
     * \f[ P = n \hat R T. \f]
     */ 
    virtual doublereal pressure() const {
      return GasConstant * molarDensity() * temperature();
    }

    
    //! Set the pressure at constant temperature and composition.
    /*!
     *  Units: Pa.
     *   This method is implemented by setting the mass density to
     * \f[
     * \rho = \frac{P \overline W}{\hat R T }.
     * \f]
     *
     * @param p Pressure (Pa) 
     */
    virtual void setPressure(doublereal p) {
      setDensity(p * meanMolecularWeight()
		 /(GasConstant * temperature()));
    }

    //! Returns  the isothermal compressibility. Units: 1/Pa.
    /**
     * The isothermal compressibility is defined as
     * \f[
     * \kappa_T = -\frac{1}{v}\left(\frac{\partial v}{\partial P}\right)_T
     * \f]
     *  For ideal gases it's equal to the negative of the inverse of the pressure
     */
    virtual doublereal isothermalCompressibility() const {
      return -1.0/pressure();
    }

    //! Return the volumetric thermal expansion coefficient. Units: 1/K.
    /*!
     * The thermal expansion coefficient is defined as
     * \f[
     * \beta = \frac{1}{v}\left(\frac{\partial v}{\partial T}\right)_P
     * \f]
     * For ideal gases, it's equal to the inverse of the temperature.
     */
    virtual doublereal thermalExpansionCoeff() const {
      return 1.0/temperature();
    }

    //@}

    /**
     * @name Chemical Potentials and Activities ------------------------------------------
     *     
     *
     * The activity \f$a_k\f$ of a species in solution is
     * related to the chemical potential by 
     * \f[
     *  \mu_k(T,P,X_k) = \mu_k^0(T,P)
     * + \hat R T \log a_k.
     *  \f] 
     * The quantity \f$\mu_k^0(T,P)\f$ is
     * the standard state chemical potential at unit activity.
     * It may depend on the pressure and the temperature. However,
     * it may not depend on the mole fractions of the species 
     * in the solution.
     *
     * The activities are related to the generalized 
     * concentrations, \f$\tilde C_k\f$, and standard 
     * concentrations, \f$C^0_k\f$, by the following formula:
     *
     *  \f[
     *  a_k = \frac{\tilde C_k}{C^0_k} 
     *  \f] 
     * The generalized concentrations are used in the kinetics classes
     * to describe the rates of progress of reactions involving the
     * species. Their formulation depends upons the specification
     * of the rate constants for reaction, especially the units used
     * in specifying the rate constants. The bridge between the
     * thermodynamic equilibrium expressions that use a_k and the
     * kinetics expressions which use the generalized concentrations
     * is provided by the multiplicative factor of the 
     * standard concentrations. 
     * @{
     */

    //! This method returns the array of generalized concentrations. 
    /*!
     *  For an ideal gas mixture, these are simply the actual concentrations.
     *
     * @param c Output array of generalized concentrations. The 
     *           units depend upon the implementation of the
     *           reaction rate expressions within the phase.
     */
    virtual void getActivityConcentrations(doublereal* c) const {
      getConcentrations(c);
    }
    
    //! Returns the standard concentration \f$ C^0_k \f$, which is used to normalize
    //! the generalized concentration.
    /*!
     * This is defined as the concentration  by which the generalized
     * concentration is normalized to produce the activity.
     * In many cases, this quantity will be the same for all species in a phase.
     * Since the activity for an ideal gas mixture is 
     * simply the mole fraction, for an ideal gas \f$ C^0_k = P/\hat R T \f$.
     *
     * @param k Optional parameter indicating the species. The default
     *          is to assume this refers to species 0.
     * @return 
     *   Returns the standard Concentration in units of m3 kmol-1.
     */
    virtual doublereal standardConcentration(int k=0) const;

    //! Returns the natural logarithm of the standard 
    //! concentration of the kth species
    /*!
     * @param k    index of the species. (defaults to zero)
     */
    virtual doublereal logStandardConc(int k=0) const;

    //! Get the array of non-dimensional activity coefficients at
    //! the current solution temperature, pressure, and solution concentration. 
    /*!
     *  For ideal gases, the activity coefficients are all equal to one.
     *
     * @param ac Output vector of activity coefficients. Length: m_kk.
     */
    virtual void getActivityCoefficients(doublereal* ac) const;


    //@}
    /// @name Partial Molar Properties of the Solution ----------------------------------
    //@{

    
    //! Get the species chemical potentials. Units: J/kmol.
    /*!
     * This function returns a vector of chemical potentials of the 
     * species in solution at the current temperature, pressure
     * and mole fraction of the solution.
     *
     * @param mu  Output vector of species chemical 
     *            potentials. Length: m_kk. Units: J/kmol
     */
    virtual void getChemPotentials(doublereal* mu) const;
 
    //!  Get the species partial molar enthalpies. Units: J/kmol.
    /*!
     * @param hbar    Output vector of species partial molar enthalpies.
     *                Length: m_kk. units are J/kmol.
     */
    virtual void getPartialMolarEnthalpies(doublereal* hbar) const;

    //! Get the species partial molar entropies. Units: J/kmol/K.
    /*!
     * @param sbar    Output vector of species partial molar entropies.
     *                Length = m_kk. units are J/kmol/K.
     */
    virtual void getPartialMolarEntropies(doublereal* sbar) const;

    //! Get the species partial molar enthalpies. Units: J/kmol.
    /*!
     * @param ubar    Output vector of speciar partial molar internal energies.
     *                Length = m_kk. units are J/kmol.
     */
    virtual void getPartialMolarIntEnergies(doublereal* ubar) const;

    //! Get the partial molar heat capacities Units: J/kmol/K
    /*!
     * @param cpbar   Output vector of species partial molar heat capacities at constant pressure.
     *                Length = m_kk. units are J/kmol/K.
     */
    virtual void getPartialMolarCp(doublereal* cpbar) const;

    //! Get the species partial molar volumes. Units: m^3/kmol.
    /*!
     *  @param vbar   Output vector of speciar partial molar volumes.
     *                Length = m_kk. units are m^3/kmol.
     */
    virtual void getPartialMolarVolumes(doublereal* vbar) const;

    //@}
    /// @name  Properties of the Standard State of the Species in the Solution ----------
    //@{

    //! Get the array of chemical potentials at unit activity for the
    //! species standard states at the current <I>T</I> and <I>P</I> of the solution.
    /*!
     * These are the standard state chemical potentials \f$ \mu^0_k(T,P)
     * \f$. The values are evaluated at the current
     * temperature and pressure of the solution
     *
     * @param mu  Output vector of chemical potentials. 
     *            Length: m_kk.
     */
    virtual void getStandardChemPotentials(doublereal* mu) const;

    //! Get the nondimensional Enthalpy functions for the species standard states
    //! at their standard states at the current <I>T</I> and <I>P</I> of the solution.
    /*!
     * @param hrt      Output vector of  nondimensional standard state enthalpies.
     *                 Length: m_kk.
     */
    virtual void getEnthalpy_RT(doublereal* hrt) const;

    //! Get the array of nondimensional Entropy functions for the
    //! species standard states at the current <I>T</I> and <I>P</I> of the solution.
    /*!
     * @param sr   Output vector of  nondimensional standard state entropies.
     *             Length: m_kk.
     */
    virtual void getEntropy_R(doublereal* sr) const;

    //! Get the nondimensional Gibbs functions for the species
    //! standard states at the current <I>T</I> and <I>P</I> of the solution.
    /*!
     * @param grt  Output vector of nondimensional standard state gibbs free energies
     *             Length: m_kk.
     */
    virtual void getGibbs_RT(doublereal* grt) const;

    //! Get the Gibbs functions for the standard
    //! state of the species at the current <I>T</I> and <I>P</I> of the solution
    /*!
     * Units are Joules/kmol
     * @param gpure  Output vector of  standard state gibbs free energies
     *               Length: m_kk.
     */
    virtual void getPureGibbs(doublereal* gpure) const;

    //!  Returns the vector of nondimensional Internal Energies  of the standard
    //!  state species at the current <I>T</I> and <I>P</I> of the solution
    /*!
     * @param urt  output vector of nondimensional standard state internal energies
     *             of the species. Length: m_kk. 
     */
    virtual void getIntEnergy_RT(doublereal *urt) const;

    //! Get the nondimensional Heat Capacities at constant
    //! pressure for the species standard states
    //! at the current <I>T</I> and <I>P</I> of the solution
    /*!
     * @param cpr   Output vector of nondimensional standard state heat capacities
     *              Length: m_kk.
     */
    virtual void getCp_R(doublereal* cpr) const;

    //!  Get the molar volumes of the species standard states at the current
    //!  <I>T</I> and <I>P</I> of the solution.
    /*!
     * units = m^3 / kmol
     *
     * @param vol     Output vector containing the standard state volumes.
     *                Length: m_kk.
     */
    virtual void getStandardVolumes(doublereal *vol) const;

    //@}
    /// @name Thermodynamic Values for the Species Reference States ---------------------
    //@{


    //!  Returns the vector of nondimensional
    //!  enthalpies of the reference state at the current temperature
    //!  of the solution and the reference pressure for the species.
    /*!
     * @param hrt     Output vector containing the nondimensional reference state 
     *                enthalpies.  Length: m_kk.
     */
    virtual void getEnthalpy_RT_ref(doublereal *hrt) const;

    //!  Returns the vector of nondimensional
    //!  Gibbs Free Energies of the reference state at the current temperature
    //!  of the solution and the reference pressure for the species.
    /*!
     * @param grt     Output vector containing the nondimensional reference state 
     *                Gibbs Free energies.  Length: m_kk.
     */
    virtual void getGibbs_RT_ref(doublereal *grt) const;

    //!  Returns the vector of the
    //!  gibbs function of the reference state at the current temperature
    //!  of the solution and the reference pressure for the species.
    /*!
     *  units = J/kmol
     *
     * @param g       Output vector containing the  reference state 
     *                Gibbs Free energies.  Length: m_kk. Units: J/kmol.
     */
    virtual void getGibbs_ref(doublereal *g) const;

    //!  Returns the vector of nondimensional
    //!  entropies of the reference state at the current temperature
    //!  of the solution and the reference pressure for each species.
    /*!
     * @param er      Output vector containing the nondimensional reference state 
     *                entropies.  Length: m_kk.
     */
    virtual void getEntropy_R_ref(doublereal *er) const;

    //! Returns the vector of nondimensional
    //!  internal Energies of the reference state at the current temperature
    //!  of the solution and the reference pressure for each species.
    /*!
     * @param urt    Output vector of nondimensional reference state
     *               internal energies of the species.
     *               Length: m_kk
     */
    virtual void getIntEnergy_RT_ref(doublereal *urt) const;
    
    //!  Returns the vector of nondimensional
    //!  constant pressure heat capacities of the reference state
    //!  at the current temperature of the solution
    //!  and reference pressure for each species.
    /*!
     * @param cprt   Output vector of nondimensional reference state
     *               heat capacities at constant pressure for the species.
     *               Length: m_kk
     */
    virtual void  getCp_R_ref(doublereal *cprt) const;

    //!  Get the molar volumes of the species standard states at the current
    //!  <I>T</I> and <I>P_ref</I> of the solution.
    /*!
     * units = m^3 / kmol
     *
     * @param vol     Output vector containing the standard state volumes.
     *                Length: m_kk.
     */
    virtual void getStandardVolumes_ref(doublereal *vol) const;

    //@}
    /// @name New Methods Defined Here  -------------------------------------------------
    //@{

    const array_fp& enthalpy_RT_ref() const {
      _updateThermo();
      return m_h0_RT;
    }

    const array_fp& gibbs_RT_ref() const {
      _updateThermo();
      return m_g0_RT;
    }

    const array_fp& expGibbs_RT_ref() const {
      _updateThermo();
      int k;
      for (k = 0; k != m_kk; k++) m_expg0_RT[k] = std::exp(m_g0_RT[k]);
      return m_expg0_RT;
    }

    const array_fp& entropy_R_ref() const {
      _updateThermo();
      return m_s0_R;
    }

    const array_fp& cp_R_ref() const {
      _updateThermo();
      return m_cp0_R;
    }

    // @}

    /**
     * @internal Initialize. This method is provided to allow
     * subclasses to perform any initialization required after all
     * species have been added. For example, it might be used to
     * resize internal work arrays that must have an entry for
     * each species.  The base class implementation does nothing,
     * and subclasses that do not require initialization do not
     * need to overload this method. 
     *
     * @see importCTML.cpp
     */
    virtual void initThermo();

    //!This method is used by the ChemEquil equilibrium solver.
    /*!
     * @internal
     * @name Chemical Equilibrium
     * @{
     * 
     * Set mixture to an equilibrium state consistent with specified 
     * element potentials and temperature.
     * It sets the state such that the chemical potentials satisfy
     * \f[ \frac{\mu_k}{\hat R T} = \sum_m A_{k,m}
     * \left(\frac{\lambda_m} {\hat R T}\right) \f] where 
     * \f$ \lambda_m \f$ is the element potential of element m. The
     * temperature is unchanged.  Any phase (ideal or not) that
     * implements this method can be equilibrated by ChemEquil.
     *
     * @param lambda_RT vector of non-dimensional element potentials
     *                  \f[ \lambda_m/RT \f].
     */
    virtual void setToEquilState(const doublereal* lambda_RT);

    //@}

  protected:

    //! Number of Elements in the phase
    /*!
     * This member is defined here, from a call to the Elements ojbect, for speed.
     */
    int m_mm;

    //! Minimum temperature for valid species standard state thermo props
    /*!
     * This is the minimum temperature at which all species have valid standard
     * state thermo props defined.
     */
    doublereal m_tmin;

    //! Maximum temperature for valid species standard state thermo props
    /*!
     * This is the maximum temperature at which all species have valid standard
     * state thermo props defined.
     */
    doublereal m_tmax;

    //! Reference state pressure
    /*!
     *  Value of the reference state pressure in Pascals. 
     *  All species must have the same reference state pressure.
     */
    doublereal m_p0;

    //! last value of the temperature processed by reference state 
    mutable doublereal    m_tlast;

    //! Temporary storage for log of p/rt
    mutable doublereal    m_logc0;

   //! Temporary storage for dimensionless reference state enthalpies
    mutable array_fp      m_h0_RT;

   //! Temporary storage for dimensionless reference state heat capacities
    mutable array_fp      m_cp0_R;

   //! Temporary storage for dimensionless reference state gibbs energies
    mutable array_fp      m_g0_RT;

    //! Temporary storage for dimensionless reference state entropies
    mutable array_fp      m_s0_R;

    //! currently unsed
    /*!
     * @deprecated
     */
    mutable array_fp      m_expg0_RT;

    //! Currently unused
    /*
     * @deprecated
     */
    mutable array_fp      m_pe;

    //! Temporary array containing internally calculated partial pressures
    mutable array_fp      m_pp;

  private:

    void _updateThermo() const;
  };
}
        
#endif
