/**
 *
 *  @file IdealGasPhase.h
 *   
 *   ThermoPhase object for the ideal gas equation of state.
 */

/*  $Author: dggoodwin $
 *  $Date: 2006/04/28 17:22:23 $
 *  $Revision: 1.8 $
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

    /**
     * Class IdealGasPhase represents low-density gases that obey the
     * ideal gas equation of state. 
     *
     * IdealGasPhase derives from class ThermoPhase,
     * and overloads the virtual methods defined there with ones that
     * use expressions appropriate for ideal gas mixtures.
     * @ingroup thermoprops
     */
    class IdealGasPhase : public ThermoPhase  {

    public:

        IdealGasPhase(): m_tlast(0.0) {}

        virtual ~IdealGasPhase() {}

        /**
         * Equation of state flag. Returns the value cIdealGas, defined 
         * in mix_defs.h.
         */
        virtual int eosType() const { return cIdealGas; }


        /**
         * @name Molar Thermodynamic Properties of the Solution ------------------------------
         * @{
         */

        /**
         * Molar enthalpy. Units: J/kmol.
         * For an ideal gas mixture,
         * \f[
         * \hat h(T) = \sum_k X_k \hat h^0_k(T),
         * \f]
         * and is a function only of temperature.
         * The standard-state pure-species enthalpies 
         * \f$ \hat h^0_k(T) \f$ are computed by the species thermodynamic 
         * property manager.
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
                sum_xlogx() - log(pressure()/m_spthermo->refPressure()));
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

        /**
         * Set the pressure at constant temperature. Units: Pa.
         * This method is implemented by setting the mass density to
         * \f[
         * \rho = \frac{P \overline W}{\hat R T }.
         * \f] 
         */
        virtual void setPressure(doublereal p) {
            setDensity(p * meanMolecularWeight()
                /(GasConstant * temperature()));
        }

        virtual doublereal isothermalCompressibility() const {
            return -1.0/pressure();
        }

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

	/**
         * This method returns the array of generalized
         * concentrations.  For an ideal gas mixture, these are simply
         * the actual concentrations.
         */
        virtual void getActivityConcentrations(doublereal* c) const {
            getConcentrations(c);
        }

	/**
         * The standard concentration. This is defined as the concentration 
         * by which the generalized concentration is normalized to produce 
         * the activity. Since the activity for an ideal gas mixture is 
         * simply the mole fraction, the standard concentration is 
         * \f$ P / R T \f$.
         */ 
         virtual doublereal standardConcentration(int k=0) const {
	     double p = pressure();
	     return p/(GasConstant * temperature());
        }

	/**
	 * Returns the natural logarithm of the standard 
	 * concentration of the kth species
	 */
         virtual doublereal logStandardConc(int k=0) const {
             _updateThermo();
	     double p = pressure();
	     double lc = log (p / (GasConstant * temperature()));
	     return lc;
        }

	/** 
         * Get the array of non-dimensional activity coefficients at
	 * the current solution temperature, pressure, and
	 * solution concentration. 
         *  For ideal gases, the activity coefficients are all equal
         *  to one.
         */
        virtual void getActivityCoefficients(doublereal* ac) const;

        /**
	 *  Get the array of chemical potentials at unit activity \f$
         * \mu^0_k \f$ at the current temperature and pressure of the
	 *  solution.
         *  These are the standard state chemical potentials.
         */
        virtual void getStandardChemPotentials(doublereal* muStar) const;

  
       //@}
        /// @name Partial Molar Properties of the Solution ----------------------------------
        //@{

        /**
         * Get the species chemical potentials. Units: J/kmol.
	 *
	 * This function returns a vector of chemical potentials of the 
	 * species in solution at the current temperature, pressure
	 * and mole fraction of the solution.
         */
        virtual void getChemPotentials(doublereal* mu) const;
 
	/**
	 * Get the array of partial molar enthalpies
	 * units = J / kmol
	 */
        virtual void getPartialMolarEnthalpies(doublereal* hbar) const;

        /**
         * Returns an array of partial molar entropies of the species in the
	 * solution. Units: J/kmol.
         */
        virtual void getPartialMolarEntropies(doublereal* sbar) const;

	/**
	 * Get the array of partial molar volumes
	 * units = m^3 / kmol
	 */
	virtual void getPartialMolarVolumes(doublereal* vbar) const;

        //@}
        /// @name  Properties of the Standard State of the Species in the Solution ----------
        //@{

        /**
         * Get the nondimensional Enthalpy functions for the species
         * at their standard states at the current
	 *  <I>T</I> and <I>P</I> of the solution.
         */
        virtual void getEnthalpy_RT(doublereal* hrt) const;

        /**
         * Get the array of nondimensional Enthalpy functions for the 
	 * standard state species
         * at the current <I>T</I> and <I>P</I> of the solution.
         */
        virtual void getEntropy_R(doublereal* sr) const;

        /**
         * Get the nondimensional gibbs function for the species
         * standard states at the current T and P of the solution.
         */
        virtual void getGibbs_RT(doublereal* grt) const;

	/**
	 * Get the Gibbs functions for the pure species
         * at the current <I>T</I> and <I>P</I> of the solution.
	 */
	virtual void getPureGibbs(doublereal* gpure) const;

	/**
	 *  Returns the vector of nondimensional
         *  internal Energies of the standard state at the current temperature
         *  and pressure of the solution for each species.
         */
        virtual void getIntEnergy_RT(doublereal *urt) const;

        /**
         * Get the nondimensional heat capacity at constant pressure
	 * function for the species
         * standard states at the current T and P of the solution.
         */
        virtual void getCp_R(doublereal* cpr) const;

	//@}
        /// @name Thermodynamic Values for the Species Reference States ---------------------
	//@{

	/**
	 *  Returns the vector of nondimensional
	 *  enthalpies of the reference state at the current temperature
	 *  and reference presssure for the species
	 */
        virtual void getEnthalpy_RT_ref(doublereal *hrt) const;
	/**
	 *  Returns the vector of nondimensional
	 *  enthalpies of the reference state at the current temperature
	 *  and reference pressure for the species.
	 */
        virtual void getGibbs_RT_ref(doublereal *grt) const;

	/**
	 *  Returns the vector of the 
	 *  gibbs function of the reference state at the current temperature
	 *  and reference pressure for the species.
	 *  units = J/kmol
	 */
        virtual void getGibbs_ref(doublereal *g) const;

	/**
	 *  Returns the vector of nondimensional
	 *  entropies of the reference state at the current temperature
	 *  and reference pressure for the species.
	 */
        virtual void getEntropy_R_ref(doublereal *er) const;

	/**
	 *  Returns the vector of nondimensional
         *  internal Energies of the reference state at the current temperature
         *  of the solution and the reference pressure for each species.
         */
        virtual void getIntEnergy_RT_ref(doublereal *urt) const;
	/**
	 *  Returns the vector of nondimensional
	 *  constant pressure heat capacities of the reference state
	 *  at the current temperature and reference pressure
	 *  for the species.
	 */
        virtual void  getCp_R_ref(doublereal *cprt) const;


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
            for (k = 0; k != m_kk; k++) m_expg0_RT[k] = exp(m_g0_RT[k]);
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

        virtual void initThermo();


	/** 
         * @internal
         * @name Chemical Equilibrium
         * @{
	 * 
	 * Set mixture to an equilibrium state consistent with specified 
	 * element potentials and temperature.
	 *
	 * @param lambda_RT vector of non-dimensional element potentials
	 * \f[ \lambda_m/RT \f].
	 * @param t temperature in K.
	 * @param work. Temporary work space. Must be dimensioned at least
	 * as large as the number of species. 
	 *
	 */
        virtual void setToEquilState(const doublereal* lambda_RT);

	// @}


    protected:

        int m_kk, m_mm;
        doublereal m_tmin, m_tmax, m_p0;

        mutable doublereal     m_tlast, m_logc0;
        mutable array_fp      m_h0_RT;
        mutable array_fp      m_cp0_R;
        mutable array_fp      m_g0_RT;
        mutable array_fp      m_s0_R;
        mutable array_fp      m_expg0_RT;
        mutable array_fp      m_pe;
        mutable array_fp      m_pp;

    private:

        void _updateThermo() const;
    };
}
        
#endif





