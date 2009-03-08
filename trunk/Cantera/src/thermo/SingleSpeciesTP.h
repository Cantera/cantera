/**
 *  @file SingleSpeciesTP.h
 *
 * Header file for class SingleSpeciesTP
 *
 */

/*
 * Copywrite (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */


/*
 *  $Author: hkmoffa $
 *  $Date: 2005/10/24 21:52:23 $
 *  $Revision: 1.2 $
 *
 */

#ifndef CT_SINGLESPECIESTP_H
#define CT_SINGLESPECIESTP_H

#include "ThermoPhase.h"


namespace Cantera {

    /**
     * @ingroup thermoprops
     *
     *  The SingleSpeciesTP class is a filter class for ThermoPhase.
     *  What it does is to simplify the construction of ThermoPhase
     *  objects by assuming that the phase consists of one and 
     *  only one type of species. In other words, it's a stoichiometric
     *  phase. However, no assumptions are made concerning the
     *  thermodynamic functions or the equation of state of the
     *  phase. Therefore it's an incomplete description of
     *  the thermodynamics. The complete description must be
     *  made in a derived class of SingleSpeciesTP.
     *  \nosubgrouping
     */
    class SingleSpeciesTP : public ThermoPhase {

    public:
        
        /// Constructor. 
        SingleSpeciesTP();

        /// Destructor
        virtual ~SingleSpeciesTP();

        /**
         *   
         * @name  Information Methods  
         * @{
         */

        /** 
         * Returns the equation of state type flag.
	 * This is a modified base class.
	 * Therefore, if not overridden in derivied classes,
	 * this call will throw an exception.
         */
        virtual int eosType() const;

        /**
         * @} 
         * @name  Molar Thermodynamic Properties of the Solution
	 *
	 *  These functions are resolved at this level, by reference
	 *  to the partial molar functions and standard state
	 *  functions for species 0. Derived classes don't need
	 *  to supply entries for these functions.
         * @{
         */

         /// Molar enthalpy. Units: J/kmol. 
	doublereal enthalpy_mole() const;

        /// Molar internal energy. Units: J/kmol. 
        doublereal intEnergy_mole() const;

        /// Molar entropy. Units: J/kmol/K. 
        doublereal entropy_mole() const;

        /// Molar Gibbs function. Units: J/kmol. 
        doublereal gibbs_mole() const;

        /// Molar heat capacity at constant pressure. Units: J/kmol/K. 
	doublereal cp_mole() const;

        /// Molar heat capacity at constant volume. Units: J/kmol/K. 
	doublereal cv_mole() const;

        /**
         * @}
         * @name Mechanical Properties
         * @{
         */

        /**
         *  Pressure. Return the thermodynamic pressure (Pa). This
	 *  method must be reimplemented in derived classes.
	 *  Since the mass density, temperature, and mass fractions 
	 *  are stored, this method should use these 
         *  values to implement the mechanical equation of state 
         *  \f$ P(T, \rho, Y_1, \dots, Y_K) \f$. 
         */
        virtual doublereal pressure() const {
            return err("pressure");
        }

        /**
         * Set the pressure. 
	 *     Sets the thermodynamic pressure -> must be reimplemented
	 *     in derived classes. Units: Pa. 
         */
        virtual void setPressure(doublereal p) {
            err("setPressure");
        }

        /**
         * The isothermal compressibility. Units: 1/Pa.
         * The isothermal compressibility is defined as
         * \f[
         * \kappa_T = -\frac{1}{v}\left(\frac{\partial v}{\partial P}\right)_T
         * \f]
         */
        virtual doublereal isothermalCompressibility() const {
            err("isothermalCompressibility"); return -1.0;
        }

        /**
         * The thermal expansion coefficient. Units: 1/K.
         * The thermal expansion coefficient is defined as
         *
         * \f[
         * \beta = \frac{1}{v}\left(\frac{\partial v}{\partial T}\right)_P
         * \f]
         */
        virtual doublereal thermalExpansionCoeff() const {
            err("thermalExpansionCoeff()"); return -1.0;
        }

	/**
         * @} 
         * @name Electric Potential
         * 
         * The phase may be at some non-zero electrical
	 * potential. These methods set or get the value of the
	 * electric potential.
	 */
         //@{

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
         * @}
         * @name Activities, Standard State, and Activity Concentrations
         *
         * The activity \f$a_k\f$ of a species in solution is
         * related to the chemical potential by \f[ \mu_k = \mu_k^0(T)
         * + \hat R T \log a_k. \f] The quantity \f$\mu_k^0(T)\f$ is
         * the chemical potential at unit activity, which depends only
         * on temperature.
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
	 *           units depend upon the implementation of the
	 *           reaction rate expressions within the phase.
         */
        virtual void getActivityConcentrations(doublereal* c) const {
            err("getActivityConcentrations");
        }

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
         virtual doublereal standardConcentration(int k=0) const {
             err("standardConcentration");
             return -1.0;
         }

        /**
	 * Returns the natural logarithm of the standard 
	 * concentration of the kth species
	 */
         virtual doublereal logStandardConc(int k=0) const {
             err("logStandardConc");
             return -1.0;
         }

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
         * Get the array of non-dimensional activities at
         * the current solution temperature, pressure, and
         * solution concentration.
	 *
	 * We redefine this function to just return 1.0 here.
         */
        virtual void getActivities(doublereal* a) {
	    a[0] = 1.0;
	}

	/**
         * Get the array of non-dimensional activity coefficients at
         * the current solution temperature, pressure, and
         * solution concentration.
         */
        virtual void getActivityCoefficients(doublereal* ac) const {
            if (m_kk == 1) {
              ac[0] = 1.0;
            } else {
              err("getActivityCoefficients");
            }
        }    

        //@}
        /// @name  Partial Molar Properties of the Solution
	///
	///  These functions are resolved at this level, by reference
	///  to the partial molar functions and standard state
	///  functions for species 0. Derived classes don't need
	///  to supply entries for these functions.
        //@{

	/*
	 * These functions are all resolved here to point to the
	 * standard state functions for species 0
	 */

	/**
         * Get the array of non-dimensional species chemical potentials
	 * These are partial molar Gibbs free energies.
         * \f$ \mu_k / \hat R T \f$.
	 * Units: unitless
         */
        void getChemPotentials_RT(doublereal* mu) const;

	/**
         * Get the species chemical potentials in the solution
	 * These are partial molar Gibbs free energies.
	 * Units: J/kmol.
         */
        void getChemPotentials(doublereal* mu) const;

	/**
         * Get the species electrochemical potentials. Units: J/kmol.
         * This method adds a term \f$ Fz_k \phi_k \f$ to 
         * each chemical potential.
	 *
	 * This is resolved here. A single single species phase
	 * is not allowed to have anything other than a zero
	 * charge.
         */
        void getElectrochemPotentials(doublereal* mu) const;

        /**
         * Get the species partial molar enthalpies. Units: J/kmol.
         */
        void getPartialMolarEnthalpies(doublereal* hbar) const;

        /**
         * Get the species partial molar internal energies. Units: J/kmol.
         */
        virtual void getPartialMolarIntEnergies(doublereal* ubar) const;

        /**
         * Get the species partial molar entropies. Units: J/kmol.
         */
        void getPartialMolarEntropies(doublereal* sbar) const;

        /**
         * Get the species partial molar volumes. Units: m^3/kmol.
         */
        void getPartialMolarVolumes(doublereal* vbar) const;

        //@}
        /// @name  Properties of the Standard State of the Species in the Solution
	/// These functions are the primary way real properties are
	/// supplied to derived thermodynamics classes of SingleSpeciesTP.
	/// These functions must be supplied in derived classes. They
	/// are not resolved at the SingleSpeciesTP level.
        //@{

	/**
         * Get the array of chemical potentials at unit activity.
         * These are the standard state chemical potentials. 
         * \f$ \mu^0_k(T,P) \f$. The values are evaluated at the current
         * temperature and pressure.
         */
        virtual void getStandardChemPotentials(doublereal* mu) const {
            err("getStandardChemPotentials");
        }

        /**
	 * Get the nondimensional Enthalpy functions for the species
         * at their standard states at the current
         * <I>T</I> and <I>P</I> of the solution.
         */
        virtual void getEnthalpy_RT(doublereal* hrt) const {
            err("getEnthalpy_RT");
        }

       /**
	* Get the nondimensional Enthalpy functions for the species
	* at their standard states at the current
	* <I>T</I> and <I>P</I> of the solution.
	*/
        virtual void getIntEnergy_RT(doublereal* urt) const {
            err("getIntEnergy_RT");
        }

        /**
	 * Get the array of nondimensional Enthalpy functions for the
         * standard state species
         * at the current <I>T</I> and <I>P</I> of the solution.
         */
        virtual void getEntropy_R(doublereal* sr) const {
            err("getEntropy_R");
        }

        /**
         * Get the nondimensional Gibbs functions for the species
         * at their standard states of solution at the current T and P
         * of the solution
         */
        virtual void getGibbs_RT(doublereal* grt) const {
            err("getGibbs_RT");
        }

        /**
         * Get the dimensional Gibbs functions for the standard
         * state of the species at the current T and P.
         */
        void getPureGibbs(doublereal* gpure) const;

        /**
         * Get the nondimensional Gibbs functions for the standard
         * state of the species at the current T and P. 
         */
        virtual void getCp_R(doublereal* cpr) const {
            err("getCp_RT");
        }

	/**
         * Get the molar volumes of each species in their standard
         * states at the current
         * <I>T</I> and <I>P</I> of the solution.
         * units = m^3 / kmol
	 *
	 * We resolve this function at this level, by assigning 
	 * the molec weight divided by the phase density
         */
	void getStandardVolumes(doublereal *vol) const;


       //@}
        /// @name Thermodynamic Values for the Species Reference State
	///
	/// Almost all functions in this group are resolved by this
	/// class. It is assumed that the m_spthermo species thermo
	/// pointer is populated and yields the reference state.
	/// The internal energy function is not given by this
	/// class, since it would involve a specification of the
	/// equation of state.
        //@{

	/**
         *  Returns the vector of nondimensional
         *  enthalpies of the reference state at the current temperature
         *  of the solution and the reference pressure for the species.
         */
        virtual void getEnthalpy_RT_ref(doublereal *hrt) const;
     
        /**
         *  Returns the vector of nondimensional
         *  enthalpies of the reference state at the current temperature
         *  of the solution and the reference pressure for the species.
         */
        virtual void getGibbs_RT_ref(doublereal *grt) const;
                   
        /**
         *  Returns the vector of the
         *  gibbs function of the reference state at the current temperature
         *  of the solution and the reference pressure for the species.
         *  units = J/kmol
         */
        virtual void  getGibbs_ref(doublereal *g) const;
      
        /**
         *  Returns the vector of nondimensional
         *  entropies of the reference state at the current temperature
         *  of the solution and the reference pressure for the species.
         */
        virtual void getEntropy_R_ref(doublereal *er) const;
                 
        /**
         *  Returns the vector of nondimensional
         *  constant pressure heat capacities of the reference state
         *  at the current temperature of the solution
         *  and reference pressure for the species.
         */
        virtual void getCp_R_ref(doublereal *cprt) const;

        /**
         * @name Setting the State
         *
         * These methods set all or part of the thermodynamic
         * state.
         * @{
         */
        /** Set the temperature (K), pressure (Pa), and mole fractions.  */
        void setState_TPX(doublereal t, doublereal p, const doublereal* x);

        /** Set the temperature (K), pressure (Pa), and mole fractions.  */
        void setState_TPX(doublereal t, doublereal p, compositionMap& x);

        /** Set the temperature (K), pressure (Pa), and mole fractions.  */
        void setState_TPX(doublereal t, doublereal p, const string& x);

        /** Set the temperature (K), pressure (Pa), and mass fractions. */
        void setState_TPY(doublereal t, doublereal p, const doublereal* y);

        /** Set the temperature (K), pressure (Pa), and mass fractions. */
        void setState_TPY(doublereal t, doublereal p, compositionMap& y);
        
        /** Set the temperature (K), pressure (Pa), and mass fractions.  */
        void setState_TPY(doublereal t, doublereal p, const string& y);

        /** Set the pressure (Pa) and mole fractions.  */
        void setState_PX(doublereal p, doublereal* x);

        /** Set the pressure (Pa) and mass fractions.  */
        void setState_PY(doublereal p, doublereal* y);


        /** Set the specific enthalpy (J/kg) and pressure (Pa). */
        virtual void setState_HP(doublereal h, doublereal p, 
            doublereal tol = 1.e-8);

        /** Set the specific enthalpy (J/kg) and specific volume (m^3/kg). */
        virtual void setState_UV(doublereal u, doublereal v, 
            doublereal tol = 1.e-8);

        /** Set the specific entropy (J/kg/K) and pressure (Pa). */
        virtual void setState_SP(doublereal s, doublereal p, 
            doublereal tol = 1.e-8);

        /** Set the specific entropy (J/kg/K) and specific volume (m^3/kg). */
        virtual void setState_SV(doublereal s, doublereal v, 
            doublereal tol = 1.e-8);

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

        //@}



        /**
         * @internal
         * Set equation of state parameters. The number and meaning of
         * these depends on the subclass. 
         * @param n number of parameters
         * @param c array of \i n coefficients
         * 
         */
        virtual void setParameters(int n, doublereal* c) {}
	virtual void getParameters(int &n, doublereal * const c) {}

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
        virtual void setParametersFromXML(const XML_Node& eosdata) {}
	
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


   protected:


        doublereal m_tmin, m_tmax, m_press, m_p0;

	/**
	 * Last temperature used to evaluate the thermodynamic
	 * polynomial.
	 */
        mutable doublereal     m_tlast;
        mutable array_fp       m_h0_RT;
        mutable array_fp       m_cp0_R;
        mutable array_fp       m_s0_R;

    protected:

        void _updateThermo() const;

    private:
        doublereal err(string msg) const;

    };

}
        
#endif



