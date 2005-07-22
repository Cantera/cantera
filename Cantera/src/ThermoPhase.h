/**
 *  @file ThermoPhase.h
 *
 * Header file for class ThermoPhase.
 *
 */

/*
 *  $Author$
 *  $Date$
 *  $Revision$
 *
 *  Copyright 2002 California Institute of Technology
 *
 */

#ifndef CT_THERMOPHASE_H
#define CT_THERMOPHASE_H

#include "Phase.h"


namespace Cantera {

    class XML_Node;


    /**
     * @defgroup thermoprops Thermodynamic Properties
     *
     * These classes are used to compute thermodynamic properties of
     * phases of matter.
     */

    /**
     * A phase with thermodynamic properties.  
     * Class ThermoPhase is the base class for the family of classes
     * that represent phases of matter of any type. It defines a
     * common public interface, and implements a few methods. Most of
     * the methods, however, are declared virtual and are meant to be
     * overloaded in derived classes.  The standard way used
     * throughout Cantera to compute properties of phases of matter is
     * through pointers of type ThermoPhase* that point to objects of
     * subclasses of ThermoPhase.
     * 
     * Class ThermoPhase
     * extends class Phase by adding methods to compute thermodynamic
     * properties in addition to the ones (temperature, density,
     * composition) that class Phase provides. The distinction is that
     * the methods declared in ThermoPhase require knowing the
     * particular equation of state of the phase of interest, while
     * those of class Phase do not, since they only involve data values
     * stored within the object.
     *
     * Instances of subclasses of ThermoPhase should be created using
     * the factory class ThermoFactory, not by calling the constructor
     * directly. This allows new classes to be used with the various
     * Cantera language interfaces.
     * 
     * To implement a new equation of state, derive a class from
     * ThermoPhase and overload the virtual methods in
     * ThermoPhase. Methods that are not needed can be left
     * unimplimented, which will cause an exception to be thrown if it
     * is called.
     * @ingroup thermoprops
     * @ingroup phases
     */
    class ThermoPhase : public Phase {

    public:

        /// Constructor. Note that ThermoPhase is meant to be used as
        /// a base class, so this constructor should not be called
        /// explicitly.
        ThermoPhase() : Phase(), m_spthermo(0), m_speciesData(0),
                        m_index(-1), m_phi(0.0), m_hasElementPotentials(false) {}


        /// Destructor. Deletes the species thermo manager.
        virtual ~ThermoPhase() {
            delete m_spthermo;
        }


        /**
         *   
         * @name  Information Methods  
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
	 * Returns the reference pressure in Pa. This function is a wrapper
	 * that calls the species thermo refPressure function.
	 */
        doublereal refPressure() const {
            return m_spthermo->refPressure();
        }

        /**
         * Minimum temperature for which the thermodynamic data for
         * the species are valid. If no argument is supplied, the
         * value returned will be the lowest temperature at which the
         * data for \e all species are valid. Otherwise, the value
         * will be only for species \a k. This function is a wrapper
	 * that calls the species thermo minTemp function.
         */
        doublereal minTemp(int k = -1) {
            return m_spthermo->minTemp(k);
        }

        /**
         * Maximum temperature for which the thermodynamic data for
         * the species are valid. If no argument is supplied, the
         * value returned will be the highest temperature at which the
         * data for \e all species are valid. Otherwise, the value
         * will be only for species \a k. This function is a wrapper
	 * that calls the species thermo maxTemp function.
         */
        doublereal maxTemp(int k = -1) {
            return m_spthermo->maxTemp(k);
        }
        

        /**
         * @} 
         * @name  Molar Thermodynamic Properties 
         * @{
         */

         /// Molar enthalpy. Units: J/kmol. 
        virtual doublereal enthalpy_mole() const {
            return err("enthalpy_mole");
        }

        /// Molar internal energy. Units: J/kmol. 
        virtual doublereal intEnergy_mole() const {
            return err("intEnergy_mole");
        }

        /// Molar entropy. Units: J/kmol/K. 
        virtual doublereal entropy_mole() const {
            return err("entropy_mole");
        }

        /// Molar Gibbs function. Units: J/kmol. 
        virtual doublereal gibbs_mole() const {
            return err("gibbs_mole");
        }

        /// Molar heat capacity at constant pressure. Units: J/kmol/K. 
        virtual doublereal cp_mole() const {
            return err("cp_mole");
        }

        /// Molar heat capacity at constant volume. Units: J/kmol/K. 
        virtual doublereal cv_mole() const {
            return err("cv_mole");
        }


        /**
         * @}
         * @name Mechanical Properties
         * @{
         */

        /**
         *  Pressure. Return the thermodynamic pressure (Pa). This
	 *  method must be overloaded in derived classes. Since the
	 *  mass density, temperature, and mass fractions are stored,
	 *  this method should use these values to implement the
	 *  mechanical equation of state \f$ P(T, \rho, Y_1, \dots,
	 *  Y_K) \f$.
         */
        virtual doublereal pressure() const {
            return err("pressure");
        }


        /**
         * Set the pressure (Pa).  This method must be reimplemented
	 *     in derived classes.
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
         * This method may optionally be defined in derived classes.
         */
        virtual doublereal isothermalCompressibility() const {
            err("isothermalCompressibility"); return -1.0;
        }


        /**
         * The volumetric thermal expansion coefficient. Units: 1/K.
         * The thermal expansion coefficient is defined as
         *
         * \f[
         * \beta = \frac{1}{v}\left(\frac{\partial v}{\partial T}\right)_P
         * \f]
         */
        virtual doublereal thermalExpansionCoeff() const {
            err("thermalExpansionCoeff()"); return -1.0;
        }

        /// @deprecated
        virtual void updateDensity() {
            deprecatedMethod("ThermoPhase","updateDensity","");
        }

        /**
         * @} 
         * @name Electric Potential
         * 
         * The phase may be at some non-zero electrical
	 * potential. These methods set or get the value of the
	 * electric potential.
         //@{
         */

//         /**
//          * Set the potential energy of species k to pe.
//          * Units: J/kmol.
// 	 * This function must be reimplemented in inherited classes
// 	 * of ThermoPhase.
//          @deprecated
//          */
//         virtual void setPotentialEnergy(int k, doublereal pe) {
//             deprecatedMethod("ThermoPhase","setPotentialEnergy","none");
//             err("setPotentialEnergy");
//         }

//         /**
//          * Get the potential energy of species k.
//          * Units: J/kmol.
// 	 * This function must be reimplemented in inherited classes
// 	 * of ThermoPhase.
//          * @deprecated
//          */
//         virtual doublereal potentialEnergy(int k) const {
//             deprecatedMethod("ThermoPhase","potentialEnergy","none");
//             return err("potentialEnergy");
//         }

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
         * @name Chemical Potentials and Activities
         *
         * The activity \f$a_k\f$ of a species in solution is related
         * to the chemical potential by \f[ \mu_k = \mu_k^0(T,P) +
         * \hat R T \log a_k. \f] The quantity \f$\mu_k^0(T,P)\f$ is
         * the chemical potential at unit activity, which depends on
         * temperature and pressure, but not on composition. The
         * activity is dimensionless.
         * @{
         */


        /**
         * This method returns an array of generalized concentrations
         * \f$ C_k\f$ that are defined such that \f$ a_k = C_k /
         * C^0_k, \f$ where \f$ C^0_k \f$ is a standard concentration
         * defined below.  These generalized concentrations are used
         * by kinetics manager classes to compute the forward and
         * reverse rates of elementary reactions. Note that they may
         * or may not have units of concentration --- they might be
         * partial pressures, mole fractions, or surface coverages,
         * for example.
         *
	 * @param c Output array of generalized concentrations. The 
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
	 * Natural logarithm of the standard 
	 * concentration of the kth species.
	 */
         virtual doublereal logStandardConc(int k=0) const {
             err("logStandardConc");
             return -1.0;
         }

        /**
         * Get the array of chemical potentials at unit activity These
         * are the standard state chemical potentials \f$ \mu^0_k(T,P)
         * \f$.. The values are evaluated at the current
         * temperature and pressure.
         */
        virtual void getStandardChemPotentials(doublereal* mu) const {
            err("getStandardChemPotentials");
        }

	/**
	 * Returns the units of the standard and generalized
	 * concentrations. Note they have the same units, as their
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


        //@}
        /// @name  Partial Molar Properties of the Solution
        //@{

	/**
         * Get the array of non-dimensional species chemical potentials
	 * These are partial molar Gibbs free energies.
         * \f$ \mu_k / \hat R T \f$.
	 * Units: unitless
         */
        virtual void getChemPotentials_RT(doublereal* mu) const {
            err("getChemPotentials_RT");
        }

	/**
         * Get the species chemical potentials in the solution
	 * These are partial molar Gibbs free energies.
	 * Units: J/kmol.
         */
        virtual void getChemPotentials(doublereal* mu) const {
            err("getChemPotentials");
        }

        /**
         * Get the species electrochemical potentials.  These are
	 * partial molar quantities.  This method adds a term \f$ Fz_k
	 * \phi_k \f$ to each chemical potential.
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

        void getActivities(doublereal* a);

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


        /**
         * Get the species partial molar enthalpies. Units: J/kmol.
         */
        virtual void getPartialMolarEnthalpies(doublereal* hbar) const {
            err("getPartialMolarEnthalpies");
        }

        /**
         * Get the species partial molar entropies. Units: J/kmol.
         */
        virtual void getPartialMolarEntropies(doublereal* sbar) const {
            err("getPartialMolarEntropies");
        }

        /**
         * Get the species partial molar volumes. Units: m^3/kmol.
         */
        virtual void getPartialMolarVolumes(doublereal* vbar) const {
            err("getPartialMolarVolumes");
        }

        //@}
        /// @name  Properties of the Standard State of the Species in the Solution 
        //@{

        /**
	 * Get the nondimensional Enthalpy functions for the species
         * at their standard states at the current
         * <I>T</I> and <I>P</I> of the solution.
         */
        virtual void getEnthalpy_RT(doublereal* hrt) const {
            err("getEnthalpy_RT");
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
         * of the solution.
         */
        virtual void getGibbs_RT(doublereal* grt) const {
            err("getGibbs_RT");
        }

	/**
	 * Get the nondimensional Gibbs functions for the standard
         * state of the species at the current T and P.
	 */
        virtual void getPureGibbs(doublereal* gpure) const {
            err("getPureGibbs");
        }

        /**
         * Get the nondimensional Heat Capacities at constant
         * pressure for the standard state of the species 
         * at the current T and P. 
         */
        virtual void getCp_R(doublereal* cpr) const {
            err("getCp_RT");
        }

	/**
         * Get the molar volumes of each species in their standard
         * states at the current
         * <I>T</I> and <I>P</I> of the solution.
         * units = m^3 / kmol
         */
        virtual void getStandardVolumes(doublereal *vol) const {
            err("getStandardVolumes");
        }


       //@}
        /// @name Thermodynamic Values for the Species Reference States 
        //@{

	/**
         *  Returns the vector of nondimensional
         *  enthalpies of the reference state at the current temperature
         *  of the solution and the reference pressure for the species.
         */
        virtual void getEnthalpy_RT_ref(doublereal *hrt) const {
            err("enthalpy_RT_ref");
        }
     
        /**
         *  Returns the vector of nondimensional
         *  enthalpies of the reference state at the current temperature
         *  of the solution and the reference pressure for the species.
         */
        virtual void getGibbs_RT_ref(doublereal *grt) const {
            err("gibbs_RT_ref");
        }
                   
        /**
         *  Returns the vector of the
         *  gibbs function of the reference state at the current temperature
         *  of the solution and the reference pressure for the species.
         *  units = J/kmol
         */
        virtual void  getGibbs_ref(doublereal *g) const {
            err("gibbs_ref");
        }
      
        /**
         *  Returns the vector of nondimensional
         *  entropies of the reference state at the current temperature
         *  of the solution and the reference pressure for the species.
         */
        virtual void getEntropy_R_ref(doublereal *er) const {
            err("entropy_R_ref");
        }
                 
        /**
         *  Returns the vector of nondimensional
         *  constant pressure heat capacities of the reference state
         *  at the current temperature of the solution
         *  and reference pressure for the species.
         */
        virtual void getCp_R_ref(doublereal *cprt) const {
            err("cp_R_ref()");
        }


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
         * Specific enthalpy. Units: J/kg. 
         */
        doublereal enthalpy_mass() const {
            return enthalpy_mole()/meanMolecularWeight();
        }

        /**
         * Specific internal energy. Units: J/kg. 
         */
        doublereal intEnergy_mass() const {
            return intEnergy_mole()/meanMolecularWeight();
        }

        /**
         * Specific entropy. Units: J/kg/K. 
         */
        doublereal entropy_mass() const {
            return entropy_mole()/meanMolecularWeight();
        }

        /**
         * Specific Gibbs function. Units: J/kg. 
         */
        doublereal gibbs_mass() const {
            return gibbs_mole()/meanMolecularWeight();
        }

        /**
         * Specific heat at constant pressure. Units: J/kg/K. 
         */
        doublereal cp_mass() const {
            return cp_mole()/meanMolecularWeight();
        }

        /**
         * Specific heat at constant volume. Units: J/kg/K. 
         */
        doublereal cv_mass() const {
            return cv_mole()/meanMolecularWeight();
        }
        //@}

        doublereal _RT() const {
            return temperature() * GasConstant;
        }

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

        /** Set the temperature (K) and pressure (Pa) */
        void setState_TP(doublereal t, doublereal p);

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

        // called by function 'equilibrate' in ChemEquil.h to transfer
        // the element potentials to this object
        void setElementPotentials(const vector_fp& lambda) {
            m_lambda = lambda;
            m_hasElementPotentials = true;
        }

        void getElementPotentials(doublereal* lambda) {
            if (m_hasElementPotentials)
                copy(m_lambda.begin(), m_lambda.end(), lambda);
        }

        //@}
    

        //---------------------------------------------------------

        /// @name Critical state properties.
        /// These methods are only implemented by some subclasses, and may 
        /// be moved out of ThermoPhase at a later date.
        
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
        /// implement full liquid-vapor equations of state. They may be
        /// moved out of ThermoPhase at a later date.
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


        /// @name For Internal Use

        /// The following methods are used in the process of constructing
        /// the phase and setting its parameters from a specification in an 
        /// input file. They are not normally used in application programs.
        /// To see how they are used, see files importCTML.cpp and 
        /// ThermoFactory.cpp.
        //@{

        /// Store a reference to the XML tree containing the species
        /// data for this phase. This is used to access data needed to
        /// construct transport manager later.
        /// @internal
        void saveSpeciesData(const XML_Node* data) {
            m_speciesData = data;
        }

        /// Return a pointer to the XML tree containing the species
        /// data for this phase.
        const XML_Node* speciesData() { 
            if (m_speciesData) 
                return m_speciesData;
            else {
                throw CanteraError("ThermoPhase::speciesData",
                    "m_speciesData is NULL");
                return 0;
            }
        }


        /**
         * @internal Install a species thermodynamic property
         * manager. The species thermodynamic property manager
         * computes properties of the pure species for use in
         * constructing solution properties. It is meant for internal
         * use, and some classes derived from ThermoPhase may not use
         * any species thermodynamic property manager. This method is
         * called by function importPhase in importCTML.cpp.
         */
        void setSpeciesThermo(SpeciesThermo* spthermo) 
            { m_spthermo = spthermo; }

        /**
         * @internal Return a reference to the species thermodynamic property
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
        virtual void initThermo() { }



        // The following methods are used by the clib interface
        // library, and should not be used by application programs.

        /**
         * @internal 
         * Index number.  This method can be used to identify the
         * location of a phase object in a list, and is used by the
         * interface library (clib) routines for this purpose.
         */
        int index() { return m_index; }


        /**
         * @internal Set the index number. The Cantera interface
         * library uses this method to set the index number to the
         * location of the pointer to this object in the pointer array
         * it maintains. Using this method for any other purpose will
         * lead to unpredictable results if used in conjunction with
         * the interface library.
        */ 
        void setIndex(int m) { m_index = m; }


        /**
         * @internal
         * Set equation of state parameters. The number and meaning of
         * these depends on the subclass. 
         * @param n number of parameters
         * @param c array of \a n coefficients
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
        virtual void setStateFromXML(const XML_Node& state);


        //@}

            
    protected:

        /// Pointer to the species thermodynamic property manager
        SpeciesThermo* m_spthermo;

        const XML_Node* m_speciesData;

        /// Index number
        int m_index;
        doublereal m_phi;
        vector_fp m_lambda;
        bool m_hasElementPotentials;

    private:

        doublereal err(string msg) const;

    };

    typedef ThermoPhase thermophase_t;
    typedef ThermoPhase thermo_t;
}
        
#endif





