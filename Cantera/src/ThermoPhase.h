/**
 *
 *  @file ThermoPhase.h
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
     * These classes are used to compute thermodynamic properties.
     */

    /**
     * A phase with thermodynamic properties.
     * Extends class Phase by adding methods that compute 
     * thermodynamic properties. 
     *
     * Class ThermoPhase is the base class for the family of classes
     * that represent phases of matter with particular equations of
     * state. Instances of subclasses of ThermoPhase should be created
     * using the factory class ThermoFactory, not by calling the
     * constructor directly.  
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
        
        /// Constructor.
        ThermoPhase() : Phase(), m_spthermo(0), m_speciesData(0),
            m_index(-1), m_phi(0.0) {}

        
        virtual ~ThermoPhase() {
            delete m_spthermo;
        }


        /**
         *   
         * @name  Utilities  
         * @{
         */


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


        /// used to access data needed to construct transport manager
        /// later.
        void saveSpeciesData(const XML_Node* data) {
            m_speciesData = data;
        }

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
        virtual void initThermo() {}


        /** 
         * Equation of state type flag. The base class returns
         * zero. Subclasses should define this to return a unique
         * non-zero value. Constants defined for this purpose are
         * listed in mix_defs.h.
         */
        virtual int eosType() const { return 0; }



        /**
         * @} 
         * @name  Molar Thermodynamic Properties 
         * @{
         */


        /**
         * Molar enthalpy. Units: J/kmol. 
         */
        virtual doublereal enthalpy_mole() const {
            return err("enthalpy_mole");
        }


        /**
         * Molar internal energy. Units: J/kmol. 
         */
        virtual doublereal intEnergy_mole() const {
            return err("intEnergy_mole");
        }


        /**
         * Molar entropy. Units: J/kmol/K. 
         */
        virtual doublereal entropy_mole() const {
            return err("entropy_mole");
        }


        /**
         * Molar Gibbs function. Units: J/kmol. 
         */
        virtual doublereal gibbs_mole() const {
            return err("gibbs_mole");
        }


        /**
         * Molar heat capacity at constant pressure. Units: J/kmol/K. 
         */
        virtual doublereal cp_mole() const {
            return err("cp_mole");
        }


        /**
         * Molar heat capacity at constant volume. Units: J/kmol/K. 
         */
        virtual doublereal cv_mole() const {
            return err("cv_mole");
        }
        

        /**
         * @}
         * @name Mechanical Properties
         * @{
         */

        /**
         * Pressure. Units: Pa. 
	 *     Returns the thermodynamic pressure -> must be reimplemented
	 *     in inherited classes.
         */
        virtual doublereal pressure() const {
            return err("pressure");
        }


        /**
         * Set the pressure. Units: Pa. 
	 *     Sets the thermodynamic pressure -> must be reimplemented
	 *     in inherited classes.
         */
        virtual void setPressure(doublereal p) {
            err("setPressure");
        }

        virtual void updateDensity() {}

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

        void setElectricPotential(doublereal v) {
            m_phi = v;
        }

        doublereal electricPotential() { return m_phi; }


        /**
         * @}
         * @name Chemical Potentials and Activities
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
         * \f$ C_k \f$ that are defined such that \f$ a_k = C_k /
         * C^0_k, \f$ where \f$ C^0_k \f$ is a standard concentration
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
         * for an ideal gas \f$ C^0_k = P^0/\hat R T \f$. For this
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
	 *
	 * Returns the natural logarithm of the standard 
	 * concentration of the kth species
	 */
         virtual doublereal logStandardConc(int k=0) const {
             err("logStandardConc");
             return -1.0;
         }

        /** Get the array of chemical potentials at unit activity \f$
         * \mu^0_k \f$.
         */
        virtual void getStandardChemPotentials(doublereal* mu) const {
            err("getStandardChemPotentials");
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
         * Get the array of non-dimensional chemical potentials \f$
         * \mu_k / \hat R T \f$.
         */
        virtual void getChemPotentials_RT(doublereal* mu) const {
            err("getChemPotentials_RT");
        }
        
        /**
         * Get the species chemical potentials. Units: J/kmol.
         */
        virtual void getChemPotentials(doublereal* mu) const {
            err("getChemPotentials_RT");
        }

        //@}
        /// @name  Partial Molar Properties
        //@{

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

        /**
         * Get the nondimensional Enthalpy functions for the species
         * at their standard states at the current T and P.
         */
        virtual void getEnthalpy_RT(doublereal* hrt) const {
            err("getEnthalpy_RT");
        }


        /**
         * Get the nondimensional Entropies for the species
         * at their standard states at the current T and P.
         */
        virtual void getEntropy_R(doublereal* sr) const {
            err("getEntropy_R");
        }


        /**
         * Get the nondimensional Gibbs functions for the species
         * at their standard states at the current T and P.
         */
        virtual void getGibbs_RT(doublereal* grt) const {
            err("getGibbs_RT");
        }


	/**
	 * Get the nondimensional Gibbs functions for the pure species
	 * at the current T and P.
	 */
        virtual void getPureGibbs(doublereal* gpure) const {
            err("getPureGibbs");
        }


        /**
         * Get the nondimensional Gibbs functions for the pure species
         * at the current T and P. @deprecated
         */
        virtual void getCp_R(doublereal* cpr) const {
            err("getCp_RT");
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
            //cout << "entropy_mass. " << endl;
            //cout << "entropy_mole = " << entropy_mole() << endl;
            //cout << "meanMolecularWeight = " << meanMolecularWeight() << endl;
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

        /// @internal
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
         * @internal
         * @name Chemical Equilibrium
         * @{
         *
         * This method is used by the ChemEquil equilibrium solver.
         * It sets the state such that the chemical potentials satisfy
         * \f[ \frac{\mu_k}{\hat R T} = \sum_m A_{k,m}
         * \left(\frac{\lambda_m} {\hat R T}\right) \f] where \f$
         * \lambda_m \f$ is the element potential of element m. The
         * temperature is unchanged.  Any phase (ideal or not) that
         * implements this method can be equilibrated by ChemEquil.
         */ 
        virtual void setToEquilState(const doublereal* lambda_RT) {
            err("setToEquilState");
        }
        //@}

        void getActivities(doublereal* a);


        /**
         * @internal
         * Set equation of state parameters. The number and meaning of
         * these depends on the subclass. 
         * @param n number of parameters
         * @param c array of \i n coefficients
         */
        virtual void setParameters(int n, doublereal* c) {}
        
        
        //---------------------------------------------------------
        /// @name Critical state properties.
        
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
	 * Returns the reference pressure in Pa. This function is a wrapper
	 * that calls the species thermo refPressure function.
	 */
        doublereal refPressure() const {
            return m_spthermo->refPressure();
        }

        doublereal minTemp(int k = -1) {
            return m_spthermo->minTemp(k);
        }

        doublereal maxTemp(int k = -1) {
            return m_spthermo->maxTemp(k);
        }
        
            
    protected:

        /// Pointer to the species thermodynamic property manager
        SpeciesThermo* m_spthermo;

        const XML_Node* m_speciesData;

        /// Index number
        int m_index;
        doublereal m_phi;

    private:

        doublereal err(string msg) const;

    };

    typedef ThermoPhase thermophase_t;
    typedef ThermoPhase thermo_t;
}
        
#endif





