/**
 *
 *  @file StoichSubstance.h
 *
 * This file contains the class declarations for the StoichSubstance
 * ThermoPhase class.
 */

/*  $Author: hkmoffa $
 *  $Date: 2006/06/08 14:05:41 $
 *  $Revision: 1.8 $
 *
 *  Copyright 2001 California Institute of Technology
 *
 */


#ifndef CT_STOICHSUBSTANCE_H
#define CT_STOICHSUBSTANCE_H

#include "mix_defs.h"
#include "ThermoPhase.h"
#include "SpeciesThermo.h"

namespace Cantera {

    /**
     * @ingroup thermoprops
     *
     * Class StoichSubstance represents a stoichiometric (fixed composition) 
     * incompressible substance.
     * \nosubgrouping
     *
     */
    class StoichSubstance : public ThermoPhase {

    public:

        StoichSubstance():
	    m_kk(0),
	    m_tmin(0.0),
	    m_tmax(0.0),
	    m_press(OneAtm),
	    m_p0(OneAtm),
	    m_tlast(-1.0)  {}

        virtual ~StoichSubstance() {}

	/**
         *   
         * @name  Utilities  
         * @{
         */

        /**
         * Equation of state flag. Returns the value cStoichSubstance,
         * defined in mix_defs.h.
         */
        virtual int eosType() const { return cStoichSubstance; }


        /**
	 * @}
         * @name Molar Thermodynamic Properties of the Solution ---------
         * @{
         */

        /**
         * Molar enthalpy. Units: J/kmol.  For an incompressible,
         * stoichiometric substance, the internal energy is
         * independent of pressure, and therefore the molar enthalpy
         * is \f[ \hat h(T, P) = \hat u(T) + P \hat v \f], where the
         * molar specific volume is constant.
         */
        virtual doublereal enthalpy_mole() const {
            double hh = intEnergy_mole() + m_press / molarDensity();
            return hh;
        }

        /**
         * Molar internal energy. J/kmol.  For an incompressible,
         * stoichiometric substance, the molar internal energy is
         * independent of pressure. Since the thermodynamic properties
         * are specified by giving the standard-state enthalpy, the
         * term \f$ P_0 \hat v\f$ is subtracted from the specified molar
         * enthalpy to compute the molar internal energy.
         */
        virtual doublereal intEnergy_mole() const {
            _updateThermo();
            return GasConstant * temperature() * m_h0_RT[0]
                - m_p0 / molarDensity();
        }

        /**
         * Molar entropy. Units: J/kmol/K.  For an incompressible,
         * stoichiometric substance, the molar entropy depends only on
         * the temperature.
         */
        virtual doublereal entropy_mole() const {
            _updateThermo();
            return GasConstant * m_s0_R[0];
        }


	/**
	 * Molar gibbs Function. Units: J/kmol. This is determined
	 * from the molar enthalpy and entropy functions.
	 */
        virtual doublereal gibbs_mole() const {
            return enthalpy_mole() - temperature() * entropy_mole();
        }


        /**
         * Molar heat capacity at constant pressure. Units: J/kmol/K.
         * For an incompressible substance, \f$ \hat c_p = \hat c_v\f$.
         */
        virtual doublereal cp_mole() const {
            _updateThermo();
            return GasConstant * m_cp0_R[0];
        }

        /**
         * Molar heat capacity at constant volume. Units: J/kmol/K.
         * For an incompressible substance, \f$ \hat c_p = \hat c_v\f$.
         */
        virtual doublereal cv_mole() const {
            return cp_mole();
        }

        //@}


        /**
         * @name Mechanical Equation of State
         * @{
         */

        /**
         * Pressure. Units: Pa.
         * For an incompressible substance, the density is independent
         * of pressure. This method simply returns the stored
         * pressure value.
         */ 
        virtual doublereal pressure() const {
            return m_press;
        }

        /**
         * Set the pressure at constant temperature. Units: Pa.
         * For an incompressible substance, the density is 
         * independent of pressure. Therefore, this method only 
         * stores the specified pressure value. It does not 
         * modify the density.
         */
        virtual void setPressure(doublereal p) {
            m_press = p;
        }

        //@}

	/**
         * @name Chemical Potentials and Activities
	 *@{
         */

        /**
         * This method returns the array of generalized
         * concentrations.  For a stoichiometric substance, there is
         * only one species, and the generalized concentration is 1.0.
         */
        virtual void getActivityConcentrations(doublereal* c) const {
            c[0] = 1.0;
        }

        /**
         * The standard concentration. This is defined as the concentration 
         * by which the generalized concentration is normalized to produce 
         * the activity. 
         */ 
         virtual doublereal standardConcentration(int k=0) const {
             return 1.0;
        }

        /**
	 * Returns the natural logarithm of the standard 
	 * concentration of the kth species
	 */
         virtual doublereal logStandardConc(int k=0) const {
            return 0.0;
        }

        /**
	 * Get the array of chemical potentials at unit activity 
         * \f$ \mu^0_k \f$.
	 *
         * For a stoichiometric substance, there is no activity term in 
         * the chemical potential expression, and therefore the
         * standard chemical potential and the chemical potential
         * are both equal to the molar Gibbs function.
         */
        virtual void getStandardChemPotentials(doublereal* mu0) const {
            mu0[0] = gibbs_mole();
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
	 *  uA[0] = kmol units - default  = 0
	 *  uA[1] = m    units - default  = 0
	 *  uA[2] = kg   units - default  = 0;
	 *  uA[3] = Pa(pressure) units - default = 0;
	 *  uA[4] = Temperature units - default = 0;
	 *  uA[5] = time units - default = 0
	 */
	virtual void getUnitsStandardConc(double *uA, int k = 0,
					  int sizeUA = 6);


	//@}
        /// @name  Partial Molar Properties of the Solution ----------------------------------
        //@{


        /** 
         * Get the array of non-dimensional chemical potentials 
         * \f$ \mu_k / \hat R T \f$.
         */
        virtual void getChemPotentials_RT(doublereal* mu) const {
            mu[0] = gibbs_mole() / (GasConstant * temperature());
        }

	/**
         * For a stoichiometric substance, there is only one species. 
         * This method returns the molar gibbs function in the
         * first element of array \c mu.
         */
        virtual void getChemPotentials(doublereal* mu) const {
            mu[0] = gibbs_mole();
        }

        /**
         * Get the species electrochemical potentials. Units: J/kmol.
         * This method adds a term \f$ Fz_k \phi_k \f$ to the 
         * to each chemical potential.
         */
        void getElectrochemPotentials(doublereal* mu) const {
	    getChemPotentials(mu);
	}

	/**
	 * Returns an array of partial molar enthalpies for the species
	 * in the mixture.
	 * Units (J/kmol)
         */
        virtual void getPartialMolarEnthalpies(doublereal* hbar) const {
	    hbar[0] = enthalpy_mole();
        }

	/**
         * Returns an array of partial molar entropies of the species in the
	 * solution. Units: J/kmol/K.
         */
        virtual void getPartialMolarEntropies(doublereal* sbar) const {
	    sbar[0] = entropy_mole();
        }

	/**
         * returns an array of partial molar volumes of the species
	 * in the solution. Units: m^3 kmol-1.
         */
        virtual void getPartialMolarVolumes(doublereal* vbar) const {
	    vbar[0] = 1.0 / molarDensity();
        }


      //@}
        /// @name  Properties of the Standard State of the Species in the Solution -------------------------------------
        //@{
       /**
         * Get the nondimensional Enthalpy functions for the species
         * at their standard states at the current 
	 * <I>T</I> and <I>P</I> of the solution.
         */
        virtual void getEnthalpy_RT(doublereal* hrt) const {
            hrt[0] = enthalpy_mole() / (GasConstant * temperature());
        }


        /**
	 * Get the array of nondimensional Enthalpy functions for the
         * standard state species
         * at the current <I>T</I> and <I>P</I> of the solution.
         */
        virtual void getEntropy_R(doublereal* sr) const {
            sr[0] = entropy_mole() / GasConstant;
        }

	/**
         * Get the nondimensional Gibbs functions for the species
         * at their standard states of solution at the current T and P
         * of the solution.
         */
        virtual void getGibbs_RT(doublereal* grt) const {
          grt[0] =  gibbs_mole() / (GasConstant * temperature());
        }

	/**
         * Get the nondimensional Heat Capacities at constant
         * pressure for the standard state of the species 
         * at the current T and P. 
         */
        virtual void getCp_R(doublereal* cpr) const {
            cpr[0] = cp_mole() / GasConstant;
        }

       /**
	* Get the standard volumes for the standard state of the species
	* at the current T and P
	*/
        virtual void getStandardVolumes(doublereal*vol) const {
	  vol[0] = 1.0 / molarDensity();
        }

	//@}
        /// @name Thermodynamic Values for the Species Reference States --------------------
	//@{

	/**
	 *  Returns the vector of nondimensional
	 *  enthalpies of the reference state at the current temperature
	 *  of the solution and the reference pressure for the species. 
	 *
	 *  This function fills in its one entry in hrt[] by calling
	 *  the underlying species thermo function for the 
	 *  dimensionless enthalpy.
	 */
        virtual void getEnthalpy_RT_ref(doublereal *hrt) const {
	    _updateThermo();
	    hrt[0] = m_h0_RT[0];
	}

	/**
	 *  Returns the vector of nondimensional
	 *  enthalpies of the reference state at the current temperature
	 *  of the solution and the reference pressure for the species.
	 *
	 *  This function fills in its one entry in hrt[] by calling
	 *  the underlying species thermo function for the 
	 *  dimensionless gibbs free energy, calculated from the
	 *  dimensionless enthalpy and entropy.
	 */
        virtual void getGibbs_RT_ref(doublereal *grt) const {
	    _updateThermo();
	    grt[0] = m_h0_RT[0] - m_s0_R[0];
        }

	/**
	 *  Returns the vector of the 
	 *  gibbs function of the reference state at the current temperature
	 *  of the solution and the reference pressure for the species.
	 *  units = J/kmol
	 *
	 *  This function fills in its one entry in g[] by calling
	 *  the underlying species thermo functions for the 
	 *  gibbs free energy, calculated from enthalpy and the
	 *  entropy, and the multiplying by RT.
	 */
        virtual void  getGibbs_ref(doublereal *g) const {
	    getGibbs_RT_ref(g);
	    g[0] *= GasConstant * temperature();
        }

	/**
	 *  Returns the vector of nondimensional
	 *  entropies of the reference state at the current temperature
	 *  of the solution and the reference pressure for the species.
	 *
	 *  This function fills in its one entry in hrt[] by calling
	 *  the underlying species thermo function for the 
	 *  dimensionless entropy.
	 */
        virtual void getEntropy_R_ref(doublereal *er) const {
	    _updateThermo();
            er[0] = m_s0_R[0];
	}


        virtual void initThermo();

        virtual void setParameters(int n, double *c);
        virtual void getParameters(int &n, double * const c);

	virtual void setParametersFromXML(const XML_Node& eosdata);

    protected:

        int m_kk;
        doublereal m_tmin, m_tmax, m_press, m_p0;

        mutable doublereal     m_tlast;
        mutable array_fp      m_h0_RT;
        mutable array_fp      m_cp0_R;
        mutable array_fp      m_s0_R;

    private:

        void _updateThermo() const;
    };
    
}
        
#endif





