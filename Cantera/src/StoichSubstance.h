/**
 *
 *  @file StoichSubstance.h
 *
 */

/*  $Author$
 *  $Date$
 *  $Revision$
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
         * Equation of state flag. Returns the value cStoichSubstance,
         * defined in mix_defs.h.
         */
        virtual int eosType() const { return cStoichSubstance; }


        /**
         * @name Molar Thermodynamic Properties
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
         * term \f$ P_0 \hat v$ is subtracted from the specified molar
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


        virtual doublereal gibbs_mole() const {
            return enthalpy_mole() - temperature() * entropy_mole();
        }


        /**
         * Molar heat capacity at constant pressure. Units: J/kmol/K.
         * For an incompressible substance, \f$ \hat c_p = \hat c_v$.
         */
        virtual doublereal cp_mole() const {
            _updateThermo();
            return GasConstant * m_cp0_R[0];
        }

        /**
         * Molar heat capacity at constant volume. Units: J/kmol/K.
         * For an incompressible substance, \f$ \hat c_p = \hat c_v$.
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
         * For a stoichiometric substance, there is only one species. 
         * This method returns the molar gibbs function in the
         * first element of array \c mu.
         */
        virtual void getChemPotentials(doublereal* mu) const {
            mu[0] = gibbs_mole();
        }

        /**
         * For a stoichiometric substance, there is no activity term in 
         * the chemical potential expression, and therefore the
         * standard chemical potential and the chemical potential
         * are both equal to the molar Gibbs function.
         */
        virtual void getStandardChemPotentials(doublereal* mu0) const {
            mu0[0] = gibbs_mole();
        }

        /**
         * This method returns the array of generalized
         * concentrations.  For a stoichiomeetric substance, there is
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

         virtual doublereal logStandardConc(int k=0) const {
            return 0.0;
        }

        virtual void initThermo();


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





