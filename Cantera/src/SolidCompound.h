/**
 *
 *  @file SolidPhase.h
 *
 */

/*  $Author$
 *  $Date$
 *  $Revision$
 *
 *  Copyright 2001 California Institute of Technology
 *
 */


#ifndef CT_SOLIDPHASE_H
#define CT_SOLIDPHASE_H

//#include "ct_defs.h"
#include "mix_defs.h"
#include "ThermoPhase.h"
#include "SpeciesThermo.h"

namespace Cantera {

    /**
     * @ingroup thermoprops
     *
     * Class IdealGasPhase represents low-density gases that obey the
     * ideal gas equation of state. It derives from class ThermoPhase,
     * and overloads the virtual methods defined there with ones that
     * use expressions appropriate for ideal gas mixtures.
     *
     */
    class SolidCompound : public ThermoPhase  {

    public:

        SolidCompound(): m_tlast(-1.0), m_tmin(0.0), m_tmax(0.0),
                         m_press(OneAtm), m_p0(OneAtm) {}

        virtual ~SolidCompound() {}

        /**
         * Equation of state flag. Returns the value cIdealGas, defined 
         * in mix_defs.h.
         */
        virtual int eosType() const { return cSolidCompound; }


        /**
         * @name Molar Thermodynamic Properties
         * @{
         */

        /**
         * Molar enthalpy. Units: J/kmol.
         */
        virtual doublereal enthalpy_mole() const {
            double hh = intEnergy_mole() + m_press / molarDensity();
            return hh;
        }

        /**
         * Molar internal energy. J/kmol. 
         */
        virtual doublereal intEnergy_mole() const {
            _updateThermo();
            //            cout << "intEnergy: " << m_h0_RT[0] << "  " << m_p0/molarDensity()
            //      << endl;
            return GasConstant * temperature() * m_h0_RT[0]
                - m_p0 / molarDensity();
        }

        /**
         * Molar entropy. Units: J/kmol/K.
         */
        virtual doublereal entropy_mole() const {
            _updateThermo();
            //cout << "s/r = " << m_s0_R[0] << endl;
            return GasConstant * m_s0_R[0];
        }


        virtual doublereal gibbs_mole() const {
            return enthalpy_mole() - temperature() * entropy_mole();
        }


        /**
         * Molar heat capacity at constant pressure. Units: J/kmol/K.
         */
        virtual doublereal cp_mole() const {
            _updateThermo();
            return GasConstant * m_cp0_R[0];
        }

        /**
         * Molar heat capacity at constant volume. Units: J/kmol/K.
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
         */ 
        virtual doublereal pressure() const {
            return m_press;
        }

        /**
         * Set the pressure at constant temperature. Units: Pa.
         */
        virtual void setPressure(doublereal p) {
            m_press = p;
        }

        //@}


        virtual void getChemPotentials(doublereal* mu) const {
            mu[0] = gibbs_mole();
        }

        virtual void getStandardChemPotentials(doublereal* mu0) const {
            mu0[0] = gibbs_mole();
            //cout << m_h0_RT[0] << "   " << m_s0_R[0] << endl;
            //cout << "std chem pot = " << mu0[0] << endl;
        }

        /**
         * This method returns the array of generalized
         * concentrations.  For a solid compound, there is only one
         * species, and the generalized concentration is 1.0.
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
        doublereal m_tmin, m_tmax, m_p0, m_press;

        mutable doublereal     m_tlast;
        mutable array_fp      m_h0_RT;
        mutable array_fp      m_cp0_R;
        mutable array_fp      m_s0_R;

    private:

        void _updateThermo() const;
    };
        
}
        
#endif





