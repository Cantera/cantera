/**
 *
 *  @file IdealGasPhase.h
 *
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
         * @name Molar Thermodynamic Properties
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
         * @see SpeciesThermo
         */
        virtual doublereal enthalpy_mole() const {
            return GasConstant * temperature() * 
                mean_X(enthalpy_RT().begin());
        }

        /**
         * Molar internal energy. J/kmol. For an ideal gas mixture,
         * \f[
         * \hat u(T) = \sum_k X_k \hat h^0_k(T) - \hat R T,
         * \f]
         * and is a function only of temperature.
         * The standard-state pure-species enthalpies 
         * \f$ \hat h^0_k(T) \f$ are computed by the species thermodynamic 
         * property manager.
         * @see SpeciesThermo
         */
        virtual doublereal intEnergy_mole() const {
            return GasConstant * temperature()
                * ( mean_X(enthalpy_RT().begin()) - 1.0);
        }


        /**
         * Molar entropy. Units: J/kmol/K.
         * For an ideal gas mixture,
         * \f[
         * \hat s(T, P) = \sum_k X_k \hat s^0_k(T) - \hat R \log (P/P^0).
         * \f]
         * The standard-state pure-species entropies 
         * \f$ \hat s^0_k(T) \f$ are computed by the species thermodynamic 
         * property manager.
         * @see SpeciesThermo
         */
        virtual doublereal entropy_mole() const {
            return GasConstant * (mean_X(entropy_R().begin()) -
                sum_xlogx() - log(pressure()/m_spthermo->refPressure()));
        }


        virtual doublereal gibbs_mole() const {
            return enthalpy_mole() - temperature() * entropy_mole();
        }


        /**
         * Molar heat capacity at constant pressure. Units: J/kmol/K.
         * For an ideal gas mixture, 
         * \f[
         * \hat c_p(t) = \sum_k \hat c^0_{p,k}(T).
         * \f]
         * The standard-state pure-species heat capacities  
         * \f$ \hat c^0_{p,k}(T) \f$ are computed by the species thermodynamic 
         * property manager.
         * @see SpeciesThermo
         */
        virtual doublereal cp_mole() const {
            return GasConstant * mean_X(cp_R().begin());
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
         * @name Mechanical Equation of State
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

        //@}


        virtual void getChemPotentials(doublereal* mu) const;

        virtual void getStandardChemPotentials(doublereal* mu0) const {
            getPureGibbs(mu0);
        }

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
         * \f$ P^0/\hat R T \f$.
         */ 
         virtual doublereal standardConcentration(int k=0) const {
            return m_p0/(GasConstant * temperature());
        }

         virtual doublereal logStandardConc(int k=0) const {
             _updateThermo();
            return m_logc0;
        }

        virtual void getPartialMolarEnthalpies(doublereal* hbar) const {
            const array_fp& _h = enthalpy_RT();
            doublereal rt = GasConstant * temperature();
            scale(_h.begin(), _h.end(), hbar, rt);
        }

        virtual void getPureGibbs(doublereal* gpure) const {
            const array_fp& gibbsrt = gibbs_RT();
            scale(gibbsrt.begin(), gibbsrt.end(), gpure, _RT());
        }

        void getEnthalpy_RT(doublereal* hrt) const {
            const array_fp& _h = enthalpy_RT();
            copy(_h.begin(), _h.end(), hrt);
        }

        void getEntropy_R(doublereal* sr) const {
            const array_fp& _s = entropy_R();
            copy(_s.begin(), _s.end(), sr);
        }

        virtual void getGibbs_RT(doublereal* grt) const {
            const array_fp& gibbsrt = gibbs_RT();
            copy(gibbsrt.begin(), gibbsrt.end(), grt);
        }

        void getCp_R(doublereal* cpr) const {
            const array_fp& _cpr = cp_R();
            copy(_cpr.begin(), _cpr.end(), cpr);
        }


        // new methods defined here

        const array_fp& enthalpy_RT() const {
            _updateThermo();
            return m_h0_RT;
        }

        const array_fp& gibbs_RT() const {
            _updateThermo();
            return m_g0_RT;
        }

        const array_fp& expGibbs_RT() const {
            _updateThermo();
            int k;
            for (k = 0; k != m_kk; k++) m_expg0_RT[k] = exp(m_g0_RT[k]);
            return m_expg0_RT;
        }

        const array_fp& entropy_R() const {
            _updateThermo();
            return m_s0_R;
        }

        const array_fp& cp_R() const {
            _updateThermo();
            return m_cp0_R;
        }

        virtual void setPotentialEnergy(int k, doublereal pe) {
            m_pe[k] = pe;
            _updateThermo();
        }

        virtual doublereal potentialEnergy(int k) const {
            return m_pe[k];
        }

        virtual void initThermo();


    /** 
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





