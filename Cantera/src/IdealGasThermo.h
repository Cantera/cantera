/**
 *
 *  @file IdealGasThermo.h
 *
 *  Template for an equation of state class that implements the ideal
 *  gas equation.
 */

/*  $Author$
 *  $Date$
 *  $Revision$
 *
 *  Copyright 2001 California Institute of Technology
 *
 */


#ifndef CT_IDEALGASTHERMO_H
#define CT_IDEALGASTHERMO_H

#include "ct_defs.h"
#include "mix_defs.h"
#include "Thermo.h"
#include "SpeciesThermo.h"

namespace Cantera {

    /**
     * Overloads the virtual methods of class Thermo to implement the
     * ideal gas equation of state.
     */
    class IdealGasThermo : public Thermo  {

    public:

        IdealGasThermo(phase_t* phase=0, SpeciesThermo* sptherm = 0) 
            : Thermo(phase, sptherm), m_tlast(0.0) {}

        virtual ~IdealGasThermo() {}

        virtual int eosType() const { return cIdealGas; }

        virtual doublereal enthalpy_mole() const {
            return GasConstant * m_s->temperature() * 
                m_s->mean_X(enthalpy_RT().begin());
        }

        virtual doublereal intEnergy_mole() const {
            return GasConstant * m_s->temperature()
                * ( m_s->mean_X(enthalpy_RT().begin()) - 1.0);
        }

        virtual doublereal entropy_mole() const {
            return GasConstant * (m_s->mean_X(entropy_R().begin()) -
                m_s->sum_xlogx() - log(pressure()/m_spthermo->refPressure()));
        }

        virtual doublereal gibbs_mole() const {
            return enthalpy_mole() - m_s->temperature() * entropy_mole();
        }

        virtual doublereal cp_mole() const {
            return GasConstant * m_s->mean_X(cp_R().begin());
        }

        virtual doublereal cv_mole() const {
            return cp_mole() - GasConstant;
        }

        virtual doublereal pressure() const {
            return GasConstant * m_s->molarDensity() * m_s->temperature();
        }

        virtual void setPressure(doublereal p) {
            m_s->setDensity(p * m_s->meanMolecularWeight()
                /(GasConstant * m_s->temperature()));
        }

        virtual void getChemPotentials(doublereal* mu) const;

        virtual void getPartialMolarEnthalpies(doublereal* hbar) const {
            const array_fp& _h = enthalpy_RT();
            doublereal rt = GasConstant * _temp();
            scale(_h.begin(), _h.end(), hbar, rt);
        }

        //virtual void getPartialMolarEntropies(doublereal* sbar) const {
        //    err("getPartialMolarEntropies");
        //}

        //virtual void getPartialMolarVolumes(doublereal* vbar) const {
        //    err("getPartialMolarVolumes");
        //}

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

        virtual doublereal minTemp(int k = -1) {
            return m_spthermo->minTemp(k);
        }

        virtual doublereal maxTemp(int k = -1) {
            return m_spthermo->maxTemp(k);            
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

        void setPotentialEnergy(int k, doublereal pe) {
            m_pe[k] = pe;
        }

        doublereal potentialEnergy(int k) {
            return m_pe[k];
        }

        virtual doublereal refPressure() const {
            return m_spthermo->refPressure();
        }

        void initThermo(Phase& s);


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

        mutable doublereal     m_tlast;
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





