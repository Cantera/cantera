/**
 *
 *  @file ConstDensityThermo.h
 *
 *  Thermo manager for incompressible substances.

 *  $Author$
 *  $Date$
 *  $Revision$
 *
 *  Copyright 2002 California Institute of Technology
 *
 */

#ifndef CT_CONSTRHOTHERMO_H
#define CT_CONSTRHOTHERMO_H

#include "ct_defs.h"
#include "mix_defs.h"
#include "ThermoPhase.h"
#include "SpeciesThermo.h"

namespace Cantera {

    /**
     * Overloads the virtual methods of class Thermo to implement the
     * incompressible equation of state.
     */
    class ConstDensityThermo : public ThermoPhase  {

    public:

        ConstDensityThermo() : m_tlast(0.0) {}

        virtual ~ConstDensityThermo() {}

        virtual int eosType() const { return cIncompressible; }

        virtual doublereal enthalpy_mole() const {
            doublereal p0 = m_spthermo->refPressure();
            return GasConstant * temperature() * 
                mean_X(enthalpy_RT().begin()) 
                + (pressure() - p0)/molarDensity();
        }

        virtual doublereal intEnergy_mole() const {
            doublereal p0 = m_spthermo->refPressure();
            return GasConstant * temperature() * 
                mean_X(enthalpy_RT().begin()) 
                - p0/molarDensity();
        }

        virtual doublereal entropy_mole() const {
            return GasConstant * (mean_X(entropy_R().begin()) -
                sum_xlogx());
        }

        virtual doublereal gibbs_mole() const {
            return enthalpy_mole() - temperature() * entropy_mole();
        }

        virtual doublereal cp_mole() const {
            return GasConstant * mean_X(cp_R().begin());
        }

        virtual doublereal cv_mole() const {
            return cp_mole();
        }

        virtual doublereal pressure() const {
            return m_press;
        }

        virtual void setPressure(doublereal p) {
            m_press = p;
        }

        virtual void getActivityConcentrations(doublereal* c) const;

        virtual void getChemPotentials(doublereal* mu) const;
        virtual void getStandardChemPotentials(doublereal* mu0) const;
        virtual doublereal standardConcentration(int k=0) const;
        virtual doublereal logStandardConc(int k=0) const;

//         virtual void getPartialMolarEnthalpies(doublereal* hbar) const {
//             const array_fp& _h = enthalpy_RT();
//             doublereal rt = GasConstant * _temp();
//             scale(_h.begin(), _h.end(), hbar, rt);
//         }

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
        }

        virtual doublereal potentialEnergy(int k) const {
            return m_pe[k];
        }

        virtual void initThermo();


        virtual void setToEquilState(const doublereal* lambda_RT);

        // set the density
        virtual void setParameters(int n, doublereal* c) {
            setDensity(c[0]);
        }


    protected:

        int m_mm;
        doublereal m_tmin, m_tmax, m_p0;

        mutable doublereal     m_tlast;
        mutable array_fp      m_h0_RT;
        mutable array_fp      m_cp0_R;
        mutable array_fp      m_g0_RT;
        mutable array_fp      m_s0_R;
        mutable array_fp      m_expg0_RT;
        mutable array_fp      m_pe;
        mutable array_fp      m_pp;

        doublereal m_press;


    private:

        void _updateThermo() const;
    };
}
        
#endif





