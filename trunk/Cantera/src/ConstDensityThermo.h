/**
 *
 *  @file ConstDensityThermo.h
 *
 *  Thermo manager for incompressible substances.
 */

/*
 *  $Author: dggoodwin $
 *  $Date: 2005/06/18 17:01:08 $
 *  $Revision: 1.8 $
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
#include "utilities.h"

namespace Cantera {

    /**
     * Overloads the virtual methods of class ThermoPhase to implement the
     * incompressible equation of state.
     */
    class ConstDensityThermo : public ThermoPhase  {

    public:

        ConstDensityThermo() : m_tlast(0.0) {}
        virtual ~ConstDensityThermo() {}

        // overloaded methods of class ThermoPhase

        virtual int eosType() const;
        virtual doublereal enthalpy_mole() const;
        virtual doublereal intEnergy_mole() const;
        virtual doublereal entropy_mole() const;
        virtual doublereal gibbs_mole() const;
        virtual doublereal cp_mole() const;
        virtual doublereal cv_mole() const;
        virtual doublereal pressure() const;
        virtual void setPressure(doublereal p);
        virtual void getActivityConcentrations(doublereal* c) const;
	virtual void getActivityCoefficients(doublereal* ac) const;
        virtual void getChemPotentials(doublereal* mu) const;
        virtual void getStandardChemPotentials(doublereal* mu0) const;
        virtual doublereal standardConcentration(int k=0) const;
        virtual doublereal logStandardConc(int k=0) const;

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

        virtual void getParameters(int &n, doublereal * const c) {
            double d = density();
            c[0] = d;
            n = 1;
        }

        virtual void setParametersFromXML(const XML_Node& eosdata);

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
