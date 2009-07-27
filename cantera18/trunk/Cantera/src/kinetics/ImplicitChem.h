/**
 *  @file ImplicitChem.h
 *
 * $Author: dggoodwin $
 * $Revision: 1.1 $
 * $Date: 2007/05/04 14:27:23 $
 */

// Copyright 2001  California Institute of Technology

#ifndef CT_IMPCHEM_H
#define CT_IMPCHEM_H

#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "FuncEval.h"
#include "Integrator.h"
#include "Kinetics.h"
#include "ThermoPhase.h"

namespace Cantera {

    /**
     * Advances the composition of an associated phase object in time
     * by implicitly integrating
     * \f[
     * \dot Y_k = \frac{\omega_k}{\rho}
     * \f]
     */
    class ImplicitChem : public FuncEval {

    public:

        /**
         * Constructor.
         */
        ImplicitChem(Kinetics& kin, ThermoPhase& therm);


        /**
         * Destructor. Deletes the integrator.
         */
        virtual ~ImplicitChem(){ delete m_integ; }
        

        /**
         * Overloads the virtual function
         * declared in FuncEval. 
         */
        virtual void initialize(doublereal t0 = 0.0);

        void adiabatic() {
            m_energy = true;
        }

        void isothermal() {
            m_energy = false;
        }

        /**
         * Integrate from t0 to t1. The integrator is reinitialized
         * first.
         */
        void integrate(doublereal t0, doublereal t1) {
            m_integ->reinitialize(t0, *this);
            m_integ->setMaxStepSize(t1 - t0);
            m_rho = m_thermo->density();
            m_integ->integrate(t1);
            updateState(m_integ->solution());
        }

        /**
         * Integrate from t0 to t1 without reinitializing the
         * integrator.
         */
        void integrate0(doublereal t0, doublereal t1) {
            m_integ->integrate(t1);
            updateState(m_integ->solution());
        }

        // overloaded methods of class FuncEval
        virtual int neq() { return m_nsp; }
	virtual void eval(doublereal t, doublereal* y, doublereal* ydot,
            doublereal* p);
        virtual void getInitialConditions(doublereal t0, size_t leny, 
            doublereal* y);


    protected:
        
        /**
         * Set the mixture to a state consistent with solution
         * vector y.
         */
        void updateState(doublereal* y);

        //Kinetics::phase_t* m_mix;
        Kinetics* m_kin;
        ThermoPhase*   m_thermo;
        int m_nsp;
        Integrator* m_integ;         // pointer to integrator
        doublereal m_atol, m_rtol;   // tolerances
        doublereal m_maxstep;        // max step size
        array_fp m_wt;
        doublereal m_rho;
        bool m_energy;
        doublereal m_h0;
        doublereal m_press;

    private:

    };
}

#endif
