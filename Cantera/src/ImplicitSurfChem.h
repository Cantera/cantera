/**
 *  @file ImplicitSurfChem.h
 *
 * Implicit integration of surface site density equations.
 *
 * $Author$
 * $Revision$
 * $Date$
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_IMPSURFCHEM_H
#define CT_IMPSURFCHEM_H

#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "FuncEval.h"
#include "CVode.h"
#include "InterfaceKinetics.h"
#include "SurfPhase.h"

namespace Cantera {

    /**
     * Advances the surface coverages of an associated SurfacePhase
     * object in time by implicitly integrating \f[ \dot \theta_k =
     * \dot s_k (\sigma_k / s_0)\f]
     */
    class ImplicitSurfChem : public FuncEval {

    public:

        /**
         * Constructor.
         */
        ImplicitSurfChem(InterfaceKinetics& kin);


        /**
         * Destructor. Deletes the integrator.
         */
        virtual ~ImplicitSurfChem(){ delete m_integ; }
        

        /**
         * Overloads the virtual function
         * declared in FuncEval. 
         */
        virtual void initialize(doublereal t0 = 0.0);


        /**
         * Integrate from t0 to t1. The integrator is reinitialized
         * first.
         */
        void integrate(doublereal t0, doublereal t1) {
            m_integ->reinitialize(t0, *this);
            m_integ->setMaxStep(t1 - t0);
            m_integ->integrate(t1);
            updateState(m_integ->solution());
        }

        /**
         * Integrate from t0 to t1 without reinitializing the
         * integrator. Use when the coverages have not changed from
         * their values on return from the last call to integrate or
         * integrate0.
         */
        void integrate0(doublereal t0, doublereal t1) {
            m_integ->integrate(t1);
            updateState(m_integ->solution());
        }

        // overloaded methods of class FuncEval
        virtual int neq() { return m_nsp; }
	virtual void eval(doublereal t, doublereal* y, doublereal* ydot);
        virtual void getInitialConditions(doublereal t0, 
            size_t leny, doublereal* y);


    protected:
        
        /**
         * Set the mixture to a state consistent with solution
         * vector y.
         */
        void updateState(doublereal* y);

        SurfPhase* m_surf;
        InterfaceKinetics* m_kin;
        int m_nsp;
        Integrator* m_integ;         // pointer to integrator
        doublereal m_atol, m_rtol;   // tolerances
        doublereal m_maxstep;        // max step size
        vector_fp m_work;

    private:

    };
}

#endif

