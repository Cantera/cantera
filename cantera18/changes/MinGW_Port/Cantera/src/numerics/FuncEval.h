/**
 *  @file FuncEval.h
 *
 */

// Copyright 2001  California Institute of Technology

#ifndef CT_FUNCEVAL_H
#define CT_FUNCEVAL_H

#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "ct_defs.h"

namespace Cantera {


    /**
     *  Virtual base class for ODE right-hand-side function evaluators.
     *  Classes derived from FuncEval evaluate the right-hand-side function
     * \f$ \vec{F}(t,\vec{y})\f$ in 
     * \f[
     *  \dot{\vec{y}} = \vec{F}(t,\vec{y}).
     * \f] 
     *  @ingroup odeGroup 
     */
    class FuncEval {

    public:

        FuncEval() {}
        virtual ~FuncEval() {}

        /**
         * Evaluate the right-hand-side function. Called by the
         * integrator.
         * @param t time. (input)
         * @param y solution vector. (input)
         * @param ydot rate of change of solution vector. (output)
         * @param p parameter vector
         */
	virtual void eval(double t, double* y, double* ydot, double* p)=0;

        /**
         * Fill the solution vector with the initial conditions 
         * at initial time t0.
         */
        virtual void getInitialConditions(double t0, size_t leny, double* y)=0;

        /** 
         * Number of equations. 
         */
        virtual int neq()=0;

        /// Number of parameters.
        virtual int nparams() { return 0; }

    protected:

    private:

    };

}

#endif
