/**
 *  @file ResidEval.h
 *
 */

// Copyright 2006  California Institute of Technology

#ifndef CT_RESIDEVAL_H
#define CT_RESIDEVAL_H

#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "ct_defs.h"

namespace Cantera {

    const int c_NONE = 0;
    const int c_GE_ZERO = 1;
    const int c_GT_ZERO = 2;
    const int c_LE_ZERO = -1;
    const int c_LT_ZERO = -2;

    /**
     *  Virtual base class for DAE residual function evaluators.
     *  Classes derived from ResidEval evaluate the residual function
     * \f[
          \vec{F}(t,\vec{y}, \vec{y^\prime})
     * \f]
     * The DAE solver attempts to find a solution y(t) such that F = 0.
     *  @ingroup DAE_Group 
     */
    class ResidEval {

    public:

        ResidEval() {}
        virtual ~ResidEval() {}

        /**
         * Constrain solution component k. Possible values for 
         * 'flag' are:
         *   - c_NONE       no constraint
         *   - c_GE_ZERO    >= 0
         *   - c_GT_ZERO    >  0
         *   - c_LE_ZERO    <= 0
         *   - c_LT_ZERO    <  0
         */
        virtual void constrain(int k, int flag) { m_constrain[k] = flag; }
        int constraint(int k) { return m_constrain[k]; }

        /** 
         * Specify that solution component k is purely algebraic -
         * that is, the derivative of this component does not appear
         * in the residual function.
         */
        virtual void setAlgebraic(int k) { m_alg[k] = 1; }
        virtual bool isAlgebraic(int k) {return (m_alg[k] == 1); }
        

        /**
         * Evaluate the residual function. Called by the
         * integrator.
         * @param t time. (input)
         * @param y solution vector. (input)
         * @param ydot rate of change of solution vector. (input)
         * @param r residual vector (output)
         */
	virtual int eval(double t, const double* y, 
            const double* ydot, double* r)=0;

        /**
         * Fill the solution and derivative vectors with the initial
         * conditions at initial time t0.  If these do not satisfy the
         * residual equation, call one of the "corrrectInitial_xxx"
         * methods before calling solve.
         */
        virtual void getInitialConditions(double t0, double* y, 
            doublereal* ydot)=0;

        /** 
         * Number of equations. 
         */
        virtual int nEquations()=0;


    protected:

        std::map<int, int> m_alg;
        std::map<int, int> m_constrain;

    private:

    };

}

#endif
