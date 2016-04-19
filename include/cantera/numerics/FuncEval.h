/**
 *  @file FuncEval.h
 */

// Copyright 2001  California Institute of Technology

#ifndef CT_FUNCEVAL_H
#define CT_FUNCEVAL_H

#include "cantera/base/ct_defs.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/global.h"

namespace Cantera
{
/**
 *  Virtual base class for ODE right-hand-side function evaluators.
 *  Classes derived from FuncEval evaluate the right-hand-side function
 * \f$ \vec{F}(t,\vec{y})\f$ in
 * \f[
 *  \dot{\vec{y}} = \vec{F}(t,\vec{y}).
 * \f]
 *  @ingroup odeGroup
 */
class FuncEval
{
public:
    FuncEval() {}
    virtual ~FuncEval() {}

    /**
     * Evaluate the right-hand-side function. Called by the integrator.
     * @param[in] t time.
     * @param[in] y solution vector, length neq()
     * @param[out] ydot rate of change of solution vector, length neq()
     * @param[in] p sensitivity parameter vector, length nparams()
     */
    virtual void eval(double t, double* y, double* ydot, double* p)=0;

    /**
     * Fill the solution vector with the initial conditions
     * at initial time t0.
     * @deprecated Use getState() instead. To be removed after Cantera 2.3.
     */
    virtual void getInitialConditions(double t0, size_t leny, double* y) {
        warn_deprecated("FuncEval::getInitialConditions",
            "Use getState instead. To be removed after Cantera 2.3.");
        getState(y);
    }

    //! Fill in the vector *y* with the current state of the system
    virtual void getState(double* y) {
        throw NotImplementedError("FuncEval::getState");
    }

    //! Number of equations.
    virtual size_t neq()=0;

    //! Number of sensitivity parameters.
    virtual size_t nparams() {
        return m_sens_params.size();
    }

    //! Values for the problem parameters for which sensitivities are computed
    //! This is the array which is perturbed and passed back as the fourth
    //! argument to eval().
    vector_fp m_sens_params;
};

}

#endif
