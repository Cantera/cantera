/**
 *  @file FuncEval.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

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
    FuncEval();
    virtual ~FuncEval() {}

    /**
     * Evaluate the right-hand-side function. Called by the integrator.
     * @param[in] t time.
     * @param[in] y solution vector, length neq()
     * @param[out] ydot rate of change of solution vector, length neq()
     * @param[in] p sensitivity parameter vector, length nparams()
     */
    virtual void eval(double t, double* y, double* ydot, double* p)=0;

    //! Evaluate the right-hand side using return code to indicate status.
    /*!
     * Errors are indicated using the return value, rather than by throwing
     *  exceptions. This method is used when calling from a C-based integrator
     *  such as CVODES. Exceptions may either be stored or printed, based on the
     *  setting of suppressErrors().
     *  @returns 0 for a successful evaluation; 1 after a potentially-
     *      recoverable error; -1 after an unrecoverable error.
     */
    int eval_nothrow(double t, double* y, double* ydot);

    /*! Evaluate the setup processes for the Jacobian preconditioner.
     * @param[in] t time.
     * @param[in] y solution vector, length neq()
     * @param gamma the gamma in M=I-gamma*J
     * @warning This function is an experimental part of the %Cantera API and may be
     * changed or removed without notice.
     */
    virtual void preconditionerSetup(double t, double* y, double gamma) {
        throw NotImplementedError("FuncEval::preconditionerSetup");
    }

    /*! Evaluate the linear system Ax=b where A is the preconditioner.
     * @param[in] rhs right hand side vector used in linear system
     * @param[out] output output vector for solution
     * @warning This function is an experimental part of the %Cantera API and may be
     * changed or removed without notice.
     */
    virtual void preconditionerSolve(double* rhs, double* output) {
        throw NotImplementedError("FuncEval::preconditionerSolve");
    }

    /*! Preconditioner setup that doesn't throw an error but returns a
     * CVODES flag. It also helps as a first level of polymorphism
     * which identifies the specific FuncEval, e.g., ReactorNet.
     * @param[in] t time.
     * @param[in] y solution vector, length neq()
     * @param gamma the gamma in M=I-gamma*J
     * @warning This function is an experimental part of the %Cantera API and may be
     * changed or removed without notice.
     */
    int preconditioner_setup_nothrow(double t, double* y, double gamma);

    /*! Preconditioner solve that doesn't throw an error but returns a
     * CVODES flag. It also helps as a first level of polymorphism
     * which identifies the specific FuncEval, e.g., ReactorNet.
     * @param[in] rhs right hand side vector used in linear system
     * @param[out] output output vector for solution
     * @warning This function is an experimental part of the %Cantera API and may be
     * changed or removed without notice.
     */
    int preconditioner_solve_nothrow(double* rhs, double* output);

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

    //! Enable or disable suppression of errors when calling eval()
    void suppressErrors(bool suppress) {
        m_suppress_errors = suppress;
    }

    //! Get current state of error suppression
    bool suppressErrors() const {
        return m_suppress_errors;
    };

    //! Return a string containing the text of any suppressed errors
    std::string getErrors() const;

    //! Clear any previously-stored suppressed errors
    void clearErrors() {
        m_errors.clear();
    };

    //! Values for the problem parameters for which sensitivities are computed
    //! This is the array which is perturbed and passed back as the fourth
    //! argument to eval().
    vector_fp m_sens_params;

    //! Scaling factors for each sensitivity parameter
    vector_fp m_paramScales;

protected:
    // If true, errors are accumulated in m_errors. Otherwise, they are printed
    bool m_suppress_errors;

    //! Errors occurring during function evaluations
    std::vector<std::string> m_errors;
};

}

#endif
