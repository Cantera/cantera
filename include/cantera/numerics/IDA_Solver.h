/**
 *  @file IDA_Solver.h
 *  Header file for class IDA_Solver
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_IDA_SOLVER_H
#define CT_IDA_SOLVER_H

#include "DAE_Solver.h"

#include "sundials/sundials_nvector.h"

// These constants are defined internally in the IDA package, ida.c
#define IDA_NN 0
#define IDA_SS 1
#define IDA_SV 2
#define IDA_WF 3

#define REAL_WORKSPACE_SIZE 0

namespace Cantera
{

class ResidData;

/**
 * Wrapper for Sundials IDA solver
 *
 * @attention This class currently does not have any test cases or examples. Its
 *     implementation may be incomplete, and future changes to Cantera may
 *     unexpectedly cause this class to stop working. If you use this class,
 *     please consider contributing examples or test cases. In the absence of
 *     new tests or examples, this class may be deprecated and removed in a
 *     future version of Cantera. See
 *     https://github.com/Cantera/cantera/issues/267 for additional information.
 */
class IDA_Solver : public DAE_Solver
{
public:
    //! Constructor.
    /*!
     * Default settings: dense Jacobian, no user-supplied Jacobian function,
     * Newton iteration.
     *
     * @param f  Function that will supply the time dependent residual to be solved
     */
    IDA_Solver(ResidJacEval& f);

    virtual ~IDA_Solver();

    virtual void setTolerances(doublereal reltol,
                               doublereal* abstol);

    virtual void setTolerances(doublereal reltol, doublereal abstol);

    virtual void setLinearSolverType(int solverType);

    //! Set up the problem to use a dense linear direct solver
    virtual void setDenseLinearSolver();

    //! Set up the problem to use a band solver
    /*!
     * @param m_upper   upper band width of the matrix
     * @param m_lower   lower band width of the matrix
     */
    virtual void setBandedLinearSolver(int m_upper, int m_lower);

    virtual void setMaxOrder(int n);

    //! Set the maximum number of time steps
    /*!
     * @param n  input of maximum number of time steps
     */
    virtual void setMaxNumSteps(int n);

    //! Set the initial step size
    /*!
     * @param h0  initial step size value
     */
    virtual void setInitialStepSize(doublereal h0);

    //! Set the stop time
    /*!
     * @param tstop the independent variable value past which the solution is
     *     not to proceed.
     */
    virtual void setStopTime(doublereal tstop);

    //! Get the current step size from IDA via a call
    /*!
     * @returns the current step size.
     */
    virtual double getCurrentStepFromIDA();

    //! Set the form of the Jacobian
    /*!
     * @param formJac  Form of the Jacobian
     *                 0 numerical Jacobian
     *                 1 analytical Jacobian given by the evalJacobianDP() function
     */
    virtual void setJacobianType(int formJac);


    virtual void setMaxErrTestFailures(int n);

    //! Set the maximum number of nonlinear iterations on a timestep
    /*!
     * @param n  Set the max iterations. The default is 4, which seems awfully low to me.
     */
    virtual void setMaxNonlinIterations(int n);

    //! Set the maximum number of nonlinear solver convergence failures
    /*!
     * @param n  Value of nonlin failures. If value is exceeded, the calculation terminates.
     */
    virtual void setMaxNonlinConvFailures(int n);


    virtual void inclAlgebraicInErrorTest(bool yesno);

    /**
     * Get the value of a solver-specific output parameter.
     */
    virtual doublereal getOutputParameter(int flag) const;

    virtual void correctInitial_Y_given_Yp(doublereal* y, doublereal* yp,
                                           doublereal tout);

    virtual void correctInitial_YaYp_given_Yd(doublereal* y, doublereal* yp, doublereal tout);

    //! Step the system to a final value of the time
    /*!
     * @param tout  Final value of the time
     * @returns the IDASolve() return flag
     *
     * The return values for IDASolve are described below. (The numerical return
     * values are defined above in this file.) All unsuccessful returns give a
     * negative return value.
     *
     * IDA_SUCCESS
     *   IDASolve succeeded and no roots were found.
     *
     * IDA_ROOT_RETURN:  IDASolve succeeded, and found one or more roots.
     *   If nrtfn > 1, call IDAGetRootInfo to see which g_i were found
     *   to have a root at (*tret).
     *
     * IDA_TSTOP_RETURN:
     *   IDASolve returns computed results for the independent variable
     *   value tstop. That is, tstop was reached.
     *
     * IDA_MEM_NULL:
     *   The IDA_mem argument was NULL.
     *
     * IDA_ILL_INPUT:
     *   One of the inputs to IDASolve is illegal. This includes the
     *   situation when a component of the error weight vectors
     *   becomes < 0 during internal stepping.  It also includes the
     *   situation where a root of one of the root functions was found
     *   both at t0 and very near t0.  The ILL_INPUT flag
     *   will also be returned if the linear solver function IDA---
     *   (called by the user after calling IDACreate) failed to set one
     *   of the linear solver-related fields in ida_mem or if the linear
     *   solver's init routine failed. In any case, the user should see
     *   the printed error message for more details.
     *
     * IDA_TOO_MUCH_WORK:
     *   The solver took mxstep internal steps but could not reach tout.
     *   The default value for mxstep is MXSTEP_DEFAULT = 500.
     *
     * IDA_TOO_MUCH_ACC:
     *   The solver could not satisfy the accuracy demanded by the user
     *   for some internal step.
     *
     * IDA_ERR_FAIL:
     *   Error test failures occurred too many times (=MXETF = 10) during
     *   one internal step.
     *
     * IDA_CONV_FAIL:
     *   Convergence test failures occurred too many times (= MXNCF = 10)
     *   during one internal step.
     *
     * IDA_LSETUP_FAIL:
     *   The linear solver's setup routine failed
     *   in an unrecoverable manner.
     *
     * IDA_LSOLVE_FAIL:
     *   The linear solver's solve routine failed
     *   in an unrecoverable manner.
     *
     * IDA_CONSTR_FAIL:
     *    The inequality constraints were violated,
     *    and the solver was unable to recover.
     *
     * IDA_REP_RES_ERR:
     *    The user's residual function repeatedly returned a recoverable
     *    error flag, but the solver was unable to recover.
     *
     * IDA_RES_FAIL:
     *    The user's residual function returned a nonrecoverable error
     *    flag.
     */
    virtual int solve(doublereal tout);

    virtual doublereal step(doublereal tout);

    virtual void init(doublereal t0);

    virtual doublereal solution(int k) const;

    virtual const doublereal* solutionVector() const;

    virtual doublereal derivative(int k) const;

    virtual const doublereal* derivativeVector() const;

    void* IDAMemory() {
        return m_ida_mem;
    }

protected:
    //! Pointer to the IDA memory for the problem
    void* m_ida_mem;
    void* m_linsol; //!< Sundials linear solver object
    void* m_linsol_matrix; //!< matrix used by Sundials

    //! Initial value of the time
    doublereal m_t0;

    //! Current value of the solution vector
    N_Vector m_y;

    //! Current value of the derivative of the solution vector
    N_Vector m_ydot;
    N_Vector m_id;
    N_Vector m_constraints;
    N_Vector m_abstol;
    int m_type;

    int m_itol;
    int m_iter;
    doublereal m_reltol;
    doublereal m_abstols;
    int m_nabs;

    //! Maximum value of the timestep allowed
    doublereal m_hmax;

    //! Minimum value of the timestep allowed
    doublereal m_hmin;

    //! Value of the initial time step
    doublereal m_h0;

    //! Maximum number of time steps allowed
    int m_maxsteps;

    //!  maximum time step order of the method
    int m_maxord;

    //! Form of the Jacobian
    /*!
     *  0 numerical Jacobian created by IDA
     *  1 analytical Jacobian. Must have populated the evalJacobianDP()
     *    function in the ResidJacEval class.
     *  2 numerical Jacobian formed by the ResidJacEval class (unimplemented)
     */
    int m_formJac;

    //! maximum time
    doublereal m_tstop;

    //! Value of the previous, previous time
    doublereal m_told_old;

    //! Value of the previous time
    doublereal m_told;

    //! Value of the current time
    doublereal m_tcurrent;

    //! Value of deltaT for the current step
    doublereal m_deltat;

    //! maximum number of error test failures
    int m_maxErrTestFails;

    //! Maximum number of nonlinear solver iterations at one solution
    /*!
     *  If zero, this is the default of 4.
     */
    int m_maxNonlinIters;

    //! Maximum number of nonlinear convergence failures
    int m_maxNonlinConvFails;

    //! If true, the algebraic variables don't contribute to error tolerances
    int m_setSuppressAlg;

    std::unique_ptr<ResidData> m_fdata;
    int m_mupper;
    int m_mlower;
};

}

#endif
