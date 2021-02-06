/**
 *  @file ResidJacEval.h
 *
 * Dense, Square (not sparse) matrices.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_RESIDJACEVAL_H
#define CT_RESIDJACEVAL_H

#include "ResidEval.h"
#include "DenseMatrix.h"

namespace Cantera
{

//! Differentiates the type of residual evaluations according to functionality
enum ResidEval_Type_Enum {
    //! Base residual calculation for the time-stepping function
    Base_ResidEval = 0,
    //! Base residual calculation for the Jacobian calculation
    JacBase_ResidEval,
    //! Delta residual calculation for the Jacbobian calculation
    JacDelta_ResidEval,
    //! Base residual calculation for the showSolution routine
    /*!
     *    We calculate this when we want to display a solution
     */
    Base_ShowSolution,
    //! Base residual calculation containing any lagged components
    /*!
     * We use this to calculate residuals when doing line searches along
     * irections determined by Jacobians that are missing contributions from
     * lagged entries.
     */
    Base_LaggedSolutionComponents
};

//! Wrappers for the function evaluators for Nonlinear solvers and Time steppers
/*!
 * @attention This class currently does not have any test cases or examples. Its
 *     implementation may be incomplete, and future changes to Cantera may
 *     unexpectedly cause this class to stop working. If you use this class,
 *     please consider contributing examples or test cases. In the absence of
 *     new tests or examples, this class may be deprecated and removed in a
 *     future version of Cantera. See
 *     https://github.com/Cantera/cantera/issues/267 for additional information.
 *
 * A class for full (non-sparse dense matrices with Fortran-compatible data
 * storage. The class adds support for identifying what types of calls are made
 * to the residual evaluator by adding the ResidEval_Type_Enum class.
 */
class ResidJacEval : public ResidEval
{
public:
    //!Default constructor
    /*!
     * @param atol   Initial value of the global tolerance (defaults to 1.0E-13)
     */
    ResidJacEval(doublereal atol = 1.0e-13);

    //! Return the number of equations in the equation system
    virtual int nEquations() const;

    //! Evaluate the residual function
    /*!
     * @param t             Time                    (input)
     * @param delta_t       The current value of the time step (input)
     * @param y             Solution vector (input, do not modify)
     * @param ydot          Rate of change of solution vector. (input, do not modify)
     * @param resid         Value of the residual that is computed (output)
     * @param evalType      Type of the residual being computed (defaults to Base_ResidEval)
     * @param id_x          Index of the variable that is being numerically differenced to find
     *                      the Jacobian (defaults to -1, which indicates that no variable is being
     *                      differenced or that the residual doesn't take this issue into account)
     * @param delta_x       Value of the delta used in the numerical differencing
     *
     * @returns a flag to indicate that operation is successful.
     *            1  Means a successful operation
     *           -0 or neg value Means an unsuccessful operation
     */
    virtual int evalResidNJ(const doublereal t, const doublereal delta_t,
                            const doublereal* const y,
                            const doublereal* const ydot,
                            doublereal* const resid,
                            const ResidEval_Type_Enum evalType = Base_ResidEval,
                            const int id_x = -1,
                            const doublereal delta_x = 0.0);

    virtual int eval(const doublereal t, const doublereal* const y,
                     const doublereal* const ydot,
                     doublereal* const r);

    virtual int getInitialConditions(const doublereal t0, doublereal* const y, doublereal* const ydot);

    //! Filter the solution predictions
    /*!
     * Codes might provide a predicted step change. This routine filters the
     * predicted solution vector eliminating illegal directions.
     *
     * @param t             Time                    (input)
     * @param ybase         Solution vector (input, output)
     * @param step          Proposed step in the solution that will be cropped
     * @returns             the norm of the amount of filtering
     */
    virtual doublereal filterNewStep(const doublereal t, const doublereal* const ybase,
                                     doublereal* const step);

    //! Filter the solution predictions
    /*!
     * Codes might provide a predicted solution vector. This routine filters the
     * predicted solution vector.
     *
     * @param t             Time                    (input)
     * @param y             Solution vector (input, output)
     * @returns the norm of the amount of filtering
     */
    virtual doublereal filterSolnPrediction(const doublereal t, doublereal* const y);

    //! Set a global value of the absolute tolerance
    /*!
     *  @param atol   Value of atol
     */
    void setAtol(doublereal atol);

    //! Evaluate the time tracking equations, if any
    /*!
     * Evaluate time integrated quantities that are calculated at the end of
     * every successful time step.  This call is made once at the end of every
     * successful time step that advances the time. It's also made once at the
     * start of the time stepping.
     *
     * @param t             Time                    (input)
     * @param delta_t       The current value of the time step (input)
     * @param y             Solution vector (input, do not modify)
     * @param ydot          Rate of change of solution vector. (input, do not modify)
     * @returns a flag to indicate that operation is successful.
     *            1  Means a successful operation
     *           -0 or neg value Means an unsuccessful operation
     */
    virtual int evalTimeTrackingEqns(const doublereal t, const doublereal delta_t, const doublereal* const y,
                                     const doublereal* const ydot);

    //! Evaluate any stopping criteria other than a final time limit
    /*!
     * If we are to stop the time integration for any reason other than reaching
     * a final time limit, tout, provide a test here. This call is made at the
     * end of every successful time step iteration
     *
     * @return    If true, the the time stepping is stopped. If false, then time
     *            stepping is stopped if t >= tout Defaults to false.
     * @param t             Time                    (input)
     * @param delta_t       The current value of the time step (input)
     * @param y             Solution vector (input, do not modify)
     * @param ydot          Rate of change of solution vector. (input, do not modify)
     */
    virtual bool evalStoppingCritera(const doublereal t,
                                     const doublereal delta_t,
                                     const doublereal* const y,
                                     const doublereal* const ydot);

    //! Return a vector of delta y's for calculation of the numerical Jacobian
    /*!
     * There is a default algorithm provided.
     *
     *     delta_y[i] = atol[i] + 1.0E-6 ysoln[i]
     *     delta_y[i] = atol[i] + MAX(1.0E-6 ysoln[i] * 0.01 * solnWeights[i])
     *
     * @param t             Time                    (input)
     * @param y             Solution vector (input, do not modify)
     * @param ydot          Rate of change of solution vector. (input, do not modify)
     * @param delta_y       Value of the delta to be used in calculating the numerical Jacobian
     * @param solnWeights   Value of the solution weights that are used in determining convergence (default = 0)
     * @returns a flag to indicate that operation is successful.
     *            1  Means a successful operation
     *           -0 or neg value Means an unsuccessful operation
     */
    virtual int
    calcDeltaSolnVariables(const doublereal t,
                           const doublereal* const y,
                           const doublereal* const ydot,
                           doublereal* const delta_y,
                           const doublereal* const solnWeights = 0);

    //! Returns a vector of column scale factors that can be used to column
    //! scale Jacobians.
    /*!
     * Default to yScales[] = 1.0
     *
     * @param t             Time                    (input)
     * @param y             Solution vector (input, do not modify)
     * @param y_old         Old Solution vector (input, do not modify)
     * @param yScales       Value of the column scales
     */
    virtual void calcSolnScales(const doublereal t, const doublereal* const y,
                                const doublereal* const y_old, doublereal* const yScales);

    //! This function may be used to create output at various points in the
    //! execution of an application.
    /*!
     *  @param ifunc     identity of the call
     *                          0  Initial call
     *                          1  Called at the end of every successful time step
     *                         -1  Called at the end of every unsuccessful time step
     *                          2  Called at the end of every call to integrateRJE()
     * @param t             Time                    (input)
     * @param delta_t       The current value of the time step (input)
     * @param y             Solution vector (input, do not modify)
     * @param ydot          Rate of change of solution vector. (input)
     */
    virtual void user_out2(const int ifunc, const doublereal t,
                           const doublereal delta_t,
                           const doublereal* const y,
                           const doublereal* const ydot);

    //! This function may be used to create output at various points in the
    //! execution of an application.
    /*!
     *  This routine calls user_out2().
     *
     *  @param ifunc     identity of the call
     * @param t             Time                    (input)
     * @param y             Solution vector (input, do not modify)
     * @param ydot          Rate of change of solution vector. (input)
     */
    virtual void user_out(const int ifunc, const doublereal t,
                          const doublereal* y,
                          const doublereal* ydot);

    //! Multiply the matrix by another matrix that leads to better conditioning
    /*!
     * Provide a left sided matrix that will multiply the current Jacobian,
     * after scaling and lead to a better conditioned system. This routine is
     * called just before the matrix is factored.
     *
     * Original Problem:
     *        J delta_x = - Resid
     *
     * New problem:
     *      M (J delta_x) = - M Resid
     *
     * @param    matrix     Pointer to the current Jacobian (if zero, it's already been factored)
     * @param    nrows      offsets for the matrix
     * @param    rhs        residual vector. This also needs to be LHS multiplied by M
     * @returns a flag to indicate that operation is successful.
     *            1  Means a successful operation
     *           -0 or neg value Means an unsuccessful operation
     */
    virtual int matrixConditioning(doublereal* const matrix, const int nrows,
                                   doublereal* const rhs);

    //! Calculate an analytical Jacobian and the residual at the current time
    //! and values.
    /*!
     *  Only called if the jacFormation method is set to analytical
     *
     * @param t             Time                    (input)
     * @param delta_t       The current value of the time step (input)
     * @param cj            Coefficient of yprime used in the evaluation of the Jacobian
     * @param y             Solution vector (input, do not modify)
     * @param ydot          Rate of change of solution vector. (input, do not modify)
     * @param J             Reference to the DenseMatrix object to be calculated (output)
     * @param resid         Value of the residual that is computed (output)
     * @returns a flag to indicate that operation is successful.
     *            1  Means a successful operation
     *           -0 or neg value Means an unsuccessful operation
     */
    virtual int evalJacobian(const doublereal t, const doublereal delta_t, doublereal cj,
                             const doublereal* const y, const doublereal* const ydot,
                             DenseMatrix& J, doublereal* const resid);

    //! Calculate an analytical Jacobian and the residual at the current time and values.
    /*!
     *  Only called if the jacFormation method is set to analytical
     *
     * @param t             Time                    (input)
     * @param delta_t       The current value of the time step (input)
     * @param cj            Coefficient of yprime used in the evaluation of the Jacobian
     * @param y             Solution vector (input, do not modify)
     * @param ydot          Rate of change of solution vector. (input, do not modify)
     * @param jacobianColPts   Pointer to the vector of pts to columns of the DenseMatrix
     *                         object to be calculated (output)
     * @param resid         Value of the residual that is computed (output)
     * @returns a flag to indicate that operation is successful.
     *            1  Means a successful operation
     *           -0 or neg value Means an unsuccessful operation
     */
    virtual int evalJacobianDP(const doublereal t, const doublereal delta_t, doublereal cj,
                               const doublereal* const y,
                               const doublereal* const ydot,
                               doublereal* const* jacobianColPts,
                               doublereal* const resid);

protected:
    //! constant value of atol
    doublereal m_atol;

    //! Number of equations
    int neq_;
};
}

#endif
