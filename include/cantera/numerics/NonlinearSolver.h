/**
 *  @file NonlinearSolver.h
 *    Class that calculates the solution to a nonlinear, dense, set
 *    of equations (see \ref numerics
 *    and class \link Cantera::NonlinearSolver NonlinearSolver\endlink).
 */

/*
 *  $Date$
 *  $Revision$
 */
/*
 * Copyright 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#ifndef CT_NONLINEARSOLVER_H
#define CT_NONLINEARSOLVER_H

#include "cantera/numerics/ResidJacEval.h"
#include "cantera/numerics/SquareMatrix.h"

namespace Cantera
{

//@{
///  @name  Constant which determines the type of the nonlinear solve
/*!
 *  I think steady state is the only option I'm gunning for
 */
//!  The nonlinear problem is part of a pseudo time dependent calculation (NOT TESTED)
#define  NSOLN_TYPE_PSEUDO_TIME_DEPENDENT 2
//! The nonlinear problem is part of a time dependent calculation
#define  NSOLN_TYPE_TIME_DEPENDENT 1
//!  The nonlinear problem is part of a steady state calculation
#define  NSOLN_TYPE_STEADY_STATE   0
//@}


//@{
///  @name  Constant which determines the Return int from the nonlinear solver
/*!
 *  This int is returned from the nonlinear solver
 */
//!  The nonlinear solve is successful.
#define  NSOLN_RETN_SUCCESS 1
//! Problem isn't solved yet
#define   NSOLN_RETN_CONTINUE  0
//! The nonlinear problem started to take too small an update step. This indicates that either the
//! Jacobian is bad, or a constraint is being bumped up against.
#define  NSOLN_RETN_FAIL_STEPTOOSMALL -1
//! The nonlinear problem didn't solve the problem
#define  NSOLN_RETN_FAIL_DAMPSTEP  -2
//!  The nonlinear problem's Jacobian is singular
#define  NSOLN_RETN_MATRIXINVERSIONERROR   -3
//!  The nonlinear problem's Jacobian formation produced an error
#define  NSOLN_RETN_JACOBIANFORMATIONERROR   -4
//!  The nonlinear problem's base residual produced an error
#define  NSOLN_RETN_RESIDUALFORMATIONERROR   -5
//!  The nonlinear problem's max number of iterations has been exceeded
#define  NSOLN_RETN_MAXIMUMITERATIONSEXCEEDED   -7
//@}
//@}

//@{
///  @name  Constant which determines the type of the Jacobian
//! The Jacobian will be calculated from a numerical method
#define NSOLN_JAC_NUM 1
//! The Jacobian is calculated from an analytical function
#define NSOLN_JAC_ANAL 2
//@}


//! Class that calculates the solution to a nonlinear system
/*!
 *  This is a small nonlinear solver that can solve highly nonlinear problems that
 *  must use a dense matrix to relax the system.
 *
 *  Newton's method is used.
 *
 *  Damping is used extensively when relaxing the system.
 *
 *   The basic idea is that we predict a direction that is parameterized by an overall coordinate
 *   value, beta, from zero to one, This may or may not be the same as the value, damp,
 *   depending upon whether the direction is straight.
 *
 * TIME STEP TYPE
 *
 *   The code solves a nonlinear problem. Frequently the nonlinear problem is created from time-dependent
 *   residual. Whenever you change the solution vector, you are also changing the derivative of the
 *   solution vector. Therefore, the code has the option of altering ydot, a vector of time derivatives
 *   of the solution in tandem with the solution vector and then feeding a residual and Jacobian routine
 *   with the time derivatives as well as the solution. The code has support for a backwards euler method
 *   and a second order  Adams-Bashforth or Trapezoidal Rule.
 *
 *   In order to use these methods, the solver must be initialized with delta_t and m_y_nm1[i] to specify
 *   the conditions at the previous time step. For second order methods, the time derivative at t_nm1 must
 *   also be supplied,  m_ydot_nm1[i]. Then the solution type  NSOLN_TYPE_TIME_DEPENDENT may be used to
 *   solve the problem.
 *
 *   For steady state problem whose residual doesn't have a solution time derivative in it, you should
 *   use the NSOLN_TYPE_STEADY_STATE problem type.
 *
 *   We have a NSOLN_TYPE_PSEUDO_TIME_DEPENDENT defined. However, this is not implemented yet. This would
 *   be a pseudo time dependent calculation, where an optional time derivative could be added in order to
 *   help equilibrate a nonlinear steady state system. The time transient is not important in and of
 *   itself. Many physical systems have a time dependence to them that provides a natural way to relax
 *   the nonlinear system.
 *
 *   MATRIX SCALING
 *
 *  @code
 *  NonlinearSolver *nls = new NonlinearSolver(&r1);
 *  int solnType =       NSOLN_TYPE_STEADY_STATE ;
 *  nls->setDeltaBoundsMagnitudes(deltaBounds);
 *  nls->solve_nonlinear_problem(solnType, y_comm, ydot_comm, CJ, time_curr, jac,
 *                                num_newt_its,  num_linear_solves, numBacktracks,
 *                               loglevelInput);
 *  @endcode
 *
 *  @ingroup numerics
 *  @deprecated Unused. To be removed after Cantera 2.2.
 */
class NonlinearSolver
{
public:
    //! Default constructor
    /*!
     * @param func   Residual and Jacobian evaluator function object
     */
    NonlinearSolver(ResidJacEval* func);

    //!Copy Constructor
    NonlinearSolver(const NonlinearSolver& right);

    //! Destructor
    virtual ~NonlinearSolver();

    //! Assignment operator
    NonlinearSolver& operator=(const NonlinearSolver& right);

    //! Create solution weights for convergence criteria
    /*!
     *  We create soln weights from the following formula
     *
     *  wt[i] = rtol * abs(y[i]) + atol[i]
     *
     *  The program always assumes that atol is specific
     *  to the solution component
     *
     * @param y  vector of the current solution values
     */
    void createSolnWeights(const doublereal* const y);

    //! L2 norm of the delta of the solution vector
    /*!
     *  calculate the norm of the solution vector. This will
     *  involve the column scaling of the matrix
     *
     *  The third argument has a default of false. However, if true, then a
     *  table of the largest values is printed out to standard output.
     *
     *  @param delta_y       Vector to take the norm of
     *  @param title         Optional title to be printed out
     *  @param printLargest  int indicating how many specific lines should be printed out
     *  @param dampFactor    Current value of the damping factor. Defaults to 1.
     *                       only used for printout out a table.
     *
     *  @return Returns the L2 norm of the delta
     */
    doublereal solnErrorNorm(const doublereal* const delta_y,  const char* title = 0, int printLargest = 0,
                             const doublereal dampFactor = 1.0) const;

    //! L2 norm of the residual of the equation system
    /*!
     * Calculate the norm of the residual vector. This may
     * involve using the row sum scaling from the matrix problem.
     *
     *  The third argument has a default of false. However, if true, then a
     *  table of the largest values is printed out to standard output.
     *
     *  @param resid    Vector of the residuals
     *  @param title    Optional title to be printed out
     *  @param printLargest  Number of specific entries to be printed
     *  @param y        Current value of y - only used for printouts
     *
     *  @return Returns the L2 norm of the delta
     */
    doublereal residErrorNorm(const doublereal* const resid, const char* title = 0, const int printLargest = 0,
                              const doublereal* const y = 0) const;

    //! Compute the current residual
    /*!
     *  The current value of the residual is stored in the internal work array
     *  m_resid, which is defined as mutable
     *
     *  @param time_curr    Value of the time
     *  @param typeCalc     Type of the calculation
     *  @param y_curr       Current value of the solution vector
     *  @param ydot_curr    Current value of the time derivative of the solution vector
     *  @param evalType     Base evaluation type. Defaults to Base_ResidEval
     *
     * @return Returns a flag to indicate that operation is successful.
     *            1  Means a successful operation.
     *            0 or neg value Means an unsuccessful operation.
     */
    int doResidualCalc(const doublereal time_curr, const int typeCalc, const doublereal* const y_curr,
                       const doublereal* const ydot_curr,
                       const ResidEval_Type_Enum evalType = Base_ResidEval) const;

    //! Compute the undamped Newton step
    /*!
     *  The residual function is evaluated at the current time, t_n, at the
     *  current values of the solution vector, m_y_n, and the solution time
     *  derivative, m_ydot_n. The Jacobian is not recomputed.
     *
     *  A factored Jacobian is reused, if available. If a factored Jacobian
     *  is not available, then the Jacobian is factored. Before factoring,
     *  the Jacobian is row and column-scaled. Column scaling is not
     *  recomputed. The row scales are recomputed here, after column
     *  scaling has been implemented.
     *
     *  @param time_curr      Current value of the time
     *  @param y_curr         Current value of the solution
     *  @param ydot_curr      Current value of the solution derivative.
     *  @param delta_y        return value of the raw change in y
     *  @param jac            Jacobian
     *
     *  @return Returns the result code from LAPACK. A zero means success.
     *          Anything else indicates a failure.
     */
    int doNewtonSolve(const doublereal time_curr, const doublereal* const y_curr,
                      const doublereal* const ydot_curr, doublereal* const delta_y,
                      GeneralMatrix& jac);

    //! Compute the Newton step, either by direct Newton's or by solving a
    //! close problem that is represented by a Hessian
    /*!
     * This is algorith A.6.5.1 in Dennis / Schnabel
     *
     * Compute the QR decomposition
     *
     * Compute the undamped Newton step.  The residual function is
     * evaluated at the current time, t_n, at the current values of the
     * solution vector, m_y_n, and the solution time derivative, m_ydot_n.
     * The Jacobian is not recomputed.
     *
     *  A factored Jacobian is reused, if available. If a factored Jacobian
     *  is not available, then the Jacobian is factored. Before factoring,
     *  the Jacobian is row and column-scaled. Column scaling is not
     *  recomputed. The row scales are recomputed here, after column
     *  scaling has been implemented.
     *
     *  @param y_curr         Current value of the solution
     *  @param ydot_curr      Current value of the solution derivative.
     *  @param delta_y        return value of the raw change in y
     *  @param jac            Jacobian
     *
     *  Internal input
     * ---------------
     *  internal m_resid      Stored residual is used as input
     *
     *  @return Returns the result code from LAPACK. A zero means success. Anything
     *          else indicates a failure.
     */
    int doAffineNewtonSolve(const doublereal* const y_curr, const doublereal* const ydot_curr,
                            doublereal* const delta_y, GeneralMatrix& jac);

    //! Calculate the length of the current trust region in terms of the solution error norm
    /*!
     *  We carry out a norm of deltaX_trust_ first. Then, we multiply that value
     *  by trustDelta_
     */
    doublereal trustRegionLength() const;

    //! Set default deulta bounds amounts
    /*!
     *     Delta bounds are set to 0.01 for all unknowns arbitrarily and capriciously
     *     Then, for each call to the nonlinear solver
     *     Then, they are increased to  1000 x atol
     *     then, they are increased to 0.1 fab(y[i])
     */
    void setDefaultDeltaBoundsMagnitudes();

    //! Adjust the step minimums
    void adjustUpStepMinimums();

    //! Set the delta Bounds magnitudes by hand
    /*!
     *  @param deltaBoundsMagnitudes  set the deltaBoundsMagnitude vector
     */
    void setDeltaBoundsMagnitudes(const doublereal* const deltaBoundsMagnitudes);

protected:
    //! Readjust the trust region vectors
    /*!
     *  The trust region is made up of the trust region vector calculation and the trustDelta_ value
     *  We periodically recalculate the trustVector_ values so that they renormalize to the
     *  correct length.  We change the trustDelta_ values regularly
     *
     *    The trust region calculate is based on
     *
     *        || delta_x   dot  1/trustDeltaX_ ||   <= trustDelta_
     */
    void readjustTrustVector();

    //! Fill a dogleg solution step vector
    /*!
     *   Previously, we have filled up deltaX_Newton_[], deltaX_CP_[], and Nuu_, so that
     *   this routine is straightforward.
     *
     *  @param leg      Leg of the dog leg you are on (0, 1, or 2)
     *  @param alpha    Relative length along the dog length that you are on.
     *  @param deltaX   Vector to be filled up
     */
    void fillDogLegStep(int leg, doublereal alpha, std::vector<doublereal>  & deltaX) const;

    //! Calculate the trust distance of a step in the solution variables
    /*!
     *  The trust distance is defined as the length of the step according to the norm wrt to the trust region.
     *  We calculate the trust distance by the following method.
     *
     *      trustDist =  || delta_x   dot  1/trustDeltaX_ ||
     *
     * @param deltaX  Current value of deltaX
     */
    doublereal calcTrustDistance(std::vector<doublereal> const& deltaX) const;

public:
    //! Bound the step
    /*!
     * Return the factor by which the undamped Newton step 'step0'
     * must be multiplied in order to keep all solution components in
     * all domains between their specified lower and upper bounds.
     * Other bounds may be applied here as well.
     *
     * Currently the bounds are hard coded into this routine:
     *
     *  Minimum value for all variables: - 0.01 * m_ewt[i]
     *  Maximum value = none.
     *
     * Thus, this means that all solution components are expected
     * to be numerical greater than zero in the limit of time step
     * truncation errors going to zero.
     *
     * Delta bounds: The idea behind these is that the Jacobian
     *               couldn't possibly be representative if the
     *               variable is changed by a lot. (true for
     *               nonlinear systems, false for linear systems)
     *  Maximum increase in variable in any one Newton iteration:
     *   factor of 2
     *  Maximum decrease in variable in any one Newton iteration:
     *   factor of 5
     *
     *   @param y         Current solution value of the old step
     *   @param step0     Proposed step change in the solution
     *
     *  @return  Returns the damping factor determined by the bounds calculation
     */
    doublereal boundStep(const doublereal* const  y, const doublereal* const step0);

    //! Set bounds constraints for all variables in the problem
    /*!
     *   @param y_low_bounds  Vector of lower bounds
     *   @param y_high_bounds Vector of high bounds
     */
    void setBoundsConstraints(const doublereal* const y_low_bounds,
                              const doublereal* const y_high_bounds);

    //! Return an editable vector of the low bounds constraints
    std::vector<doublereal> & lowBoundsConstraintVector();

    //! Return an editable vector of the high bounds constraints
    std::vector<doublereal> & highBoundsConstraintVector();

    //!   Internal function to calculate the time derivative of the solution at the new step
    /*!
     *  Previously, the user must have supplied information about the previous time step for this routine to
     *  work as intended.
     *
     *   @param  order of the BDF method
     *   @param   y_curr current value of the solution
     *   @param   ydot_curr  Calculated value of the solution derivative that is consistent with y_curr
     */
    void calc_ydot(const int order, const doublereal* const y_curr, doublereal* const ydot_curr) const;

    //! Function called to evaluate the Jacobian matrix and the current
    //! residual vector at the current time step
    /*!
     *  @param J  Jacobian matrix to be filled in
     *  @param f   Right hand side. This routine returns the current
     *             value of the RHS (output), so that it does
     *             not have to be computed again.
     *  @param time_curr Current time
     *  @param CJ  inverse of the value of deltaT
     *  @param y    value of the solution vector
     *  @param ydot  value of the time derivative of the solution vector
     *  @param num_newt_its Number of Newton iterations
     *
     * @return Returns a flag to indicate that operation is successful.
     *            1  Means a successful operation
     *            0  Means an unsuccessful operation
     */
    int  beuler_jac(GeneralMatrix& J, doublereal* const f,
                    doublereal time_curr, doublereal CJ, doublereal* const y,
                    doublereal* const ydot, int num_newt_its);

    //! Apply a filtering process to the step
    /*!
     *  @param timeCurrent    Current value of the time
     *  @param ybase          current value of the solution
     *   @param step0     Proposed step change in the solution
     *
     *  @return Returns the norm of the value of the amount filtered
     */
    doublereal filterNewStep(const doublereal timeCurrent, const doublereal* const ybase, doublereal* const step0);

    //! Apply a filter to the solution
    /*!
     *  @param timeCurrent   Current value of the time
     *  @param y_current     current value of the solution
     *  @param ydot_current   Current value of the solution derivative.
     *
     *  @return Returns the norm of the value of the amount filtered
     */
    doublereal filterNewSolution(const doublereal timeCurrent, doublereal* const y_current,
                                 doublereal* const ydot_current);

    //! Return the factor by which the undamped Newton step 'step0'
    //!  must be multiplied in order to keep the update within the bounds of an accurate Jacobian.
    /*!
     *  The idea behind these is that the Jacobian couldn't possibly be representative, if the
     *  variable is changed by a lot. (true for nonlinear systems, false for linear systems)
     *  Maximum increase in variable in any one Newton iteration: factor of 1.5
     *  Maximum decrease in variable in any one Newton iteration: factor of 2
     *
     *  @param y   Initial value of the solution vector
     *  @param step0  initial proposed step size
     *
     *  @return returns the damping factor
     */
    doublereal deltaBoundStep(const doublereal* const y, const doublereal* const step0);

    //!  Find a damping coefficient through a look-ahead mechanism
    /*!
     *    On entry, step_1 must contain an undamped Newton step for the
     *    solution y_n_curr. This method attempts to find a damping coefficient
     *    such that all components stay in bounds, and  the next
     *    undamped step would have a norm smaller than
     *    that of step_1. If successful, the new solution after taking the
     *    damped step is returned in y_n_1, and the undamped step at y_n_1 is
     *    returned in step_2.
     *
     *    @param time_curr Current physical time
     *    @param y_n_curr     Base value of the solution before any steps
     *                        are taken
     *    @param ydot_n_curr  Base value of the time derivative of the
     *                        solution
     *    @param step_1       Initial step suggested.
     *    @param y_n_1        Value of y1, the suggested solution after damping
     *    @param ydot_n_1     Value of the time derivative of the solution at y_n_1
     *    @param step_2       Value of the step change from y_n_1 to y_n_2
     *    @param stepNorm_2   norm of the step change in going from y_n_1 to y_n_2
     *    @param jac          Jacobian
     *    @param writetitle   Write a title line
     *    @param num_backtracks Number of backtracks taken
     *
     * @return   1 Successful step was taken: Next step was less than previous step.
     *                                        s1 is calculated
     *           2 Successful step: Next step's norm is less than 0.8
     *           3 Success:  The final residual is less than 1.0
     *                        A predicted deltaSoln1 is not produced however. s1 is estimated.
     *           4 Success:  The final residual is less than the residual
     *                       from the previous step.
     *                        A predicted deltaSoln1 is not produced however. s1 is estimated.
     *           0 Uncertain Success: s1 is about the same as s0
     *           NSOLN_RETN_FAIL_DAMPSTEP
     *             Unsuccessful step. We can not find a damping factor that is suitable.
     */
    int dampStep(const doublereal time_curr, const doublereal* const y_n_curr,
                 const doublereal* const ydot_n_curr, doublereal* const step_1,
                 doublereal* const y_n_1, doublereal* const ydot_n_1, doublereal* step_2,
                 doublereal& stepNorm_2, GeneralMatrix& jac, bool writetitle,
                 int& num_backtracks);

    //! Find the solution to F(X) = 0 by damped Newton iteration.
    /*!
     * On entry, x0 contains an initial estimate of the solution.  On
     * successful return, x1 contains the converged solution.
     *
     * SolnType = TRANSIENT -> we will assume we are relaxing a transient
     *        equation system for now. Will make it more general later,
     *        if an application comes up.
     *
     *  @param SolnType  Solution type
     *  @param y_comm    Initial value of the solution. On return this is the converged
     *                   value of the solution
     *  @param ydot_comm  Initial value of the solution derivative. On return this is the
     *                    converged value of the solution derivative.
     *  @param CJ        Inverse of the value of deltaT
     *  @param time_curr  Current value of the time
     *  @param jac        Matrix that will be used to store the Jacobian
     *  @param num_newt_its Number of Newton iterations taken
     *  @param num_linear_solves Number of linear solves taken
     *  @param num_backtracks Number of backtracking steps taken
     *  @param loglevelInput  Input log level determines the amount of printing.
     *
     *   @return  A positive value indicates a successful convergence
     *            -1  Failed convergence
     */
    int solve_nonlinear_problem(int SolnType, doublereal* const y_comm, doublereal* const ydot_comm, doublereal CJ,
                                doublereal time_curr, GeneralMatrix& jac, int& num_newt_its,
                                int& num_linear_solves,	int& num_backtracks, int loglevelInput);

    //! Set the values for the previous time step
    /*!
     *   We set the values for the previous time step here. These are used in the nonlinear
     *   solve because they affect the calculation of ydot.
     *
     * @param y_nm1     Value of the solution vector at the previous time step
     * @param ydot_nm1  Value of the solution vector derivative at the previous time step
     */
    virtual void
    setPreviousTimeStep(const std::vector<doublereal>& y_nm1, const std::vector<doublereal>& ydot_nm1);

private:
    //! Set the column scaling vector at the current time
    void calcColumnScales();

public:

    //! Set the column scaling that are used for the inversion of the matrix
    /*!
     *  There are three ways to do this.
     *
     *  The first method is to set the bool useColScaling to true, leaving the scaling factors unset.
     *  Then, the column scales will be set to the solution error weighting factors. This has the
     *  effect of ensuring that all delta variables will have the same order of magnitude at convergence
     *  end.
     *
     *  The second way is the explicitly set the column factors in the second argument of this function call.
     *
     *  The final way to input the scales is to override the ResidJacEval member function call,
     *
     *     calcSolnScales(double time_n, const double *m_y_n_curr, const double *m_y_nm1, double *m_colScales)
     *
     *  Overriding this function call will trump all other ways to specify the column scaling factors.
     *
     *  @param useColScaling   Turn this on if you want to use column scaling in the calculations
     *  @param scaleFactors    A vector of doubles that specifies the column factors.
     */
    void setColumnScaling(bool useColScaling, const double* const scaleFactors = 0);

    //! Set the rowscaling that are used for the inversion of the matrix
    /*!
     * Row scaling is set here. Right now the row scaling is set internally in the code.
     *
     * @param useRowScaling   Turn row scaling on or off.
     */
    void setRowScaling(bool useRowScaling);

    //! Scale the matrix
    /*!
     *  @param jac              Jacobian
     *  @param y_comm           Current value of the solution vector
     *  @param ydot_comm        Current value of the time derivative of the solution vector
     *  @param time_curr        current value of the time
     *  @param num_newt_its      Current value of the number of newt its
     */
    void scaleMatrix(GeneralMatrix& jac, doublereal* const y_comm, doublereal* const ydot_comm,
                     doublereal time_curr, int num_newt_its);

    //! Print solution norm contribution
    /*!
     *  Prints out the most important entries to the update to the solution vector for the current step
     *
     *   @param step_1                 Raw update vector for the current nonlinear step
     *   @param stepNorm_1             Norm of the vector step_1
     *   @param step_2                 Raw update vector for the next solution value based on the old matrix
     *   @param stepNorm_2             Norm of the vector step_2
     *   @param title                  title of the printout
     *   @param y_n_curr               Old value of the solution
     *   @param y_n_1                  New value of the solution after damping corrections
     *   @param damp                   Value of the damping factor
     *   @param num_entries            Number of entries to print out
     */
    void
    print_solnDelta_norm_contrib(const doublereal* const step_1, const char* const stepNorm_1,
                                 const doublereal* const step_2, const char* const stepNorm_2,
                                 const char* const title, const doublereal* const y_n_curr,
                                 const doublereal* const y_n_1,  doublereal damp, size_t num_entries);

    //! Compute the Residual Weights
    /*!
     *  The residual weights are defined here to be equal to the inverse of the row scaling factors used to
     *  row scale the matrix, after column scaling is used. They are multiplied by 10-3 because the column
     *  weights are also multiplied by that same quantity.
     *
     *  The basic idea is that a change in the solution vector on the order of the convergence tolerance
     *  multiplied by  [RJC] which is of order one after row scaling should give you the relative weight
     *  of the row. Values of the residual for that row can then be normalized by the value of this weight.
     *  When the tolerance in delta x is achieved, the tolerance in the residual is also achieved.
     */
    void computeResidWts();

    //! Return the residual weights
    /*!
     *  @param residWts  Vector of length neq_
     */
    void getResidWts(doublereal* const residWts) const;

    //! Check to see if the nonlinear problem has converged
    /*!
     * @param dampCode Code from the damping routine
     * @param s1       Value of the norm of the step change
     *
     * @return integer is returned. If positive, then the problem has converged
     *           1 Successful step was taken: Next step's norm is less than 1.0.
     *                                        The final residual norm is less than 1.0.
     *           2 Successful step: Next step's norm is less than 0.8.
     *                              This step's norm is less than 1.0.
     *                              The residual norm can be anything.
     *           3 Success:  The final residual is less than 1.0
     *                        The predicted deltaSoln is below 1.0.
     *           0 Not converged yet
     */
    int convergenceCheck(int dampCode, doublereal s1);

    //! Set the absolute tolerances for the solution variables
    /*!
     *   Set the absolute tolerances used in the calculation
     *
     *  @param atol   Vector of length neq_ that contains the tolerances to be used for the solution variables
     */
    void setAtol(const doublereal* const atol);

    //! Set the relative tolerances for the solution variables
    /*!
     *   Set the relative tolerances used in the calculation
     *
     *  @param rtol  single double
     */
    void setRtol(const doublereal rtol);

    //! Set the relative and absolute tolerances for the Residual norm comparisons, if used
    /*!
     *   Residual norms are used to calculate convergence within the nonlinear solver, since
     *   these are the norms that are associated with convergence proofs, especially for ill-conditioned systems.
     *   Usually the residual weights for each row are calculated by the program such that they
     *   correlate with the convergence requirements on the solution variables input by the user using
     *   the routines setAtol() and setRtol().
     *   The residual weights are essentially calculated from the value
     *
     *        residWeightNorm[i]  = m_ScaleSolnNormToResNorm  * sum_j ( fabs(A_i,j) ewt(j))
     *
     *    The factor,  m_ScaleSolnNormToResNorm, is computed periodically to ensure that the solution norms
     *    and the residual norms are converging at the same time and thus accounts for some-illconditioning issues
     *    but not all.
     *
     *    With this routine the user can override or add to the residual weighting norm evaluation by specifying
     *    their own vector of residual absolute and relative tolerances.
     *
     *   The user specified tolerance for the residual is given by the following quantity
     *
     *   residWeightNorm[i] = residAtol[i] +  residRtol * m_rowWtScales[i] / neq
     *
     *     @param residRtol  scalar residual relative tolerance
     *     @param residAtol  vector of residual absolute tolerances
     *
     *     @param residNormHandling  Parameter that sets the default handling of the residual norms
     *                      0   The residual weighting vector is calculated to make sure that the solution
     *                          norms are roughly 1 when the residual norm is roughly 1.
     *                          This is the default if this routine is not called.
     *                      1   Use the user residual norm specified by the parameters in this routine
     *                      2   Use the minimum value of the residual weights calculated by method 1 and 2.
     *                          This is the default if this routine is called and this parameter isn't specified.
     */
    void setResidualTols(double residRtol, double* residAtol, int residNormHandling = 2);

    //! Set the value of the maximum # of Newton iterations
    /*!
     *  @param maxNewtIts   Maximum number of Newton iterations
     */
    void setMaxNewtIts(const int maxNewtIts);

    //! Calculate the scaling factor for translating residual norms into solution norms.
    /*!
     *  This routine calls computeResidWts() a couple of times in the calculation of m_ScaleSolnNormToResNorm.
     *  A more sophisticated routine may do more with signs to get a better value. Perhaps, a series of calculations
     *  with different signs attached may be in order. Then, m_ScaleSolnNormToResNorm would be calculated
     *  as the minimum of a series of calculations.
     */
    void calcSolnToResNormVector();

    //! Calculate the steepest descent direction and the Cauchy Point where the quadratic formulation
    //! of the nonlinear problem expects a minimum along the descent direction.
    /*!
     *  @param jac   Jacobian matrix: must be unfactored.
     *
     *  @return Returns the norm of the solution update
     */
    doublereal doCauchyPointSolve(GeneralMatrix& jac);

    //! This is a utility routine that can be used to print out the rates of the initial residual decline
    /*!
     *  The residual**2 decline for various directions is printed out. The rate of decline of the
     *  square of the residuals multiplied by the number of equations along each direction is printed out
     *  This quantity can be directly related to the theory, and may be calculated from derivatives at the
     *  original point.
     *
     *   ( (r)**2 * neq   - (r0)**2 * neq ) / distance
     *
     *  What's printed out:
     *
     *  The theoretical linearized residual decline
     *  The actual residual decline in the steepest descent direction determined by numerical differencing
     *  The actual residual decline in the Newton direction determined by numerical differencing
     *
     *  This routine doesn't need to be called for the solution of the nonlinear problem.
     *
     *  @param time_curr   Current time
     *  @param ydot0       INPUT    Current value of the derivative of the solution vector
     *  @param ydot1       INPUT  Time derivatives of solution at the conditions which are evaluated for success
     *  @param numTrials   OUTPUT Counter for the number of residual evaluations
     */
    void descentComparison(doublereal time_curr ,doublereal* ydot0, doublereal* ydot1, int& numTrials);

    //!  Setup the parameters for the double dog leg
    /*!
     *  The calls to the doCauchySolve() and doNewtonSolve() routines are done at the main level. This routine comes
     *  after those calls.  We calculate the point Nuu_ here, the distances of the dog-legs,
     *  and the norms of the CP and Newton points in terms of the trust vectors.
     */
    void setupDoubleDogleg();

    //! Change the global lambda coordinate into the (leg,alpha) coordinate for the double dogleg
    /*!
     * @param lambda Global value of the distance along the double dogleg
     * @param alpha  relative value along the particular leg
     *
     * @return Returns the leg number ( 0, 1, or 2).
     */
    int lambdaToLeg(const doublereal lambda, doublereal& alpha) const;

    //! Given a trust distance, this routine calculates the intersection of the this distance with the
    //! double dogleg curve
    /*!
     *   @param      trustVal    (INPUT)     Value of the trust distance
     *   @param      lambda      (OUTPUT)    Returns the internal coordinate of the double dogleg
     *   @param      alpha       (OUTPUT)    Returns the relative distance along the appropriate leg
     *   @return     leg         (OUTPUT)    Returns the leg ID (0, 1, or 2)
     */
    int calcTrustIntersection(doublereal trustVal, doublereal& lambda, doublereal& alpha) const;

    //! Initialize the size of the trust vector.
    /*!
     *  The algorithm we use is to set it equal to the length of the Distance to the Cauchy point.
     */
    void initializeTrustRegion();

    //! Set Trust region initialization strategy
    /*!
     *  The default is use method 2 with a factor of 1.
     *  Then, on subsequent invocations of solve_nonlinear_problem() the strategy flips to method 0.
     *
     * @param method Method to set the strategy
     *            0 No strategy - Use the previous strategy
     *            1 Factor of the solution error weights
     *            2 Factor of the first Cauchy Point distance
     *            3 Factor of the first Newton step distance
     *
     *  @param factor Factor to use in combination with the method
     */
    void setTrustRegionInitializationMethod(int method, doublereal factor);

    //! Damp using the dog leg approach
    /*!
     * @param time_curr  INPUT     Current value of the time
     * @param y_n_curr   INPUT    Current value of the solution vector
     * @param ydot_n_curr INPUT   Current value of the derivative of the solution vector
     * @param step_1     INPUT    First trial step for the first iteration
     * @param y_n_1   INPUT       First trial value of the solution vector
     * @param ydot_n_1 INPUT      First trial value of the derivative of the solution vector
     * @param stepNorm_1        OUTPUT   Norm of the vector step_1
     * @param stepNorm_2        OUTPUT   Estimated norm of the vector step_2
     * @param jac        INPUT    Jacobian
     * @param num_backtracks OUTPUT  number of backtracks taken in the current damping step
     *
     *  @return  1 Successful step was taken. The predicted residual norm is less than one
     *           2 Successful step: Next step's norm is less than 0.8
     *           3 Success:  The final residual is less than 1.0
     *                        A predicted deltaSoln1 is not produced however. s1 is estimated.
     *           4 Success:  The final residual is less than the residual
     *                       from the previous step.
     *                        A predicted deltaSoln1 is not produced however. s1 is estimated.
     *           0 Uncertain Success: s1 is about the same as s0
     *          -2 Unsuccessful step.
     */
    int dampDogLeg(const doublereal time_curr, const doublereal* y_n_curr,
                   const doublereal* ydot_n_curr,  std::vector<doublereal> & step_1,
                   doublereal* const y_n_1, doublereal* const ydot_n_1,
                   doublereal& stepNorm_1, doublereal& stepNorm_2, GeneralMatrix& jac, int& num_backtracks);

    //! Decide whether the current step is acceptable and adjust the trust region size
    /*!
     *  This is an extension of algorithm 6.4.5 of Dennis and Schnabel.
     *
     *  Here we decide whether to accept the current step
     *  At the end of the calculation a new estimate of the trust region is calculated
     *
     * @param time_curr  INPUT     Current value of the time
     * @param leg        INPUT    Leg of the dogleg that we are on
     * @param alpha      INPUT    Distance down that leg that we are on
     * @param y_n_curr   INPUT    Current value of the solution vector
     * @param ydot_n_curr INPUT    Current value of the derivative of the solution vector
     * @param step_1     INPUT    Trial step
     * @param y_n_1         OUTPUT   Solution values at the conditions which are evaluated for success
     * @param ydot_n_1      OUTPUT   Time derivatives of solution at the conditions which are evaluated for success
     * @param trustDeltaOld INPUT Value of the trust length at the old conditions
     *
     * @return This function returns a code which indicates whether the step will be accepted or not.
     *        3  Step passed with flying colors. Try redoing the calculation with a bigger trust region.
     *        2  Step didn't pass deltaF requirement. Decrease the size of the next  trust region for a retry and return
     *        0  The step passed.
     *       -1  The step size is now too small (||d || < 0.1). A really small step isn't decreasing the function.
     *           This is an error condition.
     *       -2  Current value of the solution vector caused a residual error in its evaluation.
     *           Step is a failure, and the step size must be reduced in order to proceed further.
     */
    int decideStep(const doublereal time_curr, int leg, doublereal alpha, const doublereal* const y_n_curr,
                   const doublereal* const ydot_n_curr,
                   const std::vector<doublereal> & step_1,
                   const doublereal* const y_n_1, const doublereal* const ydot_n_1, doublereal trustDeltaOld);

    //! Calculated the expected residual along the double dogleg curve.
    /*!
     *  @param leg 0, 1, or 2 representing the curves of the dogleg
     *  @param alpha  Relative distance along the particular curve.
     *
     *  @return Returns the expected value of the residual at that point according to the quadratic model.
     *          The residual at the Newton point will always be zero.
     */
    doublereal expectedResidLeg(int leg, doublereal alpha) const;

    //!  Here we print out the residual at various points along the double dogleg, comparing against the quadratic model
    //!  in a table format
    /*!
     *  @param time_curr     INPUT    current time
     *  @param ydot0         INPUT    Current value of the derivative of the solution vector for non-time dependent
     *                                determinations
     *  @param legBest       OUTPUT   leg of the dogleg that gives the lowest residual
     *  @param alphaBest     OUTPUT   distance along dogleg for best result.
     */
    void residualComparisonLeg(const doublereal time_curr, const doublereal* const ydot0, int& legBest,
                               doublereal& alphaBest) const;

    //! Set the print level from the nonlinear solver
    /*!
     *   0 -> absolutely nothing is printed for a single time step.
     *   1 -> One line summary per solve_nonlinear call
     *   2 -> short description, points of interest: Table of nonlinear solve - one line per iteration
     *   3 -> Table is included -> More printing per nonlinear iteration (default) that occurs during the table
     *   4 -> Summaries of the nonlinear solve iteration as they are occurring -> table no longer printed
     *   5 -> Algorithm information on the nonlinear iterates are printed out
     *   6 -> Additional info on the nonlinear iterates are printed out
     *   7 -> Additional info on the linear solve is printed out.
     *   8 -> Info on a per iterate of the linear solve is printed out.
     *
     *  @param printLvl  integer value
     */
    void setPrintLvl(int printLvl);

    //! Parameter to turn on solution solver schemes
    /*!
     *     @param doDogLeg Parameter to turn on the double dog leg scheme
     *                     Default is to always use a damping scheme in the Newton Direction.
     *                     When this is nonzero, a model trust region approach is used using a double dog leg
     *                     with the steepest descent direction used for small step sizes.
     *
     *     @param doAffineSolve  Parameter to turn on or off the solution of the system using a Hessian
     *                           if the matrix has a bad condition number.
     */
    void setSolverScheme(int doDogLeg, int doAffineSolve);

    /*
     * -----------------------------------------------------------------------------------------------------------------
     *              MEMBER DATA
     * ------------------------------------------------------------------------------------------------
     */
private:

    //! Pointer to the residual and Jacobian evaluator for the
    //! function
    /*!
     *   See ResidJacEval.h for an evaluator.
     */
    ResidJacEval* m_func;

    //! Solution type
    int solnType_;

    //! Local copy of the number of equations
    size_t neq_;

    //! Soln error weights
    std::vector<doublereal> m_ewt;

    //! Boolean indicating whether a manual delta bounds has been input.
    int m_manualDeltaStepSet;

    //! Soln Delta bounds magnitudes
    std::vector<doublereal> m_deltaStepMinimum;

    //! Value of the delta step magnitudes
    std::vector<doublereal> m_deltaStepMaximum;

    //! Vector containing the current solution vector within the nonlinear solver
    std::vector<doublereal> m_y_n_curr;

    //! Vector containing the time derivative of the current solution vector within the nonlinear solver
    //! (where applicable)
    std::vector<doublereal> m_ydot_n_curr;

    //! Vector containing the solution at the previous time step
    std::vector<doublereal> m_y_nm1;

    //! Vector containing the solution derivative at the previous time step
    std::vector<doublereal> m_ydot_nm1;

    //! Vector containing the solution at the new point that is to be considered
    std::vector<doublereal> m_y_n_trial;

    //! Value of the solution time derivative at the new point that is to be considered
    std::vector<doublereal> m_ydot_trial;

    //! Value of the step to be taken in the solution
    std::vector<doublereal> m_step_1;

    //! Vector of column scaling factors
    std::vector<doublereal> m_colScales;

    //! Weights for normalizing the values of the residuals
    /*!
     *  These are computed if row scaling, m_rowScaling, is turned on. They are calculated currently as the
     *  sum of the absolute values of the rows of the Jacobian.
     */
    std::vector<doublereal> m_rowScales;

    //! Weights for normalizing the values of the residuals
    /*!
     *  They are calculated as the  sum of the absolute values of the Jacobian
     *  multiplied by the solution weight function.
     *  This is carried out in scaleMatrix().
     */
    std::vector<doublereal> m_rowWtScales;

    //! Value of the residual for the nonlinear problem
    mutable std::vector<doublereal> m_resid;

    //! Workspace of length neq_
    mutable std::vector<doublereal> m_wksp;

    //! Workspace of length neq_
    mutable std::vector<doublereal> m_wksp_2;

    /*****************************************************************************************
     *        INTERNAL WEIGHTS FOR TAKING SOLUTION NORMS
     ******************************************************************************************/
    //! Vector of residual weights
    /*!
     *   These are used to establish useful and informative weighted norms of the residual vector.
     */
    std::vector<doublereal> m_residWts;

    //! Norm of the residual at the start of each nonlinear iteration
    doublereal m_normResid_0;

    //! Norm of the residual after it has been bounded
    doublereal m_normResid_Bound;

    //! Norm of the residual at the end of the first leg of the current iteration
    doublereal m_normResid_1;

    //! Norm of the residual at the end of the first leg of the current iteration
    doublereal m_normResid_full;

    //! Norm of the solution update created by the iteration in its raw, undamped form, using the solution norm
    doublereal m_normDeltaSoln_Newton;

    //! Norm of the distance to the Cauchy point using the solution norm
    doublereal m_normDeltaSoln_CP;

    //! Norm of the residual for a trial calculation which may or may not be used
    doublereal m_normResidTrial;

    //! Vector of the norm
    doublereal m_normResidPoints[15];

    //! Boolean indicating whether we should scale the residual
    mutable bool m_resid_scaled;

    /*****************************************************************************************
     *        INTERNAL BOUNDARY INFO FOR SOLUTIONS
     *****************************************************************************************/

    //! Bounds vector for each species
    std::vector<doublereal> m_y_high_bounds;

    //! Lower bounds vector for each species
    std::vector<doublereal> m_y_low_bounds;

    //! Damping factor imposed by hard bounds and by delta bounds
    doublereal m_dampBound;

    //! Additional damping factor due to bounds on the residual and solution norms
    doublereal m_dampRes;

    //! Delta t for the current step
    doublereal delta_t_n;

    //! Counter for the total number of function evaluations
    mutable int m_nfe;

    /***********************************************************************************************
     *             MATRIX INFORMATION
     **************************************************************************************/

    //!  The type of column scaling used in the matrix inversion of the problem
    /*!
     *    If 1 then colScaling = m_ewt[]
     *    If 2 then colScaling = user set
     *    if 0 then colScaling = 1.0
     */
    int m_colScaling;

    //! int indicating whether row scaling is turned on (1) or not (0)
    int m_rowScaling;

    //! Total number of linear solves taken by the solver object
    int m_numTotalLinearSolves;

    //! Number of local linear solves done during the current iteration
    int m_numLocalLinearSolves;

    //! Total number of Newton iterations
    int m_numTotalNewtIts;

public:
    //! Minimum number of Newton iterations to use
    int m_min_newt_its;
private:

    //! Maximum number of Newton iterations
    int maxNewtIts_;

    //! Jacobian formation method
    /*!
     *   1 = numerical  (default)
     *   2 = analytical
     */
    int m_jacFormMethod;

    //! Number of Jacobian evaluations
    int m_nJacEval;

    //! Current system time
    /*!
     *  Note, we assume even for steady state problems that the residual
     *  is a function of a system time.
     */
    doublereal time_n;

    //! Boolean indicating matrix conditioning
    int m_matrixConditioning;

    //! Order of the time step method = 1
    int m_order;

    //! value of the relative tolerance to use in solving the equation set
    doublereal rtol_;

    //! Base value of the absolute tolerance
    doublereal atolBase_;

    //! absolute tolerance in the solution unknown
    /*!
     * This is used to evaluating the weighting factor
     */
    std::vector<doublereal> atolk_;

    //! absolute tolerance in the unscaled solution unknowns
    std::vector<doublereal> userResidAtol_;

    //! absolute tolerance in the unscaled solution unknowns
    doublereal userResidRtol_;

    //! Check the residual tolerances explicitly against user input
    /*!
     *   0 Don't calculate residual weights from residual tolerance inputs
     *   1 Calculate residual weights from residual tolerance inputs only
     *   2 Calculate residual weights from a minimum of the solution error weights process and the direct residual tolerance inputs
     */
    int checkUserResidualTols_;

    //! Determines the level of printing for each time step.
    /*!
     *   0 -> absolutely nothing is printed for a single time step.
     *   1 -> One line summary per solve_nonlinear call
     *   2 -> short description, points of interest: Table of nonlinear solve - one line per iteration
     *   3 -> Table is included -> More printing per nonlinear iteration (default) that occurs during the table
     *   4 -> Summaries of the nonlinear solve iteration as they are occurring -> table no longer printed
     *        Base_ShowSolution Residual called for residual printing at the end of convergence.
     *   5 -> Algorithm information on the nonlinear iterates are printed out
     *   6 -> Additional info on the nonlinear iterates are printed out
     *        Base_ShowSolution Residual called for residual printing at the end of each step.
     *   7 -> Additional info on the linear solve is printed out.
     *   8 -> Info on a per iterate of the linear solve is printed out.
     */
    int m_print_flag;

    //! Scale factor for turning residual norms into solution norms
    doublereal m_ScaleSolnNormToResNorm;

    //! Copy of the Jacobian that doesn't get overwritten when the inverse is determined
    /*!
     *  The Jacobian stored here is the raw matrix, before any row or column scaling is carried out
     */
    Cantera::GeneralMatrix* jacCopyPtr_;

    //! Hessian
    Cantera::GeneralMatrix* HessianPtr_;

    /*********************************************************************************************
     *      VARIABLES ASSOCIATED WITH STEPS AND ASSOCIATED DOUBLE DOGLEG PARAMETERS
     *********************************************************************************************/

    //!  Steepest descent direction. This is also the distance to the Cauchy Point
    std::vector<doublereal> deltaX_CP_;

    //! Newton Step - This is the Newton step determined from the straight Jacobian
    /*
     *  Newton step for the current step only
     */
    std::vector<doublereal> deltaX_Newton_;

    //! Expected value of the residual norm at the Cauchy point if the quadratic model
    //! were valid
    doublereal residNorm2Cauchy_;

    //! Current leg
    int dogLegID_;

    //! Current Alpha param along the leg
    doublereal dogLegAlpha_;

    //! Residual dot Jd norm
    /*!
     *  This is equal to  R_hat dot  J_hat d_y_descent
     */
    doublereal RJd_norm_;

    //! Value of lambdaStar_ which is used to calculate the Cauchy point
    doublereal lambdaStar_;

    //!  Jacobian times the steepest descent direction in the normalized coordinates.
    /*!
     *  This is equal to [ Jhat d^y_{descent} ] in the notes, Eqn. 18.
     */
    std::vector<doublereal> Jd_;

    //! Vector of trust region values.
    std::vector<doublereal> deltaX_trust_;

    //! Current norm of the vector deltaX_trust_ in terms of the solution norm
    mutable  doublereal norm_deltaX_trust_;

    //! Current value of trust radius. This is used with deltaX_trust_ to
    //! calculate the max step size.
    doublereal trustDelta_;

    //! Method for handling the trust region initialization
    /*!
     *  Then, on subsequent invocations of solve_nonlinear_problem() the strategy flips to method 0.
     *
     *  method Method to set the strategy
     *            0 No strategy - Use the previous strategy
     *            1 Factor of the solution error weights
     *            2 Factor of the first Cauchy Point distance
     *            3 Factor of the first Newton step distance
     */
    int trustRegionInitializationMethod_;

    //! Factor used to set the initial trust region
    doublereal trustRegionInitializationFactor_;

    //! Relative distance down the Newton step that the second dogleg starts
    doublereal Nuu_;

    //! Distance of the zeroeth leg of the dogleg in terms of the solution norm
    doublereal dist_R0_;

    //! Distance of the first leg of the dogleg in terms of the solution norm
    doublereal dist_R1_;

    //! Distance of the second leg of the dogleg in terms of the solution norm
    doublereal dist_R2_;

    //! Distance of the sum of all legs of the doglegs in terms of the solution norm
    doublereal dist_Total_;

    //! Dot product of the Jd_ variable defined above with itself.
    doublereal JdJd_norm_;

    //! Norm of the Newton Step wrt trust region
    doublereal normTrust_Newton_;

    //! Norm of the Cauchy Step direction wrt trust region
    doublereal normTrust_CP_;

    //! General toggle for turning on dog leg damping.
    int doDogLeg_;

    //! General toggle for turning on Affine solve with Hessian
    int doAffineSolve_;

    //! Condition number of the matrix
    doublereal m_conditionNumber;

    //! Factor indicating how much trust region has been changed this iteration - output variable
    doublereal CurrentTrustFactor_;

    //! Factor indicating how much trust region has been changed next iteration - output variable
    doublereal NextTrustFactor_;

    //! Boolean indicating that the residual weights have been reevaluated this iteration - output variable
    bool ResidWtsReevaluated_;

    //! Expected DResid_dS for the steepest descent path - output variable
    doublereal ResidDecreaseSDExp_;

    //! Actual DResid_dS for the steepest descent path - output variable
    doublereal ResidDecreaseSD_;

    //! Expected DResid_dS for the Newton path - output variable
    doublereal ResidDecreaseNewtExp_;

    //! Actual DResid_dS for the Newton path - output variable
    doublereal ResidDecreaseNewt_;

    /*******************************************************************************************
     *     STATIC VARIABLES
     *****************************************************************************************/

public:
    //! Turn off printing of time
    /*!
     *  Necessary to do for test suites
     */
    static bool s_TurnOffTiming;

    //! Turn on or off printing of the Jacobian
    static bool s_print_NumJac;

    //! Turn on extra printing of dogleg information
    static bool s_print_DogLeg;

    //! Turn on solving both the Newton and Hessian systems and comparing the results
    /*!
     *  This is off by default
     */
    static bool s_doBothSolvesAndCompare;

    //! This toggle turns off the use of the Hessian when it is warranted by the condition number.
    /*!
     *   This is a debugging option.
     */
    static bool s_alwaysAssumeNewtonGood;

};

}

#endif
