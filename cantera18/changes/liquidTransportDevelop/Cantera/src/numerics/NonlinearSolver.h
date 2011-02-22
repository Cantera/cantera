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
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#ifndef CT_NONLINEARSOLVER_H
#define CT_NONLINEARSOLVER_H

#include "ResidJacEval.h"
#include "SquareMatrix.h"

namespace Cantera {
  
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
  ///  @name  Constant which determines the type of the Jacobian
  //! The jacobian will be calculated from a numerical method
#define NSOLN_JAC_NUM 1
  //! The jacobian is calculated from an analytical function
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
   *
   *   The basic idea is that we predict a direction that is parameterized by an overall coordinate
   *   value, beta, from zero to one, This may or may not be the same as the value, damp,
   *   depending upon whether the direction is straight.
   *
   *  
   *
   *  @code 
   *
   *
   *  NonlinearSolver *nls = new NonlinearSolver(&r1);
   *  
   *  int solnType =       NSOLN_TYPE_STEADY_STATE ;
   *
   *  nls->setDeltaBoundsMagnitudes(deltaBounds);
   *
   *  nls->solve_nonlinear_problem(solnType, y_comm, ydot_comm, CJ, time_curr, jac,
   *                                num_newt_its,  num_linear_solves, numBacktracks,
   *                               loglevelInput);
   *
   *  @endcode  
   *
   *
   *  @ingroup numerics
   */
  class NonlinearSolver {

  public:
    //! Default constructor
    /*!
     * @param func   Residual and jacobian evaluator function object
     */
    NonlinearSolver(ResidJacEval *func);

    //!Copy Constructor for the %ThermoPhase object.
    /*!
     * @param right Item to be copied
     */
    NonlinearSolver(const NonlinearSolver &right);

    //! Destructor
    ~NonlinearSolver();

    //! Assignment operator
    /*!
     *  This is NOT a virtual function.
     *
     * @param right    Reference to %NonlinearSolver object to be
     *                 copied into the
     *                 current one.
     */
    NonlinearSolver& operator=(const NonlinearSolver &right);

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
    void createSolnWeights(const doublereal * const y);

    //!  L2 norm of the delta of the solution vector
    /*!
     *  calculate the norm of the solution vector. This will
     *  involve the column scaling of the matrix
     *
     *    The third argument has a default of false. However,
     *    if true, then a table of the largest values is printed
     *    out to standard output.
     *
     *  @param delta_y       Vector to take the norm of
     *  @param title         Optional title to be printed out
     *  @param printLargest  int indicating how many specific lines should be printed out
     *  @param dampFactor    Current value of the damping factor. Defaults to 1.
     *                       only used for printout out a table.
     *
     *  @return Returns the L2 norm of the delta
     */
    doublereal solnErrorNorm(const doublereal * const delta_y,  const char * title = 0, int printLargest = 0, 
			 const doublereal dampFactor = 1.0) const;

    //! L2 norm of the residual of the equation system
    /*!
     * Calculate the norm of the residual vector. This may
     * involve using the row sum scaling from the matrix problem. 
     *
     *  The second argument has a default of false. However,
     *  if true, then a table of the largest values is printed
     *  out to standard output.
     *
     *  @param resid    Vector of the residuals
     *  @param title    Optional title to be printed out
     *  @param printLargest  Number of specific entries to be printed
     *  @param y        Current value of y - only used for printouts
     *
     *
     *  @return Returns the L2 norm of the delta
     */
    doublereal residErrorNorm(const doublereal * const resid, const char * title = 0, const int printLargest = 0,
			  const doublereal * const y = 0);

    //! Compute the current residual
    /*!
     *  The current value of the residual is storred in the internal work array m_resid.
     *
     *  @param time_curr    Value of the time 
     *  @param typeCalc     Type of the calculation
     *  @param y_curr       Current value of the solution vector
     *  @param ydot_curr    Current value of the time derivative of the solution vector
     *  @param evalType     Base evalulation type
     *                        Defaults to Base_ResidEval
     *
     * @return Returns a flag to indicate that operation is successful.
     *            1  Means a successful operation
     *           -0 or neg value Means an unsuccessful operation
     */
    int doResidualCalc(const doublereal time_curr, const int typeCalc, const doublereal * const y_curr, 
		       const doublereal * const ydot_curr,
		       const ResidEval_Type_Enum evalType = Base_ResidEval);

    //! Compute the undamped Newton step
    /*!
     *
     * Compute the undamped Newton step.  The residual function is
     * evaluated at the current time, t_n, at the current values of the
     * solution vector, m_y_n, and the solution time derivative, m_ydot_n. 
     * The Jacobian is not recomputed.
     *
     *  A factored jacobian is reused, if available. If a factored jacobian
     *  is not available, then the jacobian is factored. Before factoring,
     *  the jacobian is row and column-scaled. Column scaling is not 
     *  recomputed. The row scales are recomputed here, after column
     *  scaling has been implemented.
     *
     *  @param time_curr      Current value of the time
     *  @param y_curr         Current value of the solution
     *  @param ydot_curr      Current value of the solution derivative.
     *  @param delta_y        return value of the raw change in y
     *  @param jac            Jacobian
     *  @param loglevel       Log level
     *
     *  @return Returns the result code from lapack. A zero means success. Anything
     *          else indicates a failure.
     */ 
    int doNewtonSolve(const doublereal time_curr, const doublereal * const y_curr,
		      const doublereal * const ydot_curr, doublereal * const delta_y,
		      SquareMatrix& jac, int loglevel);

    //! Set default deulta bounds amounts
    /*!
     *     Delta bounds are set to 0.01 for all unknowns arbitrarily and capriciously
     *     Then, for each call to the nonlinear solver
     *     Then, they are increased to  1000 x atol
     *     then, they are increased to 0.1 fab(y[i])
     */
    void setDefaultDeltaBoundsMagnitudes();

    //! Set the delta Bounds magnitudes by hand
    /*!
     *  @param deltaBoundsMagnitudes  set the deltaBoundsMagnitude vector      
     */
    void setDeltaBoundsMagnitudes(const doublereal * const deltaBoundsMagnitudes);
  
  protected:
    //! Calculate the trust region vectors
    /*!
     *  The trust region is made up of the trust region vector calculation and the trustDelta_ value
     *  We periodically recalculate the trustVector_ values so that they renormalize to the
     *  correct length.  We change the trustDelta_ values regularly
     *
     *    The trust region calculate is based on 
     *
     *        || delta_x   dot  1/trustDeltaX_ ||   <= trustDelta_
     *
     * @param y   current value of the solution
     */
    void calcTrustVector();

    //! Fill a dogleg solution step vector
    /*!
     *   Previously, we have filled up deltaX_Newton_[], deltaX_CP_[], and Nuu_, so that
     *   this routine is straightforward.
     *
     *  @param leg      Leg of the dog leg you are on (0, 1, or 2)
     *  @param alpha    Relative length along the dog length that you are on.
     *  @param deltaX   Vector to be filled up
     */
    void fillDogLegStep(int leg, double alpha, std::vector<doublereal>  & deltaX) const;

    //! Calculate the trust distance of a step in the solution variables
    /*!
     *  The trust distance is defined as the length of the step according to the norm wrt to the trust region.
     *  We calculate the trust distance by the following method.
     *  
     *      trustDist =  || delta_x   dot  1/trustDeltaX_ || 
     *
     * @param deltaX  Current value of deltaX
     */
    doublereal calcTrustDistance(std::vector<doublereal> const & deltaX) const;



    int  calcTrustIntersection(double trustDelta, double &lambda, double &alpha) const;
  public:
    //! Bound the step
    /*!
     *
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
     *  Maximum increase in variable in any one newton iteration: 
     *   factor of 2
     *  Maximum decrease in variable in any one newton iteration:
     *   factor of 5
     *
     *   @param y         Current solution value of the old step
     *   @param step0     Proposed step change in the solution
     *   @param loglevel  Log level
     *
     *  @return  Returns the damping factor determined by the bounds calculation
     */
    doublereal boundStep(const doublereal * const  y, const doublereal * const step0,  const int loglevel);

    //! Set bounds constraints for all variables in the problem
    /*!
     *  
     *   @param y_low_bounds  Vector of lower bounds
     *   @param y_high_bounds Vector of high bounds
     */
    void setBoundsConstraints(const doublereal * const y_low_bounds,
			      const doublereal * const y_high_bounds);

    //! Return an editable vector of the low bounds constraints
    std::vector<double> & lowBoundsConstraintVector();

    //! Return an editable vector of the high bounds constraints
    std::vector<double> & highBoundsConstraintVector();
   
    //! Internal function to calculate the time derivative at the new step
    /*!
     *   @param  order of the BDF method
     *   @param   y_curr current value of the solution
     *   @param   ydot_curr  Calculated value of the solution derivative that is consistent with y_curr
     */ 
    void calc_ydot(const int order, const doublereal * const y_curr, doublereal * const ydot_curr);

    //! Function called to evaluate the jacobian matrix and the current
    //! residual vector at the current time step
    /*!
     *  
     *
     *  @param J  Jacobian matrix to be filled in
     *  @param f   Right hand side. This routine returns the current
     *             value of the rhs (output), so that it does
     *             not have to be computed again.
     *  @param time_curr Current time
     *  @param CJ  inverse of the value of deltaT
     *  @param y    value of the solution vector
     *  @param ydot  value of the time derivative of the solution vector
     *  @param num_newt_its Number of newton iterations
     *  
     * @return Returns a flag to indicate that operation is successful.
     *            1  Means a successful operation
     *            0  Means an unsuccessful operation
     */
    int  beuler_jac(SquareMatrix &J, doublereal * const f,
		    doublereal time_curr, doublereal CJ, doublereal * const y,
		    doublereal * const ydot, int num_newt_its);

    //! Apply a filtering process to the step
    /*!
     *  @param timeCurrent    Current value of the time
     *  @param ybase          current value of the solution
     *   @param step0     Proposed step change in the solution
     *
     *  @return Returns the norm of the value of the amount filtered
     */
    doublereal filterNewStep(const doublereal timeCurrent, const doublereal * const ybase, doublereal * const step0);

    //! Apply a filter to the solution
    /*!
     *  @param timeCurrent   Current value of the time
     *  @param y_current     current value of the solution
     *  @param ydot_current   Current value of the solution derivative.
     *
     *  @return Returns the norm of the value of the amount filtered
     */
    doublereal filterNewSolution(const doublereal timeCurrent, doublereal * const y_current, 
				 doublereal * const ydot_current);

    //! Return the factor by which the undamped Newton step 'step0'
    //!  must be multiplied in order to keep the update within the bounds of an accurate jacobian.
    /*!
     *
     *  The idea behind these is that the Jacobian couldn't possibly be representative, if the
     *  variable is changed by a lot. (true for nonlinear systems, false for linear systems)
     *  Maximum increase in variable in any one newton iteration:
     *      factor of 1.5
     *  Maximum decrease in variable in any one newton iteration:
     *      factor of 2
     *
     *  @param y   Initial value of the solution vector
     *  @param step0  initial proposed step size
     *  @param loglevel log level
     *
     *  @return returns the damping factor
     */
    doublereal deltaBoundStep(const doublereal * const y, const doublereal * const step0, const int loglevel);
			       
    //!  Find a damping coefficient through a look-ahead mechanism
    /*!
     *    On entry, step0 must contain an undamped Newton step for the
     *    solution x0. This method attempts to find a damping coefficient
     *    such that all components stay in bounds, and  the next
     *    undamped step would have a norm smaller than
     *    that of step0. If successful, the new solution after taking the
     *    damped step is returned in y1, and the undamped step at y1 is
     *    returned in step1.
     *
     *    @param time_curr Current physical time
     *    @param y0        Base value of the solution before any steps 
     *                     are taken
     *    @param ydot0     Base value of the time derivative of teh
     *                     solution
     *    @param step0     Initial step suggested.
     *    @param y1        Value of y1, the suggested solution after damping
     *    @param ydot1     Value of the time derivative of the solution at y1
     *    @param step1     Value of the step change from y0 to y1
     *    @param s1        norm of the step change in going from y0 to y1
     *    @param jac       Jacobian
     *    @param loglevel  Log level to be used
     *    @param writetitle  Write a title line
     *    @param num_backtracks Number of backtracks taken
     *
     *  @return returns an integer indicating what happened.
     */
    int dampStep(const doublereal time_curr, const double* y0, 
		 const doublereal *ydot0, const double* step0, 
		 double* const y1, double* const ydot1, double* step1,
		 double& s1, SquareMatrix& jac, int& loglevel, bool writetitle,
		 int& num_backtracks);

    //! Find the solution to F(X) = 0 by damped Newton iteration. 
    /*!
     *  On
     * entry, x0 contains an initial estimate of the solution.  On
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
     *  @param jac        Matrix that will be used to store the jacobian
     *  @param num_newt_its Number of newton iterations taken 
     *  @param num_linear_solves Number of linear solves taken
     *  @param num_backtracks Number of backtracking steps taken
     *  @param loglevelInput  Input log level determines the amount of printing.
     *
     *
     *   @return  A positive value indicates a successful convergence
     *            -1  Failed convergence
     */
    int solve_nonlinear_problem(int SolnType, double* y_comm,double* ydot_comm, doublereal CJ,
				doublereal time_curr, SquareMatrix& jac,int &num_newt_its,
				int &num_linear_solves,	int &num_backtracks, int loglevelInput);

    //! Set the column scales
    void setColumnScales();

    //! Scale the matrix
    /*!
     *  @param jac              Jacobian
     *  @param y_comm           Current value of the solution vector
     *  @param ydot_comm        Current value of the time derivative of the solution vector
     *  @param time_curr        current value of the time
     */
    void scaleMatrix(SquareMatrix& jac, double* y_comm, double* ydot_comm, doublereal time_curr);

    //! Print solution norm contribution
    /*!
     *  Prints out the most important entries to the update to the solution vector for the current step
     *
     *   @param solnDelta0             Raw update vector for the current nonlinear step
     *   @param s0                     Norm of the vector solnDelta0
     *   @param solnDelta1             Raw update vector for the next solution value based on the old matrix
     *   @param s1                     Norm of the vector solnDelta1
     *   @param title                  title of the printout
     *   @param y0                     Old value of the solution
     *   @param y1                     New value of the solution after damping corrections
     *   @param damp                   Value of the damping factor
     *   @param num_entries            Number of entries to print out
     */
    void
    print_solnDelta_norm_contrib(const doublereal * const solnDelta0, const char * const s0,
				 const doublereal * const solnDelta1, const char * const s1,
				 const char * const title, const doublereal * const y0, const doublereal * const y1,
				 doublereal damp, int num_entries);

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
    void  getResidWts(doublereal * const residWts) const;

    //! Check to see if the nonlinear problem has converged
    /*!
     *
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
    void setAtol(const doublereal * const atol);

    //! Set the relative tolerances for the solution variables
    /*!
     *   Set the relative tolerances used in the calculation
     *
     *  @param rtol  single double
     */
    void setRtol(const doublereal rtol);

    //! Set the value of the maximum # of newton iterations
    /*!
     *  @param maxNewtIts   Maximum number of newton iterations
     */
    void setMaxNewtIts(const int maxNewtIts);

    //! Calculate the scaling factor for translating residual norms into 
    //! solution norms.
    void calcSolnToResNormVector();

    //! Calculate the steepest descent direction and the Cauchy Point where the quadratic formulation 
    //! of the nonlinear problem expects a minimum along the descent direction.
    /*!
     *  @param jac   Jacobian matrix: must be unfactored.
     *
     *  @return Returns the norm of the solution update
     */
    doublereal doCauchyPointSolve(SquareMatrix& jac);

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
     *  The actual residual decline in the newton direction determined by numerical differencing
     *
     *  This routine doesn't need to be called for the solution of the nonlinear problem.
     */
    void descentComparison(double time_curr ,double *ydot0, double *ydot1, const double *newtDir);

    void setupDoubleDogleg(double *newtDir);

    //! Change the global lambda coordinate into the (leg,alpha) coordinate for the double dogleg
    /*!
     * @param lambda Global value of the distance along the double dogleg
     * @param alpha  relative value along the particular leg
     *
     * @return Returns the leg number ( 0, 1, or 2).
     */
    int lambdaToLeg(const double lambda, double &alpha) const;

    int calcTrustIntersection(double trustVal, const double &lambda, double &alpha) const;

    int dampDogLeg(const doublereal time_curr, const double* y0, 
				  const doublereal *ydot0,  std::vector<doublereal> & step0,
				  double* const y1, double* const ydot1, double* step1,
				  double& s1, SquareMatrix& jac, int& loglevel, bool writetitle,
				  int& num_backtracks);

    int decideStep(const doublereal time_curr, int leg, double alpha, const double* y0, const doublereal *ydot0, 
		   std::vector<doublereal> & step0,
		   double* const y1, double* const ydot1,  int& loglevel, double trustDeltaOld);

    //! Calculated the expected residual along the double dogleg curve. 
    /*!
     *  @param leg 0, 1, or 2 representing the curves of the dogleg
     *  @param alpha  Relative distance along the particular curve.
     *
     *  @return Returns the expected value of the residual at that point according to the quadratic model.
     *          The residual at the newton point will always be zero.
     */
    double expectedResidLeg(int leg, doublereal alpha) const;

    void residualComparisonLeg(const double time_curr, const double *ydot0, const double *ydot1, const double *newtDir);

    //! Set the print level from the rootfinder
    /*!
     * 
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

    /*
     * -----------------------------------------------------------------------------------------------------------------
     *              MEMBER DATA
     * ------------------------------------------------------------------------------------------------
     */
 private:

    //! Pointer to the residual and jacobian evaluator for the 
    //! function
    /*!
     *   See ResidJacEval.h for an evaluator.
     */
    ResidJacEval *m_func;

    //! Solution type
    int solnType_;

    //! Local copy of the number of equations
    int neq_;
  
    //! Soln error weights
    std::vector<doublereal> m_ewt;

    //! Boolean indicating whether a manual delta bounds has been input.
    int m_manualDeltaStepSet;

    //! Soln Delta bounds magnitudes
    std::vector<doublereal> m_deltaStepMinimum;

    //! Value of the delta step magnitudes
    std::vector<doublereal> m_deltaStepMaximum;

    //! Vector containing the current solution vector within the nonlinear solver
    std::vector<doublereal> m_y_n;

    //! Vector containing the solution at the previous time step
    std::vector<doublereal> m_y_nm1;

    //! New value of the solution time derivative
    std::vector<doublereal> ydot_new;

    //! Vector of column scaling factors
    std::vector<doublereal> m_colScales;

    //! Weights for normalizing the values of the residuals
    /*!
     *  These are computed if row scaling, m_rowScaling, is turned on. They are calculated currently as the
     *  sum of the absolute values of the rows of the jacobian.
     */
    std::vector<doublereal> m_rowScales;

    //! Weights for normalizing the values of the residuals
    /*!
     *  These are computed if row scaling, m_rowScaling, is turned on. They are calculated currently as the
     *  sum of the absolute values  jacobian multiplied by the solution weight function
     */
    std::vector<doublereal> m_rowWtScales;

    //! Value of the residual for the nonlinear problem
    std::vector<doublereal> m_resid;

    //! Workspace of length neq_
    mutable std::vector<doublereal> m_wksp;

    /*****************************************************************************************
     *        INTERNAL WEIGHTS FOR TAKING SOLUTION NORMS
     ******************************************************************************************/
    //! Vector of residual weights
    /*!
     *   These are used to establish useful and informative weighted norms of the residual vector.
     */
    std::vector<doublereal> m_residWts;

    //! Norm of the residual at the start of each nonlinear iteration
    doublereal m_normResid0;

    //! Norm of the residual before damping
    doublereal m_normResidFRaw;

    //! Norm of the solution update created by the iteration in its raw, undamped form.
    doublereal m_normDeltaSoln_Newton;

    doublereal m_normDeltaSoln_CP;

    //! Norm of the residual for a trial calculation which may or may not be used
    doublereal m_normResidTrial;

    //! Vector of the norm
    doublereal m_normResidPoints[15];

    //! Boolean indicating whether we should scale the residual
    bool m_resid_scaled;

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
    int m_nfe;

    /***********************************************************************************************
     *             MATRIX INFORMATION
     **************************************************************************************/

    //! The type of column scaled used in the solution of the problem
    /*!
     *    If true then colScaling = m_ewt[]
     *    if false then colScaling = 1.0
     *  Currently, this is not part of the interface
     */
    bool m_colScaling;

    //! int indicating whether row scaling is turned on (1) or not (0)
    int m_rowScaling;

    //! Total number of linear solves
    int m_numTotalLinearSolves;

    //! Total number of newton iterations
    int m_numTotalNewtIts;

    //! Minimum number of newton iterations to use
    int m_min_newt_its;

    //! Maximum number of newton iterations
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

    //! Pointer containing the solution derivative at the previous time step
    doublereal *m_ydot_nm1;

    //! absolute tolerance in the solution unknown
    /*!
     * This is used to evaluating the weighting factor
     */
    std::vector<doublereal> atolk_;

    //! Determines the level of printing for each time step.
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
     */
    int m_print_flag;

    //! Scale factor for turning residual norms into solution norms
    double m_ScaleSolnNormToResNorm;

    //! Copy of the jacobian that doesn't get overwritten when the inverse is determined
    Cantera::SquareMatrix jacCopy_;

    /*********************************************************************************************
     *      VARIABLES ASSOCIATED WITH STEPS AND ASSOCIATED DOUBLE DOGLEG PARAMETERS
     *********************************************************************************************/

    //!  Steepest descent direction. This is also the distance to the Cauchy Point
    std::vector<doublereal> deltaX_CP_;

    //! Newton Step - This is the newton step determined from the straight Jacobian
    /*
     *  Newton step for the current step only
     */
    std::vector<doublereal> deltaX_Newton_;

    //! Expected value of the residual norm at the Cauchy point
    doublereal residNorm2Cauchy_;

    //! Residual dot Jd norm
    doublereal RJd_norm_;

    //! Value of lambda_ which is used to calculate the Cauchy point
    doublereal lambda_;

    //!  Jacobian times the Steepest descent direction. 
    std::vector<doublereal> Jd_;

    //! Vector of trust region values.
    std::vector<doublereal> deltaX_trust_;

    //! Current value of trust radius. This is used with trustDeltaX_ to 
    //! calculate the max step size.
    doublereal trustDelta_;

    //! Relative distance down the Newton step that the second dogleg starts
    doublereal Nuu_;

    doublereal dist_R0_;
    doublereal dist_R1_;
    doublereal dist_R2_;
    doublereal dist_Total_;
    doublereal JdJd_norm_;

    //! Norm of the Newton Step wrt trust region
    doublereal normTrust_Newton_;

    //! Norm of the Cauchy Step direction wrt trust region
    doublereal normTrust_CP_;



    /*******************************************************************************************
     *      OTHER COUNTERS
     *****************************************************************************************/


  public:
    //! Turn off printing of time
    /*!
     *  Necessary to do for test suites
     */
    static bool m_TurnOffTiming;

    //! Turn on or off printing of the Jacobian
    static bool s_print_NumJac;
  };

}

#endif
