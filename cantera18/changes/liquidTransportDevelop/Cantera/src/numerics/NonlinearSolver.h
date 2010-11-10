/**
 *  @file NonlinearSolve.h
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

namespace Cantera {

  // I think steady state is the only option I'm gunning for
#define  NSOLN_TYPE_PSEUDO_TIME_DEPENDENT 2
#define  NSOLN_TYPE_TIME_DEPENDENT 1
#define  NSOLN_TYPE_STEADY_STATE   0

#define NSOLN_JAC_NUM 1
#define NSOLN_JAC_ANAL 2



  //! Class that calculates the solution to a nonlinear system
  /*!
   *  UNDER CONSTRUCTION - do not use!!!!!!!!!!!!!!!!!!!!!!!!!
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
     * param y  vector of the current solution values
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
     */
    doublereal solnErrorNorm(const doublereal * const delta_y,  const char * title = 0, int printLargest = 0, 
			 const doublereal dampFactor = 1.0);

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
     */
    doublereal residErrorNorm(const doublereal * const resid, const char * title = 0, const int printLargest = 0,
			  const doublereal * const y = 0);

    //! Compute the current Residual
    /*!
     *   Compute the time dependent residual of 
     *   the set of equations.
     */
    // void doTDResidualCalc(const double time_curr, const int typeCalc,
    //			  const double * const y_curr, const double * const ydot_curr, int loglevel);

    //! Compute the current Residual
    /*!
     *   Compute the steady state residual of 
     *   the set of equations.
     */
    // void doSteadyResidualCalc(const double time_curr, const int typeCalc,
    //		      const double * const y_curr,  int loglevel);

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
     */
    void doResidualCalc(const doublereal time_curr, const int typeCalc, const doublereal * const y_curr, 
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
     *  @param timeCurrent    Current value of the time
     *  @param y_current      Current value of the solution
     *  @param ydot_current   Current value of the solution derivative.
     *
     */ 
    void doNewtonSolve(const doublereal time_curr, const doublereal * const y_curr, 
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
     *  @param deltaboundsMagnitudes          
     */
    void setDeltaBoundsMagnitudes(const doublereal * const deltaBoundsMagnitudes);

  
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
     */
    doublereal boundStep(const double* const  y, 
		     const double* const step0,  const int loglevel);


    //! Set bounds constraints for all variables in the problem
    /*!
     *  
     *   @param y_low_bounds  Vector of lower bounds
     *   @param y_high_bounds Vector of high bounds
     */
    void setBoundsConstraints(const doublereal * const y_low_bounds,
			      const doublereal * const y_high_bounds);

   
    //!   Internal function to calculate the time derivative at the new step
    /*!
     *   @param  order of the BDF method
     *   @param   y_curr current value of the solution
     *   @param   ydot_curr  Calculated value of the solution derivative that is consistent with y_curr
     */ 
    void calc_ydot(const int order, const doublereal * const y_curr, doublereal * const ydot_curr);

    //! Function called to evaluate the jacobian matrix and the curent
    //! residual vector.
    /*!
     *
     *
     */
    void beuler_jac(SquareMatrix &J, doublereal * const f,
		    doublereal time_curr, doublereal CJ, doublereal * const y,
		    doublereal * const ydot, int num_newt_its);


    //! Apply a filtering step
    /*!
     *  @param timeCurrent   Current value of the time
     *  @param y_current     current value of the solution
     *  @param ydot_current   Current value of the solution derivative.
     *
     *  @return Returns the norm of the value of the amount filtered
     */
    doublereal filterNewStep(const doublereal timeCurrent, doublereal * const y_current, doublereal * const ydot_current);
  
  
    //! Return the factor by which the undamped Newton step 'step0'
    //!  must be multiplied in order to keep the update within the bounds of an accurate jacobian.
    /*!
     *
     *  The idea behind these is that the Jacobian couldn't possibly be representative, if the
     *  variable is changed by a lot. (true for nonlinear systems, false for linear systems)
     *  Maximum increase in variable in any one newton iteration:
     *   factor of 1.5
     *  Maximum decrease in variable in any one newton iteration:
     *   factor of 2
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
     *    @param y1        
     */
    int dampStep(const doublereal time_curr, const double* y0, 
		 const doublereal *ydot0, const double* step0, 
		 double* const y1, double* const ydot1, double* step1,
		 double& s1, SquareMatrix& jac, 
		 int& loglevel, bool writetitle,
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
     *
     *
     *   @return  A positive value indicates a successful convergence
     *            -1  Failed convergence
     */
    int solve_nonlinear_problem(int SolnType, double* y_comm,
				double* ydot_comm, doublereal CJ,
				doublereal time_curr, 
				SquareMatrix& jac,
				int &num_newt_its,
				int &num_linear_solves,
				int &num_backtracks, 
				int loglevelInput);


    //! Set the column scales
    void setColumnScales();

    //! Scale the matrix
    /*!
     *
     */
    void scaleMatrix(SquareMatrix& jac, double* y_comm, double* ydot_comm, doublereal time_curr);


    //! Print solution norm contribution
    void
    print_solnDelta_norm_contrib(const doublereal * const solnDelta0,
				 const char * const s0,
				 const doublereal * const solnDelta1,
				 const char * const s1,
				 const char * const title,
				 const doublereal * const y0,
				 const doublereal * const y1,
				 doublereal damp,
				 int num_entries);

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
    int m_manualDeltaBoundsSet;

    //! Soln Delta bounds magnitudes
    std::vector<doublereal> m_deltaBoundsMagnitudes;

    //! Boolean indicating whether a manual delta steps have been input.
    int m_manualDeltaStepSet;
    
    std::vector<doublereal> m_deltaStepMagnitudes;

    std::vector<doublereal> m_y_n;
    std::vector<doublereal> m_y_nm1;


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
    std::vector<doublereal> m_wksp;

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
    doublereal m_normSolnFRaw;

    //! Norm of the residual for a trial calculation which may or may not be used
    doublereal m_normResidTrial;

    //! Vector of the norm
    doublereal m_normResidPoints[15];

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

    //! Boolean that turns on solution filtering
    int filterNewstep;

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

    int m_matrixConditioning;

    int m_order;

    //! value of the relative tolerance to use in solving the equation set
    doublereal rtol_;

    //! Base value of the absolute tolerance
    doublereal atolBase_;

    doublereal *  m_ydot_nm1;

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

  public:
    //! Turn off printing of time
    /*!
     *  Necessary to do for test suites
     */
    static bool m_TurnOffTiming;

    // Turn on or off printing of the Jacobian
    /*!
     *
     */
    static bool s_print_NumJac;
  };

}

#endif
