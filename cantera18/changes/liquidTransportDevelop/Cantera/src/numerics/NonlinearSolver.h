/**
 *  @file NonlinearSolve.h
 *    Class that calculates the solution to a nonlinear, dense, set
 *    of equations (see \ref numerics
 *    and class \link Cantera::NonlinearSolver NonlinearSolver\endlink).
 */

/*
 *  $Date: 2009/02/14 18:02:22 $
 *  $Revision: 1.4 $
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

#define  NSOLN_TYPE_PSEUDO_TIME_DEPENDENT 2
#define  NSOLN_TYPE_TIME_DEPENDENT 1
#define  NSOLN_TYPE_STEADY_STATE   0

  //! Class that calculates the solution to a nonlinear system
  /*!
   *
   *  @ingroup numerics
   */
  class NonlinearSolver {

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
    void createSolnWeights(const double * const y);


    //!  L2 norm of the delta of the solution vector
    /*!
     *  calculate the norm of the solution vector. This will
     *  involve the column scaling of the matrix
     *
     *    The second argument has a default of false. However,
     *    if true, then a table of the largest values is printed
     *    out to standard output.
     */
    double solnErrorNorm(const double * const delta_y, 
			 bool printLargest = false);

    //! L2 norm of the residual of the equation system
    /*!
     * Calculate the norm of the residual vector. This may
     * involve using the row sum scaling from the matrix problem. 
     *
     *  The second argument has a default of false. However,
     *  if true, then a table of the largest values is printed
     *  out to standard output.
     */
    double residErrorNorm(const double * const resid, 
			    bool printLargest = false);

    //! Compute the current Residual
    /*!
     *   Compute the time dependent residual of 
     *   the set of equations.
     */
    void doTDResidualCalc(const double time_curr, const int typeCalc,
			  const double * const y_curr, 
			const double * const ydot_curr, double* const residual,
			int loglevel);

    //! Compute the current Residual
    /*!
     *   Compute the steady state residual of 
     *   the set of equations.
     */
    void doSteadyResidualCalc(const double time_curr, const int typeCalc,
			      const double * const y_curr, 
			      double* const residual, int loglevel);

    void doResidualCalc(const double time_curr, const int typeCalc, 
			const double * const y_curr, 
			const double * const ydot_curr, double* const residual,
			int loglevel);

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
     */ 
    void doNewtonSolve(const double time_curr, const double * const y_curr, 
		       const double * const ydot_curr, double* const delta_y,
		       SquareMatrix& jac, int loglevel);
  
    //!
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
    double boundStep(const double* const  y, 
		     const double* const step0,  const int loglevel);


    //! set bounds constraints for all variables in the problem
    /*!
     *  
     *   @param y_low_bounds  Vector of lower bounds
     *   @param y_high_bounds Vector of high bounds
     */
    void setBoundsConstraints(const double * const y_low_bounds,
			      const double * const y_high_bounds);

    /**
     * Internal function to calculate the predicted solution
     * at a time step.
     */
    void calc_y_pred(int);

    /**
     * Internal function to calculate the time derivative at the
     * new step
     */
    void calc_ydot(int order, double * const y_curr, double * const ydot_curr);
    
    void beuler_jac(SquareMatrix &, double * const,
		    double, double, double * const, double * const, int);


    double filterNewStep(double, double *, double *);

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
    int dampStep(const double time_curr, const double* y0, 
		 const double *ydot0, const double* step0, 
		 double* const y1, double* const ydot1, double* step1,
		 double& s1, SquareMatrix& jac, 
		 int& loglevel, bool writetitle,
		 int& num_backtracks);

   

    // Compute the weighted norm of the undamped step size step0

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
     */
    int solve_nonlinear_problem(int SolnType, double* y_comm,
				double* ydot_comm, double CJ,
				double time_curr, 
				SquareMatrix& jac,
				int &num_newt_its,
				int &num_linear_solves,
				int &num_backtracks, 
				int loglevelInput);


    void setColumnScales();

    void
    print_solnDelta_norm_contrib(const double * const solnDelta0,
				 const char * const s0,
				 const double * const solnDelta1,
				 const char * const s1,
				 const char * const title,
				 const double * const y0,
				 const double * const y1,
				 double damp,
				 int num_entries);



    //! Pointer to the residual and jacobian evaluator for the 
    //! function
    /*!
     *   See ResidJacEval.h for an evaluator.
     */
    ResidJacEval *m_func;

    //! Local copy of the number of equations
    int neq_;
  
    //! Soln error weights
    std::vector<doublereal> m_ewt;

    std::vector<doublereal> m_y_n;
    std::vector<doublereal> m_y_nm1;
    std::vector<doublereal> m_colScales;

    //! Weights for normalizing the values of the residuals
    
    std::vector<doublereal> m_rowScales;

    std::vector<doublereal> m_resid;

    //! Bounds vector for each species
    std::vector<doublereal> m_y_high_bounds;

    //! Lower bounds vector for each species
    std::vector<doublereal> m_y_low_bounds;

    double delta_t_n;

    //! Counter for the total number of function evaluations
    int m_nfe;

    //! The type of column scaled used in the solution of the problem
    bool m_colScaling;

    //! int indicating whether row scaling is turned on (1) or not (0)
    int m_rowScaling;

    int m_numTotalLinearSolves;

    int m_numTotalNewtIts;

    int m_min_newt_its;

    int filterNewstep;

    //! Current system time
    /*!
     *  Note, we assume even for steady state problems that the residual
     *  is a function of a system time. 
     */
    double time_n;

    int m_matrixConditioning;

    int m_order;

    doublereal rtol_;

    doublereal atolBase_;

    std::vector<doublereal> atolk_;
  };

}

#endif
