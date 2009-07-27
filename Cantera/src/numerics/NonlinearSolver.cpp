/**
 *
 *  @file NonlinearSolver.cpp
 *
 *  Damped Newton solver for 0D and 1D problems
 */

/*
 *  $Date: 2009/03/27 21:32:33 $
 *  $Revision: 1.5 $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include <limits>

#include "SquareMatrix.h"
#include "NonlinearSolver.h"

#include "clockWC.h"
#include "vec_functions.h"
#include <ctime>

#include "mdp_allo.h"
#include <cfloat>

extern void print_line(const char *, int);

#include <vector>
#include <cstdio>
#include <cmath>


#ifndef MAX
#define MAX(x,y)    (( (x) > (y) ) ? (x) : (y))
#define MIN(x,y)    (( (x) < (y) ) ? (x) : (y))
#endif

using namespace std;

namespace Cantera {



  //-----------------------------------------------------------
  //                 Constants
  //-----------------------------------------------------------

  const double DampFactor = 4;
  const int NDAMP = 7;

  //-----------------------------------------------------------
  //                 Static Functions
  //-----------------------------------------------------------

  static void print_line(const char *str, int n)  {
    for (int i = 0; i < n; i++) {
      printf("%s", str);
    }
    printf("\n");
  }

  // Default constructor
  /*
   * @param func   Residual and jacobian evaluator function object
   */
  NonlinearSolver::NonlinearSolver(ResidJacEval *func) :
    m_func(func),
    neq_(0),
    delta_t_n(-1.0),
    m_nfe(0),
    m_colScaling(0),
    m_rowScaling(0),
    m_numTotalLinearSolves(0),
    m_numTotalNewtIts(0),
    m_min_newt_its(0),
    filterNewstep(0),
    time_n(0.0),
    m_matrixConditioning(0),
    m_order(1),
    rtol_(1.0E-3),
    atolBase_(1.0E-10)
  {
    neq_ = m_func->nEquations();

    m_ewt.resize(neq_, rtol_);
    m_y_n.resize(neq_, 0.0);
    m_y_nm1.resize(neq_, 0.0);
    m_colScales.resize(neq_, 1.0);
    m_rowScales.resize(neq_, 1.0);
    m_resid.resize(neq_, 0.0);
    atolk_.resize(neq_, atolBase_);
	doublereal hb = std::numeric_limits<double>::max();
    m_y_high_bounds.resize(neq_, hb);
    m_y_low_bounds.resize(neq_, -hb);    

    for (int i = 0; i < neq_; i++) {
      atolk_[i] = atolBase_;
      m_ewt[i] = atolk_[i];
    }
  }

  NonlinearSolver::NonlinearSolver(const NonlinearSolver &right) {
    *this =operator=(right);
  }


  NonlinearSolver::~NonlinearSolver() {
  }

  NonlinearSolver& NonlinearSolver::operator=(const NonlinearSolver &right) {
    if (this == &right) {
      return *this;
    }
    // rely on the ResidJacEval duplMyselfAsresidJacEval() function to
    // create a deep copy
    m_func = right.m_func->duplMyselfAsResidJacEval();

    neq_                       = right.neq_;
    m_ewt                      = right.m_ewt;
    m_y_n                      = right.m_y_n;
    m_y_nm1                    = right.m_y_nm1;
    m_colScales                = right.m_colScales;
    m_rowScales                = right.m_rowScales;
    m_resid                    = right.m_resid;
    m_y_high_bounds            = right.m_y_high_bounds;
    m_y_low_bounds             = right.m_y_low_bounds;
    delta_t_n                  = right.delta_t_n;
    m_nfe                      = right.m_nfe;
    m_colScaling               = right.m_colScaling;
    m_rowScaling               = right.m_rowScaling;
    m_numTotalLinearSolves     = right.m_numTotalLinearSolves;
    m_numTotalNewtIts          = right.m_numTotalNewtIts;
    m_min_newt_its             = right.m_min_newt_its;
    filterNewstep              = right.filterNewstep;
    time_n                     = right.time_n;
    m_matrixConditioning       = right.m_matrixConditioning;
    m_order                    = right.m_order;
    rtol_                      = right.rtol_;
    atolBase_                  = right.atolBase_;
    atolk_                     = right.atolk_;

    return *this;
  }

  // Create solution weights for convergence criteria
  /*
   *  We create soln weights from the following formula
   *
   *  wt[i] = rtol * abs(y[i]) + atol[i]
   *
   *  The program always assumes that atol is specific
   *  to the solution component
   *
   * param y  vector of the current solution values
   */
  void NonlinearSolver::createSolnWeights(const double * const y) {
    for (int i = 0; i < neq_; i++) {
      m_ewt[i] = rtol_ * fabs(y[i]) + atolk_[i];
    }
  }

  // set bounds constraints for all variables in the problem
  /*
   *  
   *   @param y_low_bounds  Vector of lower bounds
   *   @param y_high_bounds Vector of high bounds
   */
  void NonlinearSolver::setBoundsConstraints(const double * const y_low_bounds,
					     const double * const y_high_bounds) {
    for (int i = 0; i < neq_; i++) {
      m_y_low_bounds[i]  = y_low_bounds[i];
      m_y_high_bounds[i] = y_high_bounds[i];
    }
  }

  /**
   * L2 Norm of a delta in the solution
   *
   *  The second argument has a default of false. However,
   *  if true, then a table of the largest values is printed
   *  out to standard output.
   */
  double NonlinearSolver::solnErrorNorm(const double * const delta_y, 
					bool printLargest)
  {
    int    i;
    double sum_norm = 0.0, error;
    for (i = 0; i < neq_; i++) {
      error     = delta_y[i] / m_ewt[i];
      sum_norm += (error * error);
    }
    sum_norm = sqrt(sum_norm / neq_); 
    if (printLargest) {
      const int num_entries = 8;
      double dmax1, normContrib;
      int j;
      int *imax = mdp::mdp_alloc_int_1(num_entries, -1);
      printf("\t\tPrintout of Largest Contributors to norm "
	     "of value (%g)\n", sum_norm);
      printf("\t\t         I    ysoln  deltaY  weightY  "
	     "Error_Norm**2\n");
      printf("\t\t   "); print_line("-", 80);
      for (int jnum = 0; jnum < num_entries; jnum++) {
	dmax1 = -1.0;
	for (i = 0; i < neq_; i++) {
	  bool used = false;
	  for (j = 0; j < jnum; j++) {
	    if (imax[j] == i) used = true;
	  }
	  if (!used) {
	    error = delta_y[i] / m_ewt[i];
	    normContrib = sqrt(error * error);
	    if (normContrib > dmax1) {
	      imax[jnum] = i;
	      dmax1 = normContrib;
	    }
	  }
	}
	i = imax[jnum];
	if (i >= 0) {
	  printf("\t\t   %4d %12.4e %12.4e %12.4e %12.4e\n",
		 i, m_y_n[i], delta_y[i], m_ewt[i], dmax1);
	}	  
      }
      printf("\t\t   "); print_line("-", 80);
      mdp::mdp_safe_free((void **) &imax);
    }
    return sum_norm;
  }

  /**
   * L2 Norm of the residual
   *
   *  The second argument has a default of false. However,
   *  if true, then a table of the largest values is printed
   *  out to standard output.
   */
  double NonlinearSolver::residErrorNorm(const double * const resid,
					 bool printLargest)
  {
    int    i;
    double sum_norm = 0.0, error;
    for (i = 0; i < neq_; i++) {
      error     = resid[i] / m_rowScales[i];
      sum_norm += (error * error);
    }
    sum_norm = sqrt(sum_norm / neq_); 
    if (printLargest) {
      const int num_entries = 8;
      double dmax1, normContrib;
      int j;
      int *imax = mdp::mdp_alloc_int_1(num_entries, -1);
      printf("\t\tPrintout of Largest Contributors to norm "
	     "of Residual (%g)\n", sum_norm);
      printf("\t\t         I    resid  rowScale  weightN  "
	     "Error_Norm**2\n");
      printf("\t\t   "); print_line("-", 80);
      for (int jnum = 0; jnum < num_entries; jnum++) {
	dmax1 = -1.0;
	for (i = 0; i < neq_; i++) {
	  bool used = false;
	  for (j = 0; j < jnum; j++) {
	    if (imax[j] == i) used = true;
	  }
	  if (!used) {
	    error = resid[i] / m_rowScales[i];
	    normContrib = sqrt(error * error);
	    if (normContrib > dmax1) {
	      imax[jnum] = i;
	      dmax1 = normContrib;
	    }
	  }
	}
	i = imax[jnum];
	if (i >= 0) {
	  printf("\t\t   %4d %12.4e %12.4e %12.4e \n",
		 i,  resid[i], m_rowScales[i], normContrib);
	}	  
      }
      printf("\t\t   "); print_line("-", 80);
      mdp::mdp_safe_free((void **) &imax);
    }
    return sum_norm;
  }


  /**
   * setColumnScales():
   *
   * Set the column scaling vector at the current time
   */
  void NonlinearSolver::setColumnScales() {
    m_func->calcSolnScales(time_n, DATA_PTR(m_y_n), DATA_PTR(m_y_nm1),
			   DATA_PTR(m_colScales));
  }


  void NonlinearSolver::doResidualCalc(const double time_curr, const int typeCalc,
				       const double * const y_curr, 
				       const double * const ydot_curr, double* const residual,
				       int loglevel)
  {


    // Calculate the current residual
    //   Put the current residual into the vector, delta_y[]
    //   We need to pull this out of this function and carry it in.
    m_func->evalResidNJ(time_curr, delta_t_n, y_curr, ydot_curr, residual);
    m_nfe++;
  }

 
  // Compute the undamped Newton step
  /*
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
   */ 
  void NonlinearSolver::doNewtonSolve(const double time_curr, const double * const y_curr, 
				      const double * const ydot_curr, double* const delta_y,
				      SquareMatrix& jac, int loglevel)
  {	 
    int irow, jcol;


    //! multiply the residual by -1
    for (int n = 0; n < neq_; n++) {
      delta_y[n] = -delta_y[n];
    }


    /*
     * Column scaling -> We scale the columns of the Jacobian
     * by the nominal important change in the solution vector
     */
    if (m_colScaling) {
      if (!jac.m_factored) {
	/*
	 * Go get new scales -> Took this out of this inner loop. 
	 * Needs to be done at a larger scale.
	 */
	// setColumnScales();

	/*
	 * Scale the new Jacobian
	 */
	double *jptr = &(*(jac.begin()));
	for (jcol = 0; jcol < neq_; jcol++) {
	  for (irow = 0; irow < neq_; irow++) {
	    *jptr *= m_colScales[jcol];
	    jptr++;
	  }
	}
      }	  
    }

    //  if (m_matrixConditioning) {
    //  if (jac.m_factored) {
    //	m_func->matrixConditioning(0, neq_, delta_y);
    // } else {
    //double *jptr = &(*(jac.begin()));
    //	m_func->matrixConditioning(jptr, neq_, delta_y);
    // }
    //}

    /*
     * row sum scaling -> Note, this is an unequivical success
     *      at keeping the small numbers well balanced and
     *      nonnegative.
     */
    if (m_rowScaling) {
      if (! jac.m_factored) {
	/*
	 * Ok, this is ugly. jac.begin() returns an vector<double> iterator
	 * to the first data location.
	 * Then &(*()) reverts it to a double *.
	 */
	double *jptr = &(*(jac.begin()));
	for (irow = 0; irow < neq_; irow++) {
	  m_rowScales[irow] = 0.0;
	}
	for (jcol = 0; jcol < neq_; jcol++) {
	  for (irow = 0; irow < neq_; irow++) {
	    m_rowScales[irow] += fabs(*jptr);
	    jptr++;
	  }
	}
	
	jptr = &(*(jac.begin()));
	for (jcol = 0; jcol < neq_; jcol++) {
	  for (irow = 0; irow < neq_; irow++) {
	    *jptr /= m_rowScales[irow];
	    jptr++;
	  }
	}
      }
      for (irow = 0; irow < neq_; irow++) {
	delta_y[irow] /= m_rowScales[irow];
      }
    }


    /*
     * Solve the system -> This also involves inverting the
     * matrix
     */
    (void) jac.solve(delta_y);


    /*
     * reverse the column scaling if there was any.
     */
    if (m_colScaling) {
      for (irow = 0; irow < neq_; irow++) {
	delta_y[irow] *= m_colScales[irow];
      }
    }
	
#ifdef DEBUG_JAC
    if (printJacContributions) {
      for (int iNum = 0; iNum < numRows; iNum++) {
	if (iNum > 0) focusRow++;
	double dsum = 0.0;
	vector_fp& Jdata = jacBack.data();
	double dRow = Jdata[neq_ * focusRow + focusRow];
	printf("\n Details on delta_Y for row %d \n", focusRow);
	printf("  Value before = %15.5e, delta = %15.5e,"
	       "value after = %15.5e\n", y_curr[focusRow], 
	       delta_y[focusRow],
	       y_curr[focusRow] +  delta_y[focusRow]);
	if (!freshJac) {
	  printf("    Old Jacobian\n");
	}
	printf("     col          delta_y            aij     "
	       "contrib   \n");
	printf("--------------------------------------------------"
	       "---------------------------------------------\n");
	printf(" Res(%d) %15.5e  %15.5e  %15.5e  (Res = %g)\n",
	       focusRow, delta_y[focusRow],
	       dRow, RRow[iNum] / dRow, RRow[iNum]);
	dsum +=  RRow[iNum] / dRow;
	for (int ii = 0; ii < neq_; ii++) {
	  if (ii != focusRow) {
	    double aij =  Jdata[neq_ * ii + focusRow];
	    double contrib = aij * delta_y[ii] * (-1.0) / dRow;
	    dsum += contrib;
	    if (fabs(contrib) > Pcutoff) {
	      printf("%6d  %15.5e  %15.5e  %15.5e\n", ii,
		     delta_y[ii]  , aij, contrib); 
	    }
	  }
	}
	printf("--------------------------------------------------"
	       "---------------------------------------------\n");
	printf("        %15.5e                   %15.5e\n",
	       delta_y[focusRow], dsum);
      }
    }

#endif
	
    m_numTotalLinearSolves++;
  }

  /**************************************************************************
   *
   * boundStep():
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
  double NonlinearSolver::boundStep(const double* const y, 
				    const double* const step0, const int loglevel) {
    int i, i_lower = -1, i_fbounds, ifbd = 0, i_fbd = 0;
    double fbound = 1.0, f_bounds = 1.0, f_delta_bounds = 1.0;
    double ff, y_new, ff_alt;
    
    for (i = 0; i < neq_; i++) {
      y_new = y[i] + step0[i];
      /*
       * Force the step to only take 80% a step towards the lower bounds
       */
      if (step0[i] < 0.0) {
	if (y_new < m_y_low_bounds[i]) {
	  double legalDelta = 0.8*(m_y_low_bounds[i] - y[i]);
	  ff = legalDelta / step0[i];
	  if (ff < f_bounds) {
	    f_bounds = ff;
	    i_lower = i;
	  }
	}
      }
      /*
       * Force the step to only take 80% a step towards the high bounds
       */
      if (step0[i] > 0.0) {
	if (y_new > m_y_high_bounds[i]) {
	  double legalDelta = 0.8*(m_y_high_bounds[i] - y[i]);
	  ff = legalDelta / step0[i];
	  if (ff < f_bounds) {
	    f_bounds = ff;
	    i_lower = i;
	  }
	}
      }
      /**
       * Now do a delta bounds
       * Increase variables by a factor of 2 only
       * decrease variables by a factor of 5 only
       */
      ff = 1.0;
      if ((fabs(y_new) > 2.0 * fabs(y[i])) && 
	  (fabs(y_new-y[i]) > m_ewt[i])) {
	ff = fabs(y[i]/(y_new - y[i]));
	ff_alt = fabs(m_ewt[i] / (y_new - y[i]));
	ff = MAX(ff, ff_alt);
	ifbd = 1;
      }
      if ((fabs(5.0 * y_new) < fabs(y[i])) &&
	  (fabs(y_new - y[i]) > m_ewt[i])) {
	ff = y[i]/(y_new-y[i]) * (1.0 - 5.0)/5.0;
	ff_alt = fabs(m_ewt[i] / (y_new - y[i]));
	ff = MAX(ff, ff_alt);
	ifbd = 0;
      }
      if (ff < f_delta_bounds) {
	f_delta_bounds = ff;
	i_fbounds = i;
	i_fbd = ifbd;
      }
      f_delta_bounds = MIN(f_delta_bounds, ff);
    }
    fbound = MIN(f_bounds, f_delta_bounds);
    /*
     * Report on any corrections
     */
    if (loglevel > 1) {
      if (fbound != 1.0) {
	if (f_bounds < f_delta_bounds) {
	  printf("\t\tboundStep: Variable %d causing bounds "
		 "damping of %g\n",
		 i_lower, f_bounds);
	} else {
	  if (ifbd) {
	    printf("\t\tboundStep: Decrease of Variable %d causing "
		   "delta damping of %g\n",
		   i_fbd, f_delta_bounds);
	  } else {
	    printf("\t\tboundStep: Increase of variable %d causing"
		   "delta damping of %g\n",
		   i_fbd, f_delta_bounds);
	  }
	}
      }
    }
    //return fbound;
    return 1.0;
  }

  /**************************************************************************
   *
   * dampStep():
   *
   * On entry, step0 must contain an undamped Newton step to the
   * current solution y0. This method attempts to find a damping coefficient
   * such that the next undamped step would have a norm smaller than
   * that of step0. If successful, the new solution after taking the
   * damped step is returned in y1, and the undamped step at y1 is
   * returned in step1.
   */
  int NonlinearSolver::dampStep(const double time_curr, const double* y0, 
				const double *ydot0, const double* step0, 
				double* const y1, double* const ydot1, double* step1,
				double& s1, SquareMatrix& jac, 
				int& loglevel, bool writetitle,
				int& num_backtracks) {
    
          
    // Compute the weighted norm of the undamped step size step0
    double s0 = solnErrorNorm(step0);

    // Compute the multiplier to keep all components in bounds
    // A value of one indicates that there is no limitation
    // on the current step size in the nonlinear method due to
    // bounds constraints (either negative values of delta
    // bounds constraints.
    double fbound = boundStep(y0, step0, loglevel);

    // if fbound is very small, then y0 is already close to the
    // boundary and step0 points out of the allowed domain. In
    // this case, the Newton algorithm fails, so return an error
    // condition.
    if (fbound < 1.e-10) {
      if (loglevel > 1) printf("\t\t\tdampStep: At limits.\n");
      return -3;
    }

    //--------------------------------------------
    //           Attempt damped step
    //-------------------------------------------- 

    // damping coefficient starts at 1.0
    double damp = 1.0;
    int j, m;
    double ff;
    num_backtracks = 0;
    for (m = 0; m < NDAMP; m++) {

      ff = fbound*damp;

      // step the solution by the damped step size
      /*
       * Whenever we update the solution, we must also always
       * update the time derivative.
       */
      for (j = 0; j < neq_; j++) {
	y1[j] = y0[j] + ff * step0[j];
      }
      calc_ydot(m_order, y1, ydot1);

      doResidualCalc(time_curr, NSOLN_TYPE_STEADY_STATE, y1, ydot1, step1, loglevel);
	  
      // compute the next undamped step, step1[], that would result 
      // if y1[] were accepted.

      doNewtonSolve(time_curr, y1, ydot1, step1, jac, loglevel);

      // compute the weighted norm of step1
      s1 = solnErrorNorm(step1);

      // write log information
      if (loglevel > 3) {
	print_solnDelta_norm_contrib((const double *) step0, 
				     "DeltaSolnTrial",
				     (const double *) step1,
				     "DeltaSolnTrialTest",
				     "dampNewt: Important Entries for "
				     "Weighted Soln Updates:",
				     y0, y1, ff, 5);
      }
      if (loglevel > 1) {
	printf("\t\t\tdampNewt: s0 = %g, s1 = %g, fbound = %g,"
	       "damp = %g\n",  s0, s1, fbound, damp);
      }


      // if the norm of s1 is less than the norm of s0, then
      // accept this damping coefficient. Also accept it if this
      // step would result in a converged solution. Otherwise,
      // decrease the damping coefficient and try again.
	  
      if (s1 < 1.0E-5 || s1 < s0) {
	if (loglevel > 2) {
	  if (s1 > s0) {
	    if (s1 > 1.0) {
	      printf("\t\t\tdampStep: current trial step and damping"
		     " coefficient accepted because test step < 1\n");
	      printf("\t\t\t          s1 = %g, s0 = %g\n", s1, s0);
	    }
	  }
	}
	break;
      } else {
	if (loglevel > 1) {
	  printf("\t\t\tdampStep: current step rejected: (s1 = %g > "
		 "s0 = %g)", s1, s0);
	  if (m < (NDAMP-1)) {
	    printf(" Decreasing damping factor and retrying");
	  } else {
	    printf(" Giving up!!!");
	  }
	  printf("\n");
	}
      }
      num_backtracks++;
      damp /= DampFactor;
    }

    // If a damping coefficient was found, return 1 if the
    // solution after stepping by the damped step would represent
    // a converged solution, and return 0 otherwise. If no damping
    // coefficient could be found, return -2.
    if (m < NDAMP) {
      if (s1 > 1.0) return 0;
      else return 1;
    } else {
      if (s1 < 0.5 && (s0 < 0.5)) return 1;
      if (s1 < 1.0) return 0;
      return -2;
    }
  }

  /**
   *
   * solve_nonlinear_problem():
   *
   * Find the solution to F(X) = 0 by damped Newton iteration. On
   * entry, x0 contains an initial estimate of the solution.  On
   * successful return, x1 contains the converged solution.
   *
   * SolnType = TRANSIENT -> we will assume we are relaxing a transient
   *        equation system for now. Will make it more general later,
   *        if an application comes up.
   * 
   */
  int NonlinearSolver::solve_nonlinear_problem(int SolnType, double* y_comm,
					       double* ydot_comm, double CJ,
					       double time_curr, 
					       SquareMatrix& jac,
					       int &num_newt_its,
					       int &num_linear_solves,
					       int &num_backtracks, 
					       int loglevelInput)
  {
    clockWC wc;

    bool m_residCurrent = false;
    int m = 0;
    bool forceNewJac = false;
    double s1=1.e30;

    std::vector<doublereal> y_curr(neq_, 0.0); 
    std::vector<doublereal> ydot_curr(neq_, 0.0);
    std::vector<doublereal> stp(neq_, 0.0);
    std::vector<doublereal> stp1(neq_, 0.0);

    std::vector<doublereal> y_new(neq_, 0.0);
    std::vector<doublereal> ydot_new(neq_, 0.0);

    mdp::mdp_copy_dbl_1(DATA_PTR(y_curr), y_comm, neq_);
    // copyn((size_t)neq_, y_comm,    y_curr);
    mdp::mdp_copy_dbl_1(DATA_PTR(ydot_curr), ydot_comm, neq_);

    

    bool frst = true;
    num_newt_its = 0;
    num_linear_solves = - m_numTotalLinearSolves;
    num_backtracks = 0;
    int i_backtracks;
    int loglevel = loglevelInput;

    while (1 > 0) {

      /*
       * Increment Newton Solve counter
       */
      m_numTotalNewtIts++;
      num_newt_its++;


      if (loglevel > 1) {
	printf("\t\tSolve_Nonlinear_Problem: iteration %d:\n",
	       num_newt_its);
      }

      // Check whether the Jacobian should be re-evaluated.
            
      forceNewJac = true;
            
      if (forceNewJac) {
	if (loglevel > 1) {
	  printf("\t\t\tGetting a new Jacobian and solving system\n");
	}
	beuler_jac(jac, DATA_PTR(m_resid), time_curr, CJ,  DATA_PTR(y_curr), DATA_PTR(ydot_curr),
		   num_newt_its);
	m_residCurrent = true;
      } else {
	if (loglevel > 1) {
	  printf("\t\t\tSolving system with old jacobian\n");
	}
	m_residCurrent = false;
      }
      /*
       * Go get new scales
       */
      setColumnScales();


      doResidualCalc(time_curr, NSOLN_TYPE_STEADY_STATE,
		     DATA_PTR(y_curr), DATA_PTR(ydot_curr), DATA_PTR(stp), loglevel);

      // compute the undamped Newton step
      doNewtonSolve(time_curr, DATA_PTR(y_curr), DATA_PTR(ydot_curr), DATA_PTR(stp),
		    jac, loglevel);
	  
      // damp the Newton step
      m = dampStep(time_curr, DATA_PTR(y_curr), DATA_PTR(ydot_curr), 
		   DATA_PTR(stp), DATA_PTR(y_new), DATA_PTR(ydot_new), 
		   DATA_PTR(stp1), s1, jac, loglevel, frst, i_backtracks);
      frst = false;
      num_backtracks += i_backtracks;

      /*
       * Impose the minimum number of newton iterations critera
       */
      if (num_newt_its < m_min_newt_its) {
	if (m == 1) m = 0;
      }
      /*
       * Impose max newton iteration
       */
      if (num_newt_its > 20) {
	m = -1;
	if (loglevel > 1) {
	  printf("\t\t\tDampnewton unsuccessful (max newts exceeded) sfinal = %g\n", s1);
	}
      }

      if (loglevel > 1) {
	if (m == 1) {
	  printf("\t\t\tDampNewton iteration successful, nonlin "
		 "converged sfinal = %g\n", s1);
	} else if (m == 0) {
	  printf("\t\t\tDampNewton iteration successful, get new"
		 "direction, sfinal = %g\n", s1);
	} else {
	  printf("\t\t\tDampnewton unsuccessful sfinal = %g\n", s1);
	}
      }

      // If we are converged, then let's use the best solution possible
      // for an end result. We did a resolve in dampStep(). Let's update
      // the solution to reflect that.
      // HKM 5/16 -> Took this out, since if the last step was a 
      //             damped step, then adding stp1[j] is undamped, and
      //             may lead to oscillations. It kind of defeats the
      //             purpose of dampStep() anyway.
      // if (m == 1) {
      //  for (int j = 0; j < neq_; j++) {
      //   y_new[j] += stp1[j];
      //                HKM setting intermediate y's to zero was a tossup.
      //                    slightly different, equivalent results
      // #ifdef DEBUG_HKM
      //	      y_new[j] = MAX(0.0, y_new[j]);
      // #endif
      //  }
      // }
	  
      bool m_filterIntermediate = false;
      if (m_filterIntermediate) {
	if (m == 0) {
	  (void) filterNewStep(time_n, DATA_PTR(y_new), DATA_PTR(ydot_new));
	}
      }
      // Exchange new for curr solutions
      if (m == 0 || m == 1) {
	mdp::mdp_copy_dbl_1(DATA_PTR(y_curr), DATA_PTR(y_new), neq_);
	calc_ydot(m_order, DATA_PTR(y_curr), DATA_PTR(ydot_curr));
      }

      // convergence
      if (m == 1) goto done;

      // If dampStep fails, first try a new Jacobian if an old
      // one was being used. If it was a new Jacobian, then
      // return -1 to signify failure.
      else if (m < 0) {
	goto done;
      }
    }

  done:
    mdp::mdp_copy_dbl_1(y_comm, DATA_PTR(y_curr), neq_);
    mdp::mdp_copy_dbl_1(ydot_comm, DATA_PTR(ydot_curr), neq_);
   
 
    num_linear_solves += m_numTotalLinearSolves;
 
    double time_elapsed =  wc.secondsWC();
    if (loglevel > 1) {
      if (m == 1) {
	printf("\t\tNonlinear problem solved successfully in "
	       "%d its, time elapsed = %g sec\n",
	       num_newt_its, time_elapsed);
      }
    }
    return m;
  }
 
  /***************************************************************8
   *
   *
   */
  void NonlinearSolver::
  print_solnDelta_norm_contrib(const double * const solnDelta0,
			       const char * const s0,
			       const double * const solnDelta1,
			       const char * const s1,
			       const char * const title,
			       const double * const y0,
			       const double * const y1,
			       double damp,
			       int num_entries) {
    int i, j, jnum;
    bool used;
    double dmax0, dmax1, error, rel_norm;
    printf("\t\t%s currentDamp = %g\n", title, damp);
    printf("\t\t         I  ysoln %10s ysolnTrial "
	   "%10s weight relSoln0 relSoln1\n", s0, s1);
    int *imax = mdp::mdp_alloc_int_1(num_entries, -1);
    printf("\t\t   "); print_line("-", 90);
    for (jnum = 0; jnum < num_entries; jnum++) {
      dmax1 = -1.0;
      for (i = 0; i < neq_; i++) {
	used = false;
	for (j = 0; j < jnum; j++) {
	  if (imax[j] == i) used = true;
	}
	if (!used) {
	  error     = solnDelta0[i] /  m_ewt[i];
	  rel_norm = sqrt(error * error);
	  error     = solnDelta1[i] /  m_ewt[i];
	  rel_norm += sqrt(error * error);
	  if (rel_norm > dmax1) {
	    imax[jnum] = i;
	    dmax1 = rel_norm;
	  }
	}
      }
      if (imax[jnum] >= 0) {
	i = imax[jnum];
	error = solnDelta0[i] /  m_ewt[i];
	dmax0 = sqrt(error * error);
	error = solnDelta1[i] /  m_ewt[i];
	dmax1 = sqrt(error * error);
	printf("\t\t   %4d %12.4e %12.4e %12.4e  %12.4e "
	       "%12.4e %12.4e %12.4e\n",
	       i, y0[i], solnDelta0[i], y1[i],
	       solnDelta1[i], m_ewt[i], dmax0, dmax1);
      }
    }
    printf("\t\t   "); print_line("-", 90);
    mdp::mdp_safe_free((void **) &imax);
  }



}

