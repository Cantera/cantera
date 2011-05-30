/**
 *
 *  @file NonlinearSolver.cpp
 *
 *  Damped Newton solver for 0D and 1D problems
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

#include <limits>

#include "SquareMatrix.h"
#include "NonlinearSolver.h"
#include "ctlapack.h"

#include "clockWC.h"
#include "vec_functions.h"
#include <ctime>

#include "mdp_allo.h"
#include <cfloat>

#include <vector>
#include <cstdio>
#include <cmath>

//@{
extern void print_line(const char *, int);

#ifndef MAX
#define MAX(x,y)    (( (x) > (y) ) ? (x) : (y))
#define MIN(x,y)    (( (x) < (y) ) ? (x) : (y))
#endif
#ifndef CONSTD_DATA_PTR
#define CONSTD_DATA_PTR(x) (( const double *) (&x[0]))
#endif
//@}
using namespace std;

namespace Cantera {


  //====================================================================================================================
  //-----------------------------------------------------------
  //                 Constants
  //-----------------------------------------------------------
  //! Dampfactor is the factor by which the damping factor is reduced by when a reduction in step length is warranted
  const doublereal DampFactor = 4.0;
  //! Number of damping steps that are carried out before the solution is deemed a failure
  const int NDAMP = 7;
  //====================================================================================================================
  //! Print a line of a single repeated character string
  /*!
   *  @param str  Character string
   *  @param n    Iteration length
   */
  static void print_line(const char *str, int n)  {
    for (int i = 0; i < n; i++) {
      printf("%s", str);
    }
    printf("\n");
  }

  bool NonlinearSolver::m_TurnOffTiming(false);

#ifdef DEBUG_NUMJAC 
  bool NonlinearSolver::s_print_NumJac(true);
#else
  bool NonlinearSolver::s_print_NumJac(false);
#endif

  // Turn off printing of dogleg information
  bool NonlinearSolver::s_print_DogLeg(false);

  // Turn off solving the system twice and comparing the answer.
  /*
   *  Turn this on if you want to compare the Hessian and Newton solve results.
   */
  bool  NonlinearSolver::s_doBothSolvesAndCompare(false);
 //====================================================================================================================
  // Default constructor
  /*
   * @param func   Residual and jacobian evaluator function object
   */
  NonlinearSolver::NonlinearSolver(ResidJacEval *func) :
    m_func(func),
    solnType_(NSOLN_TYPE_STEADY_STATE),
    neq_(0),
    m_ewt(0),
    m_manualDeltaStepSet(0),
    m_deltaStepMinimum(0),
    m_y_n(0),
    m_y_nm1(0),
    ydot_new(0),
    m_colScales(0),
    m_rowScales(0),
    m_rowWtScales(0),
    m_resid(0),
    m_wksp(0),
    m_residWts(0),
    m_normResid0(0.0),
    m_normResidFRaw(0.0),
    m_normDeltaSoln_Newton(0.0),
    m_normDeltaSoln_CP(0.0),
    m_normResidTrial(0.0),
    m_resid_scaled(false),
    m_y_high_bounds(0),
    m_y_low_bounds(0),
    m_dampBound(1.0),
    m_dampRes(1.0),
    delta_t_n(-1.0),
    m_nfe(0),
    m_colScaling(0),
    m_rowScaling(0),
    m_numTotalLinearSolves(0),
    m_numTotalNewtIts(0),
    m_min_newt_its(0),
    maxNewtIts_(100),
    m_jacFormMethod(NSOLN_JAC_NUM),
    m_nJacEval(0),
    time_n(0.0),
    m_matrixConditioning(0),
    m_order(1),
    rtol_(1.0E-3),
    atolBase_(1.0E-10),
    m_ydot_nm1(0),
    atolk_(0),
    m_print_flag(0),
    m_ScaleSolnNormToResNorm(0.001),
    jacCopy_(0),
    Hessian_(0),
    deltaX_CP_(0),
    deltaX_Newton_(0),
    residNorm2Cauchy_(0.0),
    RJd_norm_(0.0),
    lambda_(0.0),
    Jd_(0),
    deltaX_trust_(0),
    trustDelta_(1.0),
    Nuu_(0.0),
    dist_R0_(0.0),
    dist_R1_(0.0),
    dist_R2_(0.0),
    dist_Total_(0.0),
    JdJd_norm_(0.0),
    normTrust_Newton_(0.0),
    normTrust_CP_(0.0),
    doDogLeg_(0),
    doAffineSolve_(0)
  {
    neq_ = m_func->nEquations();

    m_ewt.resize(neq_, rtol_);
    m_deltaStepMinimum.resize(neq_, 0.001);
    m_deltaStepMaximum.resize(neq_, 1.0E10);
    m_y_n.resize(neq_, 0.0);
    m_y_nm1.resize(neq_, 0.0);
    ydot_new.resize(neq_, 0.0);
    m_colScales.resize(neq_, 1.0);
    m_rowScales.resize(neq_, 1.0);
    m_rowWtScales.resize(neq_, 1.0);
    m_resid.resize(neq_, 0.0);
    m_wksp.resize(neq_, 0.0);
    m_residWts.resize(neq_, 0.0);
    atolk_.resize(neq_, atolBase_);
    deltaX_Newton_.resize(neq_, 0.0);
    doublereal hb = std::numeric_limits<double>::max();
    m_y_high_bounds.resize(neq_, hb);
    m_y_low_bounds.resize(neq_, -hb);    

    for (int i = 0; i < neq_; i++) {
      atolk_[i] = atolBase_;
      m_ewt[i] = atolk_[i];
    }


    jacCopy_.resize(neq_, neq_, 0.0);
    deltaX_CP_.resize(neq_, 0.0);
    Jd_.resize(neq_, 0.0);
    deltaX_trust_.resize(neq_, 1.0);


  }
 //====================================================================================================================
  NonlinearSolver::NonlinearSolver(const NonlinearSolver &right) :
    m_func(right.m_func), 
    solnType_(NSOLN_TYPE_STEADY_STATE),
    neq_(0),
    m_ewt(0),
    m_manualDeltaStepSet(0),
    m_deltaStepMinimum(0),
    m_y_n(0),
    m_y_nm1(0),
    ydot_new(0),
    m_colScales(0),
    m_rowScales(0),
    m_rowWtScales(0),
    m_resid(0),
    m_wksp(0),
    m_residWts(0),
    m_normResid0(0.0),
    m_normResidFRaw(0.0),
    m_normDeltaSoln_Newton(0.0), 
    m_normDeltaSoln_CP(0.0),
    m_normResidTrial(0.0),
    m_resid_scaled(false),
    m_y_high_bounds(0),
    m_y_low_bounds(0),
    m_dampBound(1.0),
    m_dampRes(1.0),
    delta_t_n(-1.0),
    m_nfe(0),
    m_colScaling(0),
    m_rowScaling(0),
    m_numTotalLinearSolves(0),
    m_numTotalNewtIts(0),
    m_min_newt_its(0),
    maxNewtIts_(100),
    m_jacFormMethod(NSOLN_JAC_NUM),
    m_nJacEval(0),
    time_n(0.0),
    m_matrixConditioning(0),
    m_order(1),
    rtol_(1.0E-3),
    atolBase_(1.0E-10),
    m_ydot_nm1(0),
    atolk_(0),
    m_print_flag(0),
    m_ScaleSolnNormToResNorm(0.001),
    jacCopy_(0),
    Hessian_(0),
    deltaX_CP_(0),
    deltaX_Newton_(0),
    residNorm2Cauchy_(0.0),
    RJd_norm_(0.0),
    lambda_(0.0),
    Jd_(0),
    deltaX_trust_(0),
    trustDelta_(1.0),
    Nuu_(0.0),
    dist_R0_(0.0),
    dist_R1_(0.0),
    dist_R2_(0.0),
    dist_Total_(0.0),
    JdJd_norm_(0.0),
    normTrust_Newton_(0.0),
    normTrust_CP_(0.0),
    doDogLeg_(0),
    doAffineSolve_(0)
  {
    *this =operator=(right);
  }

 //====================================================================================================================
  NonlinearSolver::~NonlinearSolver() {
  }
 //====================================================================================================================
  NonlinearSolver& NonlinearSolver::operator=(const NonlinearSolver &right) {
    if (this == &right) {
      return *this;
    }
    // rely on the ResidJacEval duplMyselfAsresidJacEval() function to
    // create a deep copy
    m_func = right.m_func->duplMyselfAsResidJacEval();

    solnType_                  = right.solnType_;
    neq_                       = right.neq_;
    m_ewt                      = right.m_ewt;
    m_manualDeltaStepSet       = right.m_manualDeltaStepSet;
    m_deltaStepMinimum         = right.m_deltaStepMinimum;
    m_y_n                      = right.m_y_n;
    m_y_nm1                    = right.m_y_nm1;
    ydot_new                   = right.ydot_new;
    m_colScales                = right.m_colScales;
    m_rowScales                = right.m_rowScales;
    m_rowWtScales              = right.m_rowWtScales;
    m_resid                    = right.m_resid;
    m_wksp                     = right.m_wksp;
    m_residWts                 = right.m_residWts;
    m_normResid0               = right.m_normResid0;
    m_normResidFRaw            = right.m_normResidFRaw;
    m_normDeltaSoln_Newton     = right.m_normDeltaSoln_Newton;
    m_normDeltaSoln_CP         = right.m_normDeltaSoln_CP;
    m_normResidTrial           = right.m_normResidTrial;
    m_resid_scaled             = right.m_resid_scaled;
    m_y_high_bounds            = right.m_y_high_bounds;
    m_y_low_bounds             = right.m_y_low_bounds;
    m_dampBound                = right.m_dampBound;
    m_dampRes                  = right.m_dampRes;
    delta_t_n                  = right.delta_t_n;
    m_nfe                      = right.m_nfe;
    m_colScaling               = right.m_colScaling;
    m_rowScaling               = right.m_rowScaling;
    m_numTotalLinearSolves     = right.m_numTotalLinearSolves;
    m_numTotalNewtIts          = right.m_numTotalNewtIts;
    m_min_newt_its             = right.m_min_newt_its;
    maxNewtIts_                = right.maxNewtIts_;
    m_jacFormMethod            = right.m_jacFormMethod;
    m_nJacEval                 = right.m_nJacEval;
    time_n                     = right.time_n;
    m_matrixConditioning       = right.m_matrixConditioning;
    m_order                    = right.m_order;
    rtol_                      = right.rtol_;
    atolBase_                  = right.atolBase_;
    atolk_                     = right.atolk_;
    m_print_flag               = right.m_print_flag;
    m_ScaleSolnNormToResNorm   = right.m_ScaleSolnNormToResNorm;

    jacCopy_                   = right.jacCopy_;
    Hessian_                   = right.Hessian_;
    deltaX_CP_                 = right.deltaX_CP_;
    deltaX_Newton_             = right.deltaX_Newton_;
    RJd_norm_                  = right.RJd_norm_;
    lambda_                    = right.lambda_;
    Jd_                        = right.Jd_;
    deltaX_trust_              = right.deltaX_trust_;	
    trustDelta_                = right.trustDelta_;

    Nuu_                       = right.Nuu_;
    dist_R0_                   = right.dist_R0_;
    dist_R1_                   = right.dist_R1_;
    dist_R2_                   = right.dist_R2_;
    dist_Total_                = right.dist_Total_;
    JdJd_norm_                 = right.JdJd_norm_;
    normTrust_Newton_          = right.normTrust_Newton_;
    normTrust_CP_              = right.normTrust_CP_;
    doDogLeg_                  = right.doDogLeg_;
    doAffineSolve_             = right.doAffineSolve_;

    return *this;
  }
  //====================================================================================================================
  // Create solution weights for convergence criteria
  /*
   *  We create soln weights from the following formula
   *
   *  wt[i] = rtol * abs(y[i]) + atol[i]
   *
   *  The program always assumes that atol is specific
   *  to the solution component
   *
   * @param y  vector of the current solution values
   */
  void NonlinearSolver::createSolnWeights(const doublereal * const y) {
    for (int i = 0; i < neq_; i++) {
      m_ewt[i] = rtol_ * fabs(y[i]) + atolk_[i];
    }
  }
  //====================================================================================================================
  // set bounds constraints for all variables in the problem
  /*
   *  
   *   @param y_low_bounds  Vector of lower bounds
   *   @param y_high_bounds Vector of high bounds
   */
  void NonlinearSolver::setBoundsConstraints(const doublereal * const y_low_bounds,
					     const doublereal * const y_high_bounds) {
    for (int i = 0; i < neq_; i++) {
      m_y_low_bounds[i]  = y_low_bounds[i];
      m_y_high_bounds[i] = y_high_bounds[i];
    }
  }
  //====================================================================================================================
  void NonlinearSolver::setSolverScheme(int doDogLeg, int doAffineSolve) {
    doDogLeg_ = doDogLeg;
    doAffineSolve_ = doAffineSolve;
#ifdef DEBUG_DOGLEG

#else
    if (doDogLeg_) {
      throw CanteraError("NonlinearSolver::setSolverScheme",
			 "ifdef block not on");
    }
#endif
  }
  //====================================================================================================================
  std::vector<double> &  NonlinearSolver::lowBoundsConstraintVector() {
    return  m_y_low_bounds;
  }
  //====================================================================================================================
  std::vector<double> &  NonlinearSolver::highBoundsConstraintVector() {
    return  m_y_high_bounds;
  }
  //====================================================================================================================
  //  L2 norm of the delta of the solution vector
  /*
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
  doublereal NonlinearSolver::solnErrorNorm(const doublereal * const delta_y, const char * title, int printLargest,
					    const doublereal dampFactor) const
  {
    int  i;
    doublereal sum_norm = 0.0, error;
    for (i = 0; i < neq_; i++) {
      error     = delta_y[i] / m_ewt[i];
      sum_norm += (error * error);
    }
    sum_norm = sqrt(sum_norm / neq_); 
    if (printLargest) { 
      if (m_print_flag >= 4 && m_print_flag <= 5) {

	printf("\t\t   solnErrorNorm(): ");
	if (title) {
	  printf("%s", title);
	} else {
	  printf(" Delta soln norm ");
	}
	printf(" = %-11.4E\n", sum_norm);
      } else if (m_print_flag >= 6) {



	const int num_entries = printLargest;
	printf("\t\t   ");
        print_line("-", 90);
        printf("\t\t   solnErrorNorm(): ");
        if (title) {
          printf("%s", title);
        } else {
          printf(" Delta soln norm ");
        }
        printf(" = %-11.4E\n", sum_norm);

	doublereal dmax1, normContrib;
	int j; 
	int *imax = mdp::mdp_alloc_int_1(num_entries, -1);
        printf("\t\t   Printout of Largest Contributors:\n");
        printf("\t\t                                                      (damp = %g)\n", dampFactor);
        printf("\t\t      I   weightdeltaY/sqtN|     deltaY    "
	       "ysolnOld     ysolnNew   Soln_Weights\n");
        printf("\t\t   ");
        print_line("-", 90);

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
	    error = delta_y[i] / m_ewt[i];
	    normContrib = sqrt(error * error);
	    printf("\t\t   %4d %12.4e       | %12.4e %12.4e %12.4e %12.4e\n", i, normContrib/sqrt((double)neq_), 
		   delta_y[i], m_y_n[i], m_y_n[i] + dampFactor * delta_y[i], m_ewt[i]);

	  }	  
	}
	printf("\t\t   "); 
	print_line("-", 90);
	mdp::mdp_safe_free((void **) &imax);
      }
    }
    return sum_norm;
  }
  //====================================================================================================================
  /*
   * L2 Norm of the residual
   *
   *  The second argument has a default of false. However,
   *  if true, then a table of the largest values is printed
   *  out to standard output.
   */
  doublereal NonlinearSolver::residErrorNorm(const doublereal * const resid, const char * title, const int printLargest,
					     const doublereal * const y)
  {
    int    i;
    doublereal sum_norm = 0.0, error;

    for (i = 0; i < neq_; i++) {
#ifdef DEBUG_HKM
      mdp::checkFinite(resid[i]);
#endif
      error     = resid[i] / m_residWts[i];
#ifdef DEBUG_HKM
      mdp::checkFinite(error);
#endif
      sum_norm += (error * error);
    }
    sum_norm = sqrt(sum_norm / neq_); 
#ifdef DEBUG_HKM
    mdp::checkFinite(sum_norm);
#endif
    if (printLargest) {
      const int num_entries = printLargest;
      doublereal dmax1, normContrib;
      int j;
      int *imax = mdp::mdp_alloc_int_1(num_entries, -1);


      printf("\t  ");
      print_line("-", 90);
      printf("\t\t  residErrorNorm():");
      if (title) {
        printf(" %s ", title);
      } else {
        printf("  residual L2 norm ");
      }
      printf("= %12.4E\n", sum_norm);
      
      if (m_print_flag >= 6) {
        printf("\t\t   Printout of Largest Contributors to norm:\n");
        printf("\t\t      I       |Resid/ResWt|     UnsclRes          ResWt   |  y_curr\n");
        printf("\t\t   ");
        print_line("-", 80);
	for (int jnum = 0; jnum < num_entries; jnum++) {
	  dmax1 = -1.0;
	  for (i = 0; i < neq_; i++) {
	    bool used = false;
	    for (j = 0; j < jnum; j++) {
	      if (imax[j] == i) used = true;
	    }
	    if (!used) {
	      error = resid[i] / m_residWts[i];
	      normContrib = sqrt(error * error);
	      if (normContrib > dmax1) {
		imax[jnum] = i;
		dmax1 = normContrib;
	      }
	    }
	  }
	  i = imax[jnum];
	  if (i >= 0) {
	    error = resid[i] / m_residWts[i];
	    normContrib = sqrt(error * error);
	    printf("\t\t   %4d     %12.4e     %12.4e     %12.4e | %12.4e\n", i,  normContrib, resid[i], m_residWts[i], y[i]);
	  }	  
	}
      
	printf("\t\t   "); 
	print_line("-", 80);
      }
      mdp::mdp_safe_free((void **) &imax);
    }
    return sum_norm;
  }
  //====================================================================================================================
  /*
   * setColumnScales():
   *
   * Set the column scaling vector at the current time
   */
  void NonlinearSolver::setColumnScales() {
    if (m_colScaling) {
      for (int i = 0; i < neq_; i++) {
	m_colScales[i] = m_ewt[i];
      }
    } else {
      for (int i = 0; i < neq_; i++) {
	m_colScales[i] = 1.0;
      }
    }
    m_func->calcSolnScales(time_n, DATA_PTR(m_y_n), DATA_PTR(m_y_nm1), DATA_PTR(m_colScales));
  }
  //====================================================================================================================
  // Compute the current residual
  /*
   *  @param time_curr    Value of the time 
   *  @param typeCalc     Type of the calculation
   *  @param y_curr       Current value of the solution vector
   *  @param ydot_curr    Current value of the time derivative of the solution vector
   *
   * @return Returns a flag to indicate that operation is successful.
   *            1  Means a successful operation
   *           -0 or neg value Means an unsuccessful operation
   */
  int NonlinearSolver::doResidualCalc(const doublereal time_curr, const int typeCalc, const doublereal * const y_curr, 
				       const doublereal * const ydot_curr, const ResidEval_Type_Enum evalType)
  {
    int retn = m_func->evalResidNJ(time_curr, delta_t_n, y_curr, ydot_curr, DATA_PTR(m_resid), evalType);
    m_nfe++;
    m_resid_scaled = false;
    return retn;
  }
  //====================================================================================================================
  // Scale the matrix
  /*
   *  @param jac              Jacobian
   *  @param y_comm           Current value of the solution vector
   *  @param ydot_comm        Current value of the time derivative of the solution vector
   *  @param time_curr        current value of the time
   */
  void NonlinearSolver::scaleMatrix(SquareMatrix& jac, double* y_comm,	double* ydot_comm, doublereal time_curr)
  {	 
    int irow, jcol;



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
	doublereal *jptr = &(*(jac.begin()));
	for (jcol = 0; jcol < neq_; jcol++) {
	  for (irow = 0; irow < neq_; irow++) {
	    *jptr *= m_colScales[jcol];
	    jptr++;
	  }
	}
      }	  
    }

  
    /*
     * row sum scaling -> Note, this is an unequivical success
     *      at keeping the small numbers well balanced and nonnegative.
     */
   
    if (! jac.m_factored) {
      /*
       * Ok, this is ugly. jac.begin() returns an vector<double> iterator
       * to the first data location.
       * Then &(*()) reverts it to a doublereal *.
       */
      doublereal *jptr = &(*(jac.begin()));
      for (irow = 0; irow < neq_; irow++) {
	m_rowScales[irow] = 0.0;
	m_rowWtScales[irow] = 0.0;
      }
      for (jcol = 0; jcol < neq_; jcol++) {
	for (irow = 0; irow < neq_; irow++) {
	  //m_rowScales[irow] += fabs(*jptr);
	  if (m_colScaling) {
	    // This is needed in order to mitgate the change in J_ij carried out just above this loop.
	    // Alternatively, we could move this loop up to the top
	    m_rowWtScales[irow] += fabs(*jptr) * m_ewt[jcol] / m_colScales[jcol];
	  } else {
	    m_rowWtScales[irow] += fabs(*jptr) * m_ewt[jcol];
	  }
	  jptr++;
	}
      }

      for (irow = 0; irow < neq_; irow++) {
	m_rowScales[irow] = 1.0/m_rowScales[irow];
      }
      // What we have defined is a maximum value that the residual can be and still pass.
      // This isn't sufficient.

      if (m_rowScaling) {
	
	jptr = &(*(jac.begin()));
	for (jcol = 0; jcol < neq_; jcol++) {
	  for (irow = 0; irow < neq_; irow++) {
	    *jptr *= m_rowScales[irow];
	    jptr++;
	  }
	}
      }
   
      
      computeResidWts();
   
    }

  }
  //====================================================================================================================
  // Calculate the scaling factor for translating residual norms into solution norms.
  /*
   *  This routine calls computeResidWts() a couple of times in the calculation of m_ScaleSolnNormToResNorm.
   *  A more sophisticated routine may do more with signs to get a better value. Perhaps, a series of calculations
   *  with different signs attached may be in order. Then, m_ScaleSolnNormToResNorm would be calculated
   *  as the minimum of a series of calculations.
   */
  void NonlinearSolver::calcSolnToResNormVector()
  { 
    if (! jacCopy_.m_factored) { 
      m_ScaleSolnNormToResNorm = 1.0;
      computeResidWts();
      for (int n = 0; n < neq_; n++) {
        m_wksp[n] = 0.0;
      }   
      doublereal *jptr = &(*(jacCopy_.begin()));
      for (int jcol = 0; jcol < neq_; jcol++) {
	for (int irow = 0; irow < neq_; irow++) {
	  m_wksp[irow] += (*jptr) * m_ewt[jcol];
	  jptr++;
	}
      }
      double resNormOld = residErrorNorm(DATA_PTR(m_wksp));
      if (resNormOld > 0.0) {
	m_ScaleSolnNormToResNorm = m_ScaleSolnNormToResNorm * resNormOld;
      }
      if (m_ScaleSolnNormToResNorm < 1.0E-8) {
	m_ScaleSolnNormToResNorm = 1.0E-8;
      }
      // Recalculate the residual weights now that we know the value of m_ScaleSolnNormToResNorm 
      computeResidWts();
    } else {
      throw CanteraError("NonlinearSolver::calcSolnToResNormVector()" , "Logic error");
    }
  }
  //====================================================================================================================
  // Compute the undamped Newton step based on the current jacobian and an input rhs
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
  int NonlinearSolver::doNewtonSolve(const doublereal time_curr, const doublereal * const y_curr, 
				     const doublereal * const ydot_curr,  doublereal * const delta_y, 
				     SquareMatrix& jac, int loglevel)
  {
    int irow;


    // multiply the residual by -1
    if (m_rowScaling  && !m_resid_scaled) {
      for (int n = 0; n < neq_; n++) {
	delta_y[n] = -m_rowScales[n] * m_resid[n];
      }
      m_resid_scaled = true;
    } else {
      for (int n = 0; n < neq_; n++) {
        delta_y[n] = -m_resid[n];
      }
    }



    /*
     * Solve the system -> This also involves inverting the
     * matrix
     */
    int info = jac.solve(DATA_PTR(delta_y));


    /*
     * reverse the column scaling if there was any.
     */
    if (m_colScaling) {
      for (irow = 0; irow < neq_; irow++) {
	delta_y[irow] = delta_y[irow] * m_colScales[irow];
      }
    }
	
#ifdef DEBUG_JAC
    if (printJacContributions) {
      for (int iNum = 0; iNum < numRows; iNum++) {
	if (iNum > 0) focusRow++;
	doublereal dsum = 0.0;
	vector_fp& Jdata = jacBack.data();
	doublereal dRow = Jdata[neq_ * focusRow + focusRow];
	printf("\n Details on delta_Y for row %d \n", focusRow);
	printf("  Value before = %15.5e, delta = %15.5e,"
	       "value after = %15.5e\n", y_curr[focusRow], 
	       delta_y[focusRow],  y_curr[focusRow] + delta_y[focusRow]);
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
	    doublereal aij =  Jdata[neq_ * ii + focusRow];
	    doublereal contrib = aij * delta_y[ii] * (-1.0) / dRow;
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
    return info;
  }
  //====================================================================================================================
  // Compute the newton step, either by direct newton's or by solving a close problem that is represented
  // by a Hessian (
  /*
   * This is algorith A.6.5.1 in Dennis / Schnabel
   *
   * Compute the QR decomposition
   */
  int NonlinearSolver::doAffineNewtonSolve(const doublereal * const y_curr,   const doublereal * const ydot_curr, 
					   doublereal * const delta_y, SquareMatrix& jac)
  {
    bool newtonGood = true;
    int irow;
    doublereal *delyNewton = 0;
    // We can default to QR here ( or not )
    jac.useQR_ = true;
    // multiply the residual by -1
    // Scale the residual if there is row scaling. Note, the matrix has already been scaled
    if (m_rowScaling  && !m_resid_scaled) {
      for (int n = 0; n < neq_; n++) {
	delta_y[n] = -m_rowScales[n] * m_resid[n];
      }
      m_resid_scaled = true;
    } else {
      for (int n = 0; n < neq_; n++) {
        delta_y[n] = -m_resid[n];
      }
    }

    // Factor the matrix
    double cond = 1.0E300;
    int info = 0;
    if (!jac.m_factored) { 
      if (jac.useQR_) {
	info = jac.factorQR();
      } else {
	info = jac.factor();
      }
    }
    /*
     *  Find the condition number of the matrix
     *  If we have failed to factor, we will fall back to calculating and factoring a modified Hessian
     */
    if (info == 0) {  
      double rcond = 0.0;
      if (jac.useQR_) {
	rcond = jac.rcondQR();
      } else {
	rcond = jac.rcond(jac.a1norm_);
      }
      if (rcond > 0.0) {
	cond = 1.0 / rcond;
      }
    }
    bool doHessian = false;
    if (s_doBothSolvesAndCompare) {
      doHessian = true;
    }
    bool doNewton = false;
    if (cond < 1.0E7) {
      doNewton = true;
      if (m_print_flag >= 3) {
	printf("\t\t\tdoAffineNewtonSolve: Condition number = %g during regular solve\n", cond);
      }

      /*
       * Solve the system -> This also involves inverting the matrix
       */
      int info = jac.solve(DATA_PTR(delta_y));
      if (info) {
	if (m_print_flag >= 2) {
	  printf("\t\t\tNonlinearSolver::doAffineSolve QRSolve returned INFO = %d. Switching to Hessian solve\n", info);
	}
	doHessian = true;
	newtonGood = false;
      }
      /*
       * reverse the column scaling if there was any on a successful solve
       */
      if (m_colScaling) {
	for (irow = 0; irow < neq_; irow++) {
	  delta_y[irow] = delta_y[irow] * m_colScales[irow];
	}
      }
      
    } else {
      doHessian = true;
      newtonGood = false;
      if (m_print_flag >= 3) {
	printf("\t\t\tdoAffineNewtonSolve: Condition number too large, %g. Doing a Hessian solve \n", cond);
      }
    }

    if (doHessian) {
      // Store the old value for later comparison
      if (doNewton) {
	delyNewton = mdp::mdp_alloc_dbl_1(neq_, MDP_DBL_NOINIT);
	for (irow = 0; irow < neq_; irow++) {
	  delyNewton[irow] = delta_y[irow];
	}
      }
      // Get memory if not done before
      if (Hessian_.nRows() == 0) {
	Hessian_.resize(neq_, neq_);
      }

      /*
       * Calculate the symmetric Hessian
       */
      Hessian_.zero();
      if (m_rowScaling) {
	for (int i = 0; i < neq_; i++) {
	  for (int j = i; j < neq_; j++) {
	    for (int k = 0; k < neq_; k++) {
	      Hessian_(i,j) += jacCopy_(k,i) * jacCopy_(k,j) * m_rowScales[k] * m_rowScales[k];
	    }
	    Hessian_(j,i) = Hessian_(i,j);
	  }
	}
      } else {
	for (int i = 0; i < neq_; i++) {
	  for (int j = i; j < neq_; j++) {
	    for (int k = 0; k < neq_; k++) {
	      Hessian_(i,j) += jacCopy_(k,i) * jacCopy_(k,j);
	    }
	    Hessian_(j,i) = Hessian_(i,j);
	  }
	}
      }

      /*
       * Calculate the matrix norm of the Hessian
       */
      doublereal hnorm = 0.0;
      doublereal hcol = 0.0;
      if (m_colScaling) {
	for (int i = 0; i < neq_; i++) {
	  for (int j = i; j < neq_; j++) {
	    hcol += fabs(Hessian_(j,i)) * m_colScales[j];
	  }
	  for (int j = i+1; j < neq_; j++) {
	    hcol += fabs(Hessian_(i,j)) * m_colScales[j];
	  }
	  hcol *= m_colScales[i];
	  if (hcol > hnorm) {
	    hnorm = hcol;
	  }
	}
      } else {
	for (int i = 0; i < neq_; i++) {
	  for (int j = i; j < neq_; j++) {
	    hcol += fabs(Hessian_(j,i));
	  }
	  for (int j = i+1; j < neq_; j++) {
	    hcol += fabs(Hessian_(i,j));
	  }
	  if (hcol > hnorm) {
	    hnorm = hcol;
	  }
	}
      }
      /*
       * Add junk to the Hessian diagonal
       */
      hcol = sqrt(neq_) * 1.0E-7 * hnorm;
      if (m_colScaling) {
	for (int i = 0; i < neq_; i++) {
	  Hessian_(i,i) += hcol / (m_colScales[i] * m_colScales[i]);
	}
      } else {
	for (int i = 0; i < neq_; i++) {
	  Hessian_(i,i) += hcol;
	}
      }

      /*
       *  Factor the Hessian
       */
      int info;
      ct_dpotrf(ctlapack::UpperTriangular, neq_, &(*(Hessian_.begin())), neq_, info);
      if (info) {
	if (m_print_flag >= 2) {
	  printf("\t\t\tNonlinearSolver::doAffineSolve DPOTRF returned INFO = %d\n", info);
	}
	return info;
      }

      // doublereal *JTF = delta_y;
      doublereal *delyH = mdp::mdp_alloc_dbl_1(neq_, MDP_DBL_NOINIT);
      // First recalculate the scaled residual. It got wiped out doing the newton solve
      if (m_rowScaling) {
	for (int n = 0; n < neq_; n++) {
	  delyH[n] = -m_rowScales[n] * m_resid[n];
	}
      } else {
	for (int n = 0; n < neq_; n++) {
	  delyH[n] = -m_resid[n];
	}
      }

      if (m_rowScaling) {
	for (int j = 0; j < neq_; j++) {
	  delta_y[j] = 0.0;
	  for (int i = 0; i < neq_; i++) {
	    delta_y[j] += delyH[i] * jacCopy_.value(i,j) * m_rowScales[i];
	  }
	}
      } else {
	for (int j = 0; j < neq_; j++) {
	  delta_y[j] = 0.0;
	  for (int i = 0; i < neq_; i++) {
	    delta_y[j] += delyH[i] * jacCopy_.value(i,j);
	  }
	}
      }

    
      /*
       * Solve the factored Hessian System
       */
      ct_dpotrs(ctlapack::UpperTriangular, neq_, 1,&(*(Hessian_.begin())), neq_, delta_y, neq_, info);
      if (info) {
	if (m_print_flag >= 2) {
	  printf("\t\t\tNonlinearSolver::doAffineSolve DPOTRS returned INFO = %d\n", info);
	}
	return info;
      }
      /*
       * reverse the column scaling if there was any.
       */
      if (m_colScaling) {
	for (irow = 0; irow < neq_; irow++) {
	  delta_y[irow] = delta_y[irow] * m_colScales[irow];
	}
      }


      if (s_print_DogLeg || (doDogLeg_ && m_print_flag > 3)) {
	printf("\t\t\t Comparison between Hessian deltaX and newton deltaX\n");
	printf("\t\t\t   i    Hessian+Junk  Newton \n");
	printf("\t\t\t--------------------------------------------------------\n");
	for (int i =0; i < neq_; i++) {
	  printf("\t\t\t%3d  %12.5g %12.5g\n", i, delta_y[i], delyNewton[i]);
	}
	printf("\t\t\t--------------------------------------------------------\n");

      }

      if (newtonGood) {
	mdp::mdp_copy_dbl_1(DATA_PTR(delta_y), CONSTD_DATA_PTR(delyNewton), neq_);
      }
      mdp::mdp_safe_free((void **) &delyH);
      mdp::mdp_safe_free((void **) &delyNewton);
    }
	
#ifdef DEBUG_JAC
    if (printJacContributions) {
      for (int iNum = 0; iNum < numRows; iNum++) {
	if (iNum > 0) focusRow++;
	doublereal dsum = 0.0;
	vector_fp& Jdata = jacBack.data();
	doublereal dRow = Jdata[neq_ * focusRow + focusRow];
	printf("\n Details on delta_Y for row %d \n", focusRow);
	printf("  Value before = %15.5e, delta = %15.5e,"
	       "value after = %15.5e\n", y_curr[focusRow], 
	       delta_y[focusRow],  y_curr[focusRow] + delta_y[focusRow]);
	if (!freshJac) {
	  printf("    Old Jacobian\n");
	}
	printf("     col          delta_y            aij     "
	       "contrib   \n");
	printf("-----------------------------------------------------------------------------------------------\n");
	printf(" Res(%d) %15.5e  %15.5e  %15.5e  (Res = %g)\n",
	       focusRow, delta_y[focusRow],
	       dRow, RRow[iNum] / dRow, RRow[iNum]);
	dsum +=  RRow[iNum] / dRow;
	for (int ii = 0; ii < neq_; ii++) {
	  if (ii != focusRow) {
	    doublereal aij =  Jdata[neq_ * ii + focusRow];
	    doublereal contrib = aij * delta_y[ii] * (-1.0) / dRow;
	    dsum += contrib;
	    if (fabs(contrib) > Pcutoff) {
	      printf("%6d  %15.5e  %15.5e  %15.5e\n", ii,
		     delta_y[ii]  , aij, contrib); 
	    }
	  }
	}
	printf("-----------------------------------------------------------------------------------------------\n");
	printf("        %15.5e                   %15.5e\n",
	       delta_y[focusRow], dsum);
      }
    }

#endif
	
    m_numTotalLinearSolves++;
    return info;

  }
  //====================================================================================================================
  // Do a steepest descent calculation
  /*
   *  This call must be made on the unfactored jacobian!
   */
  doublereal NonlinearSolver::doCauchyPointSolve(SquareMatrix& jac)
  {
    double rowFac = 1.0;
    double normSoln;
    //  Calculate desDir = -0.5 * R dot J
    /*
     * For confirmation of the scaling factors, see Dennis and Schnabel p, 152, p, 156 and my notes
     * 
     *  The colFac and rowFac values are used to eliminate the scaling of the matrix from the 
     *  actual equation
     *
     *  Here we calculate the steepest direction. this is equation (10) in the notes. It is
     *  storred in deltaX_CP_[].The value corresponds to d_descent[].
     */
    for (int j = 0; j < neq_; j++) {
      deltaX_CP_[j] = 0.0;
      double colFac = 1.0;
      if (m_colScaling) {
	colFac = 1.0 / m_colScales[j];
      }
      for (int i = 0; i < neq_; i++) {
	if (m_rowScaling) {
	  rowFac = 1.0 / m_rowScales[i];
	}
        deltaX_CP_[j] -= m_resid[i] * jac.value(i,j) * colFac * rowFac * m_ewt[j] * m_ewt[j] 
                	  / (m_residWts[i] * m_residWts[i]);
      }
    }

    /*
     *  Calculate J_hat d_y_descent. This is formula 17 in the notes.
     */
    for (int i = 0; i < neq_; i++) {
      Jd_[i] = 0.0;
      if (m_rowScaling) {
	rowFac = 1.0 / m_rowScales[i];
      } else {
	rowFac = 1.0;
      }
      for (int j = 0; j < neq_; j++) {
	Jd_[i] += deltaX_CP_[j] * jac.value(i,j) * rowFac/ m_residWts[i];
      }
    }

    /*
     * Calculate the distance along the steepest descent until the Cauchy point
     * This is Eqn. 16 in the notes.
     */
    RJd_norm_ = 0.0;
    JdJd_norm_ = 0.0;
    for (int i = 0; i < neq_; i++) {
      RJd_norm_ += m_resid[i] * Jd_[i] / m_residWts[i];
      JdJd_norm_ += Jd_[i] * Jd_[i];
    }
    lambda_ = - RJd_norm_ / (JdJd_norm_);

    /*
     *  Now we modify the steepest descent vector such that its length is equal to the 
     *  Cauchy distance. From now on, if we want to recreate the descent vector, we have
     *  to unnormalize it by dividing by lambda_.
     */
    for (int i = 0; i < neq_; i++) {
      deltaX_CP_[i] *= lambda_;
    }
    double normResid02 = m_normResid0 * m_normResid0 * neq_;

    residNorm2Cauchy_ = normResid02 - RJd_norm_ * RJd_norm_ / (JdJd_norm_);

    if (m_print_flag > 2) {

      double residCauchy = 0.0;
      if (residNorm2Cauchy_ > 0.0) {
	residCauchy = sqrt(residNorm2Cauchy_ / neq_);
      } else {
	residCauchy =  m_normResid0 - sqrt(RJd_norm_ * RJd_norm_ / (JdJd_norm_));
      }
      
      // Compute the weighted norm of the undamped step size descentDir_[]
      if (s_print_DogLeg || (doDogLeg_ && m_print_flag > 6)) {
	normSoln = solnErrorNorm(DATA_PTR(deltaX_CP_), "SteepestDescentDir", 10);
      } else {
	normSoln = solnErrorNorm(DATA_PTR(deltaX_CP_), "SteepestDescentDir", 0);
      }
      if (s_print_DogLeg || (doDogLeg_ && m_print_flag > 3)) {
	printf("\t\t\tdoCauchyPointSolve: Steepest descent to Cauchy point: \n");
	printf("\t\t\t      R0     = %g \n", m_normResid0);
	printf("\t\t\t      Rpred  = %g\n",  residCauchy);
	printf("\t\t\t      Rjd    = %g\n",  RJd_norm_);
	printf("\t\t\t      JdJd   = %g\n",  JdJd_norm_);
	printf("\t\t\t      deltaX = %g\n",  normSoln);
	printf("\t\t\t      lambda = %g\n",  lambda_);
      }
    }
    return normSoln;
  }
  //===================================================================================================================
  void  NonlinearSolver::descentComparison(double time_curr, double *ydot0, double *ydot1, const double *newtDir)
  {
    int info;
    double ff = 1.0E-5;
    double *y1 = DATA_PTR(m_wksp);
    double cauchyDistanceNorm = solnErrorNorm(DATA_PTR(deltaX_CP_));
    for (int i = 0; i < neq_; i++) {
      y1[i] = m_y_n[i] + ff * deltaX_CP_[i];
    }
    /*
     *  Calculate the residual that would result if y1[] were the new solution vector
     *  -> m_resid[] contains the result of the residual calculation
     */
    if (solnType_ != NSOLN_TYPE_STEADY_STATE) {
      info = doResidualCalc(time_curr, solnType_, y1, ydot1, Base_LaggedSolutionComponents);
    } else {
      info = doResidualCalc(time_curr, solnType_, y1, ydot0, Base_LaggedSolutionComponents);
    }  

    double normResid02 = m_normResid0 * m_normResid0 * neq_;
    double residSteep = residErrorNorm(DATA_PTR(m_resid));
    double residSteep2 = residSteep * residSteep * neq_;
    double funcDecrease2 = 0.5 * (residSteep2 - normResid02) / ( ff * cauchyDistanceNorm);

    double sNewt = solnErrorNorm(DATA_PTR(newtDir));
    for (int i = 0; i < neq_; i++) {
      y1[i] = m_y_n[i] + ff * newtDir[i];
    }
    /*
     *  Calculate the residual that would result if y1[] were the new solution vector
     *  -> m_resid[] contains the result of the residual calculation
     */
    if (solnType_ != NSOLN_TYPE_STEADY_STATE) {
      info = doResidualCalc(time_curr, solnType_, y1, ydot1, Base_LaggedSolutionComponents);
    } else {
      info = doResidualCalc(time_curr, solnType_, y1, ydot0, Base_LaggedSolutionComponents);
    }  
    double residNewt = residErrorNorm(DATA_PTR(m_resid));
    double residNewt2 = residNewt * residNewt * neq_;

    double funcDecreaseNewt2 = 0.5 * (residNewt2 - normResid02) / ( ff * sNewt);

    // This is the expected inital rate of decrease in the cauchy direction.
    //   -> This is Eqn. 29 = Rhat dot Jhat dy / || d ||
    double funcDecreaseSDExp = RJd_norm_ / cauchyDistanceNorm * lambda_;

    double funcDecreaseNewtExp2 = -  normResid02 / sNewt;

    /*
     * HKM These have been shown to exactly match up.
     *   The steepest direction is always largest even when there are variable solution weights
     */
    if (s_print_DogLeg || (doDogLeg_ && m_print_flag > 3)) {
      printf("descentComparison: initial rate of decrease in cauchy dir (expected) = %g\n", funcDecreaseSDExp);
      printf("descentComparison: initial rate of decrease in cauchy dir            = %g\n", funcDecrease2);
      printf("descentComparison: initial rate of decrease in newton dir (expected) = %g\n", funcDecreaseNewtExp2);
      printf("descentComparison: initial rate of decrease in newton dir            = %g\n", funcDecreaseNewt2);
    }
  }

  //====================================================================================================================
  //  Setup the line search along the double dog leg
  /*
   *  the calls the doCauchySolve() and doNewtonSolve() are done at the main level
   */
  void NonlinearSolver::setupDoubleDogleg(double * newtDir)
  {
    /*
     *  Gamma =                  ||grad f ||**4
     *         ---------------------------------------------
     *           (grad f)T H (grad f)  (grad f)T H-1 (grad f)
     */
    // doublereal sumG = 0.0;
    // doublereal sumH = 0.0;
    //    for (int i = 0; i < neq_; i++) {
    //  sumG =  deltax_cp_[i] *  deltax_cp_[i];
    //   sumH =  deltax_cp_[i] *  newtDir[i];
    //  }
    // double  fac1 = sumG / lambda_;
    // double  fac2 = sumH / lambda_;
    // double  gamma = fac1 / fac2;
    // double   gamma =  m_normDeltaSoln_CP / m_normDeltaSoln_Newton;
    /*
     * This hasn't worked. so will do it heuristically. One issue is that the newton
     * direction is not the inverse of the Hessian times the gradient. The Hession
     * is the matrix squared. Until I have the inverse of the Hessian from QR factorization
     * I may not be able to do it this way.
     */

    /*
     *  Heuristic algorithm - Find out where on the Newton line the residual is the same
     *                        as the residual at the cauchy point. Then, go halfway to
     *                        the newton point and call that Nuu.
     *                        Maybe we need to check that the linearized residual is
     *                        monotonic along that line. However, we haven't needed to yet.
     */
    double residSteepLin = expectedResidLeg(0, 1.0);
    double Nres2CP = residSteepLin * residSteepLin * neq_;
    double Nres2_o = m_normResid0 * m_normResid0 * neq_;
    double a = Nres2CP /  Nres2_o;
    double betaEqual = (2.0 - sqrt(4.0 - 4 * (1.0 - a))) / 2.0;
    double  beta = (1.0 + betaEqual) / 2.0;


    Nuu_ = beta;

    // put in a loop here to test derivative.


    dist_R0_ = m_normDeltaSoln_CP;
    dist_R1_ = solnErrorNorm( DATA_PTR(m_wksp));
    dist_R2_ = (1.0 - Nuu_) * m_normDeltaSoln_Newton;
    dist_Total_ = dist_R0_ + dist_R1_ + dist_R2_;

    /*
     * Calculate the trust distances
     */
    normTrust_Newton_ = calcTrustDistance(deltaX_Newton_);

    normTrust_CP_ = calcTrustDistance(deltaX_CP_);

  } 
  //====================================================================================================================
  int NonlinearSolver::lambdaToLeg(const double lambda, double &alpha) const {

    if (lambda < dist_R0_ / dist_Total_) {
      alpha = lambda * dist_Total_ / dist_R0_;
      return 0;
    } else if (lambda < ((dist_R0_ + dist_R1_)/ dist_Total_)) {
      alpha = (lambda * dist_Total_ - dist_R0_) / dist_R1_;
      return 1;
    } 
    alpha = (lambda * dist_Total_ - dist_R0_ - dist_R1_) / dist_R2_;
    return 2;
  }
  //====================================================================================================================
  
  double NonlinearSolver::expectedResidLeg(int leg, double alpha) const {

    double resD2, res2, resNorm;
    double normResid02 = m_normResid0 * m_normResid0 * neq_;

    if (leg == 0) {
      /*
       * We are on the steepest descent line
       *  along that line 
       *   R2 = R2 + 2 lambda R dot Jd + lambda**2 Jd dot Jd
       */

      double tmp = - 2.0 * alpha + alpha * alpha;
      double tmp2 = - RJd_norm_ * lambda_;
      resD2 = tmp2 * tmp;

    } else if (leg == 1) {

      /*
       *  Same formula as above for lambda=1.
       */
      double tmp2 = - RJd_norm_ * lambda_;
      double RdotJS = - tmp2;
      double JsJs =     tmp2;


      double res0_2 = m_normResid0 * m_normResid0 * neq_;

      res2 =  res0_2 + (1.0 - alpha) * 2 * RdotJS - 2 * alpha * Nuu_ * res0_2
	+ (1.0 - alpha) * (1.0 - alpha) * JsJs
	+ alpha * alpha * Nuu_ * Nuu_ * res0_2 
	- 2 * alpha * Nuu_ * (1.0 - alpha) * RdotJS;

      resNorm = sqrt(res2 / neq_);
      return resNorm;

    } else {
      double beta = Nuu_ + alpha * (1.0 - Nuu_);
      double tmp2 =  normResid02;
      double tmp = 1.0 - 2.0 * beta + 1.0 * beta * beta - 1.0;
      resD2 = tmp * tmp2;
    }

    res2 = m_normResid0 * m_normResid0 * neq_ + resD2;
    if (res2 < 0.0) {
      resNorm =  m_normResid0 - sqrt(resD2/neq_);
    } else {
      resNorm = sqrt(res2 / neq_);
    }
      
    return resNorm;

  } 
  //====================================================================================================================
  //  Here we print out the residual at various points along the double dogleg, comparing against the quadratic model

  void NonlinearSolver::residualComparisonLeg(const double time_curr, const double *ydot0, 
					      const double *ydot1, const double *newtDir)  {
    double *y1 = DATA_PTR(m_wksp);
    double sLen;
    if (s_print_DogLeg || (doDogLeg_ && m_print_flag > 6)) {
      printf("     residualComparisonLeg() \n");
      printf("      Point               StepLen     Residual_Actual  Residual_Linear  RelativeMatch\n");
    }
    // First compare at 1/4 along SD curve
    std::vector<double> alphaT;
    alphaT.push_back(0.00);
    alphaT.push_back(0.01);
    alphaT.push_back(0.1);
    alphaT.push_back(0.25);
    alphaT.push_back(0.50);
    alphaT.push_back(0.75);
    alphaT.push_back(1.0);
    for (int iteration = 0; iteration < (int) alphaT.size(); iteration++) {
      double alpha = alphaT[iteration];
      for (int i = 0; i < neq_; i++) {
	y1[i] = m_y_n[i] + alpha * deltaX_CP_[i];
      }
      sLen = alpha * solnErrorNorm(DATA_PTR(deltaX_CP_));
      /*
       *  Calculate the residual that would result if y1[] were the new solution vector
       *  -> m_resid[] contains the result of the residual calculation
       */
      if (solnType_ != NSOLN_TYPE_STEADY_STATE) {
	doResidualCalc(time_curr, solnType_, y1, ydot1, Base_LaggedSolutionComponents);
      } else {
	doResidualCalc(time_curr, solnType_, y1, ydot0, Base_LaggedSolutionComponents);
      }  


      double residSteep = residErrorNorm(DATA_PTR(m_resid));
      double residSteepLin = expectedResidLeg(0, alpha);

      double relFit = (residSteep - residSteepLin) / (fabs(residSteepLin) + 1.0E-10);
      if (s_print_DogLeg || (doDogLeg_ && m_print_flag > 6)) {
	printf("  (%2d - % 10.3g)  % 15.8E  % 15.8E % 15.8E  % 15.8E\n", 0, alpha, sLen, residSteep, residSteepLin , relFit);
      }
    }

    for (int iteration = 0; iteration < (int) alphaT.size(); iteration++) {
      double alpha = alphaT[iteration];
      for (int i = 0; i < neq_; i++) {
	y1[i] = m_y_n[i] + (1.0 - alpha) * deltaX_CP_[i];
	y1[i] += alpha * Nuu_ * newtDir[i];
      }
      /*
       *  Calculate the residual that would result if y1[] were the new solution vector
       *  -> m_resid[] contains the result of the residual calculation
       */
      if (solnType_ != NSOLN_TYPE_STEADY_STATE) {
	doResidualCalc(time_curr, solnType_, y1, ydot1, Base_LaggedSolutionComponents);
      } else {
	doResidualCalc(time_curr, solnType_, y1, ydot0, Base_LaggedSolutionComponents);
      }  

      for (int i = 0; i < neq_; i++) {
	y1[i] -= m_y_n[i];
      }
      sLen = solnErrorNorm(DATA_PTR(y1));

      double residSteep = residErrorNorm(DATA_PTR(m_resid));
      double residSteepLin = expectedResidLeg(1, alpha);

      double relFit = (residSteep - residSteepLin) / (fabs(residSteepLin) + 1.0E-10);
      if (s_print_DogLeg || (doDogLeg_ && m_print_flag > 6)) {
	printf("  (%2d - % 10.3g) % 15.8E   % 15.8E  % 15.8E  % 15.8E\n", 1, alpha, sLen,  residSteep, residSteepLin , relFit);
      }
    }

    for (int iteration = 0; iteration < (int) alphaT.size(); iteration++) {
      double alpha = alphaT[iteration];
      for (int i = 0; i < neq_; i++) {
	y1[i] = m_y_n[i] + ( Nuu_ + alpha * (1.0 - Nuu_))* newtDir[i];
      }
      sLen = ( Nuu_ + alpha * (1.0 - Nuu_)) * solnErrorNorm(DATA_PTR(newtDir));
      /*
       *  Calculate the residual that would result if y1[] were the new solution vector
       *  -> m_resid[] contains the result of the residual calculation
       */
      if (solnType_ != NSOLN_TYPE_STEADY_STATE) {
	doResidualCalc(time_curr, solnType_, y1, ydot1, Base_LaggedSolutionComponents);
      } else {
	doResidualCalc(time_curr, solnType_, y1, ydot0, Base_LaggedSolutionComponents);
      }  



      double residSteep = residErrorNorm(DATA_PTR(m_resid));
      double residSteepLin = expectedResidLeg(2, alpha);

      double relFit = (residSteep - residSteepLin) / (fabs(residSteepLin) + 1.0E-10);
      if (s_print_DogLeg || (doDogLeg_ && m_print_flag > 6)) {
	printf("  (%2d - % 10.3g)  % 15.8E % 15.8E  % 15.8E  % 15.8E\n", 2, alpha, sLen, residSteep, residSteepLin , relFit);
      }
    }


  }
  //====================================================================================================================
  double  NonlinearSolver::trustRegionLength() const
  {
    double dlen = solnErrorNorm(DATA_PTR(deltaX_trust_));
    return (trustDelta_ * dlen);
  }
  //====================================================================================================================
  void  NonlinearSolver::setDefaultDeltaBoundsMagnitudes()
  {
    for (int i = 0; i < neq_; i++) {
      m_deltaStepMinimum[i] = 1000. * atolk_[i];
      m_deltaStepMinimum[i] = MAX(m_deltaStepMinimum[i], 0.1 * fabs(m_y_n[i]));
    }
  }
  //====================================================================================================================
  void  NonlinearSolver::setDeltaBoundsMagnitudes(const doublereal * const deltaStepMinimum)
  {
    
    for (int i = 0; i < neq_; i++) {
      m_deltaStepMinimum[i] = deltaStepMinimum[i];
    }
    m_manualDeltaStepSet = 1;
  }
  //====================================================================================================================
  /*
   *
   * Return the factor by which the undamped Newton step 'step0'
   * must be multiplied in order to keep the update within the bounds of an accurate jacobian.
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
  double
  NonlinearSolver::deltaBoundStep(const doublereal * const y, const doublereal * const step0, const int loglevel)  {
			       
    int i_fbounds = 0;
    int ifbd = 0;
    int i_fbd = 0;
    
    doublereal sameSign = 0.0;
    doublereal ff;
    doublereal f_delta_bounds = 1.0;
    doublereal ff_alt;
    for (int i = 0; i < neq_; i++) {
      doublereal y_new = y[i] + step0[i];
      sameSign = y_new * y[i];
     
      /*
       * Now do a delta bounds
       * Increase variables by a factor of 1.5 only
       * decrease variables by a factor of 2 only
       */
      ff = 1.0;


      if (sameSign >= 0.0) {
	if ((fabs(y_new) > 1.5 * fabs(y[i])) && 
	    (fabs(y_new - y[i]) > m_deltaStepMinimum[i])) {
	  ff = 0.5 * fabs(y[i]/(y_new - y[i]));
	  ff_alt = fabs(m_deltaStepMinimum[i] / (y_new - y[i]));
	  ff = MAX(ff, ff_alt);
	  ifbd = 1;
	}
	if ((fabs(2.0 * y_new) < fabs(y[i])) &&
	    (fabs(y_new - y[i]) > m_deltaStepMinimum[i])) {
	  ff = y[i]/(y_new - y[i]) * (1.0 - 2.0)/2.0;
	  ff_alt = fabs(m_deltaStepMinimum[i] / (y_new - y[i]));
	  ff = MAX(ff, ff_alt);
	  ifbd = 0;
	}
      } else {
	/*
	 *  This handles the case where the value crosses the origin.
	 *       - First we don't let it cross the origin until its shrunk to the size of m_deltaBoundsMagnitudes[i]
	 */
	if (fabs(y[i]) > m_deltaStepMinimum[i]) {
	  ff = y[i]/(y_new - y[i]) * (1.0 - 2.0)/2.0;
	  ff_alt = fabs(m_deltaStepMinimum[i] / (y_new - y[i]));
	  ff = MAX(ff, ff_alt);
	  if (y[i] >= 0.0) {
	    ifbd = 0;
	  } else {
	    ifbd = 1;
	  }
	}
	/*
	 *  Second when it does cross the origin, we make sure that its magnitude is only 50% of the previous value.
	 */
	else if (fabs(y_new) > 0.5 * fabs(y[i])) {
	  ff = y[i]/(y_new - y[i]) * (-1.5);
	  ff_alt = fabs(m_deltaStepMinimum[i] / (y_new - y[i]));
	  ff = MAX(ff, ff_alt);
	  ifbd = 0;
	}
      }

      if (ff < f_delta_bounds) {
	f_delta_bounds = ff;
	i_fbounds = i;
	i_fbd = ifbd;
      }
    
    
    }
 
  
    /*
     * Report on any corrections
     */
    if (loglevel > 1) {
      if (f_delta_bounds < 1.0) {
	if (i_fbd) {
	  printf("\t\tdeltaBoundStep: Increase of Variable %d causing "
		 "delta damping of %g\n",
		 i_fbounds, f_delta_bounds);
	} else {
	  printf("\t\tdeltaBoundStep: Decrease of variable %d causing"
		 "delta damping of %g\n",
		 i_fbounds, f_delta_bounds);
	}
      }
    }
    
    
    return f_delta_bounds;
  }
 //====================================================================================================================
  //! Calculate the trust region vectors
  /*!
   *  The trust region is made up of the trust region vector calculation and the trustDelta_ value
   *  We periodically recalculate the trustVector_ values so that they renormalize to the
   *  correct length.
   */
  void  NonlinearSolver::calcTrustVector()  
  {
    double wtSum = 0.0;
    for (int i = 0; i < neq_; i++) {
      wtSum += m_ewt[i];
    }
    wtSum /= neq_;
    double trustNorm = solnErrorNorm(DATA_PTR(deltaX_trust_));
    double trustNormGoal = trustNorm * trustDelta_;

    // This is the size of each component.
    double trustDeltaEach = trustDelta_ * trustNorm / neq_;
    double oldVal; 
    double fabsy;
    // we use the old value of the trust region as an indicator
    for (int i = 0; i < neq_; i++) {
      oldVal = deltaX_trust_[i];
      fabsy = fabs(m_y_n[i]);
      // First off make sure that each trust region vector is 1/2 the size of each variable or smaller
      // unless overridden by the deltaStepMininum value. 
      if (oldVal > 0.5 * fabsy) {
	if (fabsy > m_deltaStepMinimum[i]) {
	  deltaX_trust_[i] = 0.5 * fabsy;
	} else {
	  deltaX_trust_[i] = m_deltaStepMinimum[i];
	}
      } else {
        double newValue =  trustDeltaEach * m_ewt[i] / wtSum;
	if (newValue > 4.0 * oldVal) { 
	  newValue = 4.0 * oldVal;
	} else if (newValue < 0.25 * oldVal) {
	  newValue = 0.25 * oldVal;
	  if (deltaX_trust_[i] < m_deltaStepMinimum[i]) {
	    newValue =  m_deltaStepMinimum[i];
	  }
	}
	deltaX_trust_[i] = newValue;
	if (deltaX_trust_[i] > 0.75 * m_deltaStepMaximum[i]) {
	  deltaX_trust_[i] = 0.75 * m_deltaStepMaximum[i];
	}
      }
    }


    // Final renormalization. 
    trustNorm = solnErrorNorm(DATA_PTR(deltaX_trust_));
    double  sum = trustNormGoal / trustNorm;
    for (int i = 0; i < neq_; i++) {
      deltaX_trust_[i] = deltaX_trust_[i] * sum;
    }
    trustDelta_ = 1.0;
    if (s_print_DogLeg || (doDogLeg_ && m_print_flag > 3)) {
      printf("calcTrustVector():  Trust vector size (SolnNorm Basis) changed from %g to %g \n",
	     trustNorm, trustNormGoal);
    }
  } 
  //====================================================================================================================
  void  NonlinearSolver::initializeTrustRegion() 
  {
    double cpd = calcTrustDistance(deltaX_CP_);
    if (s_print_DogLeg || (doDogLeg_ && m_print_flag > 3)) {
      printf("Relative Distance of Cauchy Vector wrt Trust Vector = %g\n", cpd);
    }
    trustDelta_ = trustDelta_ * cpd;
    calcTrustVector();
    cpd = calcTrustDistance(deltaX_CP_);
    if (s_print_DogLeg || (doDogLeg_ && m_print_flag > 3)) {
      printf("Relative Distance of Cauchy Vector wrt Trust Vector = %g\n", cpd);
    }
    trustDelta_ = trustDelta_ * cpd;
    calcTrustVector();
    cpd = calcTrustDistance(deltaX_CP_);
    if (s_print_DogLeg || (doDogLeg_ && m_print_flag > 3)) {
      printf("Relative Distance of Cauchy Vector wrt Trust Vector = %g\n", cpd);
    }
  }

  //====================================================================================================================
  // Fill a dogleg solution step vector
  /*
   *   Previously, we have filled up deltaX_Newton_[], deltaX_CP_[], and Nuu_, so that
   *   this routine is straightforward.
   *
   *  @param leg      Leg of the dog leg you are on (0, 1, or 2)
   *  @param alpha    Relative length along the dog length that you are on.
   *  @param deltaX   Vector to be filled up
   */
  void  NonlinearSolver::fillDogLegStep(int leg, double alpha, std::vector<doublereal>  & deltaX) const {
    if (leg == 0) {
      for (int i = 0; i < neq_; i++) {
	deltaX[i] = alpha * deltaX_CP_[i];
      }
    } else if (leg == 2) {
      for (int i = 0; i < neq_; i++) {
	deltaX[i] = (alpha + (1.0 - alpha) * Nuu_) * deltaX_Newton_[i];
      }
    } else {
      for (int i = 0; i < neq_; i++) {
	deltaX[i] = deltaX_CP_[i] * (1.0 - alpha) + alpha * Nuu_ * deltaX_Newton_[i];
      }
    }
  }
  //====================================================================================================================
  // Calculate the trust distance of a step in the solution variables
  /*
   *  The trust distance is defined as the length of the step according to the norm wrt to the trust region.
   *  We calculate the trust distance by the following method
   *  
   *      trustDist =  || delta_x   dot  1/trustDeltaX_ || / trustDelta_
   *
   * @param deltaX  Current value of deltaX
   */
  doublereal  NonlinearSolver::calcTrustDistance(std::vector<doublereal> const & deltaX) const
  {
    doublereal sum = 0.0;
    doublereal tmp = 0.0;
    for (int i = 0; i < neq_; i++) {
      tmp = deltaX[i] / deltaX_trust_[i];
      sum += tmp * tmp;
    }
    sum = sqrt(sum / neq_) / trustDelta_;
    return sum;
  }
  //====================================================================================================================
  int  NonlinearSolver::calcTrustIntersection(double trustDelta, double &lambda, double &alpha) const
  {
    double dist;
    if (normTrust_Newton_ < trustDelta) {
      lambda = 1.0;
      alpha = 1.0;
      return 2;
    }
   
    if (normTrust_Newton_ * Nuu_ < trustDelta) {
      alpha =  (trustDelta - normTrust_Newton_ * Nuu_) / (normTrust_Newton_ - normTrust_Newton_ * Nuu_);
      dist = dist_R0_ + dist_R1_ + alpha * dist_R2_;
      lambda = dist / dist_Total_;
      return 2;
    }
    if (normTrust_CP_ > trustDelta) {
      lambda = 1.0;
      dist = dist_R0_ * trustDelta / normTrust_CP_;
      lambda = dist / dist_Total_;
      alpha = trustDelta / normTrust_CP_;
      return 0;
    }
    double sumv = 0.0;
    for (int i = 0; i < neq_; i++) {
      sumv += (deltaX_Newton_[i] / deltaX_trust_[i]) * (deltaX_CP_[i] / deltaX_trust_[i]);
    }

    double a = normTrust_Newton_ * normTrust_Newton_ * Nuu_ * Nuu_;
    double b = 2.0 * Nuu_ * sumv;
    double c = normTrust_CP_ *  normTrust_CP_  - trustDelta * trustDelta;

    alpha =( -b + sqrt( b * b - 4.0 * a * c)) / (2.0 * a);
   

    dist = dist_R0_ + alpha * dist_R1_;
    lambda = dist / dist_Total_;
    return 1; 
  }
  //====================================================================================================================  
  /*
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
  doublereal NonlinearSolver::boundStep(const doublereal * const y, const doublereal * const step0, 
					const int loglevel) {
    int i, i_lower = -1;
    doublereal fbound = 1.0, f_bounds = 1.0;
    doublereal ff, y_new;
    
    for (i = 0; i < neq_; i++) {
      y_new = y[i] + step0[i];
      /*
       * Force the step to only take 80% a step towards the lower bounds
       */
      if (step0[i] < 0.0) {
	if (y_new < (y[i] + 0.8 * (m_y_low_bounds[i] - y[i]))) {
	  doublereal legalDelta = 0.8*(m_y_low_bounds[i] - y[i]);
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
	if (y_new > (y[i] + 0.8 * (m_y_high_bounds[i] - y[i]))) {
	  doublereal legalDelta = 0.8*(m_y_high_bounds[i] - y[i]);
	  ff = legalDelta / step0[i];
	  if (ff < f_bounds) {
	    f_bounds = ff;
	    i_lower = i;
	  }
	}
      }
    
    }

  
    /*
     * Report on any corrections
     */
    if (loglevel > 1) {
      if (f_bounds != 1.0) {
	printf("\t\tboundStep: Variable %d causing bounds damping of %g\n", i_lower, f_bounds);
      }
    }

    doublereal f_delta_bounds = deltaBoundStep(y, step0, loglevel);
    fbound = MIN(f_bounds, f_delta_bounds);

    return fbound;
  }
  //====================================================================================================================
  /*
   *
   * dampStep():
   *
   * On entry, step0 must contain an undamped Newton step to the
   * current solution y0. This method attempts to find a damping coefficient
   * such that the next undamped step would have a norm smaller than
   * that of step0. If successful, the new solution after taking the
   * damped step is returned in y1, and the undamped step at y1 is
   * returned in step1.
   *
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
   *          -2 Unsuccessful step.
   */
  int NonlinearSolver::dampStep(const doublereal time_curr, const double* y0, 
				const doublereal *ydot0, const double* step0, 
				double* const y1, double* const ydot1, double* step1,
				double& s1, SquareMatrix& jac, int& loglevel, bool writetitle,
				int& num_backtracks) 
  {
    int info = 0;
    int retnTrial = -2;
    // Compute the weighted norm of the undamped step size step0
    doublereal s0 = solnErrorNorm(step0);

    // Compute the multiplier to keep all components in bounds
    // A value of one indicates that there is no limitation
    // on the current step size in the nonlinear method due to
    // bounds constraints (either negative values of delta
    // bounds constraints.
    m_dampBound = boundStep(y0, step0, loglevel);

    // if fbound is very small, then y0 is already close to the
    // boundary and step0 points out of the allowed domain. In
    // this case, the Newton algorithm fails, so return an error
    // condition.
    if (m_dampBound < 1.e-30) {
      if (loglevel > 1) printf("\t\t\tdampStep: At limits.\n");
      return -3;
    }

    //--------------------------------------------
    //           Attempt damped step
    //-------------------------------------------- 

    // damping coefficient starts at 1.0
    m_dampRes = 1.0;
    int j, m;
    doublereal ff =  m_dampBound;
    num_backtracks = 0;
    for (m = 0; m < NDAMP; m++) {

      ff = m_dampBound * m_dampRes;

      // step the solution by the damped step size
      /*
       * Whenever we update the solution, we must also always
       * update the time derivative.
       */
      for (j = 0; j < neq_; j++) {
	y1[j] = y0[j] + ff * step0[j];
      }
      
      if (solnType_ != NSOLN_TYPE_STEADY_STATE) {
	calc_ydot(m_order, y1, ydot1);
      }
      /*
       *  Calculate the residual that would result if y1[] were the new solution vector
       *  -> m_resid[] contains the result of the residual calculation
       */
      if (solnType_ != NSOLN_TYPE_STEADY_STATE) {
	info = doResidualCalc(time_curr, solnType_, y1, ydot1, Base_LaggedSolutionComponents);
      } else {
	info = doResidualCalc(time_curr, solnType_, y1, ydot0, Base_LaggedSolutionComponents);
      }
      if (info != 1) {
	if (loglevel > 0) {
	  printf("\t\t\tdampStep: current trial step and damping led to Residual Calc ERROR %d. Bailing\n", info);
	}
	return -1;
      }
      m_normResidTrial = residErrorNorm(DATA_PTR(m_resid));

      bool steepEnough = (m_normResidTrial <  m_normResid0 * (0.9 * (1.0 - ff) * (1.0 - ff)* (1.0 - ff) + 0.1));

      if (m_normResidTrial < 1.0 || steepEnough) {
	if (loglevel >= 5) {
	  if (m_normResidTrial < 1.0) {
	    printf("\t  dampStep(): Current trial step and damping"
		   " coefficient accepted because residTrial test step < 1:\n");
	    printf("\t              resid0 = %g, residTrial = %g\n", m_normResid0, m_normResidTrial);
	  } else if (steepEnough) {
	    printf("\t  dampStep(): Current trial step and damping"
		   " coefficient accepted because resid0 > residTrial and steep enough:\n");
	    printf("\t              resid0 = %g, residTrial = %g\n", m_normResid0, m_normResidTrial);
	  } else {
	    printf("\t  dampStep(): Current trial step and damping"
		   " coefficient accepted because residual solution damping is turned off:\n");
	    printf("\t              resid0 = %g, residTrial = %g\n", m_normResid0, m_normResidTrial);
	  }
	}
	/*
	 *  We aren't going to solve the system if we don't need to. Therefore, return an estimate
	 *  of the next solution update based on the ratio of the residual reduction.
	 */
	if (m_normResid0 > 0.0) {
	  s1 = s0 * m_normResidTrial / m_normResid0;
	}
	else {
	  s1 = 0;
	}
	if (m_normResidTrial < 1.0) {
	  retnTrial = 3;
	} else {
	  retnTrial = 4;
	}
	break;
      }
      
      // Compute the next undamped step, step1[], that would result  if y1[] were accepted.
      //   We now have two steps that we have calculated step0[] and step1[]
      if (solnType_ != NSOLN_TYPE_STEADY_STATE) {
	info = doNewtonSolve(time_curr, y1, ydot1, step1, jac, loglevel);
      } else {
	info = doNewtonSolve(time_curr, y1, ydot0, step1, jac, loglevel);
      }
      if (info) {
	if (loglevel > 0) {
	  printf("\t\t\tdampStep: current trial step and damping led to LAPACK ERROR %d. Bailing\n", info);
	}
	return -1;
      }

      // compute the weighted norm of step1
      s1 = solnErrorNorm(step1);

      // write log information
      if (loglevel > 3) {
	print_solnDelta_norm_contrib((const doublereal *) step0, 
				     "DeltaSoln",
				     (const doublereal *) step1,
				     "DeltaSolnTrial",
				     "dampNewt: Important Entries for "
				     "Weighted Soln Updates:",
				     y0, y1, ff, 5);
      }
      if (loglevel > 1) {
	printf("\t\t\tdampStep(): s0 = %g, s1 = %g, dampBound = %g,"
	       "dampRes = %g\n",  s0, s1, m_dampBound, m_dampRes);
      }


      // if the norm of s1 is less than the norm of s0, then
      // accept this damping coefficient. Also accept it if this
      // step would result in a converged solution. Otherwise,
      // decrease the damping coefficient and try again.
	  
      if (s1 < 0.8 || s1 < s0) {
	if (s1 < 1.0) {
	  if (loglevel > 2) {   
	    if (s1 < 1.0) {
	      printf("\t\t\tdampStep: current trial step and damping"
		     " coefficient accepted because test step < 1\n");
	      printf("\t\t\t          s1 = %g, s0 = %g\n", s1, s0);
	    }
	  }  
	  retnTrial = 2;
	} else {
	  retnTrial = 1;
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
      m_dampRes /= DampFactor;
    }

    // If a damping coefficient was found, return 1 if the
    // solution after stepping by the damped step would represent
    // a converged solution, and return 0 otherwise. If no damping
    // coefficient could be found, return -2.
    if (m < NDAMP) {
      if (loglevel >= 4 ) {
	printf("\t  dampStep(): current trial step accepted retnTrial = %d, its = %d, damp = %g\n", retnTrial, m+1, ff);
      }
      return retnTrial;
    } else {
      if (s1 < 0.5 && (s0 < 0.5)) {
	if (loglevel >= 4 ) {
	  printf("\t  dampStep(): current trial step accepted kindof retnTrial = %d, its = %d, damp = %g\n", 2, m+1, ff);
	}
	return 2;
      }
      if (s1 < 1.0) {
	if (loglevel >= 4 ) {
	  printf("\t  dampStep(): current trial step accepted and soln converged retnTrial = %d, its = %d, damp = %g\n", 0, m+1, ff);
	}
	return 0;
      }
    }
    if (loglevel >= 4 ) {
      printf("\t  dampStep(): current direction is rejected! retnTrial = %d, its = %d, damp = %g\n", -2, m+1, ff);
    }
    return -2;
  }

  //====================================================================================================================

  int NonlinearSolver::dampDogLeg(const doublereal time_curr, const double* y0, 
				  const doublereal *ydot0, std::vector<doublereal> & step0,
				  double* const y_new, double* const ydot_new, double* step1,
				  double& s1, SquareMatrix& jac, int& loglevel, bool writetitle,
				  int& num_backtracks) 
  {
    double lambda;
    double alpha;
    int info;
    int leg = 2;
    bool success = false;
    int retn = 0;
    bool haveASuccess = false;
    double trustDeltaOld = trustDelta_;
    //--------------------------------------------
    //           Attempt damped step
    //-------------------------------------------- 

    // damping coefficient starts at 1.0
    m_dampRes = 1.0;
    int j, m;
    num_backtracks = 0;
    //double deltaSolnNorm = solnErrorNorm(DATA_PTR(deltaX_CP_));
    //double funcDecreaseSDExp = RJd_norm_ / deltaSolnNorm * lambda_;
    double tlen;
   

    for (m = 0; m < NDAMP; m++) {
      /*
       *  Find the initial value of lambda that satisfies the trust distance, trustDelta_
       */
      leg = calcTrustIntersection(trustDelta_, lambda, alpha);
      
      if (loglevel > 5) {
	tlen = trustRegionLength();
	printf("\tdampDogLeg: trust region with length %13.5E has intersection at leg = %d, alpha = %g\n",
	       tlen, leg, alpha);
      }
      /*
       *  Figure out the new step vector, step0, based on (leg, alpha)
       */
      fillDogLegStep(leg, alpha, step0);

      /*
       *  Bound the step 
       */
      m_dampBound = boundStep(y0, DATA_PTR(step0), loglevel);
      /*
       * Decrease the step length if we are bound
       */
      if (m_dampBound < 1.0) {
	for (j = 0; j < neq_; j++) {
	  step0[j] =  step0[j] *  m_dampBound;
	}
      }
    
      /*
       * OK, we have the step0. Now, ask the question whether it satisfies the acceptance criteria
       * as a good step. Also, make sure that it stays within bounds.
       */
      info = decideStep(time_curr, leg, alpha, y0, ydot0,  step0, y_new, ydot_new,  loglevel, trustDeltaOld);

      /*
       *  The algorithm failed to find a solution vector sufficiently different than the current point
       */
      if (info == -1) {
	if (loglevel >= 1) {
	  double stepNorm =  solnErrorNorm(DATA_PTR(step0));
	  printf("\t\t\tdampDogLeg: Current direction rejected, update became too small %g\n", stepNorm);
	  success = false;
	  retn = -1;
	  break;
	}
      }
      if (info == -2) {
	if (loglevel >= 1) {
	  printf("\t\t\tdampStep: current trial step and damping led to LAPACK ERROR %d. Bailing\n", info);
	  success = false;
	  retn = -1;
	  break;
	}
      }
      if (info == 0) {
	success = true;
	break;
      }
      if (info == 3) {
	haveASuccess = true;
	// Store the good results in step1
	mdp::mdp_copy_dbl_1(DATA_PTR(step1), CONSTD_DATA_PTR(step0), neq_);

      }
      if (info == 2) {
	if (haveASuccess) {
	  mdp::mdp_copy_dbl_1(DATA_PTR(step0), CONSTD_DATA_PTR(step1), neq_);
	  for (j = 0; j < neq_; j++) {
	    y_new[j] = y0[j] + step0[j];
	  }    
	  if (solnType_ != NSOLN_TYPE_STEADY_STATE) {
	    calc_ydot(m_order, y_new, ydot_new);
	  } 
	  success = true;
	  break;
	}
      }


    }

    /*
     * Estimate s1, the norm after the next step
     */
    double stepNorm =  solnErrorNorm(DATA_PTR(step1));
    if (  m_dampBound < 1.0) {
      stepNorm /= m_dampBound;
    }
    stepNorm *=  m_normResidTrial / m_normResid0;
    s1 = stepNorm;

    if (success) {
      if (m_normResidTrial < 1.0) {
	return 1;
      }
      return 0;
    }
    return -1;
  }  
  //====================================================================================================================
  // Decide whether the current step is acceptable and adjust the trust region size
  /*
   *  This is an extension of algorithm 6.4.5 of Dennis and Schnabel.
   *
   *  Here we decide whether to accept the current step
   *  At the end of the calculation a new estimate of the trust region is calculated
   *
   * @param time_curr  INPUT     Current value of the time
   * @param leg        INPUT    Leg of the dogleg that we are on
   * @param alpha      INPUT    Distance down that leg that we are on
   * @param y0         INPUT    Current value of the solution vector
   * @param ydot0      INPUT    Current value of the derivative of the solution vector
   * @param step0      INPUT    Trial step
   * @param y1         OUTPUT   Solution values at the conditions which are evalulated for success
   * @param ydot1      OUTPUT   Time derivates of solution at the conditions which are evalulated for success
   * @param loglevel   INPUT    Current loglevel
   * @param trustDeltaOld INPUT Value of the trust length at the old conditions
   *
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
  int NonlinearSolver::decideStep(const doublereal time_curr, int leg, double alpha, const double* y0, 
				  const doublereal *ydot0,  std::vector<doublereal> & step0,
				  double * const y1, double* const ydot1,  int& loglevel, double trustDeltaOld) 

  {
    int retn = 2;
    bool goodStep = false;
    int j;
    int info;
    double ll;
    // Calculate the solution step length
    double stepNorm =  solnErrorNorm(DATA_PTR(step0));

    // Calculate the initial (R**2 * neq) value for the old function
    double normResid0_2 = m_normResid0 * m_normResid0 * neq_;

    // Calculate the distance to the cauchy point
    double cauchyDistanceNorm = solnErrorNorm(DATA_PTR(deltaX_CP_));

    // This is the expected inital rate of decrease in the cauchy direction.
    //   -> This is Eqn. 29 = Rhat dot Jhat dy / || d ||
    double funcDecreaseSDExp = RJd_norm_ / cauchyDistanceNorm * lambda_;
    if (funcDecreaseSDExp > 0.0) {
      if (loglevel > 0) {
	printf("\t\tdecideStep(): Unexpected condition -> cauchy slope is positive\n");
      }
    }

    /*
     *  Calculate the newsolution value y1[] given the step size
     */
    for (j = 0; j < neq_; j++) {
      y1[j] = y0[j] + step0[j];
    }
    /*
     *  Calculate the new solution time derivative given the step size
     */
    if (solnType_ != NSOLN_TYPE_STEADY_STATE) {
      calc_ydot(m_order, y1, ydot1);
    } 
        
    /*
     *  Calculate the residual that would result if y1[] were the new solution vector
     *  -> m_resid[] contains the result of the residual calculation
     */
    if (solnType_ != NSOLN_TYPE_STEADY_STATE) {
      info = doResidualCalc(time_curr, solnType_, y1, ydot1, Base_LaggedSolutionComponents);
    } else {
      info = doResidualCalc(time_curr, solnType_, y1, ydot0, Base_LaggedSolutionComponents);
    }

    if (info != 1) {
      if (loglevel > 0) {
	printf("\t\tdecideStep: current trial step and damping led to Residual Calc ERROR %d. Bailing\n", info);
      }
      return -2;
    }
    /*
     *  Ok we have a successful new residual. Calculate the normalized residual value and store it in 
     *  m_normResidTrial
     */
    m_normResidTrial = residErrorNorm(DATA_PTR(m_resid));
    double normResidTrial_2 = neq_ *  m_normResidTrial *  m_normResidTrial;

    /*
     *  We have a minimal acceptance test for passage. deltaf < 1.0E-4 (CauchySlope) (deltS)
     *  This is the condition that D&S use in 6.4.5
     */
    double funcDecrease = 0.5 * (normResidTrial_2 - normResid0_2);
    double acceptableDelF =  funcDecreaseSDExp * stepNorm * 1.0E-4;
    if (funcDecrease < acceptableDelF) {
      goodStep = true;
      retn = 0;
    } else {
      trustDelta_ *= 0.33;
      retn = 2;
      // error condition if step is getting too small
      if (stepNorm * .5 < 0.2) {
	retn = -1;
      }
      return retn;
    }
    /* 
     *  Figure out the next trust region. We are here iff retn = 0
     *
     *      If we had to bounds delta the update, decrease the trust region
     */
    if (m_dampBound < 1.0) {
      trustDelta_ *= 0.5;
      ll = trustRegionLength();
      printf("\t\tdecideStep(): Trust region decreased from %g to %g due to bounds constraint\n",
	     ll*2, ll);
    } else {
      retn = 0;
      /*
       *  Calculate the expected residual from the quadratic model
       */
      double expectedNormRes = expectedResidLeg(leg, alpha);
      double expectedFuncDecrease = 0.5 * (neq_ * expectedNormRes * expectedNormRes - normResid0_2);
      if (funcDecrease > 0.1 * expectedFuncDecrease) {
	if ((m_normResidTrial > 0.5 * m_normResid0) && (m_normResidTrial > 0.1)) {
	  trustDelta_ *= 0.5;
	  ll = trustRegionLength();
	  printf("\t\tdecideStep(): Trust region decreased from %g to %g due to bad quad approximation\n",
		 ll*2, ll);
	}
      } else {
	/*
	 *  If we are doing well, consider increasing the trust region and recalculating
	 */
	if (funcDecrease < 0.8 * expectedFuncDecrease ||  (m_normResidTrial < 0.33 * m_normResid0)) {
	  if (trustDelta_ <= trustDeltaOld) {
	    trustDelta_ *= 2.0;
	    ll = trustRegionLength();
	    printf("\td\tecideStep(): Trust region increased from %g to %g due to good quad approximation\n",
		   ll*0.5, ll);
	    retn = 3;
	  } else {
	    /*
	     *  Increase the size of the trust region for the next calculation
	     */
	    if (m_normResidTrial < 0.75 * expectedNormRes) {
	      trustDelta_ *= 2.0;
	      ll = trustRegionLength();
	      printf("\t\tdecideStep(): Trust region further increased from %g to %g due to good nonlinear behavior\n",
		     ll*0.5, ll);
	    }
	  }
	}
      }
    }
    return retn;
  }
  //====================================================================================================================
  /*
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
   *  @return  A positive value indicates a successful convergence
   *           -1  Failed convergence
   */
  int NonlinearSolver::solve_nonlinear_problem(int SolnType, double* y_comm,
					       double* ydot_comm, doublereal CJ,
					       doublereal time_curr, 
					       SquareMatrix& jac,
					       int &num_newt_its,  int &num_linear_solves,
					       int &num_backtracks,  int loglevelInput)
  {
    clockWC wc;
    int convRes = 0;
    solnType_ = SolnType;
    int info = 0;

    bool m_residCurrent = false;
    int m = 0;
    bool forceNewJac = false;
    doublereal s1=1.e30;


    //  std::vector<doublereal> y_curr(neq_, 0.0); 
    std::vector<doublereal> ydot_curr(neq_, 0.0);
    std::vector<doublereal> stp(neq_, 0.0);
    std::vector<doublereal> stp1(neq_, 0.0);

    std::vector<doublereal> y_new(neq_, 0.0);
   
    mdp::mdp_copy_dbl_1(DATA_PTR(m_y_n), DATA_PTR(y_comm), neq_);
  
    if (SolnType != NSOLN_TYPE_STEADY_STATE || ydot_comm) {
      mdp::mdp_copy_dbl_1(DATA_PTR(ydot_curr), ydot_comm, neq_);
      mdp::mdp_copy_dbl_1(DATA_PTR(ydot_new), ydot_comm, neq_);
    }
    // Redo the solution weights every time we enter the function
    createSolnWeights(DATA_PTR(m_y_n));
    m_normDeltaSoln_Newton = 1.0E1;
    bool frst = true;
    num_newt_its = 0;
    num_linear_solves = - m_numTotalLinearSolves;
    num_backtracks = 0;
    int i_backtracks;
    m_print_flag = loglevelInput;
    if (m_print_flag > 1) {
      jac.m_printLevel = 1;
    } else {
      jac.m_printLevel = 0;
    }


    if (m_print_flag == 2 || m_print_flag == 3) {
      
      printf("\tsolve_nonlinear_problem():\n\n");
      printf("\t    Iter Resid NewJac |  LinearIts  Ax-b  | Fbound     Fdamp DampIts |   DeltaSolnF     ResidF\n");
      printf("\t-------------------------------------------------------------------------------------------------\n");
 
    }



    while (1 > 0) {

      /*
       * Increment Newton Solve counter
       */
      m_numTotalNewtIts++;
      num_newt_its++;

      /*
       *  If we are far enough away from the solution, redo the solution weights and the trust vectors.
       */
      if (m_normDeltaSoln_Newton > 1.0E2) {
	createSolnWeights(DATA_PTR(m_y_n));
#ifdef DEBUG_DOGLEG
	calcTrustVector();
#else
	if (doDogLeg_) {
	  calcTrustVector();
	}
#endif
      } else {
	// Do this stuff every 5 iterations
        if ( (num_newt_its % 5) == 1) {
	  createSolnWeights(DATA_PTR(m_y_n));
#ifdef DEBUG_DOGLEG
	  calcTrustVector();
#else
	  if (doDogLeg_) {
	    calcTrustVector();
	  }
#endif
	}
      }

      //mdp::mdp_copy_dbl_1(DATA_PTR(m_y_n), DATA_PTR(y_curr), neq_);
      /*
       * Set default values of Delta bounds constraints
       */
      if (!m_manualDeltaStepSet) {
	setDefaultDeltaBoundsMagnitudes();
      }

      if (m_print_flag > 3) {
	printf("\tsolve_nonlinear_problem(): iteration %d:\n",
	       num_newt_its);
      }


      // Check whether the Jacobian should be re-evaluated.
            
      forceNewJac = true;
            
      if (forceNewJac) {
	if (m_print_flag > 3) {
	  printf("\tsolve_nonlinear_problem(): Getting a new Jacobian and solving system\n");
	}
	info = beuler_jac(jac, DATA_PTR(m_resid), time_curr, CJ,  DATA_PTR(m_y_n), DATA_PTR(ydot_curr), num_newt_its);
	if (info == 0) {
	  m = -1;
	  goto done;
	}
	m_residCurrent = true;
      } else {
	if (m_print_flag > 1) {
	  printf("\tsolve_nonlinear_problem(): Solving system with old jacobian\n");
	}
	m_residCurrent = false;
      }
      /*
       * Go get new scales
       */
      setColumnScales();



      info = doResidualCalc(time_curr, NSOLN_TYPE_STEADY_STATE, DATA_PTR(m_y_n), DATA_PTR(ydot_curr));
      if (info != 1) {
	if (m_print_flag > 0) {
	  printf("\t\t\tsolve_nonlinear_problem(): Residual Calc ERROR %d. Bailing\n", info);
	}
	m = -1;
	goto done;
      }

      /*
       * Scale the matrix and the rhs, if they aren't already scaled
       * Figure out and store the residual scaling factors.
       */
      scaleMatrix(jac, DATA_PTR(m_y_n), DATA_PTR(ydot_curr), time_curr);


      /*
       *  Optional print out the initial residual
       */
      if (m_print_flag >= 6) {
	m_normResid0 = residErrorNorm(DATA_PTR(m_resid), "Initial norm of the residual", 10, DATA_PTR(m_y_n));
      } else if (m_print_flag == 4 || m_print_flag == 5) {
	m_normResid0 = residErrorNorm(DATA_PTR(m_resid), "Initial norm of the residual", 0, DATA_PTR(m_y_n));
      } else {
	m_normResid0 = residErrorNorm(DATA_PTR(m_resid), "Initial norm of the residual", 0, DATA_PTR(m_y_n));
      }
      
#ifdef DEBUG_DOGLEG
      m_normDeltaSoln_CP = doCauchyPointSolve(jac);
      if (m_numTotalNewtIts == 1) {
	initializeTrustRegion();
      }
#else
      if (doDogLeg_) { 
	m_normDeltaSoln_CP = doCauchyPointSolve(jac);
	if (m_numTotalNewtIts == 1) {
	  initializeTrustRegion();
	}
      }
#endif

      // compute the undamped Newton step
      if (doAffineSolve_) {
	info = doAffineNewtonSolve(DATA_PTR(m_y_n), DATA_PTR(ydot_curr), DATA_PTR(deltaX_Newton_), jac);
      } else {
	info = doNewtonSolve(time_curr, DATA_PTR(m_y_n), DATA_PTR(ydot_curr), DATA_PTR(deltaX_Newton_), jac, m_print_flag);
      }

      if (info) {
	m = -1;
	goto done;
      }
      mdp::mdp_copy_dbl_1(DATA_PTR(stp), CONSTD_DATA_PTR(deltaX_Newton_), neq_);

      if (m_print_flag > 3) {
	m_normDeltaSoln_Newton = solnErrorNorm(DATA_PTR(stp),  "Initial Undamped Step of the iteration", 10);
      } else {
	m_normDeltaSoln_Newton = solnErrorNorm(DATA_PTR(stp), "Initial Undamped Step of the iteration", 0);
      } 
   


      if (doDogLeg_) { 
#ifdef DEBUG_DOGLEG
	double trustD = calcTrustDistance(stp);
	if (s_print_DogLeg || m_print_flag > 3) {
	  if (trustD > trustDelta_) {
	    printf("newton's method trustD, %g, larger than trust region, %g\n", trustD, trustDelta_);
	  } else {
	    printf("newton's method trustD, %g, smaller than trust region, %g\n", trustD, trustDelta_);
	  }
	}
#endif
      }

      /*
       * Filter out bad directions
       */
      filterNewStep(time_curr, DATA_PTR(m_y_n), DATA_PTR(stp));


      
  
#ifdef DEBUG_DOGLEG
      descentComparison(time_curr, DATA_PTR(ydot_curr), DATA_PTR(ydot_new), DATA_PTR(stp));
#endif
    
 
      if (doDogLeg_) { 
	setupDoubleDogleg(DATA_PTR(stp));  
#ifdef DEBUG_DOGLEG
	residualComparisonLeg(time_curr, DATA_PTR(ydot_curr), DATA_PTR(ydot_new), DATA_PTR(stp));
#endif
	m =  dampDogLeg(time_curr, DATA_PTR(m_y_n), DATA_PTR(ydot_curr), 
			stp, DATA_PTR(y_new), DATA_PTR(ydot_new), 
			DATA_PTR(stp1), s1, jac, m_print_flag, frst, i_backtracks);
      }
#ifdef DEBUG_DOGLEG
      else {
	residualComparisonLeg(time_curr, DATA_PTR(ydot_curr), DATA_PTR(ydot_new), DATA_PTR(stp));
      }
#endif 

      // Damp the Newton step
      /*
       *  On return the recommended new solution and derivatisve is located in:
       *          y_new
       *          y_dot_new
       *  The update delta vector is located in
       *          stp1
       *  The estimate of the solution update norm for the next step is located in
       *          s1
       */
      if (!doDogLeg_) {
	m = dampStep(time_curr, DATA_PTR(m_y_n), DATA_PTR(ydot_curr), 
		     DATA_PTR(stp), DATA_PTR(y_new), DATA_PTR(ydot_new), 
		     DATA_PTR(stp1), s1, jac, m_print_flag, frst, i_backtracks);
	frst = false;
	num_backtracks += i_backtracks;
      }
      /*
       * Impose the minimum number of newton iterations critera
       */
      if (num_newt_its < m_min_newt_its) {
	if (m > 0) {
	  if (m_print_flag > 2) {
	    printf("\t  Damped Newton successful (m=%d) but minimum newton iterations not attained. Resolving ...\n", m);
	  }
	  m = 0;
	}
      }

      /*
       * Impose max newton iteration
       */
      if (num_newt_its > maxNewtIts_) {
	m = -1;
	if (m_print_flag > 1) {
	  printf("\t\tsolve_nonlinear_problem(): Damped newton unsuccessful (max newts exceeded) sfinal = %g\n", s1);
	}
      }


      info = doResidualCalc(time_curr, NSOLN_TYPE_STEADY_STATE, DATA_PTR(y_new), DATA_PTR(ydot_new));
      if (info != 1) {
	if (m_print_flag > 0) {
	  printf("\t\t\tdampStep: current trial step and damping led to Residual Calc ERROR %d. Bailing\n", info);
	}
	m = -1;
	goto done;
      }

      if (m_print_flag > 3) {
	residErrorNorm(DATA_PTR(m_resid), "Resulting Residual Norm", 10, DATA_PTR(y_new));
      }

      convRes = 0;
      if (m > 0) {
	convRes = convergenceCheck(m, s1);
      }
      
      if (m_print_flag >= 4) {
	if (convRes > 0) {
	  printf("\t  Damped Newton iteration successful, nonlin "
		 "converged, final estimate of the next solution update norm = %-12.4E\n", s1);
	} else if (m >= 0) {
	  printf("\t  Damped Newton iteration successful, "
		 "final estimate of the next solution update norm = %-12.4E\n", s1);
	} else {
	  printf("\t  Damped Newton unsuccessful, final estimate of the next solution update norm = %-12.4E\n", s1);
	}
      }



      bool m_filterIntermediate = false;
      if (m_filterIntermediate) {
	if (m == 0) {
	  (void) filterNewSolution(time_n, DATA_PTR(y_new), DATA_PTR(ydot_new));
	}
      }

      // Exchange new for curr solutions
      if (m >= 0) {
	mdp::mdp_copy_dbl_1(DATA_PTR(m_y_n), CONSTD_DATA_PTR(y_new), neq_);

	if (solnType_ != NSOLN_TYPE_STEADY_STATE) {
	  calc_ydot(m_order, DATA_PTR(m_y_n), DATA_PTR(ydot_curr));
	}
      }

      if (m_print_flag == 2 || m_print_flag == 3) {
	//   printf("\t    Iter Resid NewJac |  LinearIts  Ax-b  | Fbound     Fdamp DampIts |   DeltaSolnF     ResidF\n");

	printf("\t%4d %11.3E", num_newt_its, m_normResid0);
	bool m_jacAge = false;
	if (!m_jacAge) {
	  printf("   Y  |");
	} else {
	  printf("   N  |");
	}
	printf("%5d %11.3E  | %10.2E %10.2E %2d | %11.3E %11.3E ", m_numTotalLinearSolves, 0.0, m_dampBound, m_dampRes,
	       i_backtracks, m_normDeltaSoln_Newton, m_normResidFRaw);
	printf("\n");
      
      }
      // convergence
      if (convRes) {
	goto done;
      }

      // If dampStep fails, first try a new Jacobian if an old
      // one was being used. If it was a new Jacobian, then
      // return -1 to signify failure.
      else if (m < 0) {
	goto done;
      }
    }

  done:


    if (m_print_flag == 2 || m_print_flag == 3) {
      if (convRes > 0) {
	if (convRes == 3) {
	  printf("\t                      |                   |           converged = 3  |(%11.3E) \n", s1);
	} else {
	  printf("\t                      |                   |           converged = %1d  | %11.3E %11.3E \n", convRes,
		 s1, m_normResidTrial);
	}
      }
   
      printf("\t  --------------------------------------------------------------------------------------------\n");
      
    }
    

    mdp::mdp_copy_dbl_1(y_comm, CONSTD_DATA_PTR(m_y_n), neq_);
    if (solnType_ != NSOLN_TYPE_STEADY_STATE) {
      mdp::mdp_copy_dbl_1(ydot_comm, CONSTD_DATA_PTR(ydot_curr), neq_);
    }
 
    num_linear_solves += m_numTotalLinearSolves;
 
    doublereal time_elapsed =  wc.secondsWC();
    if (m_print_flag > 1) {
      if (m > 0) {
	if (NonlinearSolver::m_TurnOffTiming) {
	  printf("\t\tNonlinear problem solved successfully in %d its\n",
		 num_newt_its);
	} else {
	  printf("\t\tNonlinear problem solved successfully in %d its, time elapsed = %g sec\n",
		 num_newt_its, time_elapsed);
	}
      } else {
	printf("\t\tNonlinear problem failed to solve after %d its\n", num_newt_its);
      }
    }
    return m;
  }
  //====================================================================================================================
  // Print solution norm contribution
  /*
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
  void NonlinearSolver::
  print_solnDelta_norm_contrib(const doublereal * const solnDelta0,
			       const char * const s0,
			       const doublereal * const solnDelta1,
			       const char * const s1,
			       const char * const title,
			       const doublereal * const y0,
			       const doublereal * const y1,
			       doublereal damp,
			       int num_entries) {
    int i, j, jnum;
    bool used;
    doublereal dmax0, dmax1, error, rel_norm;
    printf("\t\t%s currentDamp = %g\n", title, damp);
    printf("\t\t     I     ysolnOld %13s ysolnNewRaw | ysolnNewTrial "
	   "%10s ysolnNewTrialRaw | solnWeight  wtDelSoln wtDelSolnTrial\n", s0, s1);
    int *imax = mdp::mdp_alloc_int_1(num_entries, -1);
    printf("\t\t   "); print_line("-", 125);
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
	printf("\t\t  %4d %12.4e %12.4e %12.4e |  %12.4e  %12.4e   %12.4e    |%12.4e %12.4e %12.4e\n",
	       i, y0[i], solnDelta0[i],  y0[i] + solnDelta0[i], y1[i],
	       solnDelta1[i], y1[i]+ solnDelta1[i], m_ewt[i], dmax0, dmax1);
      }
    }
    printf("\t\t   "); print_line("-", 125);
    mdp::mdp_safe_free((void **) &imax);
  }
  //====================================================================================================================
  //!  This routine subtracts two numbers for one another
  /*!
   *   This routine subtracts 2 numbers. If the difference is less
   *   than 1.0E-14 times the magnitude of the smallest number, then diff returns an exact zero. 
   *   It also returns an exact zero if the difference is less than
   *   1.0E-300.
   *
   *   returns:  a - b
   *
   *   This routine is used in numerical differencing schemes in order
   *   to avoid roundoff errors resulting in creating Jacobian terms.
   *   Note: This is a slow routine. However, jacobian errors may cause
   *         loss of convergence. Therefore, in practice this routine has proved cost-effective.
   *
   * @param a   Value of a
   * @param b   value of b
   *
   * @return    returns the difference between a and b
   */
  static inline doublereal subtractRD(doublereal a, doublereal b) {
    doublereal diff = a - b;
    doublereal d = MIN(fabs(a), fabs(b));
    d *= 1.0E-14;
    doublereal ad = fabs(diff);
    if (ad < 1.0E-300) {
      diff = 0.0;
    }
    if (ad < d) {
      diff = 0.0;
    }
    return diff;
  }
  //====================================================================================================================
  /*
   *
   *  Function called by BEuler to evaluate the Jacobian matrix and the
   *  current residual at the current time step.
   *  @param N = The size of the equation system
   *  @param J = Jacobian matrix to be filled in
   *  @param f = Right hand side. This routine returns the current
   *             value of the rhs (output), so that it does
   *             not have to be computed again.
   *  
   * @return Returns a flag to indicate that operation is successful.
   *            1  Means a successful operation
   *            0  Means an unsuccessful operation
   */
  int NonlinearSolver::beuler_jac(SquareMatrix &J, doublereal * const f,
				   doublereal time_curr, doublereal CJ,
				   doublereal * const y,
				   doublereal * const ydot,
				   int num_newt_its)
  {
    int i, j;
    double* col_j;
    int info;
    doublereal ysave, ydotsave, dy;
    int retn = 1;
    /*
     * Clear the factor flag
     */
    J.clearFactorFlag();
    if (m_jacFormMethod == NSOLN_JAC_ANAL) {
      /********************************************************************
       * Call the function to get a jacobian.
       */
      info = m_func->evalJacobian(time_curr, delta_t_n, y, ydot, J, f);
      m_nJacEval++;
      m_nfe++;
      if (info != 1) {
	return info;
      }
    }  else {
      /*******************************************************************
       * Generic algorithm to calculate a numerical Jacobian
       */
      /*
       * Calculate the current value of the rhs given the
       * current conditions.
       */

      info = m_func->evalResidNJ(time_curr, delta_t_n, y, ydot, f, JacBase_ResidEval);
      m_nfe++;
      if (info != 1) {
	return info;
      }
      m_nJacEval++;

      /*
       * Malloc a vector and call the function object to return a set of
       * deltaY's that are appropriate for calculating the numerical
       * derivative.
       */
      doublereal *dyVector = mdp::mdp_alloc_dbl_1(neq_, MDP_DBL_NOINIT);
      retn = m_func->calcDeltaSolnVariables(time_curr, y, ydot, dyVector, DATA_PTR(m_ewt));
      


      if (s_print_NumJac) {
	if (m_print_flag >= 7) {
	  if (neq_ < 20) {
	    printf("\t\tUnk            m_ewt              y                dyVector            ResN\n");
	    for (int iii = 0; iii < neq_; iii++){
	      printf("\t\t %4d       %16.8e   %16.8e   %16.8e  %16.8e \n",
		      iii,   m_ewt[iii],  y[iii], dyVector[iii], f[iii]);
	    }
	  }
	}
      }

      /*
       * Loop over the variables, formulating a numerical derivative
       * of the dense matrix.
       * For the delta in the variable, we will use a variety of approaches
       * The original approach was to use the error tolerance amount.
       * This may not be the best approach, as it could be overly large in
       * some instances and overly small in others.
       * We will first protect from being overly small, by using the usual
       * sqrt of machine precision approach, i.e., 1.0E-7,
       * to bound the lower limit of the delta.
       */
      for (j = 0; j < neq_; j++) {


        /*
         * Get a pointer into the column of the matrix
         */


        col_j = (doublereal *) J.ptrColumn(j);
        ysave = y[j];
        dy = dyVector[j];
        //dy = fmaxx(1.0E-6 * m_ewt[j], fabs(ysave)*1.0E-7);

        y[j] = ysave + dy;
        dy = y[j] - ysave;
	if (solnType_ != NSOLN_TYPE_STEADY_STATE) {
	  ydotsave = ydot[j];
	  ydot[j] += dy * CJ;
	}
        /*
         * Call the function
         */


        info =  m_func->evalResidNJ(time_curr, delta_t_n, y, ydot, DATA_PTR(m_wksp),
				    JacDelta_ResidEval, j, dy);
        m_nfe++;
	if (info != 1) {
	  mdp::mdp_safe_free((void **) &dyVector);
	  return info;
	}

        doublereal diff;
        for (i = 0; i < neq_; i++) {
          diff = subtractRD(m_wksp[i], f[i]);
          col_j[i] = diff / dy;
        }
	y[j] = ysave;
	if (solnType_ != NSOLN_TYPE_STEADY_STATE) {
	  ydot[j] = ydotsave;
	}

      }
      /*
       * Release memory
       */
      mdp::mdp_safe_free((void **) &dyVector);
    }

    if (m_print_flag >= 7 && s_print_NumJac) {
      if (neq_ < 30) {
	printf("\t\tCurrent Matrix and Residual:\n");
	printf("\t\t    I,J | ");
	for (j = 0; j < neq_; j++) {
	  printf("  %5d     ", j);
	}
	printf("|   Residual  \n");
	printf("\t\t        --");
	for (j = 0; j < neq_; j++) {
	  printf("------------");
	}
	printf("|  -----------\n");


	for (i = 0; i < neq_; i++) {
	  printf("\t\t   %4d |", i);
	  for (j = 0; j < neq_; j++) {
	    printf(" % 11.4E", J(i,j) );
	  }
	  printf(" |  % 11.4E\n", f[i]);
	}

	printf("\t\t        --");
	for (j = 0; j < neq_; j++) {
	  printf("------------");
	}
	printf("--------------\n");
      }
    }
    /*
     *  Make a copy of the data. Note, this jacobian copy occurs before any matrix scaling operations.
     *  It's the raw matrix producted by this routine.
     */
    jacCopy_.copyData(J);

    return retn;
  }
  //====================================================================================================================
  //   Internal function to calculate the time derivative at the new step
  /*
   *   @param  order of the BDF method
   *   @param   y_curr current value of the solution
   *   @param   ydot_curr  Calculated value of the solution derivative that is consistent with y_curr
   */ 
  void NonlinearSolver::
  calc_ydot(const int order, const doublereal * const y_curr, doublereal * const ydot_curr)
  {
    if (!ydot_curr) {
      return;
    }
    int    i;
    doublereal c1;
    switch (order) {
    case 0:
    case 1:             /* First order forward Euler/backward Euler */
      c1 = 1.0 / delta_t_n;
      for (i = 0; i < neq_; i++) {
        ydot_curr[i] = c1 * (y_curr[i] - m_y_nm1[i]);
      }
      return;
    case 2:             /* Second order Adams-Bashforth / Trapezoidal Rule */
      c1 = 2.0 / delta_t_n;
      for (i = 0; i < neq_; i++) {
        ydot_curr[i] = c1 * (y_curr[i] - m_y_nm1[i])  - m_ydot_nm1[i];
      }
      throw CanteraError("", "not implemented");
      return;
    }
  } 
  //====================================================================================================================
  // Apply a filtering process to the new step
  /*
   *  @param timeCurrent   Current value of the time
   *  @param y_current     current value of the solution
   *  @param ydot_current   Current value of the solution derivative. 
   *
   *  @return Returns the norm of the value of the amount filtered
   */
  doublereal NonlinearSolver::filterNewStep(const doublereal timeCurrent,
					    const doublereal * const ybase, doublereal * const step0) {
    doublereal tmp = m_func->filterNewStep(timeCurrent, ybase, step0);
    return tmp;
  }
  //====================================================================================================================
  // Apply a filtering process to the new solution
  /*
   *  @param timeCurrent   Current value of the time
   *  @param y_current     current value of the solution
   *  @param ydot_current   Current value of the solution derivative. 
   *
   *  @return Returns the norm of the value of the amount filtered
   */
  doublereal NonlinearSolver::filterNewSolution(const doublereal timeCurrent,
						doublereal * const y_current, doublereal *const ydot_current) {
   doublereal tmp = m_func->filterSolnPrediction(timeCurrent, y_current);
    return tmp;
  }
  //====================================================================================================================
  // Compute the Residual Weights
  /*
   *  The residual weights are defined here to be equal to the inverse of the row scaling factors used to
   *  row scale the matrix, after column scaling is used. They are multiplied by rtol and an atol factor
   *  is added as well so that if the residual is less than 1, then the calculation is deemed to be comverged.
   *
   *  The basic idea is that a change in the solution vector on the order of the convergence tolerance
   *  multiplied by  [RJC] which is of order one after row scaling should give you the relative weight
   *  of the row. Values of the residual for that row can then be normalized by the value of this weight.
   *  When the tolerance in delta x is achieved, the tolerance in the residual should also be achieved
   *  and should be checked
   */
  void
  NonlinearSolver::computeResidWts()
  {
    doublereal sum = 0.0;  
    for (int i = 0; i < neq_; i++) {
      m_residWts[i] = m_rowWtScales[i] / neq_;
      sum += m_residWts[i];
    }
    sum /= neq_;
    for (int i = 0; i < neq_; i++) {
      m_residWts[i] = m_ScaleSolnNormToResNorm * (m_residWts[i] + atolBase_ * atolBase_ * sum);
    }
  }
  //=====================================================================================================================
  // return the residual weights
  /*
   *  @param residWts  Vector of length neq_
   */
  void
  NonlinearSolver::getResidWts(doublereal * const residWts) const
  {
    for (int i = 0; i < neq_; i++) {
      residWts[i] = (m_residWts)[i];
    }
  }
  //=====================================================================================================================
  // Check to see if the nonlinear problem has converged
  /*
   *
   * @return integer is returned. If positive, then the problem has converged
   *           1 Successful step was taken: Next step's norm is less than 1.0.
   *                                        The final residual norm is less than 1.0.
   *           2 Successful step: Next step's norm is less than 0.8.
   *                              This step's norm is less than 1.0.
   *                              The residual norm can be anything.
   *           3 Success:  The final residual is less than 1.0E-2
   *                        The predicted deltaSoln is below 1.0E-2.
   *           0 Not converged yet
   */
  int
  NonlinearSolver::convergenceCheck(int dampCode, doublereal s1)
  {
    int retn = 0;
    if (m_dampBound < 0.9999) {
      return retn;
    }
    if (m_dampRes < 0.9999) {
      return retn;
    }
    if (dampCode <= 0) {
      return retn;
    }
    if (dampCode == 3) {
      if (s1 < 1.0E-2) {
	if (m_normResidTrial < 1.0E-6) {
	  return 3;
	}
      }
      if (s1 < 0.8) {
	if (m_normDeltaSoln_Newton < 1.0) {
	  return 2;
	}
      }
    }
    if (dampCode == 4) {
      if (s1 < 1.0E-2) {
	if (m_normResidTrial < 1.0E-6) {
	  return 3;
	}
      }
    }

    if (s1 < 0.8) {
      if (m_normDeltaSoln_Newton < 1.0) {
	return 2;
      }
    }
    if (dampCode == 1 || dampCode == 2) {
      if (s1 < 1.0) {
	if (m_normResidTrial < 1.0) {
	  return 1;
	}
      }
    }
    return retn;
  }
  //=====================================================================================================================
  // Set the absolute tolerances for the solution variables
  /*
   *   Set the absolute tolerances used in the calculation
   *
   *  @param atol   Vector of length neq_ that contains the tolerances to be used for the solution variables
   */
  void NonlinearSolver::setAtol(const doublereal * const atol)
  {
    for (int i = 0; i < neq_; i++) {
      atolk_[i]= atol[i];
    }
  }

  //=====================================================================================================================
  // Set the relative tolerances for the solution variables
  /*
   *   Set the relative tolerances used in the calculation
   *
   *  @param rtol  single double
   */
  void NonlinearSolver::setRtol(const doublereal rtol)
  {
    rtol_ = rtol;
  }
  //=====================================================================================================================
  
  void NonlinearSolver::setPrintLvl( int printLvl) 
  {
    m_print_flag = printLvl;
  }
  //=====================================================================================================================
}

