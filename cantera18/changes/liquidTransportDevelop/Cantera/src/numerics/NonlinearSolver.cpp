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
#include "GeneralMatrix.h"
#include "NonlinearSolver.h"
#include "ctlapack.h"

#include "clockWC.h"
#include "vec_functions.h"
#include "mdp_allo.h"

#include <cfloat>
#include <ctime>
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
#define CONSTD_DATA_PTR(x) (( const doublereal *) (&x[0]))
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

  bool NonlinearSolver::s_TurnOffTiming(false);

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

  // This toggle turns off the use of the Hessian when it is warranted by the condition number.
  /*
   *   This is a debugging option.
   */
  bool NonlinearSolver::s_alwaysAssumeNewtonGood(false);

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
    m_y_n_curr(0),
    m_ydot_n_curr(0),
    m_y_nm1(0),
    m_y_n_1(0),
    m_ydot_n_1(0),
    m_colScales(0),
    m_rowScales(0),
    m_rowWtScales(0),
    m_resid(0),
    m_wksp(0),
    m_wksp_2(0),
    m_residWts(0),
    m_normResid_0(0.0),
    m_normResid_Bound(0.0),
    m_normResid_1(0.0),
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
    userResidAtol_(0),
    userResidRtol_(1.0E-3),
    checkUserResidualTols_(0),
    m_print_flag(0),
    m_ScaleSolnNormToResNorm(0.001),
    jacCopyPtr_(0),
    HessianPtr_(0),
    deltaX_CP_(0),
    deltaX_Newton_(0),
    residNorm2Cauchy_(0.0),
    dogLegID_(0),
    dogLegAlpha_(1.0),
    RJd_norm_(0.0),
    lambdaStar_(0.0),
    Jd_(0),
    deltaX_trust_(0),
    norm_deltaX_trust_(0.0),
    trustDelta_(1.0),
    trustRegionInitializationMethod_(2),
    trustRegionInitializationFactor_(1.0),
    Nuu_(0.0),
    dist_R0_(0.0),
    dist_R1_(0.0),
    dist_R2_(0.0),
    dist_Total_(0.0),
    JdJd_norm_(0.0),
    normTrust_Newton_(0.0),
    normTrust_CP_(0.0),
    doDogLeg_(0),
    doAffineSolve_(0) ,
    CurrentTrustFactor_(1.0),
    NextTrustFactor_ (1.0),
    ResidWtsReevaluated_(false),
    ResidDecreaseSDExp_(0.0),
    ResidDecreaseSD_(0.0),
    ResidDecreaseNewtExp_(0.0),
    ResidDecreaseNewt_(0.0)
  {
    neq_ = m_func->nEquations();

    m_ewt.resize(neq_, rtol_);
    m_deltaStepMinimum.resize(neq_, 0.001);
    m_deltaStepMaximum.resize(neq_, 1.0E10);
    m_y_n_curr.resize(neq_, 0.0);
    m_ydot_n_curr.resize(neq_, 0.0);
    m_y_nm1.resize(neq_, 0.0);
    m_y_n_1.resize(neq_, 0.0);
    m_ydot_n_1.resize(neq_, 0.0);
    m_colScales.resize(neq_, 1.0);
    m_rowScales.resize(neq_, 1.0);
    m_rowWtScales.resize(neq_, 1.0);
    m_resid.resize(neq_, 0.0);
    m_wksp.resize(neq_, 0.0);
    m_wksp_2.resize(neq_, 0.0);
    m_residWts.resize(neq_, 0.0);
    atolk_.resize(neq_, atolBase_);
    deltaX_Newton_.resize(neq_, 0.0);
    m_step_1.resize(neq_, 0.0);
    m_y_n_1.resize(neq_, 0.0);
    doublereal hb = std::numeric_limits<double>::max();
    m_y_high_bounds.resize(neq_, hb);
    m_y_low_bounds.resize(neq_, -hb);    

    for (int i = 0; i < neq_; i++) {
      atolk_[i] = atolBase_;
      m_ewt[i] = atolk_[i];
    }


    // jacCopyPtr_->resize(neq_, 0.0);
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
    m_y_n_curr(0),
    m_ydot_n_curr(0),
    m_y_nm1(0),
    m_y_n_1(0),
    m_ydot_n_1(0),
    m_step_1(0),
    m_colScales(0),
    m_rowScales(0),
    m_rowWtScales(0),
    m_resid(0),
    m_wksp(0),
    m_wksp_2(0),
    m_residWts(0),
    m_normResid_0(0.0),
    m_normResid_Bound(0.0),
    m_normResid_1(0.0),
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
    userResidAtol_(0),
    userResidRtol_(1.0E-3),
    checkUserResidualTols_(0),
    m_print_flag(0),
    m_ScaleSolnNormToResNorm(0.001),
    jacCopyPtr_(0),
    HessianPtr_(0),
    deltaX_CP_(0),
    deltaX_Newton_(0),
    residNorm2Cauchy_(0.0),
    dogLegID_(0),
    dogLegAlpha_(1.0),
    RJd_norm_(0.0),
    lambdaStar_(0.0),
    Jd_(0),
    deltaX_trust_(0),
    norm_deltaX_trust_(0.0),
    trustDelta_(1.0),
    trustRegionInitializationMethod_(2),
    trustRegionInitializationFactor_(1.0),
    Nuu_(0.0),
    dist_R0_(0.0),
    dist_R1_(0.0),
    dist_R2_(0.0),
    dist_Total_(0.0),
    JdJd_norm_(0.0),
    normTrust_Newton_(0.0),
    normTrust_CP_(0.0),
    doDogLeg_(0),
    doAffineSolve_(0),
    CurrentTrustFactor_(1.0),
    NextTrustFactor_ (1.0),
    ResidWtsReevaluated_(false),
    ResidDecreaseSDExp_(0.0),
    ResidDecreaseSD_(0.0),
    ResidDecreaseNewtExp_(0.0),
    ResidDecreaseNewt_(0.0)
  {
    *this =operator=(right);
  }

  //====================================================================================================================
  NonlinearSolver::~NonlinearSolver() {
    if (jacCopyPtr_) {
      delete jacCopyPtr_;
    }
    if (HessianPtr_) {
      delete HessianPtr_;
    }
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
    m_y_n_curr                 = right.m_y_n_curr;
    m_ydot_n_curr              = right.m_ydot_n_curr;
    m_y_nm1                    = right.m_y_nm1;
    m_y_n_1                    = right.m_y_n_1;
    m_ydot_n_1                 = right.m_ydot_n_1;
    m_step_1                   = right.m_step_1;
    m_colScales                = right.m_colScales;
    m_rowScales                = right.m_rowScales;
    m_rowWtScales              = right.m_rowWtScales;
    m_resid                    = right.m_resid;
    m_wksp                     = right.m_wksp;
    m_wksp_2                   = right.m_wksp_2;
    m_residWts                 = right.m_residWts;
    m_normResid_0              = right.m_normResid_0;
    m_normResid_Bound          = right.m_normResid_Bound;
    m_normResid_1              = right.m_normResid_1;
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
    userResidAtol_             = right.userResidAtol_;
    userResidRtol_             = right.userResidRtol_;
    checkUserResidualTols_     = right.checkUserResidualTols_;
    m_print_flag               = right.m_print_flag;
    m_ScaleSolnNormToResNorm   = right.m_ScaleSolnNormToResNorm;

    if (jacCopyPtr_) {
      delete (jacCopyPtr_);
    }
    jacCopyPtr_                = (right.jacCopyPtr_)->duplMyselfAsGeneralMatrix();
    if (HessianPtr_) {
      delete (HessianPtr_);
    }
    HessianPtr_                = (right.HessianPtr_)->duplMyselfAsGeneralMatrix();

    deltaX_CP_                 = right.deltaX_CP_;
    deltaX_Newton_             = right.deltaX_Newton_;
    residNorm2Cauchy_          = right.residNorm2Cauchy_;
    dogLegID_                  = right.dogLegID_;
    dogLegAlpha_               = right.dogLegAlpha_;
    RJd_norm_                  = right.RJd_norm_;
    lambdaStar_                = right.lambdaStar_;
    Jd_                        = right.Jd_;
    deltaX_trust_              = right.deltaX_trust_;	
    norm_deltaX_trust_         = right.norm_deltaX_trust_;
    trustDelta_                = right.trustDelta_;
    trustRegionInitializationMethod_ = right.trustRegionInitializationMethod_;
    trustRegionInitializationFactor_ = right.trustRegionInitializationFactor_;
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
    CurrentTrustFactor_        = right.CurrentTrustFactor_;
    NextTrustFactor_           = right.NextTrustFactor_;

    ResidWtsReevaluated_       = right.ResidWtsReevaluated_;
    ResidDecreaseSDExp_        = right.ResidDecreaseSDExp_;
    ResidDecreaseSD_           = right.ResidDecreaseSD_;
    ResidDecreaseNewtExp_      = right.ResidDecreaseNewtExp_;
    ResidDecreaseNewt_         = right.ResidDecreaseNewt_;

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
#ifdef DEBUG_MODE
      if (m_ewt[i] <= 0.0) {
	throw CanteraError(" NonlinearSolver::createSolnWeights()", "ewts <= 0.0");
      }
#endif
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
  }
  //====================================================================================================================
  std::vector<doublereal> &  NonlinearSolver::lowBoundsConstraintVector() {
    return  m_y_low_bounds;
  }
  //====================================================================================================================
  std::vector<doublereal> &  NonlinearSolver::highBoundsConstraintVector() {
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
      if ((printLargest == 1) || (m_print_flag >= 4 && m_print_flag <= 5)) {

	printf("\t\t   solnErrorNorm(): ");
	if (title) {
	  printf("%s", title);
	} else {
	  printf(" Delta soln norm ");
	}
	printf(" = %-11.4E\n", sum_norm);
      } else if (m_print_flag >= 6) {



	const int num_entries = printLargest;
	printf("\t\t   "); print_line("-", 90);
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
        printf("\t\t        Printout of Largest Contributors:                     (damp = %g)\n", dampFactor);
        printf("\t\t        I   weightdeltaY/sqtN|     deltaY    "
	       "ysolnOld     ysolnNew   Soln_Weights\n");
        printf("\t\t     "); print_line("-", 88);

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
	    printf("\t\t     %4d %12.4e       | %12.4e %12.4e %12.4e %12.4e\n", i, normContrib/sqrt((double)neq_), 
		   delta_y[i], m_y_n_curr[i], m_y_n_curr[i] + dampFactor * delta_y[i], m_ewt[i]);

	  }	  
	}
	printf("\t\t   "); print_line("-", 90);
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
					     const doublereal * const y) const
  {
    int    i;
    doublereal sum_norm = 0.0, error;

    for (i = 0; i < neq_; i++) {
#ifdef DEBUG_MODE
      mdp::checkFinite(resid[i]);
#endif
      error     = resid[i] / m_residWts[i];
#ifdef DEBUG_MODE
      mdp::checkFinite(error);
#endif
      sum_norm += (error * error);
    }
    sum_norm = sqrt(sum_norm / neq_); 
#ifdef DEBUG_MODE
    mdp::checkFinite(sum_norm);
#endif
    if (printLargest) {
      const int num_entries = printLargest;
      doublereal dmax1, normContrib;
      int j;
      int *imax = mdp::mdp_alloc_int_1(num_entries, -1);

      if (m_print_flag >= 4 && m_print_flag <= 5) {
	printf("\t\t   residErrorNorm():");
	if (title) {
	  printf(" %s ", title);
	} else {
	  printf("  residual L2 norm ");
	}
	printf("= %12.4E\n", sum_norm);
      }
      if (m_print_flag >= 6) {
	printf("\t\t   "); print_line("-", 90);
        printf("\t\t   residErrorNorm(): ");
	if (title) {
	  printf(" %s ", title);
	} else {
	  printf("  residual L2 norm ");
	}
	printf("= %12.4E\n", sum_norm);
        printf("\t\t        Printout of Largest Contributors to norm:\n");
        printf("\t\t        I       |Resid/ResWt|     UnsclRes          ResWt    |  y_curr\n");
        printf("\t\t     ");
        print_line("-", 88);
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
	    printf("\t\t     %4d     %12.4e     %12.4e     %12.4e | %12.4e\n", i,  normContrib, resid[i], m_residWts[i], y[i]);
	  }	  
	}
      
	printf("\t\t   "); 
	print_line("-", 90);
      }
      mdp::mdp_safe_free((void **) &imax);
    }
    return sum_norm;
  }
  //====================================================================================================================
  // Set the column scaling that are used for the inversion of the matrix
  /*
   *  There are three ways to do this.
   *
   *  The first method is to set the bool useColScaling to true, leaving the scaling factors unset.
   *  Then, the column scales will be set to the solution error weighting factors. This has the
   *  effect of ensuring that all delta variables will have the same order of magnitude at convergence
   *  end.
   *  
   *  The second way is the explicity set the column factors in the second argument of this function call.
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
  void NonlinearSolver::setColumnScaling(bool useColScaling, const double * const scaleFactors) {
    if (useColScaling) {
      if (scaleFactors) {
        m_colScaling = 2;
        for (int i = 0; i < neq_; i++) {
          m_colScales[i] = scaleFactors[i];
          if (m_colScales[i] <= 1.0E-200) {
            throw CanteraError("NonlinearSolver::setColumnScaling() ERROR", "Bad column scale factor");
          }
        }
      } else {
        m_colScaling = 1; 
      }
    } else {
      m_colScaling = 0;
    }
  }
  //====================================================================================================================
  // Set the rowscaling that are used for the inversion of the matrix
  /*
   * Row scaling is set here. Right now the row scaling is set internally in the code.
   *
   * @param useRowScaling   Turn row scaling on or off.
   */
  void NonlinearSolver::setRowScaling(bool useRowScaling) {
     m_rowScaling = useRowScaling;
  }
  //====================================================================================================================
  /*
   * calcColumnScales():
   *
   * Set the column scaling vector at the current time
   */
  void NonlinearSolver::calcColumnScales() {
    if (m_colScaling == 1) {
      for (int i = 0; i < neq_; i++) {
	m_colScales[i] = m_ewt[i];
      }
    } else {
      for (int i = 0; i < neq_; i++) {
	m_colScales[i] = 1.0;
      }
    }
    if (m_colScaling) {
      m_func->calcSolnScales(time_n, DATA_PTR(m_y_n_curr), DATA_PTR(m_y_nm1), DATA_PTR(m_colScales));
    }
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
				      const doublereal * const ydot_curr, const ResidEval_Type_Enum evalType) const
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
  void NonlinearSolver::scaleMatrix(GeneralMatrix& jac, doublereal * const y_comm, doublereal * const ydot_comm,
				    doublereal time_curr, int num_newt_its)
  {	 
    int irow, jcol;
    int ku, kl;
    int ivec[2];
    jac.nRowsAndStruct(ivec);
    double *colP_j;

    /*
     * Column scaling -> We scale the columns of the Jacobian
     * by the nominal important change in the solution vector
     */
    if (m_colScaling) {
      if (!jac.factored()) {
	if (jac.matrixType_ == 0) {
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
	} else if (jac.matrixType_ == 1) {
	  kl = ivec[0];
	  ku = ivec[1];
	  for (jcol = 0; jcol < neq_; jcol++) {
	    colP_j = (doublereal *) jac.ptrColumn(jcol);
	    for (irow = jcol - ku; irow <= jcol + kl; irow++) {
	      if (irow >= 0 && irow < neq_) {
		colP_j[kl + ku + irow - jcol] *= m_colScales[jcol];
	      }
	    }
	  }
	}
      }
    }
    /*
     * row sum scaling -> Note, this is an unequivical success
     *      at keeping the small numbers well balanced and nonnegative.
     */
    if (! jac.factored()) {
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
      if (jac.matrixType_ == 0) {
	for (jcol = 0; jcol < neq_; jcol++) {
	  for (irow = 0; irow < neq_; irow++) {
	    if (m_rowScaling) {
	      m_rowScales[irow] += fabs(*jptr);
	    }
	    if (m_colScaling) {
	      // This is needed in order to mitgate the change in J_ij carried out just above this loop.
	      // Alternatively, we could move this loop up to the top
	      m_rowWtScales[irow] += fabs(*jptr) * m_ewt[jcol] / m_colScales[jcol];
	    } else {
	      m_rowWtScales[irow] += fabs(*jptr) * m_ewt[jcol];
	    }
#ifdef DEBUG_MODE
            mdp::checkFinite(m_rowWtScales[irow]);
#endif
	    jptr++;
	  }
	}	
      } else if (jac.matrixType_ == 1) {
	kl = ivec[0];
	ku = ivec[1];
	for (jcol = 0; jcol < neq_; jcol++) {
	  colP_j = (doublereal *) jac.ptrColumn(jcol);
	  for (irow = jcol - ku; irow <= jcol + kl; irow++) {
	    if (irow >= 0 && irow < neq_) {
	      double vv = fabs(colP_j[kl + ku + irow - jcol]);
	      if (m_rowScaling) {
		m_rowScales[irow] += vv;
	      }
	      if (m_colScaling) {
		// This is needed in order to mitgate the change in J_ij carried out just above this loop.
		// Alternatively, we could move this loop up to the top
		m_rowWtScales[irow] += vv * m_ewt[jcol] / m_colScales[jcol];
	      } else {
		m_rowWtScales[irow] += vv * m_ewt[jcol];
	      }	      
#ifdef DEBUG_MODE
              mdp::checkFinite(m_rowWtScales[irow]);
#endif
	    }
	  }
	}
      }
      if (m_rowScaling) {
	for (irow = 0; irow < neq_; irow++) {
	  m_rowScales[irow] = 1.0/m_rowScales[irow];
	}
      } else {
	for (irow = 0; irow < neq_; irow++) {
	  m_rowScales[irow] = 1.0;
	}
      }
      // What we have defined is a maximum value that the residual can be and still pass.
      // This isn't sufficient.
      
      if (m_rowScaling) {
	if (jac.matrixType_ == 0) {
	  jptr = &(*(jac.begin()));
	  for (jcol = 0; jcol < neq_; jcol++) {
	    for (irow = 0; irow < neq_; irow++) {
	      *jptr *= m_rowScales[irow];
	      jptr++;
	    }
	  }
	} else if (jac.matrixType_ == 1) {
	  kl = ivec[0];
	  ku = ivec[1];
	  for (jcol = 0; jcol < neq_; jcol++) {
	    colP_j = (doublereal *) jac.ptrColumn(jcol);
	    for (irow = jcol - ku; irow <= jcol + kl; irow++) {
	      if (irow >= 0 && irow < neq_) {
		colP_j[kl + ku + irow - jcol] *= m_rowScales[irow];
	      }
	    }
	  }
	}
      }
   
      if (num_newt_its % 5 == 1) {    
	computeResidWts();
      }
  
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
    if (! jacCopyPtr_->factored()) { 
   
      if (checkUserResidualTols_ != 1) {
	doublereal sum = 0.0;
	for (int irow = 0; irow < neq_; irow++) {
	  m_residWts[irow] = m_rowWtScales[irow] / neq_;
	  sum += m_residWts[irow];
	}
	sum /= neq_;
	for (int irow = 0; irow < neq_; irow++) {
	  m_residWts[irow] = (m_residWts[irow] + atolBase_ * atolBase_ * sum);
	}
	if (checkUserResidualTols_ == 2) {
	  for (int irow = 0; irow < neq_; irow++) {
	    m_residWts[irow] = MIN(m_residWts[irow], userResidAtol_[irow] + userResidRtol_ * m_rowWtScales[irow] / neq_);
	  }
	}
      } else {
	for (int irow = 0; irow < neq_; irow++) {
	  m_residWts[irow] = userResidAtol_[irow] + userResidRtol_ * m_rowWtScales[irow] / neq_;
	}
      }
      
      
      for (int irow = 0; irow < neq_; irow++) {
        m_wksp[irow] = 0.0;
      }   
      doublereal *jptr = &(jacCopyPtr_->operator()(0,0));
      for (int jcol = 0; jcol < neq_; jcol++) {
	for (int irow = 0; irow < neq_; irow++) {
	  m_wksp[irow] += (*jptr) * m_ewt[jcol];
	  jptr++;
	}
      }
      doublereal resNormOld = 0.0;
      doublereal error;
     
      for (int irow = 0; irow < neq_; irow++) {
	error = m_wksp[irow] /  m_residWts[irow];
        resNormOld  += error * error;
      }
      resNormOld = sqrt(resNormOld / neq_);

      if (resNormOld > 0.0) {
	m_ScaleSolnNormToResNorm = resNormOld;
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
   * solution vector, m_y_n_curr, and the solution time derivative, m_ydot_n. 
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
				     GeneralMatrix& jac)
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
    m_numLocalLinearSolves++;
    return info;
  }
  //====================================================================================================================
  // Compute the newton step, either by direct newton's or by solving a close problem that is represented
  // by a Hessian (
  /*
   * This is algorith A.6.5.1 in Dennis / Schnabel
   *
   * Compute the QR decomposition
   *
   *  Notes on banded Hessian solve:
   *   The matrix for jT j  has a larger band width. Both the top and bottom band widths
   *   are doubled, going from KU to KU+KL  and KL to KU+KL in size. This is not an impossible increase in cost, but
   *   has to be considered.
   */
  int NonlinearSolver::doAffineNewtonSolve(const doublereal * const y_curr,   const doublereal * const ydot_curr, 
					   doublereal * const delta_y, GeneralMatrix& jac)
  {
    bool newtonGood = true;
    int irow;
    doublereal *delyNewton = 0;
    // We can default to QR here ( or not )
    jac.useFactorAlgorithm(1);
    int useQR = jac.factorAlgorithm();
    // multiplyl the residual by -1
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

    // Factor the matrix using a standard Newton solve
    m_conditionNumber = 1.0E300;
    int info = 0;
    if (!jac.factored()) { 
      if (useQR) {
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
      doublereal rcond = 0.0;
      if (useQR) {
	rcond = jac.rcondQR();
      } else {
	doublereal a1norm = jac.oneNorm();
	rcond = jac.rcond(a1norm);
      }
      if (rcond > 0.0) {
	m_conditionNumber = 1.0 / rcond;
      }
    } else {
      m_conditionNumber = 1.0E300;
      newtonGood = false;
      if (m_print_flag >= 1) {
	printf("\t\t   doAffineNewtonSolve: ");
	if (useQR) {
	  printf("factorQR()");
	} else {
	  printf("factor()");
	}
	printf(" returned with info = %d, indicating a zero row or column\n", info);
      }
    }
    bool doHessian = false;
    if (s_doBothSolvesAndCompare) {
      doHessian = true;
    }
    bool useNewton = false;
    if (m_conditionNumber < 1.0E7) {
      useNewton = true;
      if (m_print_flag >= 4) {
	printf("\t\t   doAffineNewtonSolve: Condition number = %g during regular solve\n", m_conditionNumber);
      }

      /*
       * Solve the system -> This also involves inverting the matrix
       */
      int info = jac.solve(DATA_PTR(delta_y));
      if (info) {
	if (m_print_flag >= 2) {
	  printf("\t\t   doAffineNewtonSolve() ERROR: QRSolve returned INFO = %d. Switching to Hessian solve\n", info);
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
      if (jac.matrixType_ == 1) {
	useNewton = true;
	newtonGood = true;
	if (m_print_flag >= 3) {
	  printf("\t\t   doAffineNewtonSolve() WARNING: Condition number too large, %g, But Banded Hessian solve "
		 "not implemented yet \n", m_conditionNumber);
	}
      } else {
	doHessian = true;
	newtonGood = false;
	if (m_print_flag >= 3) {
	  printf("\t\t   doAffineNewtonSolve() WARNING: Condition number too large, %g. Doing a Hessian solve \n", m_conditionNumber);
	}
      }
    }

    if (doHessian) {
      // Store the old value for later comparison
 
      delyNewton = mdp::mdp_alloc_dbl_1(neq_, MDP_DBL_NOINIT);
      for (irow = 0; irow < neq_; irow++) {
	delyNewton[irow] = delta_y[irow];
      }
      
      // Get memory if not done before
      if (HessianPtr_ == 0) {
	HessianPtr_ = jac.duplMyselfAsGeneralMatrix();
      }

      /*
       * Calculate the symmetric Hessian
       */
      GeneralMatrix &hessian = *HessianPtr_;
      GeneralMatrix &jacCopy = *jacCopyPtr_;
      hessian.zero();
      if (m_rowScaling) {
	for (int i = 0; i < neq_; i++) {
	  for (int j = i; j < neq_; j++) {
	    for (int k = 0; k < neq_; k++) {
	      hessian(i,j) += jacCopy(k,i) * jacCopy(k,j) * m_rowScales[k] * m_rowScales[k];
	    }
	    hessian(j,i) = hessian(i,j);
	  }
	}
      } else {
	for (int i = 0; i < neq_; i++) {
	  for (int j = i; j < neq_; j++) {
	    for (int k = 0; k < neq_; k++) {
	      hessian(i,j) += jacCopy(k,i) * jacCopy(k,j);
	    }
	    hessian(j,i) = hessian(i,j);
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
	    hcol += fabs(hessian(j,i)) * m_colScales[j];
	  }
	  for (int j = i+1; j < neq_; j++) {
	    hcol += fabs(hessian(i,j)) * m_colScales[j];
	  }
	  hcol *= m_colScales[i];
	  if (hcol > hnorm) {
	    hnorm = hcol;
	  }
	}
      } else {
	for (int i = 0; i < neq_; i++) {
	  for (int j = i; j < neq_; j++) {
	    hcol += fabs(hessian(j,i));
	  }
	  for (int j = i+1; j < neq_; j++) {
	    hcol += fabs(hessian(i,j));
	  }
	  if (hcol > hnorm) {
	    hnorm = hcol;
	  }
	}
      }
      /*
       * Add junk to the Hessian diagonal
       *  -> Note, testing indicates that this will get too big for ill-conditioned systems.
       */
      hcol = sqrt(neq_) * 1.0E-7 * hnorm;
#ifdef DEBUG_HKM_NOT
      if (hcol > 1.0) {
	hcol = 1.0E1;
      }
#endif
      if (m_colScaling) {
	for (int i = 0; i < neq_; i++) {
	  hessian(i,i) += hcol / (m_colScales[i] * m_colScales[i]);
	}
      } else {
	for (int i = 0; i < neq_; i++) {
	  hessian(i,i) += hcol;
	}
      }

      /*
       *  Factor the Hessian
       */
      int info;
      ct_dpotrf(ctlapack::UpperTriangular, neq_, &(*(HessianPtr_->begin())), neq_, info);
      if (info) {
	if (m_print_flag >= 2) {
	  printf("\t\t    doAffineNewtonSolve() ERROR: Hessian isn't positive definate DPOTRF returned INFO = %d\n", info);
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
	    delta_y[j] += delyH[i] * jacCopy(i,j) * m_rowScales[i];
	  }
	}
      } else {
	for (int j = 0; j < neq_; j++) {
	  delta_y[j] = 0.0;
	  for (int i = 0; i < neq_; i++) {
	    delta_y[j] += delyH[i] * jacCopy(i,j);
	  }
	}
      }

    
      /*
       * Solve the factored Hessian System
       */
      ct_dpotrs(ctlapack::UpperTriangular, neq_, 1,&(*(hessian.begin())), neq_, delta_y, neq_, info);
      if (info) {
	if (m_print_flag >= 2) {
	  printf("\t\t   NonlinearSolver::doAffineNewtonSolve() ERROR: DPOTRS returned INFO = %d\n", info);
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


      if (doDogLeg_ && m_print_flag > 7) {
        double normNewt = solnErrorNorm(CONSTD_DATA_PTR(delyNewton));
	double normHess = solnErrorNorm(CONSTD_DATA_PTR(delta_y));
	printf("\t\t          doAffineNewtonSolve(): Printout Comparison between Hessian deltaX and Newton deltaX\n");
	
	printf("\t\t               I    Hessian+Junk     Newton");
	if (newtonGood || s_alwaysAssumeNewtonGood) {
	  printf(" (USING NEWTON DIRECTION)\n");
	} else {
	  printf(" (USING HESSIAN DIRECTION)\n");
	}
	printf("\t\t            Norm: %12.4E %12.4E\n", normHess, normNewt);
       
	printf("\t\t          --------------------------------------------------------\n");
	for (int i =0; i < neq_; i++) {
	  printf("\t\t             %3d  %13.5E %13.5E\n", i, delta_y[i], delyNewton[i]);
	}
	printf("\t\t          --------------------------------------------------------\n");
      } else if (doDogLeg_ && m_print_flag >= 4) {
	double normNewt = solnErrorNorm(CONSTD_DATA_PTR(delyNewton));
	double normHess = solnErrorNorm(CONSTD_DATA_PTR(delta_y));
	printf("\t\t          doAffineNewtonSolve():  Hessian update norm = %12.4E \n"
	       "\t\t                                  Newton  update norm = %12.4E \n", normHess, normNewt);
	if (newtonGood || s_alwaysAssumeNewtonGood) {
	  printf("\t\t                                 (USING NEWTON DIRECTION)\n");
	} else {
	  printf("\t\t                                 (USING HESSIAN DIRECTION)\n");
	}
      }

      /*
       *  Choose the delta_y to use
       */
      if (newtonGood || s_alwaysAssumeNewtonGood) {
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
    m_numLocalLinearSolves++;
    return info;

  }
  //====================================================================================================================
  // Do a steepest descent calculation
  /*
   *  This call must be made on the unfactored jacobian!
   */
  doublereal NonlinearSolver::doCauchyPointSolve(GeneralMatrix& jac)
  {
    doublereal rowFac = 1.0;
    doublereal colFac = 1.0;
    doublereal normSoln;
    //  Calculate the descent direction
    /*
     * For confirmation of the scaling factors, see Dennis and Schnabel p, 152, p, 156 and my notes
     * 
     *  The colFac and rowFac values are used to eliminate the scaling of the matrix from the 
     *  actual equation
     *
     *  Here we calculate the steepest descent direction. This is equation (11) in the notes. It is
     *  storred in deltaX_CP_[].The value corresponds to d_descent[].
     */
    for (int j = 0; j < neq_; j++) {
      deltaX_CP_[j] = 0.0;
      if (m_colScaling) {
	colFac = 1.0 / m_colScales[j];
      }
      for (int i = 0; i < neq_; i++) {
	if (m_rowScaling) {
	  rowFac = 1.0 / m_rowScales[i];
	}
        deltaX_CP_[j] -= m_resid[i] * jac(i,j) * colFac * rowFac * m_ewt[j] * m_ewt[j] 
	  / (m_residWts[i] * m_residWts[i]);
#ifdef DEBUG_MODE
	mdp::checkFinite(deltaX_CP_[j]);
#endif
      }
    }

    /*
     *  Calculate J_hat d_y_descent. This is formula 18 in the notes.
     */
    for (int i = 0; i < neq_; i++) {
      Jd_[i] = 0.0;
      if (m_rowScaling) {
	rowFac = 1.0 / m_rowScales[i];
      } else {
	rowFac = 1.0;
      }
      for (int j = 0; j < neq_; j++) {
        if (m_colScaling) {
          colFac = 1.0 / m_colScales[j];
        }
	Jd_[i] += deltaX_CP_[j] * jac(i,j) * rowFac * colFac / m_residWts[i];
      }
    }

    /*
     * Calculate the distance along the steepest descent until the Cauchy point
     * This is Eqn. 17 in the notes.
     */
    RJd_norm_ = 0.0;
    JdJd_norm_ = 0.0;
    for (int i = 0; i < neq_; i++) {
      RJd_norm_ += m_resid[i] * Jd_[i] / m_residWts[i];
      JdJd_norm_ += Jd_[i] * Jd_[i];
    }
    //if (RJd_norm_ > -1.0E-300) {
    //  printf("we are here: zero residual\n");
    //}
    if (fabs(JdJd_norm_) < 1.0E-290) {
      if (fabs(RJd_norm_) < 1.0E-300) {
	lambdaStar_ = 0.0;
      } else {
	throw CanteraError("NonlinearSolver::doCauchyPointSolve()", "Unexpected condition: norms are zero");
      }
    } else {
      lambdaStar_ = - RJd_norm_ / (JdJd_norm_);
    }

    /*
     *  Now we modify the steepest descent vector such that its length is equal to the 
     *  Cauchy distance. From now on, if we want to recreate the descent vector, we have
     *  to unnormalize it by dividing by lambdaStar_.
     */
    for (int i = 0; i < neq_; i++) {
      deltaX_CP_[i] *= lambdaStar_;
    }


    doublereal normResid02 = m_normResid_0 * m_normResid_0 * neq_;

    /*
     *  Calculate the expected square of the risdual at the Cauchy point if the linear model is correct
     */
    if (fabs(JdJd_norm_) < 1.0E-290) {
      residNorm2Cauchy_ = normResid02;
    } else {
      residNorm2Cauchy_ = normResid02 - RJd_norm_ * RJd_norm_ / (JdJd_norm_);
    }


    // Extra printout section
    if (m_print_flag > 2) {
      // Calculate the expected residual at the Cauchy point if the linear model is correct
      doublereal residCauchy = 0.0;
      if (residNorm2Cauchy_ > 0.0) {
	residCauchy = sqrt(residNorm2Cauchy_ / neq_);
      } else {
	if (fabs(JdJd_norm_) < 1.0E-290) {
	  residCauchy =  m_normResid_0;
	} else {
	  residCauchy =  m_normResid_0 - sqrt(RJd_norm_ * RJd_norm_ / (JdJd_norm_));
	}
      }  
      // Compute the weighted norm of the undamped step size descentDir_[]
      if ((s_print_DogLeg || doDogLeg_) && m_print_flag >= 6) {
	normSoln = solnErrorNorm(DATA_PTR(deltaX_CP_), "SteepestDescentDir", 10);
      } else {
	normSoln = solnErrorNorm(DATA_PTR(deltaX_CP_), "SteepestDescentDir", 0);
      }
      if ((s_print_DogLeg || doDogLeg_) && m_print_flag >= 5) {
	printf("\t\t   doCauchyPointSolve: Steepest descent to Cauchy point: \n");
	printf("\t\t\t      R0     = %g \n", m_normResid_0);
	printf("\t\t\t      Rpred  = %g\n",  residCauchy);
	printf("\t\t\t      Rjd    = %g\n",  RJd_norm_);
	printf("\t\t\t      JdJd   = %g\n",  JdJd_norm_);
	printf("\t\t\t      deltaX = %g\n",  normSoln);
	printf("\t\t\t      lambda = %g\n",  lambdaStar_);
      }
    } else {
      // Calculate the norm of the Cauchy solution update in any case
      normSoln = solnErrorNorm(DATA_PTR(deltaX_CP_), "SteepestDescentDir", 0);
    }
    return normSoln;
  }
  //===================================================================================================================
  void  NonlinearSolver::descentComparison(doublereal time_curr, doublereal *ydot0, doublereal *ydot1, int &numTrials)
  {
    int info;
    doublereal ff = 1.0E-5;
    doublereal ffNewt = 1.0E-5;
    doublereal *y_n_1 = DATA_PTR(m_wksp);
    doublereal cauchyDistanceNorm = solnErrorNorm(DATA_PTR(deltaX_CP_));
    if (cauchyDistanceNorm < 1.0E-2) {
      ff = 1.0E-9 / cauchyDistanceNorm;
      if (ff > 1.0E-2) {
	ff = 1.0E-2;
      }
    }
    for (int i = 0; i < neq_; i++) {
      y_n_1[i] = m_y_n_curr[i] + ff * deltaX_CP_[i];
    }
    /*
     *  Calculate the residual that would result if y1[] were the new solution vector
     *  -> m_resid[] contains the result of the residual calculation
     */
    if (solnType_ != NSOLN_TYPE_STEADY_STATE) {
      info = doResidualCalc(time_curr, solnType_, y_n_1, ydot1, Base_LaggedSolutionComponents);
    } else {
      info = doResidualCalc(time_curr, solnType_, y_n_1, ydot0, Base_LaggedSolutionComponents);
    }  

    doublereal normResid02 = m_normResid_0 * m_normResid_0 * neq_;
    doublereal residSteep = residErrorNorm(DATA_PTR(m_resid));
    doublereal residSteep2 = residSteep * residSteep * neq_;
    doublereal funcDecreaseSD = 0.5 * (residSteep2 - normResid02) / ( ff * cauchyDistanceNorm);

    doublereal sNewt = solnErrorNorm(DATA_PTR(deltaX_Newton_));
    if (sNewt > 1.0) {
      ffNewt = ffNewt / sNewt;
    }
    for (int i = 0; i < neq_; i++) {
      y_n_1[i] = m_y_n_curr[i] + ffNewt * deltaX_Newton_[i];
    }
    /*
     *  Calculate the residual that would result if y1[] were the new solution vector.
     *  Here we use the lagged solution components in the residual calculation as well. We are
     *  interested in the linear model and its agreement with the nonlinear model.
     *
     *  -> m_resid[] contains the result of the residual calculation
     */
    if (solnType_ != NSOLN_TYPE_STEADY_STATE) {
      info = doResidualCalc(time_curr, solnType_, y_n_1, ydot1, Base_LaggedSolutionComponents);
    } else {
      info = doResidualCalc(time_curr, solnType_, y_n_1, ydot0, Base_LaggedSolutionComponents);
    }  
    doublereal residNewt = residErrorNorm(DATA_PTR(m_resid));
    doublereal residNewt2 = residNewt * residNewt * neq_;

    doublereal funcDecreaseNewt2 = 0.5 * (residNewt2 - normResid02) / ( ffNewt * sNewt);

    // This is the expected inital rate of decrease in the Cauchy direction.
    //   -> This is Eqn. 29 = Rhat dot Jhat dy / || d ||
    doublereal funcDecreaseSDExp = RJd_norm_ / cauchyDistanceNorm * lambdaStar_;

    doublereal funcDecreaseNewtExp2 = -  normResid02 / sNewt;

    if ( m_normResid_0 > 1.0E-100) {
      ResidDecreaseSDExp_   = funcDecreaseSDExp    / neq_ / m_normResid_0;
      ResidDecreaseSD_      = funcDecreaseSD       / neq_ / m_normResid_0;
      ResidDecreaseNewtExp_ = funcDecreaseNewtExp2 / neq_ / m_normResid_0;
      ResidDecreaseNewt_    = funcDecreaseNewt2    / neq_ / m_normResid_0;
    } else {
      ResidDecreaseSDExp_   = 0.0;
      ResidDecreaseSD_      = funcDecreaseSD       / neq_;
      ResidDecreaseNewtExp_ = 0.0;
      ResidDecreaseNewt_    = funcDecreaseNewt2    / neq_;
    }
    numTrials += 2;

    /*
     *   HKM These have been shown to exactly match up.
     *   The steepest direction is always largest even when there are variable solution weights
     *
     *   HKM When a hessian is used with junk on the diagonal,  funcDecreaseNewtExp2 is no longer accurate as the
     *  direction gets signficantly shorter with increasing condition number. This suggests an algorithm where the 
     *  newton step from the Hessian should be increased so as to match funcDecreaseNewtExp2 =  funcDecreaseNewt2.
     *  This roughly equals the ratio of the norms of the hessian and newton steps. This increased Newton step can
     *  then be used with the trust region double dogleg algorithm. 
     */
    if ((s_print_DogLeg && m_print_flag >= 3) || (doDogLeg_ && m_print_flag >= 5)) {
      printf("\t\t   descentComparison: initial rate of decrease of func in cauchy dir (expected) = %g\n", funcDecreaseSDExp);
      printf("\t\t   descentComparison: initial rate of decrease of func in cauchy dir            = %g\n", funcDecreaseSD);
      printf("\t\t   descentComparison: initial rate of decrease of func in newton dir (expected) = %g\n", funcDecreaseNewtExp2);
      printf("\t\t   descentComparison: initial rate of decrease of func in newton dir            = %g\n", funcDecreaseNewt2);
    }
    if ((s_print_DogLeg && m_print_flag >= 3) || (doDogLeg_ && m_print_flag >= 4)) {
      printf("\t\t   descentComparison: initial rate of decrease of Resid in cauchy dir (expected) = %g\n", ResidDecreaseSDExp_);
      printf("\t\t   descentComparison: initial rate of decrease of Resid in cauchy dir            = %g\n", ResidDecreaseSD_);
      printf("\t\t   descentComparison: initial rate of decrease of Resid in newton dir (expected) = %g\n", ResidDecreaseNewtExp_);
      printf("\t\t   descentComparison: initial rate of decrease of Resid in newton dir            = %g\n", ResidDecreaseNewt_);
    }

    if ((s_print_DogLeg && m_print_flag >= 5) || (doDogLeg_ && m_print_flag >= 5)) {
     if (funcDecreaseNewt2 >= 0.0) {
       printf("\t\t                            %13.5E  %22.16E\n", funcDecreaseNewtExp2, m_normResid_0);
       double ff = ffNewt * 1.0E-5;
       for (int ii = 0; ii < 13; ii++) {
	 ff *= 10.;
	 if (ii == 12) {
	   ff = ffNewt;
	 }
	 for (int i = 0; i < neq_; i++) {
	   y_n_1[i] = m_y_n_curr[i] + ff * deltaX_Newton_[i];
	 } 
	 numTrials += 1;
	 if (solnType_ != NSOLN_TYPE_STEADY_STATE) {
	   info = doResidualCalc(time_curr, solnType_, y_n_1, ydot1, Base_LaggedSolutionComponents);
	 } else {
	   info = doResidualCalc(time_curr, solnType_, y_n_1, ydot0, Base_LaggedSolutionComponents);
	 }  
	 residNewt = residErrorNorm(DATA_PTR(m_resid));
	 residNewt2 = residNewt * residNewt * neq_;
	 funcDecreaseNewt2 = 0.5 * (residNewt2 - normResid02) / ( ff * sNewt);
	 printf("\t\t                 %10.3E %13.5E  %22.16E\n", ff, funcDecreaseNewt2, residNewt );
       }

     }


   }

  }

  //====================================================================================================================
  //  Setup the parameters for the double dog leg
  /*
   *  The calls to the doCauchySolve() and doNewtonSolve() routines are done at the main level. This routine comes
   *  after those calls.  We calculate the point Nuu_ here, the distances of the dog-legs,
   *  and the norms of the CP and Newton points in terms of the trust vectors.
   */
  void NonlinearSolver::setupDoubleDogleg()
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
    // double  fac1 = sumG / lambdaStar_;
    // double  fac2 = sumH / lambdaStar_;
    // double  gamma = fac1 / fac2;
    // doublereal   gamma =  m_normDeltaSoln_CP / m_normDeltaSoln_Newton;
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
    doublereal residSteepLin = expectedResidLeg(0, 1.0);
    doublereal Nres2CP = residSteepLin * residSteepLin * neq_;
    doublereal Nres2_o = m_normResid_0 * m_normResid_0 * neq_;
    doublereal a = Nres2CP /  Nres2_o;
    doublereal betaEqual = (2.0 - sqrt(4.0 - 4 * (1.0 - a))) / 2.0;
    doublereal  beta = (1.0 + betaEqual) / 2.0;


    Nuu_ = beta;

    dist_R0_ = m_normDeltaSoln_CP;
    for (int i = 0; i < neq_; i++) {
      m_wksp[i] = Nuu_ * deltaX_Newton_[i] - deltaX_CP_[i];
    }
    dist_R1_ = solnErrorNorm(DATA_PTR(m_wksp));
    dist_R2_ = (1.0 - Nuu_) * m_normDeltaSoln_Newton;
    dist_Total_ = dist_R0_ + dist_R1_ + dist_R2_;

    /*
     * Calculate the trust distances
     */
    normTrust_Newton_ = calcTrustDistance(deltaX_Newton_);
    normTrust_CP_ = calcTrustDistance(deltaX_CP_);

  } 
  //====================================================================================================================
  // Change the global lambda coordinate into the (leg,alpha) coordinate for the double dogleg
  /*
   * @param lambda Global value of the distance along the double dogleg
   * @param alpha  relative value along the particular leg
   *
   * @return Returns the leg number ( 0, 1, or 2).
   */
  int NonlinearSolver::lambdaToLeg(const doublereal lambda, doublereal &alpha) const {

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
  // Calculated the expected residual along the double dogleg curve. 
  /*
   *  @param leg 0, 1, or 2 representing the curves of the dogleg
   *  @param alpha  Relative distance along the particular curve.
   *
   *  @return Returns the expected value of the residual at that point according to the quadratic model.
   *          The residual at the newton point will always be zero.
   */
  doublereal NonlinearSolver::expectedResidLeg(int leg, doublereal alpha) const {

    doublereal resD2, res2, resNorm;
    doublereal normResid02 = m_normResid_0 * m_normResid_0 * neq_;

    if (leg == 0) {
      /*
       * We are on the steepest descent line
       *  along that line 
       *   R2 = R2 + 2 lambda R dot Jd + lambda**2 Jd dot Jd
       */

      doublereal tmp = - 2.0 * alpha + alpha * alpha;
      doublereal tmp2 = - RJd_norm_ * lambdaStar_;
      resD2 = tmp2 * tmp;

    } else if (leg == 1) {

      /*
       *  Same formula as above for lambda=1.
       */
      doublereal tmp2 = - RJd_norm_ * lambdaStar_;
      doublereal RdotJS = - tmp2;
      doublereal JsJs =     tmp2;


      doublereal res0_2 = m_normResid_0 * m_normResid_0 * neq_;

      res2 =  res0_2 + (1.0 - alpha) * 2 * RdotJS - 2 * alpha * Nuu_ * res0_2
	+ (1.0 - alpha) * (1.0 - alpha) * JsJs
	+ alpha * alpha * Nuu_ * Nuu_ * res0_2 
	- 2 * alpha * Nuu_ * (1.0 - alpha) * RdotJS;

      resNorm = sqrt(res2 / neq_);
      return resNorm;

    } else {
      doublereal beta = Nuu_ + alpha * (1.0 - Nuu_);
      doublereal tmp2 =  normResid02;
      doublereal tmp = 1.0 - 2.0 * beta + 1.0 * beta * beta - 1.0;
      resD2 = tmp * tmp2;
    }

    res2 = m_normResid_0 * m_normResid_0 * neq_ + resD2;
    if (res2 < 0.0) {
      resNorm =  m_normResid_0 - sqrt(resD2/neq_);
    } else {
      resNorm = sqrt(res2 / neq_);
    }
      
    return resNorm;

  } 
  //====================================================================================================================
  //  Here we print out the residual at various points along the double dogleg, comparing against the quadratic model
  //  in a table format
  /*
   *  @param time_curr     INPUT    current time
   *  @param ydot0         INPUT    Current value of the derivative of the solution vector for non-time dependent
   *                                determinations
   *  @param legBest       OUTPUT   leg of the dogleg that gives the lowest residual
   *  @param alphaBest     OUTPUT   distance along dogleg for best result.
   */
  void NonlinearSolver::residualComparisonLeg(const doublereal time_curr, const doublereal * const ydot0, int &legBest,
					      doublereal &alphaBest) const {
    doublereal *y1 = DATA_PTR(m_wksp);
    doublereal *ydot1 = DATA_PTR(m_wksp_2);
    doublereal sLen;
    doublereal alpha;

    doublereal residSteepBest = 1.0E300;
    doublereal residSteepLinBest = 0.0;
    if (s_print_DogLeg || (doDogLeg_ && m_print_flag > 6)) {
      printf("\t\t   residualComparisonLeg() \n");
      printf("\t\t          Point               StepLen     Residual_Actual  Residual_Linear  RelativeMatch\n");
    }
    // First compare at 1/4 along SD curve
    std::vector<doublereal> alphaT;
    alphaT.push_back(0.00);
    alphaT.push_back(0.01);
    alphaT.push_back(0.1);
    alphaT.push_back(0.25);
    alphaT.push_back(0.50);
    alphaT.push_back(0.75);
    alphaT.push_back(1.0);
    for (int iteration = 0; iteration < (int) alphaT.size(); iteration++) {
      alpha = alphaT[iteration];
      for (int i = 0; i < neq_; i++) {
	y1[i] = m_y_n_curr[i] + alpha * deltaX_CP_[i];
      }
      if (solnType_ != NSOLN_TYPE_STEADY_STATE) {
	calc_ydot(m_order, y1, ydot1);
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


      doublereal residSteep = residErrorNorm(DATA_PTR(m_resid));
      doublereal residSteepLin = expectedResidLeg(0, alpha);
      if (residSteep < residSteepBest) {
	legBest = 0;
        alphaBest = alpha;
	residSteepBest = residSteep;
	residSteepLinBest = residSteepLin;
      }

      doublereal relFit = (residSteep - residSteepLin) / (fabs(residSteepLin) + 1.0E-10);
      if (s_print_DogLeg || (doDogLeg_ && m_print_flag > 6)) {
	printf("\t\t      (%2d - % 10.3g)  % 15.8E  % 15.8E % 15.8E  % 15.8E\n", 0, alpha, sLen, residSteep, residSteepLin , relFit);
      }
    }

    for (int iteration = 0; iteration < (int) alphaT.size(); iteration++) {
      doublereal alpha = alphaT[iteration];
      for (int i = 0; i < neq_; i++) {
	y1[i] = m_y_n_curr[i] + (1.0 - alpha) * deltaX_CP_[i];
	y1[i] += alpha * Nuu_ * deltaX_Newton_[i];
      }
      if (solnType_ != NSOLN_TYPE_STEADY_STATE) {
	calc_ydot(m_order, y1, ydot1);
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
	y1[i] -= m_y_n_curr[i];
      }
      sLen = solnErrorNorm(DATA_PTR(y1));

      doublereal residSteep = residErrorNorm(DATA_PTR(m_resid));
      doublereal residSteepLin = expectedResidLeg(1, alpha);
      if (residSteep < residSteepBest) {
	legBest = 1;
        alphaBest = alpha;
	residSteepBest = residSteep;
	residSteepLinBest = residSteepLin;
      }

      doublereal relFit = (residSteep - residSteepLin) / (fabs(residSteepLin) + 1.0E-10);
      if (s_print_DogLeg || (doDogLeg_ && m_print_flag > 6)) {
	printf("\t\t      (%2d - % 10.3g) % 15.8E   % 15.8E  % 15.8E  % 15.8E\n", 1, alpha, sLen,  residSteep, residSteepLin , relFit);
      }
    }

    for (int iteration = 0; iteration < (int) alphaT.size(); iteration++) {
      doublereal alpha = alphaT[iteration];
      for (int i = 0; i < neq_; i++) {
	y1[i] = m_y_n_curr[i] + ( Nuu_ + alpha * (1.0 - Nuu_))* deltaX_Newton_[i];
      }
      if (solnType_ != NSOLN_TYPE_STEADY_STATE) {
	calc_ydot(m_order, y1, ydot1);
      }
      sLen = ( Nuu_ + alpha * (1.0 - Nuu_)) * solnErrorNorm(DATA_PTR(deltaX_Newton_));
      /*
       *  Calculate the residual that would result if y1[] were the new solution vector
       *  -> m_resid[] contains the result of the residual calculation
       */
      if (solnType_ != NSOLN_TYPE_STEADY_STATE) {
	doResidualCalc(time_curr, solnType_, y1, ydot1, Base_LaggedSolutionComponents);
      } else {
	doResidualCalc(time_curr, solnType_, y1, ydot0, Base_LaggedSolutionComponents);
      }  



      doublereal residSteep = residErrorNorm(DATA_PTR(m_resid));
      doublereal residSteepLin = expectedResidLeg(2, alpha);
      if (residSteep < residSteepBest) {
	legBest = 2;
        alphaBest = alpha;
	residSteepBest = residSteep;
	residSteepLinBest = residSteepLin;
      }
      doublereal relFit = (residSteep - residSteepLin) / (fabs(residSteepLin) + 1.0E-10);
      if (s_print_DogLeg || (doDogLeg_ && m_print_flag > 6)) {
	printf("\t\t      (%2d - % 10.3g)  % 15.8E % 15.8E  % 15.8E  % 15.8E\n", 2, alpha, sLen, residSteep, residSteepLin , relFit);
      }
    }
    if (s_print_DogLeg || (doDogLeg_ && m_print_flag > 6)) {
      printf("\t\t       Best Result: \n");
      doublereal relFit = (residSteepBest - residSteepLinBest) / (fabs(residSteepLinBest) + 1.0E-10);
      if (m_print_flag <= 6) {
	printf("\t\t       Leg %2d alpha %5g: NonlinResid = %g LinResid = %g, relfit = %g\n",
	       legBest, alphaBest, residSteepBest, residSteepLinBest, relFit);
      } else {
	if (legBest == 0) {
	  sLen = alpha * solnErrorNorm(DATA_PTR(deltaX_CP_));
	} else if (legBest == 1) {
	  for (int i = 0; i < neq_; i++) {
	    y1[i] = (1.0 - alphaBest) * deltaX_CP_[i];
	    y1[i] += alphaBest * Nuu_ * deltaX_Newton_[i];
	  }
	  sLen = solnErrorNorm(DATA_PTR(y1));
	} else {
	  sLen = ( Nuu_ + alpha * (1.0 - Nuu_)) * solnErrorNorm(DATA_PTR(deltaX_Newton_));
	}
	printf("\t\t      (%2d - % 10.3g)  % 15.8E % 15.8E  % 15.8E  % 15.8E\n", legBest, alphaBest, sLen, 
	       residSteepBest, residSteepLinBest , relFit);
      }
    }

  }
  //====================================================================================================================
  // Calculate the length of the current trust region in terms of the solution error norm
  /*
   *  We carry out a norm of deltaX_trust_ first. Then, we multiply that value
   *  by trustDelta_
   */
  doublereal  NonlinearSolver::trustRegionLength() const
  {
    norm_deltaX_trust_ = solnErrorNorm(DATA_PTR(deltaX_trust_));
    return (trustDelta_ * norm_deltaX_trust_);
  }
  //====================================================================================================================
  void  NonlinearSolver::setDefaultDeltaBoundsMagnitudes()
  {
    for (int i = 0; i < neq_; i++) {
      m_deltaStepMinimum[i] = 1000. * atolk_[i];
      m_deltaStepMinimum[i] = MAX(m_deltaStepMinimum[i], 0.1 * fabs(m_y_n_curr[i]));
    }
  }
  //====================================================================================================================
  void  NonlinearSolver::adjustUpStepMinimums() {
    for (int i = 0; i < neq_; i++) {
      doublereal goodVal = deltaX_trust_[i] * trustDelta_;
      if (deltaX_trust_[i] * trustDelta_ >   m_deltaStepMinimum[i]) {
	m_deltaStepMinimum[i] = 1.1 * goodVal;
      }
 
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
   *  @param y_n_curr  Initial value of the solution vector
   *  @param step_1  initial proposed step size
   *
   *  @return returns the damping factor
   */
  double
  NonlinearSolver::deltaBoundStep(const doublereal * const y_n_curr, const doublereal * const step_1) {
			       
    int i_fbounds = 0;
    int ifbd = 0;
    int i_fbd = 0;
    doublereal UPFAC = 2.0;
    
    doublereal sameSign = 0.0;
    doublereal ff;
    doublereal f_delta_bounds = 1.0;
    doublereal ff_alt;
    for (int i = 0; i < neq_; i++) {
      doublereal y_new = y_n_curr[i] + step_1[i];
      sameSign = y_new * y_n_curr[i];
     
      /*
       * Now do a delta bounds
       * Increase variables by a factor of UPFAC only
       * decrease variables by a factor of 2 only
       */
      ff = 1.0;


      if (sameSign >= 0.0) {
	if ((fabs(y_new) > UPFAC * fabs(y_n_curr[i])) && 
	    (fabs(y_new - y_n_curr[i]) > m_deltaStepMinimum[i])) {
	  ff = (UPFAC - 1.0) * fabs(y_n_curr[i]/(y_new - y_n_curr[i]));
	  ff_alt = fabs(m_deltaStepMinimum[i] / (y_new - y_n_curr[i]));
	  ff = MAX(ff, ff_alt);
	  ifbd = 1;
	}
	if ((fabs(2.0 * y_new) < fabs(y_n_curr[i])) &&
	    (fabs(y_new - y_n_curr[i]) > m_deltaStepMinimum[i])) {
	  ff = y_n_curr[i]/(y_new - y_n_curr[i]) * (1.0 - 2.0)/2.0;
	  ff_alt = fabs(m_deltaStepMinimum[i] / (y_new - y_n_curr[i]));
	  ff = MAX(ff, ff_alt);
	  ifbd = 0;
	}
      } else {
	/*
	 *  This handles the case where the value crosses the origin.
	 *       - First we don't let it cross the origin until its shrunk to the size of m_deltaStepMinimum[i]
	 */
	if (fabs(y_n_curr[i]) > m_deltaStepMinimum[i]) {
	  ff = y_n_curr[i]/(y_new - y_n_curr[i]) * (1.0 - 2.0)/2.0;
	  ff_alt = fabs(m_deltaStepMinimum[i] / (y_new - y_n_curr[i]));
	  ff = MAX(ff, ff_alt);
	  if (y_n_curr[i] >= 0.0) {
	    ifbd = 0;
	  } else {
	    ifbd = 1;
	  }
	}
	/*
	 *  Second when it does cross the origin, we make sure that its magnitude is only 50% of the previous value.
	 */
	else if (fabs(y_new) > 0.5 * fabs(y_n_curr[i])) {
	  ff = y_n_curr[i]/(y_new - y_n_curr[i]) * (-1.5);
	  ff_alt = fabs(m_deltaStepMinimum[i] / (y_new - y_n_curr[i]));
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
    if (m_print_flag >= 3) {
      if (f_delta_bounds < 1.0) {
	if (i_fbd) {
	  printf("\t\tdeltaBoundStep: Increase of Variable %d causing "
		 "delta damping of %g: origVal = %10.3g, undampedNew = %10.3g, dampedNew = %10.3g\n",
		 i_fbounds, f_delta_bounds, y_n_curr[i_fbounds], y_n_curr[i_fbounds] + step_1[i_fbounds],
		 y_n_curr[i_fbounds] +  f_delta_bounds * step_1[i_fbounds] );
	} else {
	  printf("\t\tdeltaBoundStep: Decrease of variable %d causing"
		 "delta damping of %g: origVal = %10.3g, undampedNew = %10.3g, dampedNew = %10.3g\n",
		 i_fbounds, f_delta_bounds, y_n_curr[i_fbounds], y_n_curr[i_fbounds] + step_1[i_fbounds],
		 y_n_curr[i_fbounds] +  f_delta_bounds * step_1[i_fbounds]);
	}
      }
    }
    
    
    return f_delta_bounds;
  }
  //====================================================================================================================
  // Readjust the trust region vectors
  /*
   *  The trust region is made up of the trust region vector calculation and the trustDelta_ value
   *  We periodically recalculate the trustVector_ values so that they renormalize to the
   *  correct length.
   */
  void  NonlinearSolver::readjustTrustVector()
  {
    doublereal trustDeltaOld = trustDelta_;
    doublereal wtSum = 0.0;
    for (int i = 0; i < neq_; i++) {
      wtSum += m_ewt[i];
    }
    wtSum /= neq_;
    doublereal trustNorm = solnErrorNorm(DATA_PTR(deltaX_trust_));
    doublereal deltaXSizeOld = trustNorm;
    doublereal trustNormGoal = trustNorm * trustDelta_;

    // This is the size of each component.
    //    doublereal trustDeltaEach = trustDelta_ * trustNorm / neq_;
    doublereal oldVal; 
    doublereal fabsy;
    // we use the old value of the trust region as an indicator
    for (int i = 0; i < neq_; i++) {
      oldVal = deltaX_trust_[i];
      fabsy = fabs(m_y_n_curr[i]);
      // First off make sure that each trust region vector is 1/2 the size of each variable or smaller
      // unless overridden by the deltaStepMininum value. 
      //      doublereal newValue =  trustDeltaEach * m_ewt[i] / wtSum;
      doublereal newValue =  trustNormGoal  * m_ewt[i];
      if (newValue > 0.5 * fabsy) {
	if (fabsy * 0.5 > m_deltaStepMinimum[i]) {
	  deltaX_trust_[i] = 0.5 * fabsy;
	} else {
	  deltaX_trust_[i] = m_deltaStepMinimum[i];
	}
      } else {
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
    norm_deltaX_trust_ = solnErrorNorm(DATA_PTR(deltaX_trust_));
    doublereal  sum = trustNormGoal / trustNorm;
    for (int i = 0; i < neq_; i++) {
      deltaX_trust_[i] = deltaX_trust_[i] * sum;
    } 
    norm_deltaX_trust_ = solnErrorNorm(DATA_PTR(deltaX_trust_));
    trustDelta_ = trustNormGoal /  norm_deltaX_trust_;

    if (doDogLeg_ && m_print_flag >= 4) {
      printf("\t\t   reajustTrustVector(): Trust size = %11.3E: Old deltaX size = %11.3E trustDelta_ = %11.3E\n"
	     "\t\t                                                     new deltaX size = %11.3E trustdelta_ = %11.3E\n",
	     trustNormGoal,  deltaXSizeOld,  trustDeltaOld, norm_deltaX_trust_, trustDelta_ );
    }
  } 
  //====================================================================================================================
  //! Initialize the size of the trust vector.
  /*!
   *  The algorithm we use is to set it equal to the length of the Distance to the Cauchy point.
   */
  void  NonlinearSolver::initializeTrustRegion() 
  {
    if (trustRegionInitializationMethod_ == 0) {
      return;
    }
    if (trustRegionInitializationMethod_ == 1) {
      for (int i = 0; i < neq_; i++) {
	deltaX_trust_[i] = m_ewt[i] * trustRegionInitializationFactor_;
       }
      trustDelta_ = 1.0;
    }
    if (trustRegionInitializationMethod_ == 2) {
      for (int i = 0; i < neq_; i++) {
	deltaX_trust_[i] = m_ewt[i] * m_normDeltaSoln_CP * trustRegionInitializationFactor_;
       }
      doublereal cpd = calcTrustDistance(deltaX_CP_);
      if ((doDogLeg_ && m_print_flag >= 4)) {
	printf("\t\t   initializeTrustRegion(): Relative Distance of Cauchy Vector wrt Trust Vector = %g\n", cpd);
      }
      trustDelta_ = trustDelta_ * cpd * trustRegionInitializationFactor_;
      readjustTrustVector();
      cpd = calcTrustDistance(deltaX_CP_);
      if ((doDogLeg_ && m_print_flag >= 4)) {
	printf("\t\t   initializeTrustRegion(): Relative Distance of Cauchy Vector wrt Trust Vector = %g\n", cpd);
      }
    }
    if (trustRegionInitializationMethod_ == 3) {
      for (int i = 0; i < neq_; i++) {
	deltaX_trust_[i] = m_ewt[i] *  m_normDeltaSoln_Newton * trustRegionInitializationFactor_;
       }
      doublereal cpd = calcTrustDistance(deltaX_Newton_);
      if ((doDogLeg_ && m_print_flag >= 4)) {
	printf("\t\t   initializeTrustRegion(): Relative Distance of Newton Vector wrt Trust Vector = %g\n", cpd);
      }
      trustDelta_ = trustDelta_ * cpd;
      readjustTrustVector();
      cpd = calcTrustDistance(deltaX_Newton_);
      if ((doDogLeg_ && m_print_flag >= 4)) {
	printf("\t\t   initializeTrustRegion(): Relative Distance of Newton Vector wrt Trust Vector = %g\n", cpd);
      }
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
  void NonlinearSolver::fillDogLegStep(int leg, doublereal alpha, std::vector<doublereal>  & deltaX) const {
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
  // Given a trust distance, this routine calculates the intersection of the this distance with the
  // double dogleg curve
  /*
   *   @param      trustDelta  (INPUT)     Value of the trust distance
   *   @param      lambda      (OUTPUT)    Returns the internal coordinate of the double dogleg
   *   @param      alpha       (OUTPUT)    Returns the relative distance along the appropriate leg
   *   @return     leg         (OUTPUT)    Returns the leg ID (0, 1, or 2)
   */
  int NonlinearSolver::calcTrustIntersection(doublereal trustDelta, doublereal &lambda, doublereal &alpha) const
  {
    doublereal dist;
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
    doublereal sumv = 0.0;
    for (int i = 0; i < neq_; i++) {
      sumv += (deltaX_Newton_[i] / deltaX_trust_[i]) * (deltaX_CP_[i] / deltaX_trust_[i]);
    }

    doublereal a = normTrust_Newton_ * normTrust_Newton_ * Nuu_ * Nuu_;
    doublereal b = 2.0 * Nuu_ * sumv;
    doublereal c = normTrust_CP_ *  normTrust_CP_  - trustDelta * trustDelta;

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
  doublereal NonlinearSolver::boundStep(const doublereal * const y, const doublereal * const step0)
  {
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
    if (m_print_flag >= 3) {
      if (f_bounds != 1.0) {
	printf("\t\tboundStep: Variable %d causing bounds damping of %g\n", i_lower, f_bounds);
      }
    }

    doublereal f_delta_bounds = deltaBoundStep(y, step0);
    fbound = MIN(f_bounds, f_delta_bounds);

    return fbound;
  }
  //===================================================================================================================
  // Find a damping coefficient through a look-ahead mechanism
  /* 
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
   *           NSOLN_RETN_FAIL_DAMPSTEP 
   *             Unsuccessful step. We can not find a damping factor that is suitable.
   */
  int NonlinearSolver::dampStep(const doublereal time_curr, const doublereal * const y_n_curr, 
				const doublereal * const ydot_n_curr, doublereal * const step_1, 
				doublereal * const y_n_1, doublereal * const ydot_n_1, doublereal * const step_2,
				doublereal & stepNorm_2, GeneralMatrix& jac, bool writetitle, int& num_backtracks) 
  {  
    int j, m;
    int info = 0;
    int retnTrial = NSOLN_RETN_FAIL_DAMPSTEP;
    // Compute the weighted norm of the undamped step size step_1
    doublereal stepNorm_1 = solnErrorNorm(step_1);

    doublereal * step_1_orig = DATA_PTR(m_wksp);
    for (j = 0; j < neq_; j++) {
      step_1_orig[j] = step_1[j];
    }
   

    // Compute the multiplier to keep all components in bounds.A value of one indicates that there is no limitation
    // on the current step size in the nonlinear method due to bounds constraints (either negative values of delta
    // bounds constraints.
    m_dampBound = boundStep(y_n_curr, step_1);

    // If fbound is very small, then y0 is already close to the boundary and step0 points out of the allowed domain. In
    // this case, the Newton algorithm fails, so return an error condition.
    if (m_dampBound < 1.e-30) {
      if (m_print_flag > 1) printf("\t\t\tdampStep(): At limits.\n");
      return -3;
    }

    //--------------------------------------------
    //           Attempt damped step
    //-------------------------------------------- 

    // damping coefficient starts at 1.0
    m_dampRes = 1.0;
 
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
	step_1[j] =  ff * step_1_orig[j];
	y_n_1[j] = y_n_curr[j] + step_1[j];
      }
      
      if (solnType_ != NSOLN_TYPE_STEADY_STATE) {
	calc_ydot(m_order, y_n_1, ydot_n_1);
      }
      /*
       *  Calculate the residual that would result if y1[] were the new solution vector
       *  -> m_resid[] contains the result of the residual calculation
       */
      if (solnType_ != NSOLN_TYPE_STEADY_STATE) {
	info = doResidualCalc(time_curr, solnType_, y_n_1, ydot_n_1, Base_LaggedSolutionComponents);
      } else {
	info = doResidualCalc(time_curr, solnType_, y_n_1, ydot_n_curr, Base_LaggedSolutionComponents);
      }
      if (info != 1) {
	if (m_print_flag > 0) {
	  printf("\t\t\tdampStep(): current trial step and damping led to Residual Calc ERROR %d. Bailing\n", info);
	}
	return -1;
      }
      m_normResidTrial = residErrorNorm(DATA_PTR(m_resid));
      m_normResid_1 =  m_normResidTrial;
      if (m == 0) {
	m_normResid_Bound = m_normResidTrial;
      }

      bool steepEnough = (m_normResidTrial <  m_normResid_0 * (0.9 * (1.0 - ff) * (1.0 - ff)* (1.0 - ff) + 0.1));

      if (m_normResidTrial < 1.0 || steepEnough) {
	if (m_print_flag >= 5) {
	  if (m_normResidTrial < 1.0) {
	    printf("\t  dampStep(): Current trial step and damping"
		   " coefficient accepted because residTrial test step < 1:\n");
	    printf("\t              resid0 = %g, residTrial = %g\n", m_normResid_0, m_normResidTrial);
	  } else if (steepEnough) {
	    printf("\t  dampStep(): Current trial step and damping"
		   " coefficient accepted because resid0 > residTrial and steep enough:\n");
	    printf("\t              resid0 = %g, residTrial = %g\n", m_normResid_0, m_normResidTrial);
	  } else {
	    printf("\t  dampStep(): Current trial step and damping"
		   " coefficient accepted because residual solution damping is turned off:\n");
	    printf("\t              resid0 = %g, residTrial = %g\n", m_normResid_0, m_normResidTrial);
	  }
	}
	/*
	 *  We aren't going to solve the system if we don't need to. Therefore, return an estimate
	 *  of the next solution update based on the ratio of the residual reduction.
	 */
	if (m_normResid_0 > 0.0) {
	  stepNorm_2 = stepNorm_1 * m_normResidTrial / m_normResid_0;
	}
	else {
	  stepNorm_2 = 0;
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
	info = doNewtonSolve(time_curr, y_n_1, ydot_n_1, step_2, jac);
      } else {
	info = doNewtonSolve(time_curr, y_n_1, ydot_n_curr, step_2, jac);
      }
      if (info) {
	if (m_print_flag > 0) {
	  printf("\t\t\tdampStep: current trial step and damping led to LAPACK ERROR %d. Bailing\n", info);
	}
	return -1;
      }

      // compute the weighted norm of step1
      stepNorm_2 = solnErrorNorm(step_2);

      // write log information
      if (m_print_flag >= 5) {
	print_solnDelta_norm_contrib((const doublereal *) step_1_orig,  "DeltaSoln",
				     (const doublereal *) step_2, "DeltaSolnTrial",
				     "dampNewt: Important Entries for Weighted Soln Updates:",
				     y_n_curr, y_n_1, ff, 5);
      }
      if (m_print_flag >= 4) {
	printf("\t\t\tdampStep(): s1 = %g, s2 = %g, dampBound = %g,"
	       "dampRes = %g\n",  stepNorm_1, stepNorm_2, m_dampBound, m_dampRes);
      }


      // if the norm of s1 is less than the norm of s0, then
      // accept this damping coefficient. Also accept it if this
      // step would result in a converged solution. Otherwise,
      // decrease the damping coefficient and try again.
	  
      if (stepNorm_2 < 0.8 || stepNorm_2 < stepNorm_1) {
	if (stepNorm_2 < 1.0) {
	  if (m_print_flag >= 3) {   
	    if (stepNorm_2 < 1.0) {
	      printf("\t\t\tdampStep: current trial step and damping coefficient accepted because test step < 1\n");
	      printf("\t\t\t          s2 = %g, s1 = %g\n", stepNorm_2, stepNorm_1);
	    }
	  }  
	  retnTrial = 2;
	} else {
	  retnTrial = 1;
	}
	break;
      } else {
	if (m_print_flag > 1) {
	  printf("\t\t\tdampStep: current step rejected: (s1 = %g > "
		 "s0 = %g)", stepNorm_2, stepNorm_1);
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
    // coefficient could be found, return NSOLN_RETN_FAIL_DAMPSTEP.
    if (m < NDAMP) {
      if (m_print_flag >= 4 ) {
	printf("\t  dampStep(): current trial step accepted retnTrial = %d, its = %d, damp = %g\n", retnTrial, m+1, ff);
      }
      return retnTrial;
    } else {
      if (stepNorm_2 < 0.5 && (stepNorm_1 < 0.5)) {
	if (m_print_flag >= 4 ) {
	  printf("\t  dampStep(): current trial step accepted kindof retnTrial = %d, its = %d, damp = %g\n", 2, m+1, ff);
	}
	return 2;
      }
      if (stepNorm_2 < 1.0) {
	if (m_print_flag >= 4 ) {
	  printf("\t  dampStep(): current trial step accepted and soln converged retnTrial ="
		 "%d, its = %d, damp = %g\n", 0, m+1, ff);
	}
	return 0;
      }
    }
    if (m_print_flag >= 4 ) {
      printf("\t  dampStep(): current direction is rejected! retnTrial = %d, its = %d, damp = %g\n", 
	     NSOLN_RETN_FAIL_DAMPSTEP, m+1, ff);
    }
    return NSOLN_RETN_FAIL_DAMPSTEP;
  }
  //====================================================================================================================
  // Damp using the dog leg approach
  /*
   *  
   * @param time_curr  INPUT     Current value of the time
   * @param y_n_curr   INPUT    Current value of the solution vector
   * @param ydot_n_curr INPUT   Current value of the derivative of the solution vector
   * @param step_1     INPUT    First trial step for the first iteration
   * @param y_n_1   INPUT       First trial value of the solution vector
   * @param ydot_n_1 INPUT      First trial value of the derivative of the solution vector
   * @param s1         OUTPUT   Norm of the vector step_1
   * @param jac        INPUT    jacobian
   * @param numTrials OUTPUT  number of trials taken in the current damping step
   *
   *
   *  @return  1 Success:   Good step was taken. The predicted residual norm is less than one
   *           2 Success:   Good step: Next step's norm is less than 0.8
   *           3 Success:   The final residual is less than 1.0
   *                        A predicted deltaSoln1 is not produced however. s1 is estimated.
   *           4 Success:   The final residual is less than the residual from the previous step.
   *                        A predicted deltaSoln1 is not produced however. s1 is estimated.
   *           0 Unknown    Uncertain Success: s1 is about the same as s0
   *           NSOLN_RETN_FAIL_DAMPSTEP 
   *                        Unsuccessful step. Can not find a damping coefficient that is suitable
   */
  int NonlinearSolver::dampDogLeg(const doublereal time_curr, const doublereal* y_n_curr, 
				  const doublereal *ydot_n_curr, std::vector<doublereal> & step_1,
				  doublereal* const y_n_1, doublereal* const ydot_n_1, 
				  doublereal& stepNorm_1,   doublereal& stepNorm_2, GeneralMatrix& jac, int& numTrials) 
  {
    doublereal lambda;
    int info;
   
    bool success = false;
    int retn = 0;
    bool haveASuccess = false;
    doublereal trustDeltaOld = trustDelta_;
    doublereal* stepLastGood = DATA_PTR(m_wksp);
    //--------------------------------------------
    //           Attempt damped step
    //-------------------------------------------- 

    // damping coefficient starts at 1.0
    m_dampRes = 1.0;
    int j, m;
    doublereal tlen;
   

    for (m = 0; m < NDAMP; m++) {
      numTrials++;
      /*
       *  Find the initial value of lambda that satisfies the trust distance, trustDelta_
       */
      dogLegID_ = calcTrustIntersection(trustDelta_, lambda, dogLegAlpha_);
      if (m_print_flag >= 4) {
	tlen = trustRegionLength();
	printf("\t\t   dampDogLeg: trust region with length %13.5E has intersection at leg = %d, alpha = %g, lambda = %g\n",
	       tlen, dogLegID_, dogLegAlpha_, lambda);
      }
      /*
       *  Figure out the new step vector, step_1, based on (leg, alpha). Here we are using the
       *  intersection of the trust oval with the dog-leg curve.
       */
      fillDogLegStep(dogLegID_, dogLegAlpha_, step_1);

      /*
       *  OK, now that we have step0, Bound the step 
       */
      m_dampBound = boundStep(y_n_curr, DATA_PTR(step_1));
      /*
       * Decrease the step length if we are bound
       */
      if (m_dampBound < 1.0) {
	for (j = 0; j < neq_; j++) {
	  step_1[j] =  step_1[j] * m_dampBound;
	}
      }
      /*
       *  Calculate the new solution value y1[] given the step size
       */
      for (j = 0; j < neq_; j++) {
	y_n_1[j] = y_n_curr[j] + step_1[j];
      }
      /*
       *  Calculate the new solution time derivative given the step size
       */
      if (solnType_ != NSOLN_TYPE_STEADY_STATE) {
	calc_ydot(m_order, y_n_1, ydot_n_1);
      } 
      /*
       * OK, we have the step0. Now, ask the question whether it satisfies the acceptance criteria
       * as a good step. The overall outcome is returned in the variable info.
       */
      info = decideStep(time_curr, dogLegID_, dogLegAlpha_, y_n_curr, ydot_n_curr, step_1, 
			y_n_1, ydot_n_1, trustDeltaOld);
      m_normResid_Bound = m_normResid_1;

      /*
       *  The algorithm failed to find a solution vector sufficiently different than the current point
       */
      if (info == -1) {

	if (m_print_flag >= 1) {
	  doublereal stepNorm =  solnErrorNorm(DATA_PTR(step_1));
	  printf("\t\t   dampDogLeg: Current direction rejected, update became too small %g\n", stepNorm);
	  success = false;
	  retn = NSOLN_RETN_FAIL_STEPTOOSMALL;
	  break;
	}
      }
      if (info == -2) {
	if (m_print_flag >= 1) {
	  printf("\t\t dampDogLeg: current trial step and damping led to LAPACK ERROR %d. Bailing\n", info);
	  success = false;
	  retn = NSOLN_RETN_MATRIXINVERSIONERROR;
	  break;
	}
      }
      if (info == 0) {
	success = true;
	break;
      }
      if (info == 3) { 

	haveASuccess = true;
	// Store the good results in stepLastGood
	mdp::mdp_copy_dbl_1(DATA_PTR(stepLastGood), CONSTD_DATA_PTR(step_1), neq_);
	// Within the program decideStep(), we have already increased the value of trustDelta_. We store the
	// value of step0 in step1, recalculate a larger step0 in the next fillDogLegStep(),
	// and then attempt to see if the larger step works in the next iteration
      }
      if (info == 2) {
	// Step was a failure. If we had a previous success with a smaller stepsize, haveASuccess is true
	// and we execute the next block and break. If we didn't have a previous success, trustDelta_ has
	// already been decreased in the decideStep() routine. We go back and try another iteration with
	// a smaller trust region. 
	if (haveASuccess) {
	  mdp::mdp_copy_dbl_1(DATA_PTR(step_1), CONSTD_DATA_PTR(stepLastGood), neq_);
	  for (j = 0; j < neq_; j++) {
	    y_n_1[j] = y_n_curr[j] + step_1[j];
	  }    
	  if (solnType_ != NSOLN_TYPE_STEADY_STATE) {
	    calc_ydot(m_order, y_n_1, ydot_n_1);
	  } 
	  success = true;
	  break;
	} else {

	}
      }
    }

    /*
     * Estimate s1, the norm after the next step
     */
    stepNorm_1 =  solnErrorNorm(DATA_PTR(step_1));
    stepNorm_2 = stepNorm_1;
    if (m_dampBound < 1.0) {
      stepNorm_2 /= m_dampBound;
    }
    stepNorm_2 /= lambda;
    stepNorm_2 *=  m_normResidTrial / m_normResid_0;
 

    if (success) {
      if (m_normResidTrial < 1.0) {
	if (normTrust_Newton_ < trustDelta_ && m_dampBound == 1.0) {
	  return 1;
	} else {
	  return 0;
	}
      }
      return 0;
    }
    return NSOLN_RETN_FAIL_DAMPSTEP;
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
  int NonlinearSolver::decideStep(const doublereal time_curr, int leg, doublereal alpha,
				  const doublereal * const y_n_curr, 
				  const doublereal * const ydot_n_curr, const std::vector<doublereal> & step_1,
				  const doublereal * const y_n_1, const doublereal * const ydot_n_1, 
				  doublereal trustDeltaOld) 
  {
    int retn = 2;
    bool goodStep = false;
    int info;
    doublereal ll;
    // Calculate the solution step length
    doublereal stepNorm = solnErrorNorm(DATA_PTR(step_1));

    // Calculate the initial (R**2 * neq) value for the old function
    doublereal normResid0_2 = m_normResid_0 * m_normResid_0 * neq_;

    // Calculate the distance to the cauchy point
    doublereal cauchyDistanceNorm = solnErrorNorm(DATA_PTR(deltaX_CP_));

    // This is the expected inital rate of decrease in the cauchy direction.
    //   -> This is Eqn. 29 = Rhat dot Jhat dy / || d ||
    doublereal funcDecreaseSDExp = RJd_norm_ / cauchyDistanceNorm * lambdaStar_;
    if (funcDecreaseSDExp > 0.0) {
      if (m_print_flag >= 5) {
	printf("\t\tdecideStep(): Unexpected condition -> cauchy slope is positive\n");
      }
    }
        
    /*
     *  Calculate the residual that would result if y1[] were the new solution vector.
     *  The Lagged solution components are kept lagged here. Unfortunately, it just doesn't work in some cases to use a
     *  Jacobian from a lagged state and then use a residual from an unlagged condition.  The linear model doesn't
     *  agree with the nonlinear model. 
     *  -> m_resid[] contains the result of the residual calculation
     */
    if (solnType_ != NSOLN_TYPE_STEADY_STATE) {
      info = doResidualCalc(time_curr, solnType_, y_n_1, ydot_n_1, Base_LaggedSolutionComponents);
    } else {
      info = doResidualCalc(time_curr, solnType_, y_n_1, ydot_n_curr, Base_LaggedSolutionComponents);
    }

    if (info != 1) {
      if (m_print_flag >= 2) {
	printf("\t\tdecideStep: current trial step and damping led to Residual Calc ERROR %d. Bailing\n", info);
      }
      return -2;
    }
    /*
     *  Ok we have a successful new residual. Calculate the normalized residual value and store it in 
     *  m_normResidTrial
     */
    m_normResidTrial = residErrorNorm(DATA_PTR(m_resid));
    doublereal normResidTrial_2 = neq_ *  m_normResidTrial *  m_normResidTrial;

    /*
     *  We have a minimal acceptance test for passage. deltaf < 1.0E-4 (CauchySlope) (deltS)
     *  This is the condition that D&S use in 6.4.5
     */
    doublereal funcDecrease = 0.5 * (normResidTrial_2 - normResid0_2);
    doublereal acceptableDelF =  funcDecreaseSDExp * stepNorm * 1.0E-4;
    if (funcDecrease < acceptableDelF) {
      m_normResid_1 = m_normResidTrial;
      goodStep = true;
      m_normResid_1 =  m_normResidTrial;
      retn = 0;
      if (m_print_flag >= 4) {
	printf("\t\t   decideStep: Norm Residual(leg=%1d, alpha=%10.2E) = %11.4E passes\n",
	       dogLegID_, dogLegAlpha_, m_normResidTrial);
      }
    } else {
      if (m_print_flag >= 4) {
	printf("\t\t   decideStep: Norm Residual(leg=%1d, alpha=%10.2E) = %11.4E failes\n", 
	       dogLegID_, dogLegAlpha_, m_normResidTrial);
      }
      trustDelta_ *= 0.33;
      CurrentTrustFactor_ *= 0.33;
      retn = 2;
      // error condition if step is getting too small
      if (rtol_ * stepNorm  < 1.0E-6) {
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
      //  trustDelta_ *= 0.5;
      // NextTrustFactor_ *= 0.5;
      // ll = trustRegionLength();
      // if (m_print_flag >= 5) {
      //	printf("\t\tdecideStep(): Trust region decreased from %g to %g due to bounds constraint\n", ll*2, ll);
      //}
    } else {
      retn = 0;
      /*
       *  Calculate the expected residual from the quadratic model
       */
      doublereal expectedNormRes = expectedResidLeg(leg, alpha);
      doublereal expectedFuncDecrease = 0.5 * (neq_ * expectedNormRes * expectedNormRes - normResid0_2);
      if (funcDecrease > 0.1 * expectedFuncDecrease) {
	if ((m_normResidTrial > 0.5 * m_normResid_0) && (m_normResidTrial > 0.1)) {
	  trustDelta_ *= 0.5;
	  NextTrustFactor_ *= 0.5;
	  ll = trustRegionLength();
	  if (m_print_flag >= 4) {
	    printf("\t\t   decideStep: Trust region decreased from %g to %g due to bad quad approximation\n",
		   ll*2, ll);
	  }
	}
      } else {
	/*
	 *  If we are doing well, consider increasing the trust region and recalculating
	 */
	if (funcDecrease < 0.8 * expectedFuncDecrease ||  (m_normResidTrial < 0.33 * m_normResid_0)) {
	  if (trustDelta_ <= trustDeltaOld && (leg != 2 || alpha < 0.75) ) {
	    trustDelta_ *= 2.0;
	    CurrentTrustFactor_ *= 2;
	    adjustUpStepMinimums();
	    ll = trustRegionLength();
	    if (m_print_flag >= 4) {
	      if  (m_normResidTrial < 0.33 * m_normResid_0) { 
		printf("\t\t   decideStep: Redo line search with trust region increased from %g to %g due to good nonlinear behavior\n",
		       ll*0.5, ll);
	      } else {
		printf("\t\t   decideStep: Redi line search with trust region increased from %g to %g due to good linear model approximation\n",
		       ll*0.5, ll);
	      }
	    }
	    retn = 3;
	  } else {
	    /*
	     *  Increase the size of the trust region for the next calculation
	     */
	    if (m_normResidTrial < 0.99 * expectedNormRes ||  (m_normResidTrial < 0.20 * m_normResid_0) ||
		(funcDecrease < -1.0E-50 && ( funcDecrease < 0.9 *expectedFuncDecrease)) ) {
	      if (leg == 2 && alpha == 1.0 ) {
		ll = trustRegionLength();
		if (ll < 2.0 * m_normDeltaSoln_Newton) {
		  trustDelta_ *= 2.0;
		  NextTrustFactor_ *= 2.0;
		  adjustUpStepMinimums();
		  ll = trustRegionLength();
		  if (m_print_flag >= 4) {
		    printf("\t\t   decideStep: Trust region further increased from %g to %g next step due to good linear model behavior\n",
			   ll*0.5, ll);
		  }
		}
	      } else {
		ll = trustRegionLength();
		trustDelta_ *= 2.0;
		NextTrustFactor_ *= 2.0;
		adjustUpStepMinimums();
		ll = trustRegionLength();
		if (m_print_flag >= 4) {
		  printf("\t\t   decideStep: Trust region further increased from %g to %g next step due to good linear model behavior\n",
			 ll*0.5, ll);
		}
	      }
	    }
	  }
	}
      }
    }
    return retn;
  }
  //====================================================================================================================
  /*
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
  int NonlinearSolver::solve_nonlinear_problem(int SolnType, doublereal * const y_comm, doublereal * const ydot_comm, 
					       doublereal CJ, doublereal time_curr,  GeneralMatrix& jac,
					       int &num_newt_its,  int &num_linear_solves,
					       int &num_backtracks,  int loglevelInput)
  {
    clockWC wc;
    int convRes = 0;
    solnType_ = SolnType;
    int info = 0;
    if (neq_ <= 0) {
      return 1;
    }

    bool m_residCurrent = false;
    num_linear_solves -= m_numTotalLinearSolves;
    int retnDamp = 0;
    int retnCode = 0;
    bool forceNewJac = false;

    if (jacCopyPtr_) {
      delete jacCopyPtr_;
    }
    jacCopyPtr_ = jac.duplMyselfAsGeneralMatrix();
  
    doublereal stepNorm_1;
    doublereal stepNorm_2;
#ifdef DEBUG_MODE
    int legBest;
    doublereal alphaBest;
#endif
    bool trInit = false;

  
    mdp::mdp_copy_dbl_1(DATA_PTR(m_y_n_curr), DATA_PTR(y_comm), neq_);
  
    if (SolnType != NSOLN_TYPE_STEADY_STATE || ydot_comm) {
      mdp::mdp_copy_dbl_1(DATA_PTR(m_ydot_n_curr), ydot_comm, neq_);
      mdp::mdp_copy_dbl_1(DATA_PTR(m_ydot_n_1), ydot_comm, neq_);
    }
    // Redo the solution weights every time we enter the function
    createSolnWeights(DATA_PTR(m_y_n_curr));
    m_normDeltaSoln_Newton = 1.0E1;
    bool frst = true;
    num_newt_its = 0;
    num_backtracks = 0;
    int i_numTrials;
    m_print_flag = loglevelInput;

    if (trustRegionInitializationMethod_ == 0) {
      trInit = true;
    } else if (trustRegionInitializationMethod_ == 1) {
      trInit = true;
      initializeTrustRegion();
    } else {
      mdp::mdp_init_dbl_1(DATA_PTR(deltaX_trust_), 1.0, neq_); 
      trustDelta_ = 1.0;
    }

    if (m_print_flag == 2 || m_print_flag == 3) {
      printf("\tsolve_nonlinear_problem():\n\n");
      if (doDogLeg_) {
	printf("\tWt Iter Resid NewJac log(CN)| dRdS_CDexp dRdS_CD dRdS_Newtexp dRdS_Newt |"
	       "DS_Cauchy  DS_Newton  DS_Trust | legID legAlpha Fbound  | CTF  NTF  | nTr|"
	       "DS_Final    ResidLag   ResidFull\n");
	printf("\t---------------------------------------------------------------------------------------------------"
	       "--------------------------------------------------------------------------------\n");
      } else {
	printf("\t Wt  Iter Resid NewJac |  Fbound  ResidBound  | DampIts Fdamp  DS_Step1   DS_Step2"
	       "ResidLag |  DS_Damp    DS_Newton  ResidFull\n");
	printf("\t--------------------------------------------------------------------------------------------------"
	       "----------------------------------\n");
      }
    }

    while (1 > 0) {

      CurrentTrustFactor_ = 1.0;
      NextTrustFactor_  = 1.0;
      ResidWtsReevaluated_ = false;
      i_numTrials = 0;
      /*
       * Increment Newton Solve counter
       */
      m_numTotalNewtIts++;
      num_newt_its++;
      m_numLocalLinearSolves = 0;

      if (m_print_flag > 3) {
	printf("\t");
	print_line("=", 119);
	printf("\tsolve_nonlinear_problem(): iteration %d:\n",
	       num_newt_its);
      }
      /*
       *  If we are far enough away from the solution, redo the solution weights and the trust vectors.
       */
      if (m_normDeltaSoln_Newton > 1.0E2) {
	createSolnWeights(DATA_PTR(m_y_n_curr));
#ifdef DEBUG_MODE
	if (trInit) {
	  readjustTrustVector();
	}
#else
	if (doDogLeg_ && trInit) {
	  readjustTrustVector();
	}
#endif
      } else {
	// Do this stuff every 5 iterations
        if ((num_newt_its % 5) == 1) {
	  createSolnWeights(DATA_PTR(m_y_n_curr));
#ifdef DEBUG_MODE
	  if (trInit) {
	    readjustTrustVector();
	  }
#else
	  if (doDogLeg_ && trInit) {
	    readjustTrustVector();
	  }
#endif
	}
      }
    
      /*
       * Set default values of Delta bounds constraints
       */
      if (!m_manualDeltaStepSet) {
	setDefaultDeltaBoundsMagnitudes();
      }

      // Check whether the Jacobian should be re-evaluated.
            
      forceNewJac = true;
            
      if (forceNewJac) {
	if (m_print_flag > 3) {
	  printf("\t   solve_nonlinear_problem(): Getting a new Jacobian\n");
	}
	info = beuler_jac(jac, DATA_PTR(m_resid), time_curr, CJ,  DATA_PTR(m_y_n_curr), 
			  DATA_PTR(m_ydot_n_curr), num_newt_its);
	if (info != 1) {
	  if (m_print_flag > 0) {
	    printf("\t   solve_nonlinear_problem(): Jacobian Formation Error: %d Bailing\n", info);
	  }
	  retnDamp = NSOLN_RETN_JACOBIANFORMATIONERROR ;
	  goto done;
	}
	m_residCurrent = true;
      } else {
	if (m_print_flag > 1) {
	  printf("\t   solve_nonlinear_problem(): Solving system with old jacobian\n");
	}
	m_residCurrent = false;
      }
      /*
       * Go get new scales
       */
      calcColumnScales();


      /*
       *  Calculate the base residual 
       */
      if (m_print_flag >= 6) {
	printf("\t   solve_nonlinear_problem(): Calculate the base residual\n");
      }
      info = doResidualCalc(time_curr, NSOLN_TYPE_STEADY_STATE, DATA_PTR(m_y_n_curr), DATA_PTR(m_ydot_n_curr));
      if (info != 1) {
	if (m_print_flag > 0) {
	  printf("\t   solve_nonlinear_problem(): Residual Calc ERROR %d. Bailing\n", info);
	}
	retnDamp = NSOLN_RETN_RESIDUALFORMATIONERROR;
	goto done;
      }

      /*
       * Scale the matrix and the rhs, if they aren't already scaled
       * Figure out and store the residual scaling factors.
       */
      scaleMatrix(jac, DATA_PTR(m_y_n_curr), DATA_PTR(m_ydot_n_curr), time_curr, num_newt_its);


      /*
       *  Optional print out the initial residual
       */
      if (m_print_flag >= 6) {
	m_normResid_0 = residErrorNorm(DATA_PTR(m_resid), "Initial norm of the residual", 10, DATA_PTR(m_y_n_curr));
      } else {
	m_normResid_0 = residErrorNorm(DATA_PTR(m_resid), "Initial norm of the residual", 0, DATA_PTR(m_y_n_curr));
	if (m_print_flag == 4 || m_print_flag == 5 ) {
	  printf("\t   solve_nonlinear_problem(): Initial Residual Norm = %13.4E\n", m_normResid_0);
	}
      }
    
      
#ifdef DEBUG_MODE
      if (m_print_flag > 3) {
	printf("\t   solve_nonlinear_problem(): Calculate the steepest descent direction and Cauchy Point\n");
      }
      m_normDeltaSoln_CP = doCauchyPointSolve(jac);
   
#else
      if (doDogLeg_) {
	if (m_print_flag > 3) {
	  printf("\t   solve_nonlinear_problem(): Calculate the steepest descent direction and Cauchy Point\n");
	}
	m_normDeltaSoln_CP = doCauchyPointSolve(jac);
      }
#endif

      // compute the undamped Newton step
      if (doAffineSolve_) {
	if (m_print_flag >= 4) {
	  printf("\t   solve_nonlinear_problem(): Calculate the Newton direction via an Affine solve\n");
	}
	info = doAffineNewtonSolve(DATA_PTR(m_y_n_curr), DATA_PTR(m_ydot_n_curr), DATA_PTR(deltaX_Newton_), jac);
      } else {
	if (m_print_flag >= 4) {
	  printf("\t   solve_nonlinear_problem(): Calculate the Newton direction via a Newton solve\n");
	}
	info = doNewtonSolve(time_curr, DATA_PTR(m_y_n_curr), DATA_PTR(m_ydot_n_curr), DATA_PTR(deltaX_Newton_), jac);
      }

      if (info) {
	retnDamp = NSOLN_RETN_MATRIXINVERSIONERROR;
	if (m_print_flag > 0) {
	  printf("\t   solve_nonlinear_problem(): Matrix Inversion Error: %d Bailing\n", info);
	}
	goto done;
      }
      mdp::mdp_copy_dbl_1(DATA_PTR(m_step_1), CONSTD_DATA_PTR(deltaX_Newton_), neq_);

      if (m_print_flag >= 6) {
	m_normDeltaSoln_Newton = solnErrorNorm(DATA_PTR(deltaX_Newton_),  "Initial Undamped Newton Step of the iteration", 10);
      } else {
	m_normDeltaSoln_Newton = solnErrorNorm(DATA_PTR(deltaX_Newton_), "Initial Undamped Newton Step of the iteration", 0);
      }

      if (m_numTotalNewtIts == 1) {
	if (trustRegionInitializationMethod_ == 2 || trustRegionInitializationMethod_ == 3) {
	  if (m_print_flag > 3) {
	    if (trustRegionInitializationMethod_ == 2) {
	      printf("\t   solve_nonlinear_problem(): Initialize the trust region size as the length of the Cauchy Vector times %f\n",
		     trustRegionInitializationFactor_);
	    } else {
	      printf("\t   solve_nonlinear_problem(): Initialize the trust region size as the length of the Newton Vector times %f\n",
		     trustRegionInitializationFactor_);
	    }
	  }
	  initializeTrustRegion();
	  trInit = true;
	}
      }
   

      if (doDogLeg_) { 



#ifdef DEBUG_MODE
	doublereal trustD = calcTrustDistance(m_step_1);
	if (m_print_flag >= 4) {
	  if (trustD > trustDelta_) {
	    printf("\t\t   Newton's method step size, %g trustVectorUnits, larger than trust region, %g trustVectorUnits\n",
		   trustD, trustDelta_);
	    printf("\t\t   Newton's method step size, %g trustVectorUnits, larger than trust region, %g trustVectorUnits\n",
		   trustD, trustDelta_);
	  } else {
	    printf("\t\t   Newton's method step size, %g trustVectorUnits, smaller than trust region, %g trustVectorUnits\n",
		   trustD, trustDelta_);
	  }
	}
#endif
      }

      /*
       * Filter out bad directions
       */
      filterNewStep(time_curr, DATA_PTR(m_y_n_curr), DATA_PTR(m_step_1));


      
      if (s_print_DogLeg && m_print_flag >= 4) {
	printf("\t   solve_nonlinear_problem(): Compare descent rates for Cauchy and Newton directions\n");
	descentComparison(time_curr, DATA_PTR(m_ydot_n_curr), DATA_PTR(m_ydot_n_1), i_numTrials);
      } else { 
	if (doDogLeg_) { 
	  descentComparison(time_curr, DATA_PTR(m_ydot_n_curr), DATA_PTR(m_ydot_n_1), i_numTrials);
	}
      }

    
 
      if (doDogLeg_) { 
	setupDoubleDogleg();  
#ifdef DEBUG_MODE
	if (s_print_DogLeg && m_print_flag >= 5) {
	  printf("\t   solve_nonlinear_problem(): Compare Linear and nonlinear residuals along double dog-leg path\n");
	  residualComparisonLeg(time_curr, DATA_PTR(m_ydot_n_curr), legBest, alphaBest);
	}
#endif
	if (m_print_flag >= 4) {
	  printf("\t   solve_nonlinear_problem(): Calculate damping along dog-leg path to ensure residual decrease\n");
	}
	retnDamp = dampDogLeg(time_curr, DATA_PTR(m_y_n_curr), DATA_PTR(m_ydot_n_curr), 
			      m_step_1, DATA_PTR(m_y_n_1), DATA_PTR(m_ydot_n_1), stepNorm_1, stepNorm_2, jac, i_numTrials);
      }
#ifdef DEBUG_MODE
      else {
	if (s_print_DogLeg && m_print_flag >= 5) {
	  printf("\t   solve_nonlinear_problem(): Compare Linear and nonlinear residuals along double dog-leg path\n");	  
	  residualComparisonLeg(time_curr, DATA_PTR(m_ydot_n_curr), legBest, alphaBest);
	}
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
	retnDamp = dampStep(time_curr, DATA_PTR(m_y_n_curr), DATA_PTR(m_ydot_n_curr), 
			    DATA_PTR(m_step_1), DATA_PTR(m_y_n_1), DATA_PTR(m_ydot_n_1), 
			    DATA_PTR(m_wksp_2), stepNorm_2, jac, frst, i_numTrials);
	frst = false;
	num_backtracks += i_numTrials;
	stepNorm_1 = solnErrorNorm(DATA_PTR(m_step_1)); 
      }


      /*
       * Impose the minimum number of newton iterations critera
       */
      if (num_newt_its < m_min_newt_its) {
	if (retnDamp > NSOLN_RETN_CONTINUE) {
	  if (m_print_flag > 2) {
	    printf("\t   solve_nonlinear_problem(): Damped Newton successful (m=%d) but minimum newton"
		   "iterations not attained. Resolving ...\n", retnDamp);
	  }
	  retnDamp = NSOLN_RETN_CONTINUE;
	}
      }

      /*
       * Impose max newton iteration
       */
      if (num_newt_its > maxNewtIts_) {
	retnDamp = NSOLN_RETN_MAXIMUMITERATIONSEXCEEDED;
	if (m_print_flag > 1) {
	  printf("\t   solve_nonlinear_problem(): Damped newton unsuccessful (max newts exceeded) sfinal = %g\n", 
		 stepNorm_1);
	}
      }

      /*
       *  Do a full residual calculation with the unlagged solution components.
       *  Then get the norm of the residual
       */
      info = doResidualCalc(time_curr, NSOLN_TYPE_STEADY_STATE, DATA_PTR(m_y_n_1), DATA_PTR(m_ydot_n_1));
      if (info != 1) {
	if (m_print_flag > 0) {
	  printf("\t   solve_nonlinear_problem(): current trial step and damping led to Residual Calc "
		 "ERROR %d. Bailing\n", info);
	}
        retnDamp = NSOLN_RETN_RESIDUALFORMATIONERROR;
	goto done;
      }
      if (m_print_flag >= 4) {
	m_normResid_full  = residErrorNorm(DATA_PTR(m_resid), " Resulting full residual norm", 10, DATA_PTR(m_y_n_1));
	if (fabs(m_normResid_full - m_normResid_1) >  1.0E-3 * ( m_normResid_1 + m_normResid_full + 1.0E-4)) {
	  if (m_print_flag >= 4) {
	    printf("\t    solve_nonlinear_problem(): Full residual norm changed from %g to %g due to "
		   "lagging of components\n",  m_normResid_1, m_normResid_full);
	  }
	}
      } else {
	m_normResid_full = residErrorNorm(DATA_PTR(m_resid));
      }

      /*
       *  Check the convergence criteria
       */
      convRes = 0;
      if (retnDamp > NSOLN_RETN_CONTINUE) {
	convRes = convergenceCheck(retnDamp, stepNorm_1);
      }      
   



      bool m_filterIntermediate = false;
      if (m_filterIntermediate) {
	if (retnDamp == NSOLN_RETN_CONTINUE) {
	  (void) filterNewSolution(time_n, DATA_PTR(m_y_n_1), DATA_PTR(m_ydot_n_1));
	}
      }

      // Exchange new for curr solutions
      if (retnDamp >= NSOLN_RETN_CONTINUE) {
	mdp::mdp_copy_dbl_1(DATA_PTR(m_y_n_curr), CONSTD_DATA_PTR(m_y_n_1), neq_);

	if (solnType_ != NSOLN_TYPE_STEADY_STATE) {
	  calc_ydot(m_order, DATA_PTR(m_y_n_curr), DATA_PTR(m_ydot_n_curr));
	}
      }

      if (m_print_flag == 2 || m_print_flag == 3) {
	//   printf("\t    Iter Resid NewJac  | Fbound  | ResidBound |   Fdamp DampIts |   DeltaSolnNewton   ResidFinal \n");
	if (ResidWtsReevaluated_) {
	  printf("\t*");
	} else {
	  printf("\t ");
	}
	printf(" %3d %11.3E", num_newt_its, m_normResid_0);
	bool m_jacAge = false;
	if (!m_jacAge) {
	  printf("   Y ");
	} else {
	  printf("   N ");
	}
	if (doDogLeg_) {
	  printf("%5.1F |", log10(m_conditionNumber));
	  //   printf("\t    Iter Resid NewJac  | DS_Cauchy  DS_Newton  DS_Trust |  legID legAlpha  Fbound  |     |   DS_F   ResidFinal \n");
	  printf("%10.3E %10.3E %10.3E %10.3E|", ResidDecreaseSDExp_, ResidDecreaseSD_,
		 ResidDecreaseNewtExp_, ResidDecreaseNewt_);
	  printf("%10.3E %10.3E %10.3E|", m_normDeltaSoln_CP , m_normDeltaSoln_Newton, norm_deltaX_trust_ * trustDelta_);
	  printf("%2d %10.2E %10.2E", dogLegID_ , dogLegAlpha_, m_dampBound);
	  printf("| %3.2f %3.2f |", CurrentTrustFactor_, NextTrustFactor_);
	  printf(" %2d ", i_numTrials);
	  printf("| %10.3E %10.3E %10.3E", stepNorm_1, m_normResid_1, m_normResid_full);
	} else {
	  printf(" |");
	  printf("%10.2E %10.3E |", m_dampBound, m_normResid_Bound);
	  printf("%2d %10.2E %10.3E  %10.3E %10.3E", i_numTrials + 1, m_dampRes,
		 stepNorm_1 / ( m_dampRes * m_dampBound),  stepNorm_2,  m_normResid_1);
	  printf("| %10.3E %10.3E %10.3E", stepNorm_1, m_normDeltaSoln_Newton, m_normResid_full);
	}
	printf("\n");
      
      }
      if (m_print_flag >= 4) {
	if (doDogLeg_) {
	  if (convRes > 0) {
	    printf("\t   solve_nonlinear_problem(): Problem Converged, stepNorm = %11.3E, reduction of res from %11.3E to %11.3E\n",
		   stepNorm_1, m_normResid_0, m_normResid_full);
	    printf("\t");
	    print_line("=", 119);
	  } else {
	    printf("\t   solve_nonlinear_problem(): Successfull step taken with stepNorm = %11.3E, reduction of res from %11.3E to %11.3E\n",
		   stepNorm_1, m_normResid_0, m_normResid_full);
	  }
	} else {
	  if (convRes > 0) {
	    printf("\t   solve_nonlinear_problem(): Damped Newton iteration successful, nonlin "
		   "converged, final estimate of the next solution update norm = %-12.4E\n", stepNorm_2);
	    printf("\t");
	    print_line("=", 119);
	  } else if (retnDamp >= NSOLN_RETN_CONTINUE) {
	    printf("\t   solve_nonlinear_problem(): Damped Newton iteration successful, "
		   "estimate of the next solution update norm = %-12.4E\n", stepNorm_2);
	  } else {
	    printf("\t   solve_nonlinear_problem(): Damped Newton unsuccessful, final estimate "
		   "of the next solution update norm = %-12.4E\n", stepNorm_2);
	  }
	}
      }
      // convergence
      if (convRes) {
	goto done;
      }

      // If dampStep fails, first try a new Jacobian if an old
      // one was being used. If it was a new Jacobian, then
      // return -1 to signify failure.
      else if (retnDamp < NSOLN_RETN_CONTINUE) {
	goto done;
      }
    }

  done:


    if (m_print_flag == 2 || m_print_flag == 3) {
      if (convRes > 0) {
	if (doDogLeg_) {
	  if (convRes == 3) {
	    printf("\t                            |                                           |                      "
		   "          |                        | converged = 3  |(%11.3E) \n", stepNorm_2);
	  } else {
	    printf("\t                            |                                           |                      "
		   "          |                        | converged = %1d  | %10.3E %10.3E\n", convRes,
		   stepNorm_2, m_normResidTrial);
	  }
	  printf("\t-----------------------------------------------------------------------------------------------------"
		 "------------------------------------------------------------------------------\n");
	} else {
	  if (convRes == 3) {
	    printf("\t                       |                 "
		   "     |                                converged = 3  |            (%11.3E) \n", stepNorm_2);
	  } else {
	    printf("\t                       |                 "
		   "     |                                converged = %1d  |            %10.3E %10.3E\n", convRes,
		   stepNorm_2, m_normResidTrial);
	  }
	  printf("\t------------------------------------------------------------------------------------"
		 "-----------------------------------------------\n");
	}
      }
   
    
      
    }
    
    mdp::mdp_copy_dbl_1(y_comm, CONSTD_DATA_PTR(m_y_n_curr), neq_);
    if (solnType_ != NSOLN_TYPE_STEADY_STATE) {
      mdp::mdp_copy_dbl_1(ydot_comm, CONSTD_DATA_PTR(m_ydot_n_curr), neq_);
    }
 
    num_linear_solves += m_numTotalLinearSolves;
 
    doublereal time_elapsed =  wc.secondsWC();
    if (m_print_flag > 1 ) {
      if (retnDamp > 0) {
	if (NonlinearSolver::s_TurnOffTiming) {
	  printf("\tNonlinear problem solved successfully in %d its\n",
		 num_newt_its);
	} else {
	  printf("\tNonlinear problem solved successfully in %d its, time elapsed = %g sec\n",
		 num_newt_its, time_elapsed);
	}
      } else {
	printf("\tNonlinear problem failed to solve after %d its\n", num_newt_its);
      }
    }
    retnCode = retnDamp;
    if (retnDamp > 0) {
      retnCode = NSOLN_RETN_SUCCESS;
    }


    return retnCode;
  }
  //====================================================================================================================
  // Print solution norm contribution
  /*
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
  void NonlinearSolver::
  print_solnDelta_norm_contrib(const doublereal * const step_1,
			       const char * const stepNorm_1,
			       const doublereal * const step_2,
			       const char * const stepNorm_2,
			       const char * const title,
			       const doublereal * const y_n_curr,
			       const doublereal * const y_n_1,
			       doublereal damp,
			       int num_entries) {
    int i, j, jnum;
    bool used;
    doublereal dmax0, dmax1, error, rel_norm;
    printf("\t\t%s currentDamp = %g\n", title, damp);
    printf("\t\t     I     ysolnOld %13s ysolnNewRaw | ysolnNewTrial "
	   "%10s ysolnNewTrialRaw | solnWeight  wtDelSoln wtDelSolnTrial\n", stepNorm_1, stepNorm_2);
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
	  error     = step_1[i] /  m_ewt[i];
	  rel_norm = sqrt(error * error);
	  error     = step_2[i] /  m_ewt[i];
	  rel_norm += sqrt(error * error);
	  if (rel_norm > dmax1) {
	    imax[jnum] = i;
	    dmax1 = rel_norm;
	  }
	}
      }
      if (imax[jnum] >= 0) {
	i = imax[jnum];
	error = step_1[i] /  m_ewt[i];
	dmax0 = sqrt(error * error);
	error = step_2[i] /  m_ewt[i];
	dmax1 = sqrt(error * error);
	printf("\t\t  %4d %12.4e %12.4e %12.4e |  %12.4e  %12.4e   %12.4e    |%12.4e %12.4e %12.4e\n",
	       i, y_n_curr[i], step_1[i],  y_n_curr[i] + step_1[i], y_n_1[i],
	       step_2[i], y_n_1[i]+ step_2[i], m_ewt[i], dmax0, dmax1);
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
  int NonlinearSolver::beuler_jac(GeneralMatrix &J, doublereal * const f,
				  doublereal time_curr, doublereal CJ,
				  doublereal * const y,  doublereal * const ydot,
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
      info = m_func->evalJacobian(time_curr, delta_t_n, CJ, y, ydot, J, f);
      m_nJacEval++;
      m_nfe++;
      if (info != 1) {
	return info;
      }
    }  else {
      if (J.matrixType_ == 0) {
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
#ifdef DEBUG_MODE
	for (int ii = 0; ii < neq_; ii++) {
	  mdp::checkFinite(f[ii]);
	}
#endif

	/*
	 * Malloc a vector and call the function object to return a set of
	 * deltaY's that are appropriate for calculating the numerical
	 * derivative.
	 */
	doublereal *dyVector = mdp::mdp_alloc_dbl_1(neq_, MDP_DBL_NOINIT);
	retn = m_func->calcDeltaSolnVariables(time_curr, y, ydot, dyVector, DATA_PTR(m_ewt));
	if (s_print_NumJac) {
	  if (m_print_flag >= 7) {
	    if (retn != 1) {
	      printf("\t\t  beuler_jac ERROR! calcDeltaSolnVariables() returned an error flag\n");
	      printf("\t\t                    We will bail from the nonlinear solver after calculating the jacobian");
	    }
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
	   
#ifdef DEBUG_MODE
	  if (fabs(dy) < 1.0E-300) {
	    throw CanteraError("NonlinearSolver::beuler_jac", "dy is equal to zero");
	  }
	  for (int ii = 0; ii < neq_; ii++) {
	    mdp::checkFinite(m_wksp[ii]);
	  }
#endif
      
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
      } else if (J.matrixType_ == 1) {
	int ku, kl;
	int ivec[2];
	int n = J.nRowsAndStruct(ivec);
	kl = ivec[0];
	ku = ivec[1];
	if (n != neq_) {
	  printf("we have probs\n"); exit(-1);
	}

	// --------------------------------- BANDED MATRIX BRAIN DEAD ---------------------------------------------------
	info = m_func->evalResidNJ(time_curr, delta_t_n, y, ydot, f, JacBase_ResidEval);
	m_nfe++;
	if (info != 1) {
	  return info;
	}
	m_nJacEval++;


	doublereal *dyVector = mdp::mdp_alloc_dbl_1(neq_, MDP_DBL_NOINIT);
	retn = m_func->calcDeltaSolnVariables(time_curr, y, ydot, dyVector, DATA_PTR(m_ewt));
      	if (s_print_NumJac) {
	  if (m_print_flag >= 7) {
	    if (retn != 1) {
	      printf("\t\t  beuler_jac ERROR! calcDeltaSolnVariables() returned an error flag\n");
	      printf("\t\t                    We will bail from the nonlinear solver after calculating the jacobian");
	    }
	    if (neq_ < 20) {
	      printf("\t\tUnk            m_ewt              y                dyVector            ResN\n");
	      for (int iii = 0; iii < neq_; iii++){
		printf("\t\t %4d       %16.8e   %16.8e   %16.8e  %16.8e \n",
		       iii,   m_ewt[iii],  y[iii], dyVector[iii], f[iii]);
	      }
	    }
	  }
	}


	for (j = 0; j < neq_; j++) {


	  col_j = (doublereal *) J.ptrColumn(j);
	  ysave = y[j];
	  dy = dyVector[j];
	 

	  y[j] = ysave + dy;
	  dy = y[j] - ysave;
	  if (solnType_ != NSOLN_TYPE_STEADY_STATE) {
	    ydotsave = ydot[j];
	    ydot[j] += dy * CJ;
	  }

	  info =  m_func->evalResidNJ(time_curr, delta_t_n, y, ydot, DATA_PTR(m_wksp), JacDelta_ResidEval, j, dy);
	  m_nfe++;
#ifdef DEBUG_MODE
	  if (fabs(dy) < 1.0E-300) {
	    throw CanteraError("NonlinearSolver::beuler_jac", "dy is equal to zero");
	  }
	  for (int ii = 0; ii < neq_; ii++) {
	    mdp::checkFinite(m_wksp[ii]);
	  }
#endif
	  if (info != 1) {
	    mdp::mdp_safe_free((void **) &dyVector);
	    return info;
	  }

	  doublereal diff;



	  for (int i = j - ku; i <= j + kl; i++) {
	    if (i >= 0 &&  i < neq_) {
	      diff = subtractRD(m_wksp[i], f[i]);
	      col_j[kl + ku + i - j] = diff / dy;
	    }
	  }
	  y[j] = ysave;
	  if (solnType_ != NSOLN_TYPE_STEADY_STATE) {
	    ydot[j] = ydotsave;
	  }

	}

	mdp::mdp_safe_free((void **) &dyVector);
	double vSmall;
	int ismall = J.checkRows(vSmall);
	if (vSmall < 1.0E-100) {
	  printf("WE have a zero row, %d\n", ismall);
	  exit(-1);
	}
        ismall = J.checkColumns(vSmall);
	if (vSmall < 1.0E-100) {
	  printf("WE have a zero column, %d\n", ismall);
	  exit(-1);
	}

	// ---------------------BANDED MATRIX BRAIN DEAD -----------------------
      }
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
    jacCopyPtr_->copyData(J);

    return retn;
  }
  //====================================================================================================================
  //   Internal function to calculate the time derivative of the solution at the new step
  /*
   *  Previously, the user must have supplied information about the previous time step for this routine to
   *  work as intended.
   *
   *   @param  order of the BDF method
   *   @param   y_curr current value of the solution
   *   @param   ydot_curr  Calculated value of the solution derivative that is consistent with y_curr
   */ 
  void NonlinearSolver::
  calc_ydot(const int order, const doublereal * const y_curr, doublereal * const ydot_curr) const
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
  
      return;
    default:
      throw CanteraError("calc_ydot()", "Case not covered");
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
   *  is added as well so that if the residual is less than 1, then the calculation is deemed to be converged.
   *
   *  The basic idea is that a change in the solution vector on the order of the convergence tolerance
   *  multiplied by  [RJC] which is of order one after row scaling should give you the relative weight
   *  of the row. Values of the residual for that row can then be normalized by the value of this weight.
   *  When the tolerance in delta x is achieved, the tolerance in the residual should also be achieved
   *  and should be checked.
   */
  void
  NonlinearSolver::computeResidWts()
  {
    ResidWtsReevaluated_ = true;
    if (checkUserResidualTols_ == 1) {
      for (int i = 0; i < neq_; i++) {
	m_residWts[i] = userResidAtol_[i] + userResidRtol_ * m_rowWtScales[i] / neq_;
#ifdef DEBUG_MODE
        mdp::checkFinite(m_residWts[i]);
#endif
      }
    } else {
      doublereal sum = 0.0;  
      for (int i = 0; i < neq_; i++) {
	m_residWts[i] = m_rowWtScales[i] / neq_;
#ifdef DEBUG_MODE
        mdp::checkFinite(m_residWts[i]);
#endif
	sum += m_residWts[i];
      }
      sum /= neq_;
      for (int i = 0; i < neq_; i++) {
	m_residWts[i] = m_ScaleSolnNormToResNorm * (m_residWts[i] + atolBase_ * atolBase_ * sum);
      }
      if (checkUserResidualTols_ == 2) {
	for (int i = 0; i < neq_; i++) {
	  double uR = userResidAtol_[i] + userResidRtol_ * m_rowWtScales[i] / neq_;
	  m_residWts[i] = MIN(m_residWts[i], uR);
	}
      }
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
      if (atol[i] <= 0.0) {
	throw CanteraError("NonlinearSolver::setAtol()",
			   "Atol is less than or equal to zero");
      }
      atolk_[i]= atol[i];
    }
  }
  //=====================================================================================================================
  // Set the relative tolerances for the solution variables
  /*
   *   Set the relative tolerances used in the calculation for the solution variables.
   *
   *  @param rtol  single double
   */
  void NonlinearSolver::setRtol(const doublereal rtol)
  {
    if (rtol <= 0.0) {
      throw CanteraError( "NonlinearSolver::setRtol()",
			  "Rtol is <= zero");
    }
    rtol_ = rtol;
  }
  //=====================================================================================================================
  // Set the relative and absolute tolerances for the Residual norm comparisons, if used
  /*
   * 
   *   residWeightNorm[i] = residAtol[i] +  residRtol * m_rowWtScales[i] / neq 
   *
   *     @param residNormHandling  Parameter that sets the default handling of the residual norms
   *                      0   The residual weighting vector is calculated to make sure that the solution
   *                          norms are roughly 1 when the residual norm is roughly 1.
   *                          This is the default if this routine is not called.
   *                      1   Use the user residual norm specified by the parameters in this routine
   *                      2   Use the minimum value of the residual weights calculcated by method 0 and 1.
   *                          This is the default if this routine is called and this parameter isn't specified.
   */
  void  NonlinearSolver::setResidualTols(double residRtol, double * residATol, int residNormHandling)
  {
    if (residNormHandling < 0 ||  residNormHandling > 2) {
      throw CanteraError("NonlinearSolver::setResidualTols()", 
			 "Unknown int for residNormHandling");
    }
    checkUserResidualTols_ = residNormHandling;
    userResidRtol_ = residRtol;
    if (residATol) {
      userResidAtol_.resize(neq_);
      for (int i = 0; i < neq_; i++) {
	userResidAtol_[i] = residATol[i];
      }
    } else {
      if (residNormHandling ==1  ||  residNormHandling == 2) {
	throw CanteraError("NonlinearSolver::setResidualTols()", 
			   "Must set residATol vector");
      }
    }
  }
  //=====================================================================================================================
  void NonlinearSolver::setPrintLvl(int printLvl) 
  {
    m_print_flag = printLvl;
  }
  //=====================================================================================================================
}

