/**
 *  @file ResidJacEval.cpp
 *
 */
/*
 * $Revision$
 * $Date$
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "ct_defs.h"
#include "ctlapack.h"
#include "ResidJacEval.h"

#include <iostream>
#include <vector>

using namespace std;

namespace Cantera {
  //====================================================================================================================
  
  ResidJacEval::ResidJacEval(doublereal atol) :
    ResidEval(),
    m_atol(atol)
  {
  }
  //====================================================================================================================
  //   Copy Constructor for the %ResidJacEval object
  /*
   */
  ResidJacEval::ResidJacEval(const ResidJacEval &right) :
    ResidEval()
  {
    *this = operator=(right);
  }
  //====================================================================================================================
  ResidJacEval::~ResidJacEval() 
  {
  }
  //====================================================================================================================
  ResidJacEval& ResidJacEval::operator=(const ResidJacEval &right) {
    if (this == &right) {
      return *this;
    }

    ResidEval::operator=(right);

    m_atol = right.m_atol;
    neq_   = right.neq_;

    return *this;
  }
  //====================================================================================================================
  // Duplication routine for objects which inherit from %ResidJacEval
  /*
   *  This virtual routine can be used to duplicate %ResidJacEval objects
   *  inherited from %ResidJacEval even if the application only has
   *  a pointer to %ResidJacEval to work with.
   *
   *  These routines are basically wrappers around the derived copy
   *  constructor.
   */
  ResidJacEval *ResidJacEval::duplMyselfAsResidJacEval() const {
    ResidJacEval *ff = new ResidJacEval(*this);
    return ff;
  }
  //====================================================================================================================
  int ResidJacEval::nEquations() const {
    return neq_;
  }
  //====================================================================================================================
  // Set a global value of the absolute tolerance
  /*
   *  @param atol   Value of atol
   */
  void ResidJacEval::setAtol(doublereal atol)
  {
    m_atol = atol;
    if (m_atol <= 0.0) {
      throw CanteraError("ResidJacEval::setAtol",
			 "atol must be greater than zero");
    }
  }
  //====================================================================================================================
  //! Fill in the initial conditions
  /*!
   * Values for both the solution and the value of ydot may be provided.
   *
   * @param t0            Time                    (input) 
   * @param y             Solution vector (output)
   * @param ydot          Rate of change of solution vector. (output)
   */
  void ResidJacEval::
  getInitialConditions(doublereal t0, doublereal * const y, doublereal * const ydot) {
    for (int i = 0; i < neq_; i++) {
      y[i] = 0.0;
    }
    if (ydot) {
      for (int i = 0; i < neq_; i++) {
	ydot[i] = 0.0;
      }
    }
  }
  //====================================================================================================================
  // This function may be used to create output at various points in the execution of an application.
  /*
   *
   *  @param ifunc     identity of the call
   *                          0  Initial call
   *                          1  Called at the end of every successful time step
   *                         -1  Called at the end of every unsuccessful time step
   *                          2  Called at the end of every call to integrateRJE()
   *
   * @param t             Time                    (input) 
   * @param delta_t       The current value of the time step (input)
   * @param y             Solution vector (input, do not modify)
   * @param ydot          Rate of change of solution vector. (input)
   */
  void ResidJacEval::
  user_out2(const int ifunc, const doublereal t, const doublereal deltaT,
	    const doublereal *y, const doublereal *ydot) {
  }

  //====================================================================================================================
  // This function may be used to create output at various points in the execution of an application.
  /*
   *  This routine calls user_out2().
   *
   *  @param ifunc     identity of the call
   * @param t             Time                    (input) 
   * @param y             Solution vector (input, do not modify)
   * @param ydot          Rate of change of solution vector. (input)
   */
  void ResidJacEval::
  user_out(const int ifunc, const doublereal t,
	   const doublereal *y, const doublereal *ydot) {
    user_out2(ifunc, t, 0.0, y, ydot);
  }
 //====================================================================================================================
  //! Evaluate the time tracking equations, if any
  /*!
   * Evaluate time integrated quantities that are calculated at the
   * end of every successful time step.  This call is made once at the end of every successful
   * time step that advances the time. It's also made once at the start of the time stepping.
   *
   * @param t             Time                    (input) 
   * @param delta_t       The current value of the time step (input)
   * @param y             Solution vector (input, do not modify)
   * @param ydot          Rate of change of solution vector. (input, do not modify)
   */
  void ResidJacEval::
  evalTimeTrackingEqns(const doublereal t, const doublereal delta_t, const doublereal *y, 
		       const doublereal *ydot)
  {
  }
  //====================================================================================================================
  //  Return a vector of delta y's for calculation of the numerical Jacobian 
  /*
   *   There is a default algorithm provided. 
   *
   *        delta_y[i] = atol[i] + 1.0E-6 ysoln[i]
   *        delta_y[i] = atol[i] + MAX(1.0E-6 ysoln[i] * 0.01 * solnWeights[i])
   *
   * @param t             Time                    (input) 
   * @param y             Solution vector (input, do not modify)
   * @param ydot          Rate of change of solution vector. (input, do not modify)
   * @param delta_y       Value of the delta to be used in calculating the numerical jacobian
   * @param solnWeights   Value of the solution weights that are used in determining convergence (default = 0)
   *
   * @return Returns a flag to indicate that operation is successful.
   *            1  Means a successful operation
   *            0  Means an unsuccessful operation
   */
  int ResidJacEval::
  calcDeltaSolnVariables(const doublereal t, const doublereal * const ySoln,
			 const doublereal * const ySolnDot, doublereal * const deltaYSoln,
			 const doublereal *const solnWeights)
  {
    if (!solnWeights) {
      for (int i = 0; i < neq_; i++) {
	deltaYSoln[i] = m_atol + fabs(1.0E-6 * ySoln[i]);
      }
    } else {
      for (int i = 0; i < neq_; i++) {
	deltaYSoln[i] = fmaxx(1.0E-2 * solnWeights[i], 1.0E-6 * fabs(ySoln[i]));
      }
    }
    return 1;
  }
  //====================================================================================================================
  //  Returns a vector of column scale factors that can be used to column scale Jacobians.
  /*
   *  Default to yScales[] = 1.0
   *
   * @param t             Time                    (input) 
   * @param y             Solution vector (input, do not modify)
   * @param y_old         Old Solution vector (input, do not modify)
   * @param yScales       Value of the column scales
   */
  void ResidJacEval::
  calcSolnScales(const doublereal t, const doublereal * const ysoln, const doublereal * const ysolnOld,
		 doublereal * const ysolnScales)
  {
    if (ysolnScales[0] == 0.0) {
      for (int i = 0; i < neq_; i++) {
	ysolnScales[i] = 1.0;
      }
    }
  }
  //====================================================================================================================
  // Filter the solution predictions
  /*
   * Codes might provide a predicted step change. This routine filters the predicted
   * solution vector eliminating illegal directions.
   *
   * @param t             Time                    (input) 
   * @param y             Solution vector (input, output)
   * @param step          Proposed step in the solution that will be cropped
   */
  doublereal ResidJacEval::filterNewStep(doublereal t, const doublereal * const ybase, doublereal * const step)
  {
    return 0.0;
  }
  //====================================================================================================================
  // Filter the solution predictions
  /*
   * Codes might provide a predicted solution vector. This routine filters the predicted
   * solution vector.
   *
   * @param t             Time                    (input) 
   * @param y             Solution vector (input, output)
   */
  doublereal ResidJacEval::filterSolnPrediction(doublereal t, doublereal * const y)
  {
    return 0.0;
  }
  //====================================================================================================================
  // Evalulate any stopping criteria other than a final time limit
  /*
   *  If we are to stop the time integration for any reason other than reaching a final time limit, tout,
   *  provide a test here. This call is made at the end of every succesful time step iteration
   *
   *  @return    If true, the the time stepping is stopped. If false, then time stepping is stopped if t >= tout
   *             Defaults to false.
   *
   * @param t             Time                    (input) 
   * @param delta_t       The current value of the time step (input)
   * @param y             Solution vector (input, do not modify)
   * @param ydot          Rate of change of solution vector. (input, do not modify)
   */
  bool ResidJacEval::
  evalStoppingCritera(const doublereal t,
		      const doublereal delta_t,
		      const doublereal * const y, 
		      const doublereal * const ydot)
  {
    return false;
  }
  //====================================================================================================================
  // Multiply the matrix by another matrix that leads to better conditioning
  /*
   *  Provide a left sided matrix that will multiply the current jacobian, after scaling 
   *  and lead to a better conditioned system.
   *  This routine is called just before the matrix is factored.
   * 
   *  Original Problem:
   *        J delta_x = - Resid
   *
   *  New problem:
   *      M (J delta_x) = - M Resid
   *
   *  @param    matrix     Pointer to the current jacobian (if zero, it's already been factored)
   *  @param    nrows      offsets for the matrix
   *  @param    rhs        residual vector. This also needs to be lhs multiplied by M
   */
  void ResidJacEval::
  matrixConditioning(doublereal * const matrix, const int nrows, doublereal * const rhs)
  {
  }
  //====================================================================================================================
  // Evaluate the residual function
  /*
   * @param t             Time                    (input) 
   * @param delta_t       The current value of the time step (input)
   * @param y             Solution vector (input, do not modify)
   * @param ydot          Rate of change of solution vector. (input, do not modify)
   * @param resid         Value of the residual that is computed (output)
   * @param evalType      Type of the residual being computed (defaults to Base_ResidEval)
   * @param id_x          Index of the variable that is being numerically differenced to find
   *                      the jacobian (defaults to -1, which indicates that no variable is being
   *                      differenced or that the residual doesn't take this issue into account)
   * @param delta_x       Value of the delta used in the numerical differencing
   */ 
  void ResidJacEval::
  evalResidNJ(const doublereal t, const doublereal deltaT,  const doublereal * y,
	      const doublereal * ydot, doublereal * const resid, const ResidEval_Type_Enum evalType,
	      const int id_x, const doublereal delta_x) 
  {  
    throw CanteraError("ResidJacEval::evalResidNJ()", "Not implemented\n");
  }
  //====================================================================================================================
  // Calculate an analytical jacobian and the residual at the current time and values.
  /*
   *  Only called if the jacFormation method is set to analytical
   *
   * @param t             Time                    (input) 
   * @param delta_t       The current value of the time step (input)
   * @param y             Solution vector (input, do not modify)
   * @param ydot          Rate of change of solution vector. (input, do not modify)
   * @param J             Reference to the SquareMatrix object to be calculated (output)
   * @param resid         Value of the residual that is computed (output)
   */
  void ResidJacEval::
  evalJacobian(const doublereal t, const doublereal delta_t,
	       const doublereal * const y,
	       const doublereal * const ydot,
	       SquareMatrix &J,
	       doublereal * const resid)
  {
    throw CanteraError("ResidJacEval::evalJacobian()", "Not implemented\n");
  }
  //====================================================================================================================

}
