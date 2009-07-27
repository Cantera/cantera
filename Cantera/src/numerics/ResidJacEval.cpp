/**
 *  @file ResidJacEval.cpp
 *
 */
/*
 * $Revision: 1.2 $
 * $Date: 2009/03/03 17:55:25 $
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

  /*************************************************************************
   *
   *  ResidJacEval():
   *
   *  Default constructor for the ResidJacEval class. 
   *
   *  atol has a default of 1.0E-13.
   */
  ResidJacEval::ResidJacEval(doublereal atol) :
    ResidEval(),
    m_atol(atol)
  {
  }

  //   Copy Constructor for the %ResidJacEval object
  /*
   */
  ResidJacEval::ResidJacEval(const ResidJacEval &right) :
    ResidEval()
  {
    *this = operator=(right);
  }

  /*
   *
   */
  ResidJacEval::~ResidJacEval() 
  {
  }

  ResidJacEval& ResidJacEval::operator=(const ResidJacEval &right) {
    if (this == &right) {
      return *this;
    }

    ResidEval::operator=(right);

    m_atol = right.m_atol;
    neq_   = right.neq_;

    return *this;
  }

  // Duplication routine for objects which inherit from
  // %ResidJacEval
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

  int ResidJacEval::nEquations() const {
    return neq_;
  }

  /*
   *
   * setAtol():
   *
   * Set the absolute tolerance value
   */
  void ResidJacEval::setAtol(doublereal atol)
  {
    m_atol = atol;
    if (m_atol <= 0.0) {
      throw CanteraError("ResidJacEval::setAtol",
			 "atol must be greater than zero");
    }
  }

  /**************************************************************************
   *
   *
   *
   *       Fill the solution vector with the initial conditions
   *       at initial time t0.
   */
  void ResidJacEval::
  getInitialConditionsDot(const doublereal t0, const size_t leny, 
			  doublereal * const y, doublereal * const ydot) {
    for (int i = 0; i < neq_; i++) {
      y[i] = 0.0;
    }
    if (ydot) {
      for (int i = 0; i < neq_; i++) {
	ydot[i] = 0.0;
      }
    }
  }

  /**************************************************************************
   *
   *
   *
   * Fill the solution vector with the initial conditions
   *       at initial time t0.
   *
   */
  void ResidJacEval::
  getInitialConditions(doublereal t0, 
		       doublereal * const y, doublereal * const ydot) {
    size_t leny = neq_;
    getInitialConditionsDot(t0, leny, y, 0);
  }

  /**************************************************************************
   *
   * user_out():
   *
   * This function may be used to create output at various points in the
   * execution of an application.
   *
   */
  void ResidJacEval::
  user_out2(const int ifunc, const doublereal t, const doublereal deltaT,
	    const doublereal *y, const doublereal *ydot) {

  }

  void ResidJacEval::
  user_out(const int ifunc, const doublereal t,
	   const doublereal *y, const doublereal *ydot) {
    user_out2(ifunc, t, 0.0, y, ydot);
  }

  /**************************************************************************
   *
   *
   */
  void ResidJacEval::
  evalTimeTrackingEqns(const doublereal t, const doublereal deltaT,
		       const doublereal *y, 
		       const doublereal *ydot) {

  }

  /********************************************************************
   *
   *
   *
   * Return a vector of delta y's for calculation of the
   * numerical Jacobian
   */
  void ResidJacEval::
  calcDeltaSolnVariables(const doublereal t,
			 const doublereal * const ySoln,
			 const doublereal * const ySolnDot,
			 doublereal * const deltaYSoln,
			 const doublereal *const solnWeights)
  {
    if (!solnWeights) {
      for (int i = 0; i < neq_; i++) {
	deltaYSoln[i] = m_atol + fabs(1.0E-6 * ySoln[i]);
      }
    } else {
      for (int i = 0; i < neq_; i++) {
	deltaYSoln[i] = m_atol +
	  fmaxx(1.0E-2 * solnWeights[i], 1.0E-6 * fabs(ySoln[i]));
      }
    }
  }

  /******************************************************************
   *
   * calcSolnScales():
   *
   *  Returns a vector of ysolnScales[] that can be used to column scale
   *  Jacobians.
   */
  void ResidJacEval::
  calcSolnScales(const doublereal t,
		 const doublereal * const ysoln,
		 const doublereal * const ysolnOld,
		 doublereal * const ysolnScales)
  {
    for (int i = 0; i < neq_; i++) {
      ysolnScales[i] = 1.0;
    }
  }

  void ResidJacEval::filterSolnPrediction(doublereal t,
					  doublereal * const y) {

  }

  /**************************************************************************
   *
   *  evalStoppingCriteria()
   *
   * If there is a stopping critera other than time set it here.
   *
   */
  bool ResidJacEval::
  evalStoppingCritera(doublereal &time_current,
		      doublereal &delta_t_n,
		      doublereal *y_n, 
		      doublereal *ydot_n)
  {
    return false;
  }

  /**************************************************************************
   *
   * matrixConditioning()
   *
   * Multiply the matrix by the inverse of a matrix which lead to a 
   * better conditioned system. The default, specified here, is to 
   * do nothing.
   */
  void ResidJacEval::
  matrixConditioning(doublereal * const matrix, const int nrows,
		     doublereal * const rhs)
  {
  }

  /**************************************************************************
   *
   */
  void ResidJacEval::
  evalResidNJ(doublereal t, const doublereal deltaT, 
	      const doublereal * y,
	      const doublereal * ydot,
	      doublereal * resid,
	      bool NJevaluation,
	      int id_x, 
	      doublereal delta_x) 
  {
    printf("Not implemented\n");
    std::exit(-1);
  }

  /**************************************************************************
   *
   * evalJacobian() 
   *
   *  Calculate the jacobian and the residual at the current
   *  time and values.
   *  Backwards Euler is assumed.
   */
  void ResidJacEval::
  evalJacobian(const doublereal t, const doublereal deltaT,
	       const doublereal * const y,
	       const doublereal * const ydot,
	       SquareMatrix &J,
	       doublereal * const resid)
  {
    printf("Not implemented\n");
    std::exit(-1);
  }


}

