/**
 *  @file ResidJacEval.h
 *
 * Dense, Square (not sparse) matrices.
 */

/*
 *  $Date: 2009/04/02 14:24:05 $
 *  $Revision: 1.3 $
 *
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */


#ifndef CT_RESIDJACEVAL_H
#define CT_RESIDJACEVAL_H

#include "ResidEval.h"
#include "SquareMatrix.h"

namespace Cantera { 

  /**
   *  A class for full (non-sparse) matrices with Fortran-compatible
   *  data storage. Adds matrix operations to class Array2D.
   */
  class ResidJacEval : public ResidEval {

  public:

    /**
     * Default constructor
     */
    ResidJacEval(doublereal atol = 1.0e-13);

    //!Copy Constructor for the %ResidJacEval object
    /*!
     * @param right   Item to be copied
     */
    ResidJacEval(const ResidJacEval &right);
    
    /// Destructor. Does nothing.
    virtual ~ResidJacEval();

    //! Assignment operator
    /*!
     *  This is NOT a virtual function.
     *
     * @param right    Reference to %ResidJacEval object to be copied into the
     *                 current one.
     */
    ResidJacEval& operator=(const ResidJacEval &right);

    //! Duplication routine for objects which inherit from
    //! residJacEval
     /*!
      *  This virtual routine can be used to duplicate %ResidJacEval objects
      *  inherited from %ResidJacEval even if the application only has
      *  a pointer to %ResidJacEval to work with.
      *
      *  These routines are basically wrappers around the derived copy
      *  constructor.
      */
    virtual ResidJacEval *duplMyselfAsResidJacEval() const;

    //! Return the number of equations in the equation system
    virtual int nEquations() const;

    /**
     * Evaluate the residual function.  
     * @param t time (input, do not modify) 
     * @param y solution vector (input, do not modify)
     * @param ydot rate of change of solution vector. (input, do
     * not modify)
     */
    virtual void evalResidNJ(doublereal t, const doublereal deltaT, 
			     const doublereal * const y,
			     const doublereal * const ydot,
			     doublereal * const resid,
			     bool NJevaluation = false,
			     int id_x = 0, 
			     doublereal delta_x = 0.0);

    /**
     *       Fill the solution vector with the initial conditions
     *       at initial time t0.
     */
    virtual void getInitialConditionsDot(const doublereal t0, size_t leny, 
					 doublereal * const y,
					 doublereal * const ydot);

    virtual void getInitialConditions(const doublereal t0, 
				      doublereal * const y, 
                                      doublereal * const ydot);
  
    virtual void filterSolnPrediction(doublereal t,
				      doublereal * const y);

    void setAtol(doublereal atol);

    virtual void evalTimeTrackingEqns(const doublereal t, const doublereal deltaT,
				      const doublereal * const y, 
				      const doublereal * const ydot);

    virtual bool evalStoppingCritera(doublereal &time_current,
				     doublereal &delta_t_n,
				     doublereal *y_n, 
				     doublereal *ydot_n);
    /**
     * Return a vector of delta y's for calculation of the
     *  numerical Jacobian 
     */
    virtual void 
    calcDeltaSolnVariables(const doublereal t, 
			   const doublereal * const ysoln,
			   const doublereal * const ysolnDot,
			   doublereal * const deltaYsoln,
			   const doublereal * const solnWeights=0);

    /**
     *  Returns a vector of ysolnScales[] that can be used to column
     *  scale Jacobians.
     */
    virtual void calcSolnScales(const doublereal t,
				const doublereal * const ysoln,
				const doublereal * const ysolnOld,
				doublereal * const ysolnScales);

    /**
     * This function may be used to create output at various points in the
     * execution of an application.
     *
     */
    virtual void user_out2(const int ifunc, const doublereal t, 
			   const doublereal deltaT,
			   const doublereal * const y,
			   const doublereal * const ydot);

    virtual void user_out(const int ifunc, const doublereal t, 
			  const doublereal *y,
			  const doublereal *ydot);

	
    virtual void matrixConditioning(doublereal * const matrix, const int nrows,
				    doublereal * const rhs);

    /*********************************************************************
     *
     * evalJacobian() 
     *
     *  Calculate the jacobian and the residual at the current
     *  time and values.
     *  Backwards Euler is assumed.
     */
    virtual void evalJacobian(const doublereal t, const doublereal deltaT,
                              
			      const double* const y,
			      const double* const ydot,
			      SquareMatrix &J,
			      doublereal * const resid);



  protected:

    doublereal m_atol;

    int neq_;

  };
}

#endif



