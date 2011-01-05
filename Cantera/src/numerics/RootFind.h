/**
 * @file RootFind.h
 *       Header file for implicit nonlinear solver of a one dimensional function
 *  (see \ref numerics and class \link Cantera::RootFind RootFind\endlink).
 */
/*
 * $Id$
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#ifndef CT_ROOTFIND_H
#define CT_ROOTFIND_H
/**
 * @defgroup solverGroup Solvers for Equation Systems
 */

#include <vector>
#include "ResidEval.h"

namespace Cantera {

#define ROOTFIND_SUCCESS             0
#define ROOTFIND_FAILEDCONVERGENCE  -1
#define ROOTFIND_BADINPUT           -2

  //! Root finder for 1D problems
  /*!
   *
   *
   *
   */
  class RootFind {

  public:

    //! Constructor for the object
    /*!
     */ 
    RootFind(ResidEval* resid);
   
    //! Destructor. Deletes the integrator.
    ~RootFind();

  private:

    //! Unimplemented private copy constructor
    RootFind(const RootFind &right);

    //! Unimplemented private assignment operator
    RootFind& operator=(const RootFind &right);

  
    double delXNonzero(double x1) const;
  
    double delXMeaningful(double x1) const;
 
    double deltaXControlled(double x2, double x1) const;
   
    bool theSame(double x2, double x1) const;

  public:

    //!  Using a line search method, find the root of a 1D function
    /*!
     *  This routine solves the following equation.
     * 
     *    \f[
     *       R(x) = f(x) - f_o = 0
     *    \f]
     *
     *    @param   xmin    Minimum value of x to be used.
     *    @param   xmax    Maximum value of x to be used
     *    @param   itmax   maximum number of iterations. Usually, it can be less than 50.
     *    @param   funcTargetValue   
     *                     Value of \f$ f_o \f$ in the equation.
     *                     On return, it contains the value of the function actually obtained.
     *    @param   xbest   Returns the x that satisfies the function
     *                     On input, xbest should contain the best estimate of the solution.
     *                     An attempt to find the solution near xbest is made.
     *
     *   @return:
     *    0  =  ROOTFIND_SUCCESS            Found function     
     *   -1  =  ROOTFIND_FAILEDCONVERGENCE  Failed to find the answer
     *   -2  =  ROOTFIND_BADINPUT           Bad input was detected
     */
    int solve(doublereal xmin, doublereal xmax, int itmax, doublereal &funcTargetValue, doublereal *xbest);


    //! Return the function value
    /*! 
     * This routine evaluates the following equation.
     * 
     *    \f[
     *       R(x) = f(x) - f_o = 0
     *    \f]
     *
     *  @param x  Value of the independent variable
     *
     *  @return   The routine returns the value of \f$ R(x) \f$
     */
    doublereal func(doublereal x);

    void setTol(doublereal rtol, doublereal atol);

    void setPrintLvl(int printLvl);
    void setFuncIsGenerallyIncreasing(bool value);
    void setFuncIsGenerallyDecreasing(bool value);
    void setDeltaX(doublereal deltaXNorm); 

  public:
    ResidEval *m_residFunc;
    doublereal m_funcTargetValue;
    doublereal m_atol;
    doublereal m_rtol;
    doublereal m_maxstep;
  protected:
    int    printLvl;
    doublereal DeltaXnorm_;
    bool FuncIsGenerallyIncreasing_;
    bool FuncIsGenerallyDecreasing_;
    doublereal deltaXConverged_;

  };
}
#endif
