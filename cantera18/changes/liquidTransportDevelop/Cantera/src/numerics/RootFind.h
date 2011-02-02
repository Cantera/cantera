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

  //@{
  ///  @name  Constant which determines the return integer from the routine
 
  //!  This means that the root solver was a success
#define ROOTFIND_SUCCESS             0
  //!  This means that the root solver failed to achieve convergence
#define ROOTFIND_FAILEDCONVERGENCE  -1
  //!  This means that the input to the root solver was defective
#define ROOTFIND_BADINPUT           -2
 //@{

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
     *
     * @param resid  Pointer to the residual function to be used to calculate f(x)
     */ 
    RootFind(ResidEval* resid);
   
    //! Destructor. Deletes the integrator.
    ~RootFind();

  private:

    //! Unimplemented private copy constructor
    /*!
     *  @param right object to be copied
     */
    RootFind(const RootFind &right);

    //! Unimplemented private assignment operator
    /*!
     *  @param right object to be copied
     */
    RootFind& operator=(const RootFind &right);

    //! Calculate a deltaX from an input value of x
    /*!
     *  This routine ensure that the deltaX will be greater or equal to DeltaXNorm_
     *  or 1.0E-14 x
     *
     * @param x1  input value of x
     */
    doublereal delXNonzero(doublereal x1) const;
  
    //! Calculate a deltaX from an input value of x
    /*!
     *  This routine ensure that the deltaX will be greater or equal to DeltaXNorm_
     *  or 1.0E-14 x or deltaXConverged_.
     *
     * @param x1  input value of x
     */
    doublereal delXMeaningful(doublereal x1) const;
 
    //! Calcuated a controlled, nonzero delta between two numbers
    /*!
     *  The delta is designed to be greater than or equal to delXMeaningful(x) defined above
     *  with the same sign as the original delta. Therefore if you subtract it from either
     *  of the two original numbers, you get a different number.
     *
     *  @param x2   first number
     *  @param x1   second number
     */
    doublereal deltaXControlled(doublereal x2, doublereal x1) const;
   
    //! Function to decide whether two real numbers are the same or not
    /*!
     *  A comparison is made between the two numbers to decide whether they
     *  are close to one another. This is defined as being within delXMeaningful() of each other
     *
     * @param x1  First number
     * @param x2  second number
     *
     * @return Returns a boolean indicating whether the two numbers are the same or not.
     */
    bool theSame(doublereal x2, doublereal x1) const;

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

    //! Set the tolerance parameters for the rootfinder
    /*!
     *  These tolerance parameters are used on the function value to determine convergence
     *  
     *
     * @param rtol  Relative tolerance. The default is 10^-5
     * @param atol  absolute tolerance. The default is 10^-11
     */
    void setTol(doublereal rtol, doublereal atol);

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

    //! Set the function behavior flag
    /*!
     *  If this is true, the function is generally an increasing function of x.
     *  In particular, if the algorithm is seeking a higher value of f, it will look
     *  in the positive x direction.
     *
     *  This type of function is needed because this algorithm must deal with regions of f(x) where 
     *  f is not changing with x.
     *  
     *  @param value   boolean value
     */
    void setFuncIsGenerallyIncreasing(bool value);

    //! Set the function behavior flag
    /*!
     *  If this is true, the function is generally a decreasing function of x.
     *  In particular, if the algorithm is seeking a higher value of f, it will look
     *  in the negative x direction.
     *
     *  This type of function is needed because this algorithm must deal with regions of f(x) where 
     *  f is not changing with x.
     *  
     *  @param value   boolean value
     */
    void setFuncIsGenerallyDecreasing(bool value);

    //! Set the minimum value of deltaX
    /*!
     *  This sets the value of deltaXNorm_
     *
     *  @param deltaXNorm
     */
    void setDeltaX(doublereal deltaXNorm); 

  public:

    //!   Pointer to the residual function evaluator
    ResidEval *m_residFunc;

    //!  Target value for the function.  We seek the value of f that is equal to this value
    doublereal m_funcTargetValue;

    //!  Absolute tolerance for the value of f
    doublereal m_atol;

    //!  Relative tolerance for the value of f
    doublereal m_rtol;

    //!  Maximum number of step sizes
    doublereal m_maxstep;
  protected:

    //!  Print level
    int    printLvl;

    //! Delta X norm. This is the minimum value of deltaX that will be used by the program
    doublereal DeltaXnorm_;

    //! Boolean indicating whether the function is an increasing with x
    bool FuncIsGenerallyIncreasing_;

    //!  Boolean indicating whether the function is decreasing with x
    bool FuncIsGenerallyDecreasing_;

    //! Value of delta X that is needed for convergence
    /*!
     *  X will be considered as converged if we are within deltaXConverged_ of the solution
     *  The default is zero.
     */
    doublereal deltaXConverged_;

    doublereal x_maxTried_;
    doublereal fx_maxTried_;

    doublereal x_minTried_;
    doublereal fx_minTried_;

  };
}
#endif
