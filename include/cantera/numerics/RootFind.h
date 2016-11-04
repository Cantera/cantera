/**
 * @file RootFind.h Header file for implicit nonlinear solver of a one
 *       dimensional function (see \ref numerics and class \link
 *       Cantera::RootFind RootFind\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_ROOTFIND_H
#define CT_ROOTFIND_H
/**
 * @defgroup solverGroup Solvers for Equation Systems
 */

#include "ResidEval.h"

namespace Cantera
{

//@{
///  @name  Constant which determines the return integer from the routine

//!  This means that the root solver was a success
#define ROOTFIND_SUCCESS 0
//! This return value means that the root finder resolved a solution in the x
//! coordinate, however, convergence in F was not achieved.
/*!
 * A common situation for this to happen is that f(x) is discontinuous about
 * f(x) = f_0, where we seek the x where the function is equal to f_0. f(x)
 * spans the f_0 while not being equal to f_0 anywhere.
 */
#define ROOTFIND_SUCCESS_XCONVERGENCEONLY 1
//! This means that the root solver failed to achieve convergence
#define ROOTFIND_FAILEDCONVERGENCE -1
//! This means that the input to the root solver was defective
#define ROOTFIND_BADINPUT -2
//! This means that the rootfinder believes the solution is lower than xmin
#define ROOTFIND_SOLNLOWERTHANXMIN -3
//! This means that the rootfinder believes the solution is higher than xmax
#define ROOTFIND_SOLNHIGHERTHANXMAX -4
//@}

//! Root finder for 1D problems
/*!
 * @deprecated Unused. To be removed after Cantera 2.3. See
 *     boost::math::tools::toms748_solve for an alternative.
 *
 * The root finder solves a single nonlinear equation described below.
 *
 * \f[
 *     f(x) = f_0
 * \f]
 *
 * \f$ f(x) \f$ is assumed to be single valued as a function of x.\f$ f(x) \f$
 * is not assumed to be continuous nor is its derivative assumed to be well
 * formed.
 *
 * Root finders are significantly different in the sense that do not have to
 * rely solely on Newton's method to find the answer to the problem. Instead
 * they use a method to bound the solution between high and low values and then
 * use a method to refine that bound.  The eventual solution to the problem is
 * presented as x_best and as a bound, delta_X, on the solution component.
 * Because of this, they are far more stable for functions and Jacobians that
 * have discontinuities or noise associated with them.
 *
 * The algorithm is a convolution of a local Secant method with an approach of
 * finding a straddle in x. The Jacobian is never required.
 *
 * There is a general breakdown of the algorithm into stages. The first stage
 * seeks to find a straddle of the function. The second stage seeks to reduce
 * the bounds in x and f in order to satisfy the specification of the stopping
 * criteria. In the last stage the algorithm seeks to find the base value of x
 * that satisfies the original equation given what it current knows about the
 * function.
 *
 * Globalization strategy
 *
 * Specifying the General Changes in x
 *
 * Supplying Hints with General Function Behavior Flags
 *
 * Stopping Criteria
 *
 * Specification of the Stopping Criteria
 *
 * Additional constraints
 *
 * Bounds Criteria For the Routine
 *
 * Example
 *
 * @code
 * // Define a residual. The definition of a residual involves a lot more work than is shown here.
 * ResidEval * ec;
 * // Instantiate the root finder with the residual to be solved, ec.
 * RootFind rf(&ec);
 * // Set the relative and absolute tolerances for f and x.
 * rf.setTol(1.0E-5, 1.0E-10, 1.0E-5, 1.0E-11);
 * // Give a hint about the function's dependence on x. This is needed, for example, if the function has
 * // flat regions.
 * rf.setFuncIsGenerallyIncreasing(true);
 * rf.setDeltaX(0.01);
 * // Supply an initial guess for the solution
 * double xbest = phiM;
 * double oldP = printLvl_;
 * // Set the print level for the solver. Zero produces no output. Two produces a summary table of each iteration.
 * rf.setPrintLvl(2);
 * // Define a minimum and maximum for the independent variable.
 * double phimin = 1.3;
 * double phimax = 2.2;
 * // Define a maximum iteration number
 * int itmax = 100;
 * // Define the f_0 value, and on return will contain the actual value of f(x) obtained
 * double currentObtained;
 * // Call the solver
 * status = rf.solve(phimin, phimax, 100, currentObtained, &xbest);
 * if (status == 0) {
 *   if (printLvl_ > 1) {
 *     printf("Electrode::integrateConstantCurrent(): Volts (%g amps) = %g\n", currentObtained, xbest);
 *   }
 * } else {
 *   if (printLvl_) {
 *      printf("Electrode::integrateConstantCurrent(): bad status = %d Volts (%g amps) = %g\n",
 *             status, currentObtained, xbest);
 *   }
 * }
 * @endcode
 *
 * @todo  Noise
 * @todo  General Search to be done when all else fails
 */
class RootFind
{
public:
    //! Constructor for the object
    /*!
     * @param resid  Pointer to the residual function to be used to calculate f(x)
     */
    RootFind(ResidEval* resid);

    //! @deprecated To be removed after Cantera 2.3.
    RootFind(const RootFind& r);
    ~RootFind() {}
    //! @deprecated To be removed after Cantera 2.3.
    RootFind& operator=(const RootFind& right);

private:
    //! Calculate a deltaX from an input value of x
    /*!
     * This routine ensure that the deltaX will be greater or equal to
     * DeltaXNorm_ or 1.0E-14 x
     *
     * @param x1  input value of x
     */
    doublereal delXNonzero(doublereal x1) const;

    //! Calculate a deltaX from an input value of x
    /*!
     * This routine ensure that the deltaX will be greater or equal to
     * DeltaXNorm_ or 1.0E-14 x or deltaXConverged_.
     *
     * @param x1  input value of x
     */
    doublereal delXMeaningful(doublereal x1) const;

    //! Calculate a controlled, nonzero delta between two numbers
    /*!
     * The delta is designed to be greater than or equal to delXMeaningful(x)
     * defined above with the same sign as the original delta. Therefore if you
     * subtract it from either of the two original numbers, you get a different
     * number.
     *
     * @param x2   first number
     * @param x1   second number
     */
    doublereal deltaXControlled(doublereal x2, doublereal x1) const;

    //! Function to decide whether two real numbers are the same or not
    /*!
     * A comparison is made between the two numbers to decide whether they are
     * close to one another. This is defined as being within factor *
     * delXMeaningful() of each other.
     *
     * The basic premise here is that if the two numbers are too close, the
     * noise will prevent an accurate calculation of the function and its slope.
     *
     * @param x1  First number
     * @param x2  second number
     * @param factor  Multiplicative factor to multiple deltaX with
     * @returns a boolean indicating whether the two numbers are the same or not.
     */
    bool theSame(doublereal x2, doublereal x1, doublereal factor = 1.0) const;

public:
    //! Using a line search method, find the root of a 1D function
    /*!
     * This routine solves the following equation.
     *
     * \f[
     *     R(x) = f(x) - f_o = 0
     * \f]
     *
     * @param xmin    Minimum value of x to be used.
     * @param xmax    Maximum value of x to be used
     * @param itmax   maximum number of iterations. Usually, it can be less than 50.
     * @param funcTargetValue Value of \f$ f_o \f$ in the equation. On return,
     *                it contains the value of the function actually obtained.
     * @param xbest   Returns the x that satisfies the function On input, xbest
     *                should contain the best estimate of the solution. An
     *                attempt to find the solution near xbest is made.
     * @return:
     *    0  =  ROOTFIND_SUCCESS            Found function
     *   -1  =  ROOTFIND_FAILEDCONVERGENCE  Failed to find the answer
     *   -2  =  ROOTFIND_BADINPUT           Bad input was detected
     */
    int solve(doublereal xmin, doublereal xmax, int itmax, doublereal& funcTargetValue, doublereal* xbest);

    //! Return the function value
    /*!
     * This routine evaluates the following equation.
     *
     *    \f[
     *       R(x) = f(x) - f_o = 0
     *    \f]
     *
     * @param x  Value of the independent variable
     *
     * @return   The routine returns the value of \f$ R(x) \f$
     */
    doublereal func(doublereal x);

    //! Set the tolerance parameters for the rootfinder
    /*!
     * These tolerance parameters are used on the function value and the
     * independent value to determine convergence
     *
     * @param rtolf  Relative tolerance. The default is 10^-5
     * @param atolf  absolute tolerance. The default is 10^-11
     * @param rtolx  Relative tolerance. The default is 10^-5. Default parameter
     *                   is 0.0, in which case rtolx is set equal to rtolf
     * @param atolx  absolute tolerance. The default is 10^-11. Default
     *                   parameter is 0.0, in which case atolx is set equal to
     *                   atolf
     */
    void setTol(doublereal rtolf, doublereal atolf, doublereal rtolx = 0.0, doublereal atolx = 0.0);

    //! Set the print level from the rootfinder
    /*!
     * - 0: No printing of any kind
     * - 1: Single print line indicating success or failure of the routine.
     * - 2: Summary table printed at the end of the routine, with a convergence
     *      history
     * - 3: Printouts during the iteration are added. Summary table is printed
     *      out at the end. if writeLogAllowed_ is turned on, a file is written
     *      out with the convergence history.
     *
     *  @param printLvl  integer value
     */
    void setPrintLvl(int printLvl);

    //! Set the function behavior flag
    /*!
     * If this is true, the function is generally an increasing function of x.
     * In particular, if the algorithm is seeking a higher value of f, it will
     * look in the positive x direction.
     *
     * This type of function is needed because this algorithm must deal with
     * regions of f(x) where f is not changing with x.
     *
     *  @param value   boolean value
     */
    void setFuncIsGenerallyIncreasing(bool value);

    //! Set the function behavior flag
    /*!
     * If this is true, the function is generally a decreasing function of x. In
     * particular, if the algorithm is seeking a higher value of f, it will look
     * in the negative x direction.
     *
     * This type of function is needed because this algorithm must deal with
     * regions of f(x) where f is not changing with x.
     *
     *  @param value   boolean value
     */
    void setFuncIsGenerallyDecreasing(bool value);

    //! Set the minimum value of deltaX
    /*!
     * @param deltaXNorm
     */
    void setDeltaX(doublereal deltaXNorm);

    //! Set the maximum value of deltaX
    /*!
     *  @param deltaX
     */
    void setDeltaXMax(doublereal deltaX);

    //! Print the iteration history table
    void printTable();

public:
    //! Pointer to the residual function evaluator
    ResidEval* m_residFunc;

    //! Target value for the function. We seek the value of f that is equal to
    //! this value
    doublereal m_funcTargetValue;

    //! Absolute tolerance for the value of f
    doublereal m_atolf;

    //! Absolute tolerance for the value of x
    doublereal m_atolx;

    //! Relative tolerance for the value of f and x
    doublereal m_rtolf;

    //! Relative tolerance for the value of x
    doublereal m_rtolx;

    //! Maximum number of step sizes
    doublereal m_maxstep;

protected:
    //! Print level. @see setPrintLvl
    int printLvl;

public:
    //! Boolean to turn on the possibility of writing a log file.
    bool writeLogAllowed_;

protected:
    //! Delta X norm. This is the nominal value of deltaX that will be used by
    //! the program
    doublereal DeltaXnorm_;

    //! Boolean indicating whether DeltaXnorm_ has been specified by the user or
    //! not
    int specifiedDeltaXnorm_;

    //! Delta X Max.
    /*!
     * This is the maximum value of deltaX that will be used by the program.
     * Sometimes a large change in x causes problems.
     */
    doublereal DeltaXMax_;

    //! Boolean indicating whether DeltaXMax_ has been specified by the user or
    //! not
    int specifiedDeltaXMax_;

    //! Boolean indicating whether the function is an increasing with x
    bool FuncIsGenerallyIncreasing_;

    //! Boolean indicating whether the function is decreasing with x
    bool FuncIsGenerallyDecreasing_;

    //! Value of delta X that is needed for convergence
    /*!
     *  X will be considered as converged if we are within deltaXConverged_ of
     *  the solution The default is zero.
     */
    doublereal deltaXConverged_;

    //! Internal variable tracking largest x tried.
    doublereal x_maxTried_;

    //! Internal variable tracking f(x) of largest x tried.
    doublereal fx_maxTried_;

    //! Internal variable tracking smallest x tried.
    doublereal x_minTried_;

    //! Internal variable tracking f(x) of smallest x tried.
    doublereal fx_minTried_;

    //! Structure containing the iteration history
    struct rfTable {
        //@{
        int its;
        int TP_its;
        double slope;
        double xval;
        double fval;
        int foundPos;
        int foundNeg;
        double deltaXConverged;
        double deltaFConverged;
        double delX;

        std::string reasoning;

        void clear() {
            its = 0;
            TP_its = 0;
            slope = -1.0E300;
            xval = -1.0E300;
            fval = -1.0E300;
            reasoning = "";
        };

        rfTable() :
            its(-2),
            TP_its(0),
            slope(-1.0E300),
            xval(-1.0E300),
            fval(-1.0E300),
            foundPos(0),
            foundNeg(0),
            deltaXConverged(-1.0E300),
            deltaFConverged(-1.0E300),
            delX(-1.0E300) {
        };
        //@}
    };

    //! Vector of iteration histories
    std::vector<struct rfTable> rfHistory_;
};
}
#endif
