/**
 *  @file vcs_internal.h Internal declarations for the VCSnonideal package
 */

/*
 * Copyright (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#ifndef _VCS_INTERNAL_H
#define _VCS_INTERNAL_H

#include "cantera/base/global.h"
namespace Cantera
{
//! Points to the data in a std::vector<> object
#define VCS_DATA_PTR(vvv) (&(vvv[0]))

//! define this Cantera function to replace printf
/*!
 * We can replace this with printf easily
 */
#define plogf writelogf

//! define this Cantera function to replace cout << endl;
/*!
 * We use this to place an endl in the log file, and
 * ensure that the IO buffers are flushed.
 */
#define plogendl() writelogendl()

//! Global hook for turning on and off time printing.
/*!
 * Default is to allow printing. But, you can assign this to zero globally to
 * turn off all time printing. This is helpful for test suite purposes where
 * you are interested in differences in text files.
 */
extern int vcs_timing_print_lvl;

// Forward references
class VCS_SPECIES_THERMO;
class VCS_PROB;

//! Class to keep track of time and iterations
/*!
 * class keeps all of the counters together.
 */
class VCS_COUNTERS
{
public:
    //! Total number of iterations in the main loop
    //! of vcs_TP() to solve for thermo equilibrium
    int T_Its;

    //!  Current number of iterations in the main loop
    //!  of vcs_TP() to solve for thermo equilibrium
    int Its;

    //! Total number of optimizations of the components basis set done
    int T_Basis_Opts;

    //! number of optimizations of the components basis set done
    int Basis_Opts;

    //! Current number of times the initial thermo
    //! equilibrium estimator has been called
    int T_Calls_Inest;

    //! Current number of calls to vcs_TP
    int T_Calls_vcs_TP;

    //! Current time spent in vcs_TP
    double T_Time_vcs_TP;

    //! Current time spent in vcs_TP
    double Time_vcs_TP;

    //! Total Time spent in basopt
    double T_Time_basopt;

    //! Current Time spent in basopt
    double Time_basopt;

    //! Time spent in initial estimator
    double T_Time_inest;

    //! Time spent in the vcs suite of programs
    double T_Time_vcs;
};

//! Returns the value of the gas constant in the units specified by parameter
/*!
 *  @param mu_units Specifies the units.
 *           -  VCS_UNITS_KCALMOL: kcal gmol-1 K-1
 *           -  VCS_UNITS_UNITLESS:  1.0 K-1
 *           -  VCS_UNITS_KJMOL:   kJ gmol-1 K-1
 *           -  VCS_UNITS_KELVIN:    1.0 K-1
 *           -  VCS_UNITS_MKS:   joules  kmol-1 K-1 =  kg m2 s-2 kmol-1 K-1
 */
double vcsUtil_gasConstant(int mu_units);

//! Definition of the function pointer for the root finder
/*!
 *  see vcsUtil_root1d for a definition of how to use this.
 */
typedef double(*VCS_FUNC_PTR)(double xval, double Vtarget,
                              int varID, void* fptrPassthrough,
                              int* err);

//! One dimensional root finder
/*!
 *  This root finder will find the root of a one dimensional equation
 *  \f[
 *     f(x) = 0
 *  \f]
 *  where x is a bounded quantity: \f$ x_{min} < x < x_max \f$
 *
 *  The function to be minimized must have the following call structure:
 *
 *  @code
 *  typedef double (*VCS_FUNC_PTR)(double xval, double Vtarget,
 *                                 int varID, void *fptrPassthrough,
 *                                 int *err);  @endcode
 *
 *  xval is the current value of the x variable. Vtarget is the requested
 *  value of f(x), usually 0. varID is an integer that is passed through.
 *  fptrPassthrough is a void pointer that is passed through. err is a return
 *  error indicator. err = 0 is the norm. anything else is considered a fatal
 *  error. The return value of the function is the current value of f(xval).
 *
 *  @param xmin  Minimum permissible value of the x variable
 *  @param xmax  Maximum permissible value of the x parameter
 *  @param itmax Maximum number of iterations
 *  @param func  function pointer, pointing to the function to be
 *               minimized
 *  @param fptrPassthrough Pointer to void that gets passed through
 *                         the rootfinder, unchanged, to the func.
 *  @param FuncTargVal Target value of the function. This is usually set
 *                     to zero.
 *  @param varID       Variable ID. This is usually set to zero.
 *  @param xbest Pointer to the initial value of x on input. On output
 *               This contains the root value.
 *  @param printLvl Print level of the routine.
 *
 * Following is a nontrial example for vcs_root1d() in which the position of a
 * cylinder floating on the water is calculated.
 *
 * @code
 * #include <cmath>
 * #include <cstdlib>
 *
 * #include "equil/vcs_internal.h"
 *
 * const double g_cgs = 980.;
 * const double mass_cyl = 0.066;
 * const double diam_cyl = 0.048;
 * const double rad_cyl = diam_cyl / 2.0;
 * const double len_cyl  = 5.46;
 * const double vol_cyl  = Pi * diam_cyl * diam_cyl / 4 * len_cyl;
 * const double rho_cyl = mass_cyl / vol_cyl;
 * const double rho_gas = 0.0;
 * const double rho_liq = 1.0;
 * const double sigma = 72.88;
 * // Contact angle in radians
 * const double alpha1 = 40.0 / 180. * Pi;
 *
 * double func_vert(double theta1, double h_2, double rho_c) {
 *   double f_grav = - Pi * rad_cyl * rad_cyl * rho_c * g_cgs;
 *   double tmp = rad_cyl * rad_cyl * g_cgs;
 *   double tmp1 = theta1 + sin(theta1) * cos(theta1) - 2.0 * h_2 / rad_cyl * sin(theta1);
 *   double f_buoy = tmp * (Pi * rho_gas + (rho_liq - rho_gas) * tmp1);
 *   double f_sten = 2 * sigma * sin(theta1 + alpha1 - Pi);
 *   return f_grav +  f_buoy +  f_sten;
 * }
 * double calc_h2_farfield(double theta1) {
 *   double rhs = sigma * (1.0 + cos(alpha1 + theta1));
 *   rhs *= 2.0;
 *   rhs = rhs / (rho_liq - rho_gas) / g_cgs;
 *   double sign = -1.0;
 *   if (alpha1 + theta1 < Pi) sign = 1.0;
 *   double res = sign * sqrt(rhs);
 *   return res + rad_cyl * cos(theta1);
 * }
 * double funcZero(double xval, double Vtarget, int varID, void *fptrPassthrough, int *err) {
 *   double theta = xval;
 *   double h2 = calc_h2_farfield(theta);
 *   return func_vert(theta, h2, rho_cyl);
 * }
 * int main () {
 *   double thetamax = Pi;
 *   double thetamin = 0.0;
 *   int maxit = 1000;
 *   int iconv;
 *   double thetaR = Pi/2.0;
 *   int printLvl = 4;
 *
 *   iconv =  vcsUtil_root1d(thetamin, thetamax, maxit,
 *                           funcZero,
 *                           (void *) 0, 0.0, 0,
 *                           &thetaR, printLvl);
 *   printf("theta = %g\n", thetaR);
 *   double h2Final = calc_h2_farfield(thetaR);
 *   printf("h2Final = %g\n", h2Final);
 *   return 0;
 * }
 * @endcode
 * @deprecated Unused. To be removed after Cantera 2.2.
 */
int vcsUtil_root1d(double xmin, double xmax, size_t itmax, VCS_FUNC_PTR func,
                   void* fptrPassthrough,
                   double FuncTargVal, int varID, double* xbest,
                   int printLvl = 0);

//! determine the l2 norm of a vector of doubles
/*!
 * @param vec vector of doubles
 *
 * @return  Returns the l2 norm of the vector
 */
double vcs_l2norm(const std::vector<double> vec);

//! Finds the location of the maximum component in a double vector
/*!
 * @param x pointer to a vector of doubles
 * @param xSize pointer to a vector of doubles used as a multiplier to x[]
 *              before making the decision. Ignored if set to NULL.
 * @param j lowest index to search from
 * @param n highest index to search from
 * @return  Return index of the greatest value on X(i) searched
 *             j <= i < n
 */
size_t vcs_optMax(const double* x, const double* xSize, size_t j, size_t n);

//! Returns the maximum integer in a list
/*!
 *  @param vector pointer to a vector of ints
 *  @param length length of the integer vector
 *
 * @return returns the max integer value in the list
 * @deprecated Unused. To be removed after Cantera 2.2.
 */
int vcs_max_int(const int* vector, int length);

//! Returns a const char string representing the type of the
//! species given by the first argument
/*!
 * @param speciesStatus  Species status integer representing the type
 *                       of the species.
 * @param length         Maximum length of the string to be returned.
 *                       Shorter values will yield abbreviated strings.
 *                       Defaults to a value of 100.
 */
const char* vcs_speciesType_string(int speciesStatus, int length = 100);

//! Print a string within a given space limit
/*!
 *   This routine limits the amount of the string that will be printed to a
 *   maximum of "space" characters. Printing is done to
 *   to Cantera's writelog() function.
 *
 *     @param str  String, which must be null terminated.
 *     @param space   space limit for the printing.
 *     @param alignment Alignment of string within the space:
 *                     -  0 centered
 *                     -  1 right aligned
 *                     -  2 left aligned
 */
void vcs_print_stringTrunc(const char* str, size_t space, int alignment);

//! Simple routine to check whether two doubles are equal up to
//! roundoff error
/*!
 *  Currently it's set to check for 10 digits of relative accuracy.
 *
 * @param d1 first double
 * @param d2 second double
 *
 * @return returns true if the doubles are "equal" and false otherwise
 */
bool vcs_doubleEqual(double d1, double d2);
}

#endif
