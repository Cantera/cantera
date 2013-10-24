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

#include <cstring>

#include "cantera/equil/vcs_defs.h"
#include "cantera/base/global.h"

namespace VCSnonideal
{
using Cantera::npos;

//! Points to the data in a std::vector<> object
#define VCS_DATA_PTR(vvv) (&(vvv[0]))

//! define this Cantera function to replace printf
/*!
 * We can replace this with printf easily
 */
#define plogf  Cantera::writelogf

//! define this Cantera function to replace cout << endl;
/*!
 * We use this to place an endl in the log file, and
 * ensure that the IO buffers are flushed.
 */
#define plogendl() Cantera::writelogendl()

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

//! Invert an n x n matrix and solve m rhs's
/*!
 * Solve a square matrix with multiple right hand sides
 *
 * \f[
 *     C X + B = 0;
 * \f]
 *
 * This routine uses Gauss elimination and is optimized for the solution
 * of lots of rhs's. A crude form of row pivoting is used here.
 * The matrix C is destroyed during the solve.
 *
 * @return  The solution x[] is returned in the matrix <I>B</I>.
 *          Routine returns an integer representing success:
 *     -   1 : Matrix is singular
 *     -   0 : solution is OK
 *
 *  @param c  Matrix to be inverted. c is in fortran format, i.e., rows
 *            are the inner loop. Row  numbers equal to idem.
 *            c[i+j*idem] = c_i_j = Matrix to be inverted:
 *                   -  i = row number
 *                   -  j = column number
 *
 *  @param idem number of row dimensions in c
 *  @param n  Number of rows and columns in c
 *  @param b  Multiple RHS. Note, b is actually the negative of
 *            most formulations.  Row  numbers equal to idem.
 *             b[i+j*idem] = b_i_j = vectors of rhs's:
 *                   -  i = row number
 *                   -  j = column number
 *            (each column is a new rhs)
 *  @param m  number of rhs's
 */
int  vcsUtil_mlequ(double* c, size_t idem, size_t n, double* b, size_t m);

//! Invert an n x n matrix and solve m rhs's
/*!
 * Solve a square matrix with multiple right hand sides
 *
 * \f[
 *     C X + B = 0;
 * \f]
 *
 * This routine uses Gauss-Jordan elimination and is optimized for the solution
 * of lots of rhs's. Full row and column pivoting is used here. It's been
 * shown to be necessary in at least one case.
 * The matrix C is destroyed during the solve.
 *
 * @return  The solution x[] is returned in the matrix <I>B</I>.
 *          Routine returns an integer representing success:
 *     -   1 : Matrix is singular
 *     -   0 : solution is OK
 *
 *  @param c  Matrix to be inverted. c is in fortran format, i.e., rows
 *            are the inner loop. Row  numbers equal to idem.
 *            c[i+j*idem] = c_i_j = Matrix to be inverted:
 *                   -  i = row number
 *                   -  j = column number
 *
 *  @param idem number of row dimensions in c
 *  @param n  Number of rows and columns in c
 *  @param b  Multiple RHS. Note, b is actually the negative of
 *            most formulations.  Row  numbers equal to idem.
 *             b[i+j*idem] = b_i_j = vectors of rhs's:
 *                   -  i = row number
 *                   -  j = column number
 *            (each column is a new rhs)
 *  @param m  number of rhs's
 */
int vcsUtil_gaussj(double* c, size_t idem, size_t n, double* b, size_t m);


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
 *   iconv =  VCSnonideal::vcsUtil_root1d(thetamin, thetamax, maxit,
 *                                        funcZero,
 *                                        (void *) 0, 0.0, 0,
 *                                        &thetaR, printLvl);
 *   printf("theta = %g\n", thetaR);
 *   double h2Final = calc_h2_farfield(thetaR);
 *   printf("h2Final = %g\n", h2Final);
 *   return 0;
 * }
 * @endcode
 */
int vcsUtil_root1d(double xmin, double xmax, size_t itmax, VCS_FUNC_PTR func,
                   void* fptrPassthrough,
                   double FuncTargVal, int varID, double* xbest,
                   int printLvl = 0);

//! Returns the system wall clock time in seconds
/*!
 * @return time in seconds.
 */
double vcs_second();

//! This define turns on using memset and memcpy. I have not run into
//! any systems where this is a problem. It's the fastest way to do
//! low lvl operations where applicable. There are alternative routines
//! available if this ever fails.
#define USE_MEMSET
#ifdef USE_MEMSET

//! Zero a double vector
/*!
 * @param vec_to vector of doubles
 * @param length length of the vector to zero.
 */
inline void vcs_dzero(double* const vec_to, const size_t length)
{
    (void) memset((void*) vec_to, 0, length * sizeof(double));
}

//! Zero an int vector
/*!
 * @param vec_to vector of ints
 * @param length length of the vector to zero.
 */
inline void vcs_izero(int* const vec_to, const size_t length)
{
    (void) memset((void*) vec_to, 0, length * sizeof(int));
}

//! Copy a double vector
/*!
 * @param vec_to    Vector to copy into. This vector must be dimensioned
 *                  at least as large as the vec_from vector.
 * @param vec_from  Vector to copy from
 * @param length    Number of doubles to copy.
 */
inline void vcs_dcopy(double* const vec_to,
                      const double* const vec_from, const size_t length)
{
    (void) memcpy((void*) vec_to, (const void*) vec_from,
                  (length) * sizeof(double));
}

//! Copy an int vector
/*!
 * @param vec_to    Vector to copy into. This vector must be dimensioned
 *                  at least as large as the vec_from vector.
 * @param vec_from  Vector to copy from
 * @param length    Number of int to copy.
 */
inline void vcs_icopy(int* const vec_to,
                      const int* const vec_from, const size_t length)
{
    (void) memcpy((void*) vec_to, (const void*) vec_from,
                  (length) * sizeof(int));
}

//! Zero a std double vector
/*!
 * @param vec_to vector of doubles
 * @param length length of the vector to zero.
 */
inline void vcs_vdzero(std::vector<double> &vec_to, const size_t length)
{
    (void) memset((void*)VCS_DATA_PTR(vec_to), 0, (length) * sizeof(double));
}

//! Zero a std int vector
/*!
 * @param vec_to vector of ints
 * @param length length of the vector to zero.
 */
inline void vcs_vizero(std::vector<int> &vec_to, const size_t length)
{
    (void) memset((void*)VCS_DATA_PTR(vec_to), 0, (length) * sizeof(int));
}

//! Copy one std double vector into another
/*!
 * This is an inlined function that uses memcpy. memcpy is probably
 * the fastest way to do this. This routine requires the vectors to be
 * previously dimensioned appropriately. No error checking is done.
 *
 * @param vec_to  Vector to copy into. This vector must be dimensioned
 *                at least as large as the vec_from vector.
 * @param vec_from Vector to copy from
 * @param length  Number of doubles to copy.
 */
inline void vcs_vdcopy(std::vector<double> & vec_to,
                       const std::vector<double> & vec_from, size_t length)
{
    (void) memcpy((void*)&(vec_to[0]), (const void*) &(vec_from[0]),
                  (length) * sizeof(double));
}

//! Copy one std integer vector into another
/*!
 * This is an inlined function that uses memcpy. memcpy is probably
 * the fastest way to do this.
 *
 * @param vec_to  Vector to copy into. This vector must be dimensioned
 *                at least as large as the vec_from vector.
 * @param vec_from Vector to copy from
 * @param length  Number of integers to copy.
 */
inline void vcs_vicopy(std::vector<int> & vec_to,
                       const std::vector<int> & vec_from, const int length)
{
    (void) memcpy((void*)&(vec_to[0]), (const void*) &(vec_from[0]),
                  (length) * sizeof(int));
}
#else
extern void vcs_dzero(double* const, const int);
extern void vcs_izero(int* const , const int);
extern void vcs_dcopy(double* const, const double* const, const int);
extern void vcs_icopy(int* const, const int* const, const int);
extern void vcs_vdzero(std::vector<double> &vvv, const int len = -1);
extern void vcs_vizero(std::vector<double> &vvv, const int len = -1);
void vcs_vdcopy(std::vector<double> &vec_to,
                const std::vector<double> vec_from, const int len = -1);
void vcs_vicopy(std::vector<int> &vec_to,
                const std::vector<int> vec_from, const int len = -1);
#endif

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
 */
int vcs_max_int(const int* vector, int length);

//! Prints a line consisting of multiple occurrences of the same string
/*!
 *  This prints a string num times, and then terminate with a
 *  end of line character
 *
 * @param str C string that is null terminated
 * @param num number of times the string is to be printed
 */
void vcs_print_line(const char* str, int num);

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

//! Sorts a vector of ints in place from lowest to the highest values
/*!
 *  The vector is returned sorted from lowest to highest.
 *
 * @param x Reference to a vector of ints.
 * @deprecated
 */
void vcs_heapsort(std::vector<int> &x);

//! Sorts a vector of ints and eliminates duplicates from the resulting list
/*!
 * @param xOrderedUnique       Ordered vector of unique ints that were part of the original list
 * @param x                    Reference to a constant vector of ints.
 * @deprecated
 */
void vcs_orderedUnique(std::vector<int> & xOrderedUnique, const std::vector<int> & x);

}

#endif
