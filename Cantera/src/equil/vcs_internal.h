/**
 *  @file vcs_internal.h
 *      Internal declarations for the VCSnonideal package
 */

/*
 * $Id$
 */

/*
 * Copywrite (2005) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#ifndef _VCS_INTERNAL_H
#define _VCS_INTERNAL_H

#include "vcs_defs.h"
#include "vcs_DoubleStarStar.h"
#include "vcs_Exception.h"

#include "global.h"

namespace VCSnonideal {

  //! Points to the data in a std::vector<> object
#define VCS_DATA_PTR(vvv) (&(vvv[0]))

  //! define this Cantera function to replace printf
  /*!
   * We can replace this with printf easily
   */
#define plogf  Cantera::writelogf


  /*****************************************************************************/
  /*****************************************************************************/
  /*****************************************************************************/

  /************ ERROR HANDLING *****************/

  struct VCS_ERR {
    int Flag;
    int Species1;
    int Species2;
    int PrintLevel;
    double Value1;
    double Value2;
    char Mess[120];
  };
  typedef struct VCS_ERR VCS_ERR_STRUCT;


  extern int vcsUtil_err_check(VCS_ERR_STRUCT &vcsE, char *string1, int ival);
  extern void vcsUtil_err_reset(VCS_ERR_STRUCT &vcsE);


  /*********************************/
  /* Function Pointer Typedefs */
  /*********************************/

  typedef double (*VCS_FUNC_PTR)(double, double, int, void *, int *);

  /*
   * Forward references
   */
  class VCS_SPECIES_THERMO;
  class VCS_PROB;


  /****************************************************************************/
  /****************************************************************************/
  /****************************************************************************/

  //!  Amount of extra printing that is done while in debug mode.
  /*!
   *                     0 -> none
   *                     1 -> some    
   *                     2 -> alot      (default)
   *                     3 -> everything
   */

  class VCS_COUNTERS {
  public:
    int T_Its;         /* Total number of iterations in the main loop 
			  of vcs_TP() to solve for thermo equilibrium */
    int Its;               /* Current number of iterations in the main loop 
			      of vcs_TP() to solve for thermo equilibrium */
    int T_Basis_Opts;  /* Total number of optimizations of the 
			  components basis set done */
    int Basis_Opts;        /* Total number of optimizations of the 
			      components basis set done */
    int T_Calls_Inest;       /* Current number of times the initial thermo
				equilibrium estimator has been called */
    int T_Calls_vcs_TP;      /* Current number of calls to vcs_TP */
    double T_Time_vcs_TP;  /* Total time spent in vcs_TP */
    double Time_vcs_TP;    /* Current time spent in vcs_TP */
    double T_Time_basopt;  /* Total Time spent in basopt () */
    double Time_basopt;    /* Current Time spent in basopt () */
    double T_Time_inest;     /* Time spent in initial estimator */
    double T_Time_vcs;       /* Time spent in the vcs suite of programs */
  };

  /*****************************************************************************/
  /**************** Prototypes *************************************************/
  /*****************************************************************************/

  /* Externals for vcs_funcVtot.c */

  extern double vcs_funcVtot(double, double, int, void *, int *);

  /* Externals defined in vcs_nondim.c */

  extern double vcsUtil_gasConstant(int);

  /* Externals in vcs_solve_TP.c */

  extern int    vcsUtil_mlequ(double *, int, int, double *, int);
  extern void   vcsUtil_dsw(double *, int ,int);
  extern void   vcsUtil_isw(int [], int, int);
  extern void   vcsUtil_ssw(char **, int, int); 
  extern void   vcsUtil_stsw(std::vector<std::string> & vecStrings, int, int); 

  /* Externals for vcs_root1d.c */

  //! One dimensional root finder
  /*!
   *
   *  This root finder will find the root of a one dimensional
   *  equation
   *
   *  \f[
   *     f(x) = 0
   *  \f]
   *  where x is a bounded quantity: \f$ x_{min} < x < x_max \f$ 
   *
   *  The functional to be minimized must have the following call
   *  structure:
   *
   *    @verbatim
           typedef double (*VCS_FUNC_PTR)(double xval, double Vtarget,
                                          int varID, void *fptrPassthrough, 
                                          int *err);  @endverbatim
   *
   *    xval is the current value of the x variable. Vtarget is the
   *    requested value of f(x), usually 0. varID is an integer
   *    that is passed through. fptrPassthrough is a void pointer
   *    that is passed through. err is a return error indicator.
   *    err = 0 is the norm. anything else is considered a fatal
   *    error.
   *    The return value of the function is the current value of
   *    f(xval).
   *
   *  @param xmin  Minimum permissible value of the x variable
   *  @param xmax  Maximum permissible value of the x paramerer
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
   *
   * Following is a nontrial example for vcs_root1d() in which the position of a 
   * cylinder floating on the water is calculated.
   *
   *    @verbatim
     #include "math.h"
     #include "stdlib.h"
    
     #include "Cantera.h"
     #include "kernel/vcs_internal.h"
     using namespace Cantera;
     using namespace VCSnonideal;

     const double g_cgs = 980.;
     const double mass_cyl = 0.066;
     const double diam_cyl = 0.048;
     const double rad_cyl = diam_cyl / 2.0; 
     const double len_cyl  = 5.46;
     const double vol_cyl  = Pi * diam_cyl * diam_cyl / 4 * len_cyl;
     const double rho_cyl = mass_cyl / vol_cyl; 
     const double rho_gas = 0.0;
     const double rho_liq = 1.0;
     const double sigma = 72.88;
     // Contact angle in radians
     const double alpha1 = 40.0 / 180. * Pi;
      
     double func_vert(double theta1, double h_2, double rho_c) {
       double f_grav = - Pi * rad_cyl * rad_cyl * rho_c * g_cgs;
       double tmp = rad_cyl * rad_cyl * g_cgs;
       double tmp1 = theta1 + sin(theta1) * cos(theta1) - 2.0 * h_2 / rad_cyl * sin(theta1);
       double f_buoy = tmp * (Pi * rho_gas + (rho_liq - rho_gas) * tmp1);
       double f_sten = 2 * sigma * sin(theta1 + alpha1 - Pi);
       double f_net =  f_grav +  f_buoy +  f_sten;
       return f_net;
     }
     double calc_h2_farfield(double theta1) {
       double rhs = sigma * (1.0 + cos(alpha1 + theta1));
       rhs *= 2.0;
       rhs = rhs / (rho_liq - rho_gas) / g_cgs;
       double sign = -1.0;
       if (alpha1 + theta1 < Pi) sign = 1.0;
       double res = sign * sqrt(rhs);
       double h2 = res + rad_cyl * cos(theta1);
       return h2;
     }
     double funcZero(double xval, double Vtarget, int varID, void *fptrPassthrough, int *err) {
       double theta = xval;
       double h2 = calc_h2_farfield(theta);
       double fv = func_vert(theta, h2, rho_cyl);
       return fv;
     }
     int main () {
       double thetamax = Pi;
       double thetamin = 0.0;
       int maxit = 1000;
       int iconv;
       double thetaR = Pi/2.0;
       int printLvl = 4;
   
       iconv =  VCSnonideal::vcsUtil_root1d(thetamin, thetamax, maxit, funcZero,
                                            (void *) 0, 0.0, 0, &thetaR, printLvl);
       printf("theta = %g\n", thetaR);
       double h2Final = calc_h2_farfield(thetaR);
       printf("h2Final = %g\n", h2Final);
       return 0;
     }   @endverbatim
   *
   */
  int vcsUtil_root1d(double xmin, double xmax, int itmax, VCS_FUNC_PTR func,
		     void *fptrPassthrough, 
		     double FuncTargVal, int varID, double *xbest, int printLvl = 0);

  /* Externals defined in vcs_timer_generic.c */

  extern double vcs_second(void);


  /* Externals defined in vcs_util.c */

#define USE_MEMSET
#ifdef USE_MEMSET
#include <string.h>
#  define vcs_dzero(vector, length) (void) memset((void *) (vector), 0, \
                                    (length) * sizeof(double))
#  define vcs_izero(vector, length) (void)  memset((void *) (vector), 0, \
                                    (length) * sizeof(int))
#  define vcs_dcopy(vec_to, vec_from, length) \
                     (void) memcpy((void *) (vec_to), (const void *) (vec_from), \
                                   (length) * sizeof(double))
#  define vcs_icopy(vec_to, vec_from, length) \
                     (void) memcpy((void *) (vec_to), (const void *) (vec_from), \
                                   (length) * sizeof(int))
#  define vcs_vdzero(vector, length) (void) memset(VCS_DATA_PTR(vector), 0, \
                                    (length) * sizeof(double))
#  define vcs_vizero(vector, length) (void)  memset(VCS_DATA_PTR(vector), 0, \
                                    (length) * sizeof(int))

  inline void  vcs_vdcopy(std::vector<double> & vec_to, 
			  const std::vector<double> & vec_from, int length) {
    (void) memcpy((void *)&(vec_to[0]), (const void *) &(vec_from[0]), 
		  (length) * sizeof(double));
  }
  inline void  vcs_vicopy(std::vector<int> & vec_to, 
			  const std::vector<int> & vec_from, int length) {
    (void) memcpy((void *)&(vec_to[0]), (const void *) &(vec_from[0]), 
		  (length) * sizeof(int));
  }
#else
  extern void vcs_dzero(double *, int);
  extern void vcs_izero(int *, int);
  extern void vcs_dcopy(double *, double *, int);
  extern void vcs_icopy(int *, int *, int);
  extern void vcs_vdzero(std::vector<double> &vvv, int len = -1);
  extern void vcs_vizero(std::vector<double> &vvv, int len = -1);
  void vcs_vdcopy(std::vector<double> &vec_to, 
		  const std::vector<double> vec_from, int len = -1);
  void vcs_vicopy(std::vector<int> &vec_to, 
		  const std::vector<int> vec_from, int len = -1);
#endif
  extern int  vcs_amax(double *, int, int);
  extern int  vcs_max_int(int *, int);
  extern void vcs_print_line(const char *, int);
  extern void vcs_print_stringTrunc(const char *, int, int);
  extern bool vcs_doubleEqual(double, double);

  /****************************************************************************/
}

#endif
/****************************************************************************/
