/******************************************************************
 *                                                                *
 * File          : cvbandpre.h                                    *
 * Programmers   : Michael Wittman and Alan C. Hindmarsh @ LLNL   *
 * Version of    : 23 March 2000                                  *
 *----------------------------------------------------------------*
 * This is the header file for the CVBANDPRE module, which        *
 * provides a banded difference quotient Jacobian-based           *
 * preconditioner and solver routines for use with CVSPGMR.       *
 *                                                                *
 * Summary:                                                       *
 * These routines provide a band matrix preconditioner based on   *
 * difference quotients of the ODE right-hand side function f.    *
 * The user supplies parameters                                   *
 *   mu = upper half-bandwidth (number of super-diagonals)        *
 *   ml = lower half-bandwidth (number of sub-diagonals)          *
 * The routines generate a band matrix of bandwidth ml + mu + 1   *
 * and use this to form a preconditioner for use with the Krylov  *
 * linear solver in CVSPGMR.  Although this matrix is intended    *
 * to approximate the Jacobian df/dy, it may be a very crude      *
 * approximation.  The true Jacobian need not be banded, or its   *
 * true bandwith may be larger than ml + mu + 1, as long as the   *
 * banded approximation generated here is sufficiently accurate   *
 * to speed convergence as a preconditioner.                      *
 *                                                                *
 * Usage:                                                         *
 *   The following is a summary of the usage of this module.      *
 *   Details of the calls to CVodeMalloc, CVSpgmr, and CVode are  *
 *   available in the CVODE User Guide.                           *
 *   To use these routines, the sequence of calls in the user     *
 *   main program should be as follows:                           *
 *                                                                *
 *   CVBandPreData bp_data;                                       *
 *   ...                                                          *
 *   bp_data = CVBandPreAlloc(N, f, f_data, mu, ml);              *
 *   ...                                                          *
 *   cvode_mem = CVodeMalloc(...);                                *
 *   ...                                                          *
 *   CVSpgmr(cvode_mem, pretype, gstype, maxl, delt,              *
 *           CVBandPrecond, CVBandPSolve, bp_data);               *
 *   ...                                                          *
 *   flag = CVode(...);                                           *
 *   ...                                                          *
 *   CVBandPreFree(bp_data);                                      *
 *   ...                                                          *
 *   CVodeFree(cvode_mem);                                        *
 *                                                                *
 * Notes:                                                         *
 * (1) Include this file for the CVBandPreData type definition.   *
 * (2) In the CVBandPreAlloc call, the arguments N, f, and f_data *
 *     are the same as in the call to CVodeMalloc.                *
 * (3) In the CVSpgmr call, the user is free to specify the inputs*
 *     pretype and gstype, and the optional inputs maxl and delt. *
 *     But the last three arguments must be as shown, with the    *
 *     last argument being the pointer returned by CVBandPreAlloc.*
 * (4) The CVBandPrecond and CVBandPSolve functions are never     *
 *     called by the user explicitly; they are simply passed to   *
 *     the CVSpgmr function.                                      *
 *                                                                *
 ******************************************************************/

 
#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _cvbandpre_h
#define _cvbandpre_h

#include "cvode.h"
#include "llnltyps.h"
#include "nvector.h"
#include "band.h"


/************* CVBandPreData type definition ************/

typedef struct {
  /* Data set by user in CVBandPreAlloc: */
  RhsFn f;
  void *f_data;
  integer ml, mu;

  /* Data set by CVBandPrecond: */
  BandMat savedJ;
  BandMat savedP;
  integer *pivots;
} *CVBandPreData;


/******************************************************************
 *                                                                *
 * Function : CVBandPreAlloc                                      *
 *----------------------------------------------------------------*
 * CVBandPreAlloc allocates and initializes a CVBandPreData       *
 * structure to be passed to CVSpgmr (and subsequently used by    *
 * CVBandPrecond and CVBandPSolve).                               *
 *                                                                *
 * The parameters of CVBandPreAlloc are as follows:               *
 *                                                                *
 * N       is the length of all vector arguments.                 *
 *                                                                *
 * f       is the right hand side function.                       *
 *                                                                *
 * f_data  is a pointer to the optional extra data for f.         *
 *                                                                *
 * mu      is the upper half bandwidth.                           *
 *                                                                *
 * ml      is the lower half bandwidth.                           *
 *                                                                *
 * CVBandPreAlloc returns the storage pointer (type CVBandPreData)*
 * or NULL if the request for storage cannot be satisfied.        *
 *                                                                *
 ******************************************************************/

CVBandPreData CVBandPreAlloc(integer N, RhsFn f, void *f_data,
                             integer mu, integer ml);


/******************************************************************
 *                                                                *
 * Function : CVBandPreFree                                       *
 *----------------------------------------------------------------*
 * CVBandPreFree frees the memory allocated by CVBandPreAlloc in  *
 * the argument pdata.                                            *
 *                                                                *
 ******************************************************************/

void CVBandPreFree(CVBandPreData pdata);



/* Prototypes of CVBandPrecond and CVBandPSolve */

  
int CVBandPrecond(integer N, real t, N_Vector y, N_Vector fy, boole jok,
                  boole *jcurPtr, real gamma, N_Vector ewt, real h,
                  real uround, long int *nfePtr, void *bp_data,
                  N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);


int CVBandPSolve(integer N, real t, N_Vector y, N_Vector fy, N_Vector vtemp,
                 real gamma, N_Vector ewt, real delta, long int *nfePtr,
                 N_Vector r, int lr, void *bp_data, N_Vector z);


#endif

#ifdef __cplusplus
}
#endif
