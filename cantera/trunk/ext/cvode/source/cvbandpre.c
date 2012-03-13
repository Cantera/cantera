/******************************************************************
 *                                                                *
 * File          : cvbandpre.c                                    *
 * Programmers   : Michael Wittman and Alan C. Hindmarsh @ LLNL   *
 * Version of    : 23 March 2000                                  *
 *----------------------------------------------------------------*
 * This file contains implementations of the banded difference    *
 * quotient Jacobian-based preconditioner and solver routines for *
 * use with CVSpgmr.                                              *
 *                                                                *
 ******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "cvbandpre.h"
#include "cvode.h"
#include "llnltyps.h"
#include "nvector.h"
#include "llnlmath.h"
#include "band.h"

#define MIN_INC_MULT RCONST(1000.0)
#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)

/* Prototype for difference quotient Jacobian calculation routine */

static void CVBandPDQJac(integer N, integer mupper, integer mlower, BandMat J,
		         RhsFn f, void *f_data, real tn, N_Vector y,
		         N_Vector fy, N_Vector ewt, real h, real uround,
		         N_Vector ftemp, N_Vector ytemp);


/********************** Malloc and Free Functions **********************/

CVBandPreData CVBandPreAlloc(integer N, RhsFn f, void *f_data,
                             integer mu, integer ml)
{
  CVBandPreData pdata;
  integer mup, mlp, storagemu;

  pdata = (CVBandPreData) malloc(sizeof *pdata);  /* Allocate data memory */
  if (pdata == NULL) return(NULL);

  /* Load pointers and bandwidths into pdata block. */
  pdata->f = f;
  pdata->f_data = f_data;
  pdata->mu = mup = MIN( N-1, MAX(0,mu) );
  pdata->ml = mlp = MIN( N-1, MAX(0,ml) );

  /* Allocate memory for saved banded Jacobian approximation. */
  pdata->savedJ = BandAllocMat(N, mup, mlp, mup);
  if (pdata->savedJ == NULL) {
    free(pdata);
    return(NULL);
  }

  /* Allocate memory for banded preconditioner. */
  storagemu = MIN( N-1, mup + mlp);
  pdata->savedP = BandAllocMat(N, mup, mlp, storagemu);
  if (pdata->savedP == NULL) {
    BandFreeMat(pdata->savedJ);
    free(pdata);
    return(NULL);
  }

  /* Allocate memory for pivot array. */
  pdata->pivots = BandAllocPiv(N);
  if (pdata->savedJ == NULL) {
    BandFreeMat(pdata->savedP);
    BandFreeMat(pdata->savedJ);
    free(pdata);
    return(NULL);
  }

  return(pdata);
}

void CVBandPreFree(CVBandPreData pdata)
{
  BandFreeMat(pdata->savedJ);
  BandFreeMat(pdata->savedP);
  BandFreePiv(pdata->pivots);
  free(pdata);
}


/***************** Preconditioner setup and solve functions *******/


/* Readability Replacements */

#define f         (pdata->f)
#define f_data    (pdata->f_data)
#define mu        (pdata->mu)
#define ml        (pdata->ml)
#define pivots    (pdata->pivots)
#define savedJ    (pdata->savedJ)
#define savedP    (pdata->savedP)


/* Preconditioner setup routine CVBandPrecond. */

/******************************************************************
 * Together CVBandPrecond and CVBandPSolve use a banded           *
 * difference quotient Jacobian to create a preconditioner.       *
 * CVBandPrecond calculates a new J, if necessary, then           *
 * calculates P = I - gamma*J, and does an LU factorization of P. *
 *                                                                *
 * The parameters of CVBandPrecond are as follows:                *
 *                                                                *
 * N       is the length of all vector arguments.                 *
 *                                                                *
 * t       is the current value of the independent variable.      *
 *                                                                *
 * y       is the current value of the dependent variable vector, *
 *           namely the predicted value of y(t).                  *
 *                                                                *
 * fy      is the vector f(t,y).                                  *
 *                                                                *
 * jok     is an input flag indicating whether Jacobian-related   *
 *         data needs to be recomputed, as follows:               *
 *           jok == FALSE means recompute Jacobian-related data   *
 *                  from scratch.                                 *
 *           jok == TRUE  means that Jacobian data from the       *
 *                  previous Precond call will be reused          *
 *                  (with the current value of gamma).            *
 *         A CVBandPrecond call with jok == TRUE should only      *
 *         occur after a call with jok == FALSE.                  *
 *                                                                *
 * jcurPtr is a pointer to an output integer flag which is        *
 *         set by CVBandPrecond as follows:                       *
 *           *jcurPtr = TRUE if Jacobian data was recomputed.     *
 *           *jcurPtr = FALSE if Jacobian data was not recomputed,*
 *                      but saved data was reused.                *
 *                                                                *
 * gamma   is the scalar appearing in the Newton matrix.          *
 *                                                                *
 * ewt     is the error weight vector.                            *
 *                                                                *
 * h       is a tentative step size in t.                         *
 *                                                                *
 * uround  is the machine unit roundoff.                          *
 *                                                                *
 * nfePtr  is a pointer to the memory location containing the     *
 *           CVODE problem data nfe = number of calls to f.       *
 *           The routine calls f a total of ml+mu+1 times, so     *
 *           it increments *nfePtr by ml+mu+1.                    *
 *                                                                *
 * bp_data is a pointer to preconditoner data - the same as the   *
 *            bp_data parameter passed to CVSpgmr.                *
 *                                                                *
 * vtemp1, vtemp2, and vtemp3 are pointers to memory allocated    *
 *           for vectors of length N for work space.  This        *
 *           routine uses only vtemp1 and vtemp2.                 *
 *                                                                *
 *                                                                *
 * The value to be returned by the CVBandPrecond function is      *
 *   0  if successful, or                                         *
 *   1  if the band factorization failed.                         *
 ******************************************************************/

int CVBandPrecond(integer N, real t, N_Vector y, N_Vector fy,
                  boole jok, boole *jcurPtr, real gamma,
                  N_Vector ewt, real h, real uround,
                  long int *nfePtr, void *bp_data,
                  N_Vector vtemp1, N_Vector vtemp2,
                  N_Vector vtemp3)
{
  integer ier;
  CVBandPreData pdata;

  /* Assume matrix and pivots have already been allocated. */
  pdata = (CVBandPreData) bp_data;

  if (jok) {
    /* If jok = TRUE, use saved copy of J. */
    *jcurPtr = FALSE;
    BandCopy(savedJ, savedP, mu, ml);
  } else {
    /* If jok = FALSE, call CVBandPDQJac for new J value. */
    *jcurPtr = TRUE;
    BandZero(savedJ);
    CVBandPDQJac(N, mu, ml, savedJ, f, f_data, t, y, fy, ewt,
	         h, uround, vtemp1, vtemp2);
    BandCopy(savedJ, savedP, mu, ml);
    *nfePtr += MIN( N, ml + mu + 1 );
  }
  
  /* Scale and add I to get savedP = I - gamma*J. */
  BandScale(-gamma, savedP);
  BandAddI(savedP);
 
  /* Do LU factorization of matrix. */
  ier = BandFactor(savedP, pivots);
 
  /* Return 0 if the LU was complete; otherwise return 1. */
  if (ier > 0) return(1);
  return(0);
}


/* Preconditioner solve routine CVBandPSolve */

/******************************************************************
 * CVBandPSolve solves a linear system P z = r, where P is the    *
 * matrix computed by CVBandPrecond.                              *
 *                                                                *
 * The parameters of CVBandPSolve used here are as follows:       *
 *                                                                *
 * r       is the right-hand side vector of the linear system.    *
 *                                                                *
 * bp_data is a pointer to preconditioner data - the same as the  *
 *          bp_data parameter passed to CVSpgmr.                  *
 *                                                                *
 * z       is the output vector computed by CVBandPSolve.         *
 *                                                                *
 * The value returned by the CVBandPSolve function is always 0,   *
 * indicating success.                                            *
 *                                                                *
 ******************************************************************/

int CVBandPSolve(integer N, real t, N_Vector y, N_Vector fy,
                 N_Vector vtemp, real gamma, N_Vector ewt,
                 real delta, long int *nfePtr, N_Vector r,
                 int lr, void *bp_data, N_Vector z)
{
  CVBandPreData pdata;

  /* Assume matrix and pivots have already been allocated. */
  pdata = (CVBandPreData) bp_data;

  /* Copy r to z. */
  N_VScale(ONE, r, z);

  /* Do band backsolve on the vector z. */
  BandBacksolve(savedP, pivots, z);
 
  return(0);
}


#undef f
#undef f_data
#undef mu
#undef ml
#undef pivots
#undef savedJ
#undef savedP




/*************** CVBandPDQJac ****************************************

 This routine generates a banded difference quotient approximation to
 the Jacobian of f(t,y).  It assumes that a band matrix of type
 BandMat is stored column-wise, and that elements within each column
 are contiguous. This makes it possible to get the address of a column
 of J via the macro BAND_COL and to write a simple for loop to set
 each of the elements of a column in succession.

**********************************************************************/

static void CVBandPDQJac(integer N, integer mupper, integer mlower, BandMat J,
		         RhsFn f, void *f_data, real tn, N_Vector y,
		         N_Vector fy, N_Vector ewt, real h, real uround,
		         N_Vector ftemp, N_Vector ytemp)
{
  real    fnorm, minInc, inc, inc_inv, srur;
  integer group, i, j, width, ngroups, i1, i2;
  real *col_j, *ewt_data, *fy_data, *ftemp_data, *y_data, *ytemp_data;

  /* Obtain pointers to the data for ewt, fy, ftemp, y, ytemp. */
  ewt_data   = N_VDATA(ewt);
  fy_data    = N_VDATA(fy);
  ftemp_data = N_VDATA(ftemp);
  y_data     = N_VDATA(y);
  ytemp_data = N_VDATA(ytemp);

  /* Load ytemp with y = predicted y vector. */
  N_VScale(ONE, y, ytemp);

  /* Set minimum increment based on uround and norm of f. */
  srur = RSqrt(uround);
  fnorm = N_VWrmsNorm(fy, ewt);
  minInc = (fnorm != ZERO) ?
           (MIN_INC_MULT * ABS(h) * uround * N * fnorm) : ONE;

  /* Set bandwidth and number of column groups for band differencing. */
  width = mlower + mupper + 1;
  ngroups = MIN(width, N);
  
  for (group = 1; group <= ngroups; group++) {
    
    /* Increment all y_j in group. */
    for(j = group-1; j < N; j += width) {
      inc = MAX(srur*ABS(y_data[j]), minInc/ewt_data[j]);
      ytemp_data[j] += inc;
    }

    /* Evaluate f with incremented y. */
    f(N, tn, ytemp, ftemp, f_data);

    /* Restore ytemp, then form and load difference quotients. */
    for (j = group-1; j < N; j += width) {
      ytemp_data[j] = y_data[j];
      col_j = BAND_COL(J,j);
      inc = MAX(srur*ABS(y_data[j]), minInc/ewt_data[j]);
      inc_inv = ONE/inc;
      i1 = MAX(0, j-mupper);
      i2 = MIN(j+mlower, N-1);
      for (i=i1; i <= i2; i++)
	BAND_COL_ELEM(col_j,i,j) =
	  inc_inv * (ftemp_data[i] - fy_data[i]);
    }
  }
}
