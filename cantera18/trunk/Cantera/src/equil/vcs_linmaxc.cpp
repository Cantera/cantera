/*  $Author: hkmoffa $
 *  $Date: 2008/12/17 16:34:18 $
 *  $Revision: 1.5 $
 */

/*
 * Copywrite (2005) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */


#ifdef hpux
#define dbocls_ dbocls
#endif
#ifdef DEBUG_MODE
//extern int vcs_debug_print_lvl;
#endif

#include "vcs_internal.h"
#include "vcs_solve.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>

extern "C" void dbocls_(double *W, int *MDW, int *MCON, int *MROWS,
			int *NCOLS,
			double *BL, double *BU, int *IND, int *IOPT, 
			double *X, double *RNORMC, double *RNORM, 
			int *MODE, double *RW, int *IW);

/*****************************************************************************
* This whole program is a wrapper for the slatec routine, DBOCLS()
* DBOCLS solves a bounded and constrained least squares problem.
******************************************************************************/
namespace VCSnonideal {

#ifdef ALTLINPROG
#else
int linprogmax(double *XMOLES, double *CC, double *AX, double *BB, 
	       int NE, int M, int NE0)
   
   /*-----------------------------------------------------------------------
   * Find XMOLES(I), i = 1, M such that
   *   Maximize CC dot W,  subject to the NE constraints:
   *
   *             [AX] [XMOLES] = [BB]
   *             and XMOLES(i) > 0
   *
   * Input
   * ---------
   *   AX(NE, M) - matrix of constraints  AX(I,J) = ax(i + j*ne0)
   *   BB(NE)    - contraint values
   *   CC(M)     -  Vector of "Good Values" to maximize 
   *
   * Output
   * ---------
   *   XMOLES(M) - optimal value of XMOLES()
   *----------------------------------------------------------------------*/
{
   int MROWS, MCON, NCOLS, NX, NI, MDW, i, j, MODE;
   double sum, F[1], RNORMC, RNORM, *W, *BL, *BU, *RW, *X;
   int    *IND, *IW, *IOPT;

   MROWS = 1;
   MCON = NE;
   NCOLS = M;
   MDW = MCON + NCOLS;
   NX = 0;
   NI = 0;
   
   sum = 0.0;
   for (i = 0; i < M; i++) {
      sum += fabs(CC[i]);
   }
   F[0] = sum * 1000.;
   if (F[0] <= 0.0) F[0] = 1000.;
   
   BL = (double *) malloc(2*(NCOLS+MCON)            * sizeof(double));
   BU = BL + (NCOLS+MCON);
   IND =   (int *) malloc((NCOLS+MCON)              * sizeof(int));
   RW = (double *) malloc((6*NCOLS + 5*MCON)        * sizeof(double));
   IW =    (int *) malloc((2*NCOLS + 2*MCON)        * sizeof(int));
   IOPT =  (int *) malloc((17 + NI)                 * sizeof(int));
   X =  (double *) malloc((2*(NCOLS+MCON) + 2 + NX) * sizeof(double));
   W =  (double *) malloc((MDW*(NCOLS+MCON+1))      * sizeof(double));
   if (W == NULL) {
     plogf("linproxmax ERROR: can not malloc memory of size %d bytes\n",
	    (int) ((MDW*(NCOLS+MCON+1)) * sizeof(double)));
     if (BL != NULL)   free((void *) BL);
     if (IND != NULL)  free((void *) IND);
     if (RW != NULL)   free((void *) RW);
     if (IW != NULL)   free((void *) IW);
     if (IOPT != NULL) free((void *) IOPT);
     if (W != NULL)    free((void *) W);
     return -1;
   }
   for (j = 0; j < MCON; j++) {
      for (i = 0; i < NCOLS; i++) {
	 W[j + i*MDW] = AX[j + i*NE0];
      } 
   }
   for (i = 0; i < NCOLS; i++) {
      W[MCON + i*MDW] = CC[i];
   }
   W[MCON + (NCOLS)*MDW] = F[0];
   IOPT[0] = 99;
   
   for (j = 0; j < NCOLS; j++) {
      IND[j] = 1;
      BL[j]  = 0.0;
      BU[j]  = 1.0e200;
   }
   for (j = 0; j < MCON; j++) {
      IND[j + NCOLS] = 3;
      BL[j + NCOLS]  = BB[j];
      BU[j + NCOLS]  = BL[j + NCOLS];
   }
   
   
   dbocls_(W, &MDW, &MCON, &MROWS, &NCOLS, BL, BU, IND, IOPT, 
	   X, &RNORMC, &RNORM, &MODE, RW, IW);
   if (MODE != 0) {
      plogf("Return from DBOCLS was not normal, MODE = %d\n", MODE);
      plogf("       refer to subroutine DBOCLS for resolution\n");
      plogf("       RNORMC = %g\n", RNORMC);
   }

   for (j = 0; j < NCOLS; j++) {
      XMOLES[j] = X[j];
   }
#ifdef DEBUG_MODE
   //sum = 0.0;
   //for (j = 0; j < NCOLS; j++) {
   //   sum += XMOLES[j] * CC[j];
   //}
   //if (vcs_debug_print_lvl >= 2) {
   //  plogf(" -- linmaxc: Final Maximized Value = %g\n", sum);
   //}
#endif
   
   free((void *)W);
   free((void *)BL);
   free((void *)IND);
   free((void *)RW);
   free((void *)IW);
   free((void *)IOPT);
   free((void *)X);
   
   return 0;
}
#endif
}
