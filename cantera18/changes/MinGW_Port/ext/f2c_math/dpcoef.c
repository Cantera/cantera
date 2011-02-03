/* dpcoef.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* DECK DPCOEF */
/* Subroutine */ int dpcoef_(integer *l, doublereal *c__, doublereal *tc, 
	doublereal *a)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__, ll, nr;
    doublereal fac;
    integer new__, llp1, llp2;
    doublereal save;
    extern /* Subroutine */ int dp1vlu_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);

/* ***BEGIN PROLOGUE  DPCOEF */
/* ***PURPOSE  Convert the DPOLFT coefficients to Taylor series form. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  K1A1A2 */
/* ***TYPE      DOUBLE PRECISION (PCOEF-S, DPCOEF-D) */
/* ***KEYWORDS  CURVE FITTING, DATA FITTING, LEAST SQUARES, POLYNOMIAL FIT */
/* ***AUTHOR  Shampine, L. F., (SNLA) */
/*           Davenport, S. M., (SNLA) */
/* ***DESCRIPTION */

/*     Abstract */

/*     DPOLFT  computes the least squares polynomial fit of degree  L  as */
/*     a sum of orthogonal polynomials.  DPCOEF  changes this fit to its */
/*     Taylor expansion about any point  C , i.e. writes the polynomial */
/*     as a sum of powers of (X-C).  Taking  C=0.  gives the polynomial */
/*     in powers of X, but a suitable non-zero  C  often leads to */
/*     polynomials which are better scaled and more accurately evaluated. */

/*     The parameters for  DPCOEF  are */

/*     INPUT -- All TYPE REAL variables are DOUBLE PRECISION */
/*         L -      Indicates the degree of polynomial to be changed to */
/*                  its Taylor expansion.  To obtain the Taylor */
/*                  coefficients in reverse order, input  L  as the */
/*                  negative of the degree desired.  The absolute value */
/*                  of L  must be less than or equal to NDEG, the highest */
/*                  degree polynomial fitted by  DPOLFT . */
/*         C -      The point about which the Taylor expansion is to be */
/*                  made. */
/*         A -      Work and output array containing values from last */
/*                  call to  DPOLFT . */

/*     OUTPUT -- All TYPE REAL variables are DOUBLE PRECISION */
/*         TC -     Vector containing the first LL+1 Taylor coefficients */
/*                  where LL=ABS(L).  If  L.GT.0 , the coefficients are */
/*                  in the usual Taylor series order, i.e. */
/*                    P(X) = TC(1) + TC(2)*(X-C) + ... + TC(N+1)*(X-C)**N */
/*                  If L .LT. 0, the coefficients are in reverse order, */
/*                  i.e. */
/*                    P(X) = TC(1)*(X-C)**N + ... + TC(N)*(X-C) + TC(N+1) */

/* ***REFERENCES  L. F. Shampine, S. M. Davenport and R. E. Huddleston, */
/*                 Curve fitting by polynomials in one variable, Report */
/*                 SLA-74-0270, Sandia Laboratories, June 1974. */
/* ***ROUTINES CALLED  DP1VLU */
/* ***REVISION HISTORY  (YYMMDD) */
/*   740601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891006  Cosmetic changes to prologue.  (WRB) */
/*   891006  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DPCOEF */

/* ***FIRST EXECUTABLE STATEMENT  DPCOEF */
    /* Parameter adjustments */
    --a;
    --tc;

    /* Function Body */
    ll = abs(*l);
    llp1 = ll + 1;
    dp1vlu_(&ll, &ll, c__, &tc[1], &tc[2], &a[1]);
    if (ll < 2) {
	goto L2;
    }
    fac = 1.;
    i__1 = llp1;
    for (i__ = 3; i__ <= i__1; ++i__) {
	fac *= i__ - 1;
/* L1: */
	tc[i__] /= fac;
    }
L2:
    if (*l >= 0) {
	goto L4;
    }
    nr = llp1 / 2;
    llp2 = ll + 2;
    i__1 = nr;
    for (i__ = 1; i__ <= i__1; ++i__) {
	save = tc[i__];
	new__ = llp2 - i__;
	tc[i__] = tc[new__];
/* L3: */
	tc[new__] = save;
    }
L4:
    return 0;
} /* dpcoef_ */

