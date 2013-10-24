/* xgetua.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__5 = 5;
static integer c__0 = 0;
static logical c_false = FALSE_;

/* DECK XGETUA */
/* Subroutine */ int xgetua_(integer *iunita, integer *n)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__, index;
    extern integer j4save_(integer *, integer *, logical *);

/* ***BEGIN PROLOGUE  XGETUA */
/* ***PURPOSE  Return unit number(s) to which error messages are being */
/*            sent. */
/* ***LIBRARY   SLATEC (XERROR) */
/* ***CATEGORY  R3C */
/* ***TYPE      ALL (XGETUA-A) */
/* ***KEYWORDS  ERROR, XERROR */
/* ***AUTHOR  Jones, R. E., (SNLA) */
/* ***DESCRIPTION */

/*     Abstract */
/*        XGETUA may be called to determine the unit number or numbers */
/*        to which error messages are being sent. */
/*        These unit numbers may have been set by a call to XSETUN, */
/*        or a call to XSETUA, or may be a default value. */

/*     Description of Parameters */
/*      --Output-- */
/*        IUNIT - an array of one to five unit numbers, depending */
/*                on the value of N.  A value of zero refers to the */
/*                default unit, as defined by the I1MACH machine */
/*                constant routine.  Only IUNIT(1),...,IUNIT(N) are */
/*                defined by XGETUA.  The values of IUNIT(N+1),..., */
/*                IUNIT(5) are not defined (for N .LT. 5) or altered */
/*                in any way by XGETUA. */
/*        N     - the number of units to which copies of the */
/*                error messages are being sent.  N will be in the */
/*                range from 1 to 5. */

/* ***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC */
/*                 Error-handling Package, SAND82-0800, Sandia */
/*                 Laboratories, 1982. */
/* ***ROUTINES CALLED  J4SAVE */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790801  DATE WRITTEN */
/*   861211  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  XGETUA */
/* ***FIRST EXECUTABLE STATEMENT  XGETUA */
    /* Parameter adjustments */
    --iunita;

    /* Function Body */
    *n = j4save_(&c__5, &c__0, &c_false);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	index = i__ + 4;
	if (i__ == 1) {
	    index = 3;
	}
	iunita[i__] = j4save_(&index, &c__0, &c_false);
/* L30: */
    }
    return 0;
} /* xgetua_ */

