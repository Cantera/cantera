/* fdump.f -- translated by f2c (version 20031025).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#ifdef _cpluscplus
extern "C" {
#endif
#include "f2c.h"

/* DECK FDUMP */
/* Subroutine */ int fdump_(void)
{
/* ***BEGIN PROLOGUE  FDUMP */
/* ***PURPOSE  Symbolic dump (should be locally written). */
/* ***LIBRARY   SLATEC (XERMSG) */
/* ***CATEGORY  R3 */
/* ***TYPE      ALL (FDUMP-A) */
/* ***KEYWORDS  ERROR, XERMSG */
/* ***AUTHOR  Jones, R. E., (SNLA) */
/* ***DESCRIPTION */

/*        ***Note*** Machine Dependent Routine */
/*        FDUMP is intended to be replaced by a locally written */
/*        version which produces a symbolic dump.  Failing this, */
/*        it should be replaced by a version which prints the */
/*        subprogram nesting list.  Note that this dump must be */
/*        printed on each of up to five files, as indicated by the */
/*        XGETUA routine.  See XSETUA and XGETUA for details. */

/*     Written by Ron Jones, with SLATEC Common Math Library Subcommittee */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790801  DATE WRITTEN */
/*   861211  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  FDUMP */
/* ***FIRST EXECUTABLE STATEMENT  FDUMP */
    return 0;
} /* fdump_ */
    
    
/* integer isamax_(integer *n, real *sx, integer *incx) */
/* { */
/*     /\* System generated locals *\/ */
/*     integer ret_val, i__1; */
/*     real r__1; */

/*     /\* Local variables *\/ */
/*     static integer i__, ix; */
/*     static real smax; */


/* /\*     finds the index of element having max. absolute value. *\/ */
/* /\*     jack dongarra, linpack, 3/11/78. *\/ */
/* /\*     modified 3/93 to return if incx .le. 0. *\/ */


/*     /\* Parameter adjustments *\/ */
/*     --sx; */

/*     /\* Function Body *\/ */
/*     ret_val = 0; */
/*     if (*n < 1 || *incx <= 0) { */
/* 	return ret_val; */
/*     } */
/*     ret_val = 1; */
/*     if (*n == 1) { */
/* 	return ret_val; */
/*     } */
/*     if (*incx == 1) { */
/* 	goto L20; */
/*     } */

/* /\*        code for increment not equal to 1 *\/ */

/*     ix = 1; */
/*     smax = dabs(sx[1]); */
/*     ix += *incx; */
/*     i__1 = *n; */
/*     for (i__ = 2; i__ <= i__1; ++i__) { */
/* 	if ((r__1 = sx[ix], dabs(r__1)) <= smax) { */
/* 	    goto L5; */
/* 	} */
/* 	ret_val = i__; */
/* 	smax = (r__1 = sx[ix], dabs(r__1)); */
/* L5: */
/* 	ix += *incx; */
/* /\* L10: *\/ */
/*     } */
/*     return ret_val; */

/* /\*        code for increment equal to 1 *\/ */

/* L20: */
/*     smax = dabs(sx[1]); */
/*     i__1 = *n; */
/*     for (i__ = 2; i__ <= i__1; ++i__) { */
/* 	if ((r__1 = sx[i__], dabs(r__1)) <= smax) { */
/* 	    goto L30; */
/* 	} */
/* 	ret_val = i__; */
/* 	smax = (r__1 = sx[i__], dabs(r__1)); */
/* L30: */
/* 	; */
/*     } */
/*     return ret_val; */
/* } /\* isamax_ *\/ */
 
#ifdef _cpluscplus
}
#endif
