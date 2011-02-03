/* simp1.f -- translated by f2c (version 20031025).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#ifdef __cplusplus
extern "C" {
#endif
#include "f2c.h"

/* Subroutine */ int simp1_(doublereal *a, integer *mp, integer *np, integer *
	mm, integer *ll, integer *nll, integer *iabf, integer *kp, doublereal 
	*bmax)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;
    doublereal d__1;

    /* Local variables */
    static integer k;
    static doublereal test;

    /* Parameter adjustments */
    --ll;
    a_dim1 = *mp;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    *kp = ll[1];
    *bmax = a[*mm + 1 + (*kp + 1) * a_dim1];
    if (*nll < 2) {
	return 0;
    }
    i__1 = *nll;
    for (k = 2; k <= i__1; ++k) {
	if (*iabf == 0) {
	    test = a[*mm + 1 + (ll[k] + 1) * a_dim1] - *bmax;
	} else {
	    test = (d__1 = a[*mm + 1 + (ll[k] + 1) * a_dim1], abs(d__1)) - 
		    abs(*bmax);
	}
	if (test > (float)0.) {
	    *bmax = a[*mm + 1 + (ll[k] + 1) * a_dim1];
	    *kp = ll[k];
	}
/* L11: */
    }
    return 0;
} /* simp1_ */

#ifdef __cplusplus
	}
#endif
