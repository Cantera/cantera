/* simp2.f -- translated by f2c (version 20031025).
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

/* Subroutine */ int simp2_(doublereal *a, integer *m, integer *n, integer *
	mp, integer *np, integer *l2, integer *nl2, integer *ip, integer *kp, 
	doublereal *q1)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, k;
    static doublereal q, q0;
    static integer ii;
    static doublereal qp;

    /* Parameter adjustments */
    --l2;
    a_dim1 = *mp;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    *ip = 0;
    if (*nl2 < 1) {
	return 0;
    }
    i__1 = *nl2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (a[l2[i__] + 1 + (*kp + 1) * a_dim1] < -1e-6) {
	    goto L2;
	}
/* L11: */
    }
    return 0;
L2:
    *q1 = -a[l2[i__] + 1 + a_dim1] / a[l2[i__] + 1 + (*kp + 1) * a_dim1];
    *ip = l2[i__];
    if (i__ + 1 > *nl2) {
	return 0;
    }
    i__1 = *nl2;
    for (++i__; i__ <= i__1; ++i__) {
	ii = l2[i__];
	if (a[ii + 1 + (*kp + 1) * a_dim1] < -1e-6) {
	    q = -a[ii + 1 + a_dim1] / a[ii + 1 + (*kp + 1) * a_dim1];
	    if (q < *q1) {
		*ip = ii;
		*q1 = q;
	    } else if (q == *q1) {
		i__2 = *n;
		for (k = 1; k <= i__2; ++k) {
		    qp = -a[*ip + 1 + (k + 1) * a_dim1] / a[*ip + 1 + (*kp + 
			    1) * a_dim1];
		    q0 = -a[ii + 1 + (k + 1) * a_dim1] / a[ii + 1 + (*kp + 1) 
			    * a_dim1];
		    if (q0 != qp) {
			goto L6;
		    }
/* L12: */
		}
L6:
		if (q0 < qp) {
		    *ip = ii;
		}
	    }
	}
/* L13: */
    }
    return 0;
} /* simp2_ */

#ifdef __cplusplus
	}
#endif
