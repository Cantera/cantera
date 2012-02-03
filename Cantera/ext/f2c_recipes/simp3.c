/* simp3.f -- translated by f2c (version 20031025).
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

/* Subroutine */ int simp3_(doublereal *a, integer *mp, integer *np, integer *
	i1, integer *k1, integer *ip, integer *kp)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer ii, kk;
    static doublereal piv;

    /* Parameter adjustments */
    a_dim1 = *mp;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    piv = 1. / a[*ip + 1 + (*kp + 1) * a_dim1];
    if (*i1 >= 0) {
	i__1 = *i1 + 1;
	for (ii = 1; ii <= i__1; ++ii) {
	    if (ii - 1 != *ip) {
		a[ii + (*kp + 1) * a_dim1] *= piv;
		i__2 = *k1 + 1;
		for (kk = 1; kk <= i__2; ++kk) {
		    if (kk - 1 != *kp) {
			a[ii + kk * a_dim1] -= a[*ip + 1 + kk * a_dim1] * a[
				ii + (*kp + 1) * a_dim1];
		    }
/* L11: */
		}
	    }
/* L12: */
	}
    }
    i__1 = *k1 + 1;
    for (kk = 1; kk <= i__1; ++kk) {
	if (kk - 1 != *kp) {
	    a[*ip + 1 + kk * a_dim1] = -a[*ip + 1 + kk * a_dim1] * piv;
	}
/* L13: */
    }
    a[*ip + 1 + (*kp + 1) * a_dim1] = piv;
    return 0;
} /* simp3_ */

#ifdef __cplusplus
	}
#endif
