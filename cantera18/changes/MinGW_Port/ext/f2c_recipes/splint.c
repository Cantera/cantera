/* splint.f -- translated by f2c (version 20031025).
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

/* Subroutine */ int splint_(doublereal *xa, doublereal *ya, doublereal *y2a, 
	integer *n, doublereal *x, doublereal *y)
{
    /* System generated locals */
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    /* Subroutine */ int s_paus(char *, ftnlen);

    /* Local variables */
    static doublereal a, b, h__;
    static integer k, khi, klo;

    /* Parameter adjustments */
    --y2a;
    --ya;
    --xa;

    /* Function Body */
    klo = 1;
    khi = *n;
L1:
    if (khi - klo > 1) {
	k = (khi + klo) / 2;
	if (xa[k] > *x) {
	    khi = k;
	} else {
	    klo = k;
	}
	goto L1;
    }
    h__ = xa[khi] - xa[klo];
    if (h__ == (float)0.) {
	s_paus("Bad XA input.", (ftnlen)13);
    }
    a = (xa[khi] - *x) / h__;
    b = (*x - xa[klo]) / h__;
/* Computing 3rd power */
    d__1 = a;
/* Computing 3rd power */
    d__2 = b;
/* Computing 2nd power */
    d__3 = h__;
    *y = a * ya[klo] + b * ya[khi] + ((d__1 * (d__1 * d__1) - a) * y2a[klo] + 
	    (d__2 * (d__2 * d__2) - b) * y2a[khi]) * (d__3 * d__3) / 6.;
    return 0;
} /* splint_ */

#ifdef __cplusplus
	}
#endif
