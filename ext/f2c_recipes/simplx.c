/* simplx.f -- translated by f2c (version 20031025).
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
#include <stdio.h>
    
/* Table of constant values */

static integer c__0 = 0;
static integer c__1 = 1;

/* Subroutine */ int simplx_(doublereal *a, integer *m, integer *n, integer *
	mp, integer *np, integer *m1, integer *m2, integer *m3, integer *
	icase, integer *izrov, integer *iposv)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Builtin functions */
    /* Subroutine */ int s_paus(char *, ftnlen);

    /* Local variables */
    static integer i__, k, l1[1000], l2[1000], l3[1000];
    static doublereal q1;
    static integer m12, kh, ip, ir, kp, is, nl1, nl2;
    static doublereal bmax;
    extern /* Subroutine */ int simp1_(doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *), simp2_(doublereal *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, doublereal *), 
	    simp3_(doublereal *, integer *, integer *, integer *, integer *, 
	    integer *, integer *);

    /* Parameter adjustments */
    --iposv;
    --izrov;
    a_dim1 = *mp;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    if (*m != *m1 + *m2 + *m3) {
        printf(" %d %d %d %d ", *m, *m1, *m2, *m3);
	s_paus("Bad input constraint counts.", (ftnlen)28);
    }
    nl1 = *n;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	l1[k - 1] = k;
	izrov[k] = k;
/* L11: */
    }
    nl2 = *m;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (a[i__ + 1 + a_dim1] < 0.) {
/*          write(*,*) 'The A matrix input to SIMPLX is invalid.' */
	    s_paus("Bad input tableau.", (ftnlen)18);
	}
	l2[i__ - 1] = i__;
	iposv[i__] = *n + i__;
/* L12: */
    }
    i__1 = *m2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l3[i__ - 1] = 1;
/* L13: */
    }
    ir = 0;
    if (*m2 + *m3 == 0) {
	goto L30;
    }
    ir = 1;
    i__1 = *n + 1;
    for (k = 1; k <= i__1; ++k) {
	q1 = (float)0.;
	i__2 = *m;
	for (i__ = *m1 + 1; i__ <= i__2; ++i__) {
	    q1 += a[i__ + 1 + k * a_dim1];
/* L14: */
	}
	a[*m + 2 + k * a_dim1] = -q1;
/* L15: */
    }
L10:
    i__1 = *m + 1;
    simp1_(&a[a_offset], mp, np, &i__1, l1, &nl1, &c__0, &kp, &bmax);
    if (bmax <= 1e-6 && a[*m + 2 + a_dim1] < -1e-6) {
	*icase = -1;
	return 0;
    } else if (bmax <= 1e-6 && a[*m + 2 + a_dim1] <= 1e-6) {
	m12 = *m1 + *m2 + 1;
	if (m12 <= *m) {
	    i__1 = *m;
	    for (ip = m12; ip <= i__1; ++ip) {
		if (iposv[ip] == ip + *n) {
		    simp1_(&a[a_offset], mp, np, &ip, l1, &nl1, &c__1, &kp, &
			    bmax);
		    if (bmax > (float)0.) {
			goto L1;
		    }
		}
/* L16: */
	    }
	}
	ir = 0;
	--m12;
	if (*m1 + 1 > m12) {
	    goto L30;
	}
	i__1 = m12;
	for (i__ = *m1 + 1; i__ <= i__1; ++i__) {
	    if (l3[i__ - *m1 - 1] == 1) {
		i__2 = *n + 1;
		for (k = 1; k <= i__2; ++k) {
		    a[i__ + 1 + k * a_dim1] = -a[i__ + 1 + k * a_dim1];
/* L17: */
		}
	    }
/* L18: */
	}
	goto L30;
    }
    simp2_(&a[a_offset], m, n, mp, np, l2, &nl2, &ip, &kp, &q1);
    if (ip == 0) {
	*icase = -1;
	return 0;
    }
L1:
    i__1 = *m + 1;
    simp3_(&a[a_offset], mp, np, &i__1, n, &ip, &kp);
    if (iposv[ip] >= *n + *m1 + *m2 + 1) {
	i__1 = nl1;
	for (k = 1; k <= i__1; ++k) {
	    if (l1[k - 1] == kp) {
		goto L2;
	    }
/* L19: */
	}
L2:
	--nl1;
	i__1 = nl1;
	for (is = k; is <= i__1; ++is) {
	    l1[is - 1] = l1[is];
/* L21: */
	}
    } else {
	if (iposv[ip] < *n + *m1 + 1) {
	    goto L20;
	}
	kh = iposv[ip] - *m1 - *n;
	if (l3[kh - 1] == 0) {
	    goto L20;
	}
	l3[kh - 1] = 0;
    }
    a[*m + 2 + (kp + 1) * a_dim1] += (float)1.;
    i__1 = *m + 2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	a[i__ + (kp + 1) * a_dim1] = -a[i__ + (kp + 1) * a_dim1];
/* L22: */
    }
L20:
    is = izrov[kp];
    izrov[kp] = iposv[ip];
    iposv[ip] = is;
    if (ir != 0) {
	goto L10;
    }
L30:
    simp1_(&a[a_offset], mp, np, &c__0, l1, &nl1, &c__0, &kp, &bmax);
    if (bmax <= (float)0.) {
	*icase = 0;
	return 0;
    }
    simp2_(&a[a_offset], m, n, mp, np, l2, &nl2, &ip, &kp, &q1);
    if (ip == 0) {
	*icase = 1;
	return 0;
    }
    simp3_(&a[a_offset], mp, np, m, n, &ip, &kp);
    goto L20;
} /* simplx_ */

#ifdef __cplusplus
	}
#endif
