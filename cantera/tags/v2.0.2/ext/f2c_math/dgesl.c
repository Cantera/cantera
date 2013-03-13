/* dgesl.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int dgesl_(doublereal *a, integer *lda, integer *n, integer *
	ipvt, doublereal *b, integer *job)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    integer k, l;
    doublereal t;
    integer kb, nm1;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);


/*     dgesl solves the double precision system */
/*     a * x = b  or  trans(a) * x = b */
/*     using the factors computed by dgeco or dgefa. */

/*     on entry */

/*        a       double precision(lda, n) */
/*                the output from dgeco or dgefa. */

/*        lda     integer */
/*                the leading dimension of the array  a . */

/*        n       integer */
/*                the order of the matrix  a . */

/*        ipvt    integer(n) */
/*                the pivot vector from dgeco or dgefa. */

/*        b       double precision(n) */
/*                the right hand side vector. */

/*        job     integer */
/*                = 0         to solve  a*x = b , */
/*                = nonzero   to solve  trans(a)*x = b  where */
/*                            trans(a)  is the transpose. */

/*     on return */

/*        b       the solution vector  x . */

/*     error condition */

/*        a division by zero will occur if the input factor contains a */
/*        zero on the diagonal.  technically this indicates singularity */
/*        but it is often caused by improper arguments or improper */
/*        setting of lda .  it will not occur if the subroutines are */
/*        called correctly and if dgeco has set rcond .gt. 0.0 */
/*        or dgefa has set info .eq. 0 . */

/*     to compute  inverse(a) * c  where  c  is a matrix */
/*     with  p  columns */
/*           call dgeco(a,lda,n,ipvt,rcond,z) */
/*           if (rcond is too small) go to ... */
/*           do 10 j = 1, p */
/*              call dgesl(a,lda,n,ipvt,c(1,j),0) */
/*        10 continue */

/*     linpack. this version dated 08/14/78 . */
/*     cleve moler, university of new mexico, argonne national lab. */

/*     subroutines and functions */

/*     blas daxpy,ddot */

/*     internal variables */


    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipvt;
    --b;

    /* Function Body */
    nm1 = *n - 1;
    if (*job != 0) {
	goto L50;
    }

/*        job = 0 , solve  a * x = b */
/*        first solve  l*y = b */

    if (nm1 < 1) {
	goto L30;
    }
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
	l = ipvt[k];
	t = b[l];
	if (l == k) {
	    goto L10;
	}
	b[l] = b[k];
	b[k] = t;
L10:
	i__2 = *n - k;
	daxpy_(&i__2, &t, &a[k + 1 + k * a_dim1], &c__1, &b[k + 1], &c__1);
/* L20: */
    }
L30:

/*        now solve  u*x = y */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	b[k] /= a[k + k * a_dim1];
	t = -b[k];
	i__2 = k - 1;
	daxpy_(&i__2, &t, &a[k * a_dim1 + 1], &c__1, &b[1], &c__1);
/* L40: */
    }
    goto L100;
L50:

/*        job = nonzero, solve  trans(a) * x = b */
/*        first solve  trans(u)*y = b */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k - 1;
	t = ddot_(&i__2, &a[k * a_dim1 + 1], &c__1, &b[1], &c__1);
	b[k] = (b[k] - t) / a[k + k * a_dim1];
/* L60: */
    }

/*        now solve trans(l)*x = y */

    if (nm1 < 1) {
	goto L90;
    }
    i__1 = nm1;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n - kb;
	i__2 = *n - k;
	b[k] += ddot_(&i__2, &a[k + 1 + k * a_dim1], &c__1, &b[k + 1], &c__1);
	l = ipvt[k];
	if (l == k) {
	    goto L70;
	}
	t = b[l];
	b[l] = b[k];
	b[k] = t;
L70:
/* L80: */
	;
    }
L90:
L100:
    return 0;
} /* dgesl_ */

