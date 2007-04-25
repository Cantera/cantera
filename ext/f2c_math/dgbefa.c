/* dgbefa.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int dgbfa_(doublereal *abd, integer *lda, integer *n, 
	integer *ml, integer *mu, integer *ipvt, integer *info)
{
    /* System generated locals */
    integer abd_dim1, abd_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    integer i__, j, k, l, m;
    doublereal t;
    integer i0, j0, j1, lm, mm, ju, jz, kp1, nm1;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), daxpy_(integer *, doublereal *, doublereal *, integer 
	    *, doublereal *, integer *);
    extern integer idamax_(integer *, doublereal *, integer *);


/*     dgbfa factors a double precision band matrix by elimination. */

/*     dgbfa is usually called by dgbco, but it can be called */
/*     directly with a saving in time if  rcond  is not needed. */

/*     on entry */

/*        abd     double precision(lda, n) */
/*                contains the matrix in band storage.  the columns */
/*                of the matrix are stored in the columns of  abd  and */
/*                the diagonals of the matrix are stored in rows */
/*                ml+1 through 2*ml+mu+1 of  abd . */
/*                see the comments below for details. */

/*        lda     integer */
/*                the leading dimension of the array  abd . */
/*                lda must be .ge. 2*ml + mu + 1 . */

/*        n       integer */
/*                the order of the original matrix. */

/*        ml      integer */
/*                number of diagonals below the main diagonal. */
/*                0 .le. ml .lt. n . */

/*        mu      integer */
/*                number of diagonals above the main diagonal. */
/*                0 .le. mu .lt. n . */
/*                more efficient if  ml .le. mu . */
/*     on return */

/*        abd     an upper triangular matrix in band storage and */
/*                the multipliers which were used to obtain it. */
/*                the factorization can be written  a = l*u  where */
/*                l  is a product of permutation and unit lower */
/*                triangular matrices and  u  is upper triangular. */

/*        ipvt    integer(n) */
/*                an integer vector of pivot indices. */

/*        info    integer */
/*                = 0  normal value. */
/*                = k  if  u(k,k) .eq. 0.0 .  this is not an error */
/*                     condition for this subroutine, but it does */
/*                     indicate that dgbsl will divide by zero if */
/*                     called.  use  rcond  in dgbco for a reliable */
/*                     indication of singularity. */

/*     band storage */

/*           if  a  is a band matrix, the following program segment */
/*           will set up the input. */

/*                   ml = (band width below the diagonal) */
/*                   mu = (band width above the diagonal) */
/*                   m = ml + mu + 1 */
/*                   do 20 j = 1, n */
/*                      i1 = max0(1, j-mu) */
/*                      i2 = min0(n, j+ml) */
/*                      do 10 i = i1, i2 */
/*                         k = i - j + m */
/*                         abd(k,j) = a(i,j) */
/*                10    continue */
/*                20 continue */

/*           this uses rows  ml+1  through  2*ml+mu+1  of  abd . */
/*           in addition, the first  ml  rows in  abd  are used for */
/*           elements generated during the triangularization. */
/*           the total number of rows needed in  abd  is  2*ml+mu+1 . */
/*           the  ml+mu by ml+mu  upper left triangle and the */
/*           ml by ml  lower right triangle are not referenced. */

/*     linpack. this version dated 08/14/78 . */
/*     cleve moler, university of new mexico, argonne national lab. */

/*     subroutines and functions */

/*     blas daxpy,dscal,idamax */
/*     fortran max0,min0 */

/*     internal variables */



    /* Parameter adjustments */
    abd_dim1 = *lda;
    abd_offset = 1 + abd_dim1;
    abd -= abd_offset;
    --ipvt;

    /* Function Body */
    m = *ml + *mu + 1;
    *info = 0;

/*     zero initial fill-in columns */

    j0 = *mu + 2;
    j1 = min(*n,m) - 1;
    if (j1 < j0) {
	goto L30;
    }
    i__1 = j1;
    for (jz = j0; jz <= i__1; ++jz) {
	i0 = m + 1 - jz;
	i__2 = *ml;
	for (i__ = i0; i__ <= i__2; ++i__) {
	    abd[i__ + jz * abd_dim1] = 0.;
/* L10: */
	}
/* L20: */
    }
L30:
    jz = j1;
    ju = 0;

/*     gaussian elimination with partial pivoting */

    nm1 = *n - 1;
    if (nm1 < 1) {
	goto L130;
    }
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
	kp1 = k + 1;

/*        zero next fill-in column */

	++jz;
	if (jz > *n) {
	    goto L50;
	}
	if (*ml < 1) {
	    goto L50;
	}
	i__2 = *ml;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    abd[i__ + jz * abd_dim1] = 0.;
/* L40: */
	}
L50:

/*        find l = pivot index */

/* Computing MIN */
	i__2 = *ml, i__3 = *n - k;
	lm = min(i__2,i__3);
	i__2 = lm + 1;
	l = idamax_(&i__2, &abd[m + k * abd_dim1], &c__1) + m - 1;
	ipvt[k] = l + k - m;

/*        zero pivot implies this column already triangularized */

	if (abd[l + k * abd_dim1] == 0.) {
	    goto L100;
	}

/*           interchange if necessary */

	if (l == m) {
	    goto L60;
	}
	t = abd[l + k * abd_dim1];
	abd[l + k * abd_dim1] = abd[m + k * abd_dim1];
	abd[m + k * abd_dim1] = t;
L60:

/*           compute multipliers */

	t = -1. / abd[m + k * abd_dim1];
	dscal_(&lm, &t, &abd[m + 1 + k * abd_dim1], &c__1);

/*           row elimination with column indexing */

/* Computing MIN */
/* Computing MAX */
	i__3 = ju, i__4 = *mu + ipvt[k];
	i__2 = max(i__3,i__4);
	ju = min(i__2,*n);
	mm = m;
	if (ju < kp1) {
	    goto L90;
	}
	i__2 = ju;
	for (j = kp1; j <= i__2; ++j) {
	    --l;
	    --mm;
	    t = abd[l + j * abd_dim1];
	    if (l == mm) {
		goto L70;
	    }
	    abd[l + j * abd_dim1] = abd[mm + j * abd_dim1];
	    abd[mm + j * abd_dim1] = t;
L70:
	    daxpy_(&lm, &t, &abd[m + 1 + k * abd_dim1], &c__1, &abd[mm + 1 + 
		    j * abd_dim1], &c__1);
/* L80: */
	}
L90:
	goto L110;
L100:
	*info = k;
L110:
/* L120: */
	;
    }
L130:
    ipvt[*n] = *n;
    if (abd[m + *n * abd_dim1] == 0.) {
	*info = *n;
    }
    return 0;
} /* dgbfa_ */

