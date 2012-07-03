/* pcoef.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;

/* DECK PCOEF */
/* Subroutine */ int pcoef_(integer *l, real *c__, real *tc, real *a)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__, ll, nr;
    real fac;
    integer new__, llp1, llp2;
    real save;
    extern /* Subroutine */ int pvalue_(integer *, integer *, real *, real *, 
	    real *, real *);

/* ***BEGIN PROLOGUE  PCOEF */
/* ***PURPOSE  Convert the POLFIT coefficients to Taylor series form. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  K1A1A2 */
/* ***TYPE      SINGLE PRECISION (PCOEF-S, DPCOEF-D) */
/* ***KEYWORDS  CURVE FITTING, DATA FITTING, LEAST SQUARES, POLYNOMIAL FIT */
/* ***AUTHOR  Shampine, L. F., (SNLA) */
/*           Davenport, S. M., (SNLA) */
/* ***DESCRIPTION */

/*     Written BY L. F. Shampine and S. M. Davenport. */

/*     Abstract */

/*     POLFIT  computes the least squares polynomial fit of degree  L  as */
/*     a sum of orthogonal polynomials.  PCOEF  changes this fit to its */
/*     Taylor expansion about any point  C , i.e. writes the polynomial */
/*     as a sum of powers of (X-C).  Taking  C=0.  gives the polynomial */
/*     in powers of X, but a suitable non-zero  C  often leads to */
/*     polynomials which are better scaled and more accurately evaluated. */

/*     The parameters for  PCOEF  are */

/*     INPUT -- */
/*         L -      Indicates the degree of polynomial to be changed to */
/*                  its Taylor expansion.  To obtain the Taylor */
/*                  coefficients in reverse order, input  L  as the */
/*                  negative of the degree desired.  The absolute value */
/*                  of L  must be less than or equal to NDEG, the highest */
/*                  degree polynomial fitted by  POLFIT . */
/*         C -      The point about which the Taylor expansion is to be */
/*                  made. */
/*         A -      Work and output array containing values from last */
/*                  call to  POLFIT . */

/*     OUTPUT -- */
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
/* ***ROUTINES CALLED  PVALUE */
/* ***REVISION HISTORY  (YYMMDD) */
/*   740601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  PCOEF */

/* ***FIRST EXECUTABLE STATEMENT  PCOEF */
    /* Parameter adjustments */
    --a;
    --tc;

    /* Function Body */
    ll = abs(*l);
    llp1 = ll + 1;
    pvalue_(&ll, &ll, c__, &tc[1], &tc[2], &a[1]);
    if (ll < 2) {
	goto L2;
    }
    fac = 1.f;
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
} /* pcoef_ */

/* $$$ */
/* $$$      subroutine  dscal(n,da,dx,incx) */
/* $$$c */
/* $$$c     scales a vector by a constant. */
/* $$$c     uses unrolled loops for increment equal to one. */
/* $$$c     jack dongarra, linpack, 3/11/78. */
/* $$$c     modified 3/93 to return if incx .le. 0. */
/* $$$c */
/* $$$      double precision da,dx(1) */
/* $$$      integer i,incx,m,mp1,n,nincx */
/* $$$c */
/* $$$      if( n.le.0 .or. incx.le.0 )return */
/* $$$      if(incx.eq.1)go to 20 */
/* $$$c */
/* $$$c        code for increment not equal to 1 */
/* $$$c */
/* $$$      nincx = n*incx */
/* $$$      do 10 i = 1,nincx,incx */
/* $$$        dx(i) = da*dx(i) */
/* $$$   10 continue */
/* $$$      return */
/* $$$c */
/* $$$c        code for increment equal to 1 */
/* $$$c */
/* $$$c */
/* $$$c        clean-up loop */
/* $$$c */
/* $$$   20 m = mod(n,5) */
/* $$$      if( m .eq. 0 ) go to 40 */
/* $$$      do 30 i = 1,m */
/* $$$        dx(i) = da*dx(i) */
/* $$$   30 continue */
/* $$$      if( n .lt. 5 ) return */
/* $$$   40 mp1 = m + 1 */
/* $$$      do 50 i = mp1,n,5 */
/* $$$        dx(i) = da*dx(i) */
/* $$$        dx(i + 1) = da*dx(i + 1) */
/* $$$        dx(i + 2) = da*dx(i + 2) */
/* $$$        dx(i + 3) = da*dx(i + 3) */
/* $$$        dx(i + 4) = da*dx(i + 4) */
/* $$$   50 continue */
/* $$$      return */
/* $$$      end */
/* Subroutine */ int dgbco_(doublereal *abd, integer *lda, integer *n, 
	integer *ml, integer *mu, integer *ipvt, doublereal *rcond, 
	doublereal *z__)
{
    /* System generated locals */
    integer abd_dim1, abd_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);

    /* Local variables */
    integer j, k, l, m;
    doublereal s, t;
    integer kb, la;
    doublereal ek;
    integer lm, mm, is, ju;
    doublereal sm, wk;
    integer lz, kp1;
    doublereal wkm;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    integer info;
    extern /* Subroutine */ int dgbfa_(doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, integer *), dscal_(integer *, 
	    doublereal *, doublereal *, integer *);
    extern doublereal dasum_(integer *, doublereal *, integer *);
    doublereal anorm;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    doublereal ynorm;


/*     dgbco factors a double precision band matrix by gaussian */
/*     elimination and estimates the condition of the matrix. */

/*     if  rcond  is not needed, dgbfa is slightly faster. */
/*     to solve  a*x = b , follow dgbco by dgbsl. */
/*     to compute  inverse(a)*c , follow dgbco by dgbsl. */
/*     to compute  determinant(a) , follow dgbco by dgbdi. */

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

/*        rcond   double precision */
/*                an estimate of the reciprocal condition of  a . */
/*                for the system  a*x = b , relative perturbations */
/*                in  a  and  b  of size  epsilon  may cause */
/*                relative perturbations in  x  of size  epsilon/rcond . */
/*                if  rcond  is so small that the logical expression */
/*                           1.0 + rcond .eq. 1.0 */
/*                is true, then  a  may be singular to working */
/*                precision.  in particular,  rcond  is zero  if */
/*                exact singularity is detected or the estimate */
/*                underflows. */

/*        z       double precision(n) */
/*                a work vector whose contents are usually unimportant. */
/*                if  a  is close to a singular matrix, then  z  is */
/*                an approximate null vector in the sense that */
/*                norm(a*z) = rcond*norm(a)*norm(z) . */

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

/*     example..  if the original matrix is */

/*           11 12 13  0  0  0 */
/*           21 22 23 24  0  0 */
/*            0 32 33 34 35  0 */
/*            0  0 43 44 45 46 */
/*            0  0  0 54 55 56 */
/*            0  0  0  0 65 66 */

/*      then  n = 6, ml = 1, mu = 2, lda .ge. 5  and abd should contain */

/*            *  *  *  +  +  +  , * = not used */
/*            *  * 13 24 35 46  , + = used for pivoting */
/*            * 12 23 34 45 56 */
/*           11 22 33 44 55 66 */
/*           21 32 43 54 65  * */

/*     linpack. this version dated 08/14/78 . */
/*     cleve moler, university of new mexico, argonne national lab. */

/*     subroutines and functions */

/*     linpack dgbfa */
/*     blas daxpy,ddot,dscal,dasum */
/*     fortran dabs,dmax1,max0,min0,dsign */

/*     internal variables */



/*     compute 1-norm of a */

    /* Parameter adjustments */
    abd_dim1 = *lda;
    abd_offset = 1 + abd_dim1;
    abd -= abd_offset;
    --ipvt;
    --z__;

    /* Function Body */
    anorm = 0.;
    l = *ml + 1;
    is = l + *mu;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	d__1 = anorm, d__2 = dasum_(&l, &abd[is + j * abd_dim1], &c__1);
	anorm = max(d__1,d__2);
	if (is > *ml + 1) {
	    --is;
	}
	if (j <= *mu) {
	    ++l;
	}
	if (j >= *n - *ml) {
	    --l;
	}
/* L10: */
    }

/*     factor */

    dgbfa_(&abd[abd_offset], lda, n, ml, mu, &ipvt[1], &info);

/*     rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) . */
/*     estimate = norm(z)/norm(y) where  a*z = y  and  trans(a)*y = e . */
/*     trans(a)  is the transpose of a .  the components of  e  are */
/*     chosen to cause maximum local growth in the elements of w  where */
/*     trans(u)*w = e .  the vectors are frequently rescaled to avoid */
/*     overflow. */

/*     solve trans(u)*w = e */

    ek = 1.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	z__[j] = 0.;
/* L20: */
    }
    m = *ml + *mu + 1;
    ju = 0;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	if (z__[k] != 0.) {
	    d__1 = -z__[k];
	    ek = d_sign(&ek, &d__1);
	}
	if ((d__1 = ek - z__[k], abs(d__1)) <= (d__2 = abd[m + k * abd_dim1], 
		abs(d__2))) {
	    goto L30;
	}
	s = (d__1 = abd[m + k * abd_dim1], abs(d__1)) / (d__2 = ek - z__[k], 
		abs(d__2));
	dscal_(n, &s, &z__[1], &c__1);
	ek = s * ek;
L30:
	wk = ek - z__[k];
	wkm = -ek - z__[k];
	s = abs(wk);
	sm = abs(wkm);
	if (abd[m + k * abd_dim1] == 0.) {
	    goto L40;
	}
	wk /= abd[m + k * abd_dim1];
	wkm /= abd[m + k * abd_dim1];
	goto L50;
L40:
	wk = 1.;
	wkm = 1.;
L50:
	kp1 = k + 1;
/* Computing MIN */
/* Computing MAX */
	i__3 = ju, i__4 = *mu + ipvt[k];
	i__2 = max(i__3,i__4);
	ju = min(i__2,*n);
	mm = m;
	if (kp1 > ju) {
	    goto L90;
	}
	i__2 = ju;
	for (j = kp1; j <= i__2; ++j) {
	    --mm;
	    sm += (d__1 = z__[j] + wkm * abd[mm + j * abd_dim1], abs(d__1));
	    z__[j] += wk * abd[mm + j * abd_dim1];
	    s += (d__1 = z__[j], abs(d__1));
/* L60: */
	}
	if (s >= sm) {
	    goto L80;
	}
	t = wkm - wk;
	wk = wkm;
	mm = m;
	i__2 = ju;
	for (j = kp1; j <= i__2; ++j) {
	    --mm;
	    z__[j] += t * abd[mm + j * abd_dim1];
/* L70: */
	}
L80:
L90:
	z__[k] = wk;
/* L100: */
    }
    s = 1. / dasum_(n, &z__[1], &c__1);
    dscal_(n, &s, &z__[1], &c__1);

/*     solve trans(l)*y = w */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
/* Computing MIN */
	i__2 = *ml, i__3 = *n - k;
	lm = min(i__2,i__3);
	if (k < *n) {
	    z__[k] += ddot_(&lm, &abd[m + 1 + k * abd_dim1], &c__1, &z__[k + 
		    1], &c__1);
	}
	if ((d__1 = z__[k], abs(d__1)) <= 1.) {
	    goto L110;
	}
	s = 1. / (d__1 = z__[k], abs(d__1));
	dscal_(n, &s, &z__[1], &c__1);
L110:
	l = ipvt[k];
	t = z__[l];
	z__[l] = z__[k];
	z__[k] = t;
/* L120: */
    }
    s = 1. / dasum_(n, &z__[1], &c__1);
    dscal_(n, &s, &z__[1], &c__1);

    ynorm = 1.;

/*     solve l*v = y */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	l = ipvt[k];
	t = z__[l];
	z__[l] = z__[k];
	z__[k] = t;
/* Computing MIN */
	i__2 = *ml, i__3 = *n - k;
	lm = min(i__2,i__3);
	if (k < *n) {
	    daxpy_(&lm, &t, &abd[m + 1 + k * abd_dim1], &c__1, &z__[k + 1], &
		    c__1);
	}
	if ((d__1 = z__[k], abs(d__1)) <= 1.) {
	    goto L130;
	}
	s = 1. / (d__1 = z__[k], abs(d__1));
	dscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L130:
/* L140: */
	;
    }
    s = 1. / dasum_(n, &z__[1], &c__1);
    dscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

/*     solve  u*z = w */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	if ((d__1 = z__[k], abs(d__1)) <= (d__2 = abd[m + k * abd_dim1], abs(
		d__2))) {
	    goto L150;
	}
	s = (d__1 = abd[m + k * abd_dim1], abs(d__1)) / (d__2 = z__[k], abs(
		d__2));
	dscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L150:
	if (abd[m + k * abd_dim1] != 0.) {
	    z__[k] /= abd[m + k * abd_dim1];
	}
	if (abd[m + k * abd_dim1] == 0.) {
	    z__[k] = 1.;
	}
	lm = min(k,m) - 1;
	la = m - lm;
	lz = k - lm;
	t = -z__[k];
	daxpy_(&lm, &t, &abd[la + k * abd_dim1], &c__1, &z__[lz], &c__1);
/* L160: */
    }
/*     make znorm = 1.0 */
    s = 1. / dasum_(n, &z__[1], &c__1);
    dscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

    if (anorm != 0.) {
	*rcond = ynorm / anorm;
    }
    if (anorm == 0.) {
	*rcond = 0.;
    }
    return 0;
} /* dgbco_ */

/* Subroutine */ int dgeco_(doublereal *a, integer *lda, integer *n, integer *
	ipvt, doublereal *rcond, doublereal *z__)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);

    /* Local variables */
    integer j, k, l;
    doublereal s, t;
    integer kb;
    doublereal ek, sm, wk;
    integer kp1;
    doublereal wkm;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    integer info;
    extern /* Subroutine */ int dgefa_(doublereal *, integer *, integer *, 
	    integer *, integer *), dscal_(integer *, doublereal *, doublereal 
	    *, integer *);
    extern doublereal dasum_(integer *, doublereal *, integer *);
    doublereal anorm;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    doublereal ynorm;


/*     dgeco factors a double precision matrix by gaussian elimination */
/*     and estimates the condition of the matrix. */

/*     if  rcond  is not needed, dgefa is slightly faster. */
/*     to solve  a*x = b , follow dgeco by dgesl. */
/*     to compute  inverse(a)*c , follow dgeco by dgesl. */
/*     to compute  determinant(a) , follow dgeco by dgedi. */
/*     to compute  inverse(a) , follow dgeco by dgedi. */

/*     on entry */

/*        a       double precision(lda, n) */
/*                the matrix to be factored. */

/*        lda     integer */
/*                the leading dimension of the array  a . */

/*        n       integer */
/*                the order of the matrix  a . */

/*     on return */

/*        a       an upper triangular matrix and the multipliers */
/*                which were used to obtain it. */
/*                the factorization can be written  a = l*u  where */
/*                l  is a product of permutation and unit lower */
/*                triangular matrices and  u  is upper triangular. */

/*        ipvt    integer(n) */
/*                an integer vector of pivot indices. */

/*        rcond   double precision */
/*                an estimate of the reciprocal condition of  a . */
/*                for the system  a*x = b , relative perturbations */
/*                in  a  and  b  of size  epsilon  may cause */
/*                relative perturbations in  x  of size  epsilon/rcond . */
/*                if  rcond  is so small that the logical expression */
/*                           1.0 + rcond .eq. 1.0 */
/*                is true, then  a  may be singular to working */
/*                precision.  in particular,  rcond  is zero  if */
/*                exact singularity is detected or the estimate */
/*                underflows. */

/*        z       double precision(n) */
/*                a work vector whose contents are usually unimportant. */
/*                if  a  is close to a singular matrix, then  z  is */
/*                an approximate null vector in the sense that */
/*                norm(a*z) = rcond*norm(a)*norm(z) . */

/*     linpack. this version dated 08/14/78 . */
/*     cleve moler, university of new mexico, argonne national lab. */

/*     subroutines and functions */

/*     linpack dgefa */
/*     blas daxpy,ddot,dscal,dasum */
/*     fortran dabs,dmax1,dsign */

/*     internal variables */



/*     compute 1-norm of a */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipvt;
    --z__;

    /* Function Body */
    anorm = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	d__1 = anorm, d__2 = dasum_(n, &a[j * a_dim1 + 1], &c__1);
	anorm = max(d__1,d__2);
/* L10: */
    }

/*     factor */

    dgefa_(&a[a_offset], lda, n, &ipvt[1], &info);

/*     rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) . */
/*     estimate = norm(z)/norm(y) where  a*z = y  and  trans(a)*y = e . */
/*     trans(a)  is the transpose of a .  the components of  e  are */
/*     chosen to cause maximum local growth in the elements of w  where */
/*     trans(u)*w = e .  the vectors are frequently rescaled to avoid */
/*     overflow. */

/*     solve trans(u)*w = e */

    ek = 1.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	z__[j] = 0.;
/* L20: */
    }
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	if (z__[k] != 0.) {
	    d__1 = -z__[k];
	    ek = d_sign(&ek, &d__1);
	}
	if ((d__1 = ek - z__[k], abs(d__1)) <= (d__2 = a[k + k * a_dim1], abs(
		d__2))) {
	    goto L30;
	}
	s = (d__1 = a[k + k * a_dim1], abs(d__1)) / (d__2 = ek - z__[k], abs(
		d__2));
	dscal_(n, &s, &z__[1], &c__1);
	ek = s * ek;
L30:
	wk = ek - z__[k];
	wkm = -ek - z__[k];
	s = abs(wk);
	sm = abs(wkm);
	if (a[k + k * a_dim1] == 0.) {
	    goto L40;
	}
	wk /= a[k + k * a_dim1];
	wkm /= a[k + k * a_dim1];
	goto L50;
L40:
	wk = 1.;
	wkm = 1.;
L50:
	kp1 = k + 1;
	if (kp1 > *n) {
	    goto L90;
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    sm += (d__1 = z__[j] + wkm * a[k + j * a_dim1], abs(d__1));
	    z__[j] += wk * a[k + j * a_dim1];
	    s += (d__1 = z__[j], abs(d__1));
/* L60: */
	}
	if (s >= sm) {
	    goto L80;
	}
	t = wkm - wk;
	wk = wkm;
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    z__[j] += t * a[k + j * a_dim1];
/* L70: */
	}
L80:
L90:
	z__[k] = wk;
/* L100: */
    }
    s = 1. / dasum_(n, &z__[1], &c__1);
    dscal_(n, &s, &z__[1], &c__1);

/*     solve trans(l)*y = w */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	if (k < *n) {
	    i__2 = *n - k;
	    z__[k] += ddot_(&i__2, &a[k + 1 + k * a_dim1], &c__1, &z__[k + 1],
		     &c__1);
	}
	if ((d__1 = z__[k], abs(d__1)) <= 1.) {
	    goto L110;
	}
	s = 1. / (d__1 = z__[k], abs(d__1));
	dscal_(n, &s, &z__[1], &c__1);
L110:
	l = ipvt[k];
	t = z__[l];
	z__[l] = z__[k];
	z__[k] = t;
/* L120: */
    }
    s = 1. / dasum_(n, &z__[1], &c__1);
    dscal_(n, &s, &z__[1], &c__1);

    ynorm = 1.;

/*     solve l*v = y */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	l = ipvt[k];
	t = z__[l];
	z__[l] = z__[k];
	z__[k] = t;
	if (k < *n) {
	    i__2 = *n - k;
	    daxpy_(&i__2, &t, &a[k + 1 + k * a_dim1], &c__1, &z__[k + 1], &
		    c__1);
	}
	if ((d__1 = z__[k], abs(d__1)) <= 1.) {
	    goto L130;
	}
	s = 1. / (d__1 = z__[k], abs(d__1));
	dscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L130:
/* L140: */
	;
    }
    s = 1. / dasum_(n, &z__[1], &c__1);
    dscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

/*     solve  u*z = v */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	if ((d__1 = z__[k], abs(d__1)) <= (d__2 = a[k + k * a_dim1], abs(d__2)
		)) {
	    goto L150;
	}
	s = (d__1 = a[k + k * a_dim1], abs(d__1)) / (d__2 = z__[k], abs(d__2))
		;
	dscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L150:
	if (a[k + k * a_dim1] != 0.) {
	    z__[k] /= a[k + k * a_dim1];
	}
	if (a[k + k * a_dim1] == 0.) {
	    z__[k] = 1.;
	}
	t = -z__[k];
	i__2 = k - 1;
	daxpy_(&i__2, &t, &a[k * a_dim1 + 1], &c__1, &z__[1], &c__1);
/* L160: */
    }
/*     make znorm = 1.0 */
    s = 1. / dasum_(n, &z__[1], &c__1);
    dscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

    if (anorm != 0.) {
	*rcond = ynorm / anorm;
    }
    if (anorm == 0.) {
	*rcond = 0.;
    }
    return 0;
} /* dgeco_ */

/* Subroutine */ int dgedi_(doublereal *a, integer *lda, integer *n, integer *
	ipvt, doublereal *det, doublereal *work, integer *job)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    integer i__, j, k, l;
    doublereal t;
    integer kb, kp1, nm1;
    doublereal ten;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dswap_(integer *, doublereal *, integer *, doublereal 
	    *, integer *), daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);


/*     dgedi computes the determinant and inverse of a matrix */
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

/*        work    double precision(n) */
/*                work vector.  contents destroyed. */

/*        job     integer */
/*                = 11   both determinant and inverse. */
/*                = 01   inverse only. */
/*                = 10   determinant only. */

/*     on return */

/*        a       inverse of original matrix if requested. */
/*                otherwise unchanged. */

/*        det     double precision(2) */
/*                determinant of original matrix if requested. */
/*                otherwise not referenced. */
/*                determinant = det(1) * 10.0**det(2) */
/*                with  1.0 .le. dabs(det(1)) .lt. 10.0 */
/*                or  det(1) .eq. 0.0 . */

/*     error condition */

/*        a division by zero will occur if the input factor contains */
/*        a zero on the diagonal and the inverse is requested. */
/*        it will not occur if the subroutines are called correctly */
/*        and if dgeco has set rcond .gt. 0.0 or dgefa has set */
/*        info .eq. 0 . */

/*     linpack. this version dated 08/14/78 . */
/*     cleve moler, university of new mexico, argonne national lab. */

/*     subroutines and functions */

/*     blas daxpy,dscal,dswap */
/*     fortran dabs,mod */

/*     internal variables */



/*     compute determinant */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipvt;
    --det;
    --work;

    /* Function Body */
    if (*job / 10 == 0) {
	goto L70;
    }
    det[1] = 1.;
    det[2] = 0.;
    ten = 10.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (ipvt[i__] != i__) {
	    det[1] = -det[1];
	}
	det[1] = a[i__ + i__ * a_dim1] * det[1];
/*        ...exit */
	if (det[1] == 0.) {
	    goto L60;
	}
L10:
	if (abs(det[1]) >= 1.) {
	    goto L20;
	}
	det[1] = ten * det[1];
	det[2] += -1.;
	goto L10;
L20:
L30:
	if (abs(det[1]) < ten) {
	    goto L40;
	}
	det[1] /= ten;
	det[2] += 1.;
	goto L30;
L40:
/* L50: */
	;
    }
L60:
L70:

/*     compute inverse(u) */

    if (*job % 10 == 0) {
	goto L150;
    }
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	a[k + k * a_dim1] = 1. / a[k + k * a_dim1];
	t = -a[k + k * a_dim1];
	i__2 = k - 1;
	dscal_(&i__2, &t, &a[k * a_dim1 + 1], &c__1);
	kp1 = k + 1;
	if (*n < kp1) {
	    goto L90;
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    t = a[k + j * a_dim1];
	    a[k + j * a_dim1] = 0.;
	    daxpy_(&k, &t, &a[k * a_dim1 + 1], &c__1, &a[j * a_dim1 + 1], &
		    c__1);
/* L80: */
	}
L90:
/* L100: */
	;
    }

/*        form inverse(u)*inverse(l) */

    nm1 = *n - 1;
    if (nm1 < 1) {
	goto L140;
    }
    i__1 = nm1;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n - kb;
	kp1 = k + 1;
	i__2 = *n;
	for (i__ = kp1; i__ <= i__2; ++i__) {
	    work[i__] = a[i__ + k * a_dim1];
	    a[i__ + k * a_dim1] = 0.;
/* L110: */
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    t = work[j];
	    daxpy_(n, &t, &a[j * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &
		    c__1);
/* L120: */
	}
	l = ipvt[k];
	if (l != k) {
	    dswap_(n, &a[k * a_dim1 + 1], &c__1, &a[l * a_dim1 + 1], &c__1);
	}
/* L130: */
    }
L140:
L150:
    return 0;
} /* dgedi_ */

