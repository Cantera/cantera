#include "blaswrap.h"
#ifdef _cpluscplus
extern "C" {
#endif
#include "f2c.h"

/* Subroutine */ int dlacpy_(char *uplo, integer *m, integer *n, doublereal *
	a, integer *lda, doublereal *b, integer *ldb)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DLACPY copies all or part of a two-dimensional matrix A to another   
    matrix B.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            Specifies the part of the matrix A to be copied to B.   
            = 'U':      Upper triangular part   
            = 'L':      Lower triangular part   
            Otherwise:  All of the matrix A   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.   

    A       (input) DOUBLE PRECISION array, dimension (LDA,N)   
            The m by n matrix A.  If UPLO = 'U', only the upper triangle   
            or trapezoid is accessed; if UPLO = 'L', only the lower   
            triangle or trapezoid is accessed.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,M).   

    B       (output) DOUBLE PRECISION array, dimension (LDB,N)   
            On exit, B = A in the locations specified by UPLO.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= max(1,M).   

    =====================================================================   


       Parameter adjustments */
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;
    /* Local variables */
    static integer i__, j;
    extern logical lsame_(char *, char *);
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]
#define b_ref(a_1,a_2) b[(a_2)*b_dim1 + a_1]

    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1 * 1;
    b -= b_offset;

    /* Function Body */
    if (lsame_(uplo, "U")) {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = min(j,*m);
	    for (i__ = 1; i__ <= i__2; ++i__) {
		b_ref(i__, j) = a_ref(i__, j);
/* L10: */
	    }
/* L20: */
	}
    } else if (lsame_(uplo, "L")) {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *m;
	    for (i__ = j; i__ <= i__2; ++i__) {
		b_ref(i__, j) = a_ref(i__, j);
/* L30: */
	    }
/* L40: */
	}
    } else {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *m;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		b_ref(i__, j) = a_ref(i__, j);
/* L50: */
	    }
/* L60: */
	}
    }
    return 0;

/*     End of DLACPY */

} /* dlacpy_ */

#undef b_ref
#undef a_ref


#ifdef _cpluscplus
}
#endif
