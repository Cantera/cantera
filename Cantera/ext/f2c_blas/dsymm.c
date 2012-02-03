#include "blaswrap.h"
#ifdef _cpluscplus
extern "C" {
#endif
#include "f2c.h"

/* Subroutine */ int dsymm_(char *side, char *uplo, integer *m, integer *n, 
	doublereal *alpha, doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *beta, doublereal *c__, integer *ldc)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3;
    /* Local variables */
    static integer info;
    static doublereal temp1, temp2;
    static integer i__, j, k;
    extern logical lsame_(char *, char *);
    static integer nrowa;
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *);
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]
#define b_ref(a_1,a_2) b[(a_2)*b_dim1 + a_1]
#define c___ref(a_1,a_2) c__[(a_2)*c_dim1 + a_1]
/*  Purpose   
    =======   
    DSYMM  performs one of the matrix-matrix operations   
       C := alpha*A*B + beta*C,   
    or   
       C := alpha*B*A + beta*C,   
    where alpha and beta are scalars,  A is a symmetric matrix and  B and   
    C are  m by n matrices.   
    Parameters   
    ==========   
    SIDE   - CHARACTER*1.   
             On entry,  SIDE  specifies whether  the  symmetric matrix  A   
             appears on the  left or right  in the  operation as follows:   
                SIDE = 'L' or 'l'   C := alpha*A*B + beta*C,   
                SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,   
             Unchanged on exit.   
    UPLO   - CHARACTER*1.   
             On  entry,   UPLO  specifies  whether  the  upper  or  lower   
             triangular  part  of  the  symmetric  matrix   A  is  to  be   
             referenced as follows:   
                UPLO = 'U' or 'u'   Only the upper triangular part of the   
                                    symmetric matrix is to be referenced.   
                UPLO = 'L' or 'l'   Only the lower triangular part of the   
                                    symmetric matrix is to be referenced.   
             Unchanged on exit.   
    M      - INTEGER.   
             On entry,  M  specifies the number of rows of the matrix  C.   
             M  must be at least zero.   
             Unchanged on exit.   
    N      - INTEGER.   
             On entry, N specifies the number of columns of the matrix C.   
             N  must be at least zero.   
             Unchanged on exit.   
    ALPHA  - DOUBLE PRECISION.   
             On entry, ALPHA specifies the scalar alpha.   
             Unchanged on exit.   
    A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is   
             m  when  SIDE = 'L' or 'l'  and is  n otherwise.   
             Before entry  with  SIDE = 'L' or 'l',  the  m by m  part of   
             the array  A  must contain the  symmetric matrix,  such that   
             when  UPLO = 'U' or 'u', the leading m by m upper triangular   
             part of the array  A  must contain the upper triangular part   
             of the  symmetric matrix and the  strictly  lower triangular   
             part of  A  is not referenced,  and when  UPLO = 'L' or 'l',   
             the leading  m by m  lower triangular part  of the  array  A   
             must  contain  the  lower triangular part  of the  symmetric   
             matrix and the  strictly upper triangular part of  A  is not   
             referenced.   
             Before entry  with  SIDE = 'R' or 'r',  the  n by n  part of   
             the array  A  must contain the  symmetric matrix,  such that   
             when  UPLO = 'U' or 'u', the leading n by n upper triangular   
             part of the array  A  must contain the upper triangular part   
             of the  symmetric matrix and the  strictly  lower triangular   
             part of  A  is not referenced,  and when  UPLO = 'L' or 'l',   
             the leading  n by n  lower triangular part  of the  array  A   
             must  contain  the  lower triangular part  of the  symmetric   
             matrix and the  strictly upper triangular part of  A  is not   
             referenced.   
             Unchanged on exit.   
    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared   
             in the calling (sub) program.  When  SIDE = 'L' or 'l'  then   
             LDA must be at least  max( 1, m ), otherwise  LDA must be at   
             least  max( 1, n ).   
             Unchanged on exit.   
    B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).   
             Before entry, the leading  m by n part of the array  B  must   
             contain the matrix B.   
             Unchanged on exit.   
    LDB    - INTEGER.   
             On entry, LDB specifies the first dimension of B as declared   
             in  the  calling  (sub)  program.   LDB  must  be  at  least   
             max( 1, m ).   
             Unchanged on exit.   
    BETA   - DOUBLE PRECISION.   
             On entry,  BETA  specifies the scalar  beta.  When  BETA  is   
             supplied as zero then C need not be set on input.   
             Unchanged on exit.   
    C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).   
             Before entry, the leading  m by n  part of the array  C must   
             contain the matrix  C,  except when  beta  is zero, in which   
             case C need not be set on entry.   
             On exit, the array  C  is overwritten by the  m by n updated   
             matrix.   
    LDC    - INTEGER.   
             On entry, LDC specifies the first dimension of C as declared   
             in  the  calling  (sub)  program.   LDC  must  be  at  least   
             max( 1, m ).   
             Unchanged on exit.   
    Level 3 Blas routine.   
    -- Written on 8-February-1989.   
       Jack Dongarra, Argonne National Laboratory.   
       Iain Duff, AERE Harwell.   
       Jeremy Du Croz, Numerical Algorithms Group Ltd.   
       Sven Hammarling, Numerical Algorithms Group Ltd.   
       Set NROWA as the number of rows of A.   
       Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1 * 1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1 * 1;
    c__ -= c_offset;
    /* Function Body */
    if (lsame_(side, "L")) {
	nrowa = *m;
    } else {
	nrowa = *n;
    }
    upper = lsame_(uplo, "U");
/*     Test the input parameters. */
    info = 0;
    if (! lsame_(side, "L") && ! lsame_(side, "R")) {
	info = 1;
    } else if (! upper && ! lsame_(uplo, "L")) {
	info = 2;
    } else if (*m < 0) {
	info = 3;
    } else if (*n < 0) {
	info = 4;
    } else if (*lda < max(1,nrowa)) {
	info = 7;
    } else if (*ldb < max(1,*m)) {
	info = 9;
    } else if (*ldc < max(1,*m)) {
	info = 12;
    }
    if (info != 0) {
	xerbla_("DSYMM ", &info);
	return 0;
    }
/*     Quick return if possible. */
    if (*m == 0 || *n == 0 || *alpha == 0. && *beta == 1.) {
	return 0;
    }
/*     And when  alpha.eq.zero. */
    if (*alpha == 0.) {
	if (*beta == 0.) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    c___ref(i__, j) = 0.;
/* L10: */
		}
/* L20: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    c___ref(i__, j) = *beta * c___ref(i__, j);
/* L30: */
		}
/* L40: */
	    }
	}
	return 0;
    }
/*     Start the operations. */
    if (lsame_(side, "L")) {
/*        Form  C := alpha*A*B + beta*C. */
	if (upper) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    temp1 = *alpha * b_ref(i__, j);
		    temp2 = 0.;
		    i__3 = i__ - 1;
		    for (k = 1; k <= i__3; ++k) {
			c___ref(k, j) = c___ref(k, j) + temp1 * a_ref(k, i__);
			temp2 += b_ref(k, j) * a_ref(k, i__);
/* L50: */
		    }
		    if (*beta == 0.) {
			c___ref(i__, j) = temp1 * a_ref(i__, i__) + *alpha * 
				temp2;
		    } else {
			c___ref(i__, j) = *beta * c___ref(i__, j) + temp1 * 
				a_ref(i__, i__) + *alpha * temp2;
		    }
/* L60: */
		}
/* L70: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		for (i__ = *m; i__ >= 1; --i__) {
		    temp1 = *alpha * b_ref(i__, j);
		    temp2 = 0.;
		    i__2 = *m;
		    for (k = i__ + 1; k <= i__2; ++k) {
			c___ref(k, j) = c___ref(k, j) + temp1 * a_ref(k, i__);
			temp2 += b_ref(k, j) * a_ref(k, i__);
/* L80: */
		    }
		    if (*beta == 0.) {
			c___ref(i__, j) = temp1 * a_ref(i__, i__) + *alpha * 
				temp2;
		    } else {
			c___ref(i__, j) = *beta * c___ref(i__, j) + temp1 * 
				a_ref(i__, i__) + *alpha * temp2;
		    }
/* L90: */
		}
/* L100: */
	    }
	}
    } else {
/*        Form  C := alpha*B*A + beta*C. */
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    temp1 = *alpha * a_ref(j, j);
	    if (*beta == 0.) {
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    c___ref(i__, j) = temp1 * b_ref(i__, j);
/* L110: */
		}
	    } else {
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    c___ref(i__, j) = *beta * c___ref(i__, j) + temp1 * b_ref(
			    i__, j);
/* L120: */
		}
	    }
	    i__2 = j - 1;
	    for (k = 1; k <= i__2; ++k) {
		if (upper) {
		    temp1 = *alpha * a_ref(k, j);
		} else {
		    temp1 = *alpha * a_ref(j, k);
		}
		i__3 = *m;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    c___ref(i__, j) = c___ref(i__, j) + temp1 * b_ref(i__, k);
/* L130: */
		}
/* L140: */
	    }
	    i__2 = *n;
	    for (k = j + 1; k <= i__2; ++k) {
		if (upper) {
		    temp1 = *alpha * a_ref(j, k);
		} else {
		    temp1 = *alpha * a_ref(k, j);
		}
		i__3 = *m;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    c___ref(i__, j) = c___ref(i__, j) + temp1 * b_ref(i__, k);
/* L150: */
		}
/* L160: */
	    }
/* L170: */
	}
    }
    return 0;
/*     End of DSYMM . */
} /* dsymm_ */
#undef c___ref
#undef b_ref
#undef a_ref

#ifdef _cpluscplus
}
#endif
