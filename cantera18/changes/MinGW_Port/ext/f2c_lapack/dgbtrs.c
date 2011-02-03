#include "blaswrap.h"
#ifdef _cpluscplus
extern "C" {
#endif
#include "f2c.h"

/* Subroutine */ int dgbtrs_(char *trans, integer *n, integer *kl, integer *
	ku, integer *nrhs, doublereal *ab, integer *ldab, integer *ipiv, 
	doublereal *b, integer *ldb, integer *info)
{
/*  -- LAPACK routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DGBTRS solves a system of linear equations   
       A * X = B  or  A' * X = B   
    with a general band matrix A using the LU factorization computed   
    by DGBTRF.   

    Arguments   
    =========   

    TRANS   (input) CHARACTER*1   
            Specifies the form of the system of equations.   
            = 'N':  A * X = B  (No transpose)   
            = 'T':  A'* X = B  (Transpose)   
            = 'C':  A'* X = B  (Conjugate transpose = Transpose)   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    KL      (input) INTEGER   
            The number of subdiagonals within the band of A.  KL >= 0.   

    KU      (input) INTEGER   
            The number of superdiagonals within the band of A.  KU >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrix B.  NRHS >= 0.   

    AB      (input) DOUBLE PRECISION array, dimension (LDAB,N)   
            Details of the LU factorization of the band matrix A, as   
            computed by DGBTRF.  U is stored as an upper triangular band   
            matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and   
            the multipliers used during the factorization are stored in   
            rows KL+KU+2 to 2*KL+KU+1.   

    LDAB    (input) INTEGER   
            The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.   

    IPIV    (input) INTEGER array, dimension (N)   
            The pivot indices; for 1 <= i <= N, row i of the matrix was   
            interchanged with row IPIV(i).   

    B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)   
            On entry, the right hand side matrix B.   
            On exit, the solution matrix X.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= max(1,N).   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0: if INFO = -i, the i-th argument had an illegal value   

    =====================================================================   


       Test the input parameters.   

       Parameter adjustments */
    /* Table of constant values */
    static doublereal c_b7 = -1.;
    static integer c__1 = 1;
    static doublereal c_b23 = 1.;
    
    /* System generated locals */
    integer ab_dim1, ab_offset, b_dim1, b_offset, i__1, i__2, i__3;
    /* Local variables */
    extern /* Subroutine */ int dger_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer i__, j, l;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *), dswap_(integer *, 
	    doublereal *, integer *, doublereal *, integer *), dtbsv_(char *, 
	    char *, char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical lnoti;
    static integer kd, lm;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static logical notran;
#define b_ref(a_1,a_2) b[(a_2)*b_dim1 + a_1]
#define ab_ref(a_1,a_2) ab[(a_2)*ab_dim1 + a_1]


    ab_dim1 = *ldab;
    ab_offset = 1 + ab_dim1 * 1;
    ab -= ab_offset;
    --ipiv;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1 * 1;
    b -= b_offset;

    /* Function Body */
    *info = 0;
    notran = lsame_(trans, "N");
    if (! notran && ! lsame_(trans, "T") && ! lsame_(
	    trans, "C")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*kl < 0) {
	*info = -3;
    } else if (*ku < 0) {
	*info = -4;
    } else if (*nrhs < 0) {
	*info = -5;
    } else if (*ldab < (*kl << 1) + *ku + 1) {
	*info = -7;
    } else if (*ldb < max(1,*n)) {
	*info = -10;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGBTRS", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0 || *nrhs == 0) {
	return 0;
    }

    kd = *ku + *kl + 1;
    lnoti = *kl > 0;

    if (notran) {

/*        Solve  A*X = B.   

          Solve L*X = B, overwriting B with X.   

          L is represented as a product of permutations and unit lower   
          triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1),   
          where each transformation L(i) is a rank-one modification of   
          the identity matrix. */

	if (lnoti) {
	    i__1 = *n - 1;
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
		i__2 = *kl, i__3 = *n - j;
		lm = min(i__2,i__3);
		l = ipiv[j];
		if (l != j) {
		    dswap_(nrhs, &b_ref(l, 1), ldb, &b_ref(j, 1), ldb);
		}
		dger_(&lm, nrhs, &c_b7, &ab_ref(kd + 1, j), &c__1, &b_ref(j, 
			1), ldb, &b_ref(j + 1, 1), ldb);
/* L10: */
	    }
	}

	i__1 = *nrhs;
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Solve U*X = B, overwriting B with X. */

	    i__2 = *kl + *ku;
	    dtbsv_("Upper", "No transpose", "Non-unit", n, &i__2, &ab[
		    ab_offset], ldab, &b_ref(1, i__), &c__1);
/* L20: */
	}

    } else {

/*        Solve A'*X = B. */

	i__1 = *nrhs;
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Solve U'*X = B, overwriting B with X. */

	    i__2 = *kl + *ku;
	    dtbsv_("Upper", "Transpose", "Non-unit", n, &i__2, &ab[ab_offset],
		     ldab, &b_ref(1, i__), &c__1);
/* L30: */
	}

/*        Solve L'*X = B, overwriting B with X. */

	if (lnoti) {
	    for (j = *n - 1; j >= 1; --j) {
/* Computing MIN */
		i__1 = *kl, i__2 = *n - j;
		lm = min(i__1,i__2);
		dgemv_("Transpose", &lm, nrhs, &c_b7, &b_ref(j + 1, 1), ldb, &
			ab_ref(kd + 1, j), &c__1, &c_b23, &b_ref(j, 1), ldb);
		l = ipiv[j];
		if (l != j) {
		    dswap_(nrhs, &b_ref(l, 1), ldb, &b_ref(j, 1), ldb);
		}
/* L40: */
	    }
	}
    }
    return 0;

/*     End of DGBTRS */

} /* dgbtrs_ */

#undef ab_ref
#undef b_ref


#ifdef _cpluscplus
}
#endif
