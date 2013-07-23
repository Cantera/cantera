#include "blaswrap.h"
#ifdef _cpluscplus
extern "C" {
#endif
#include "f2c.h"

/* Subroutine */ int dlaswp_(integer *n, doublereal *a, integer *lda, integer 
	*k1, integer *k2, integer *ipiv, integer *incx)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       June 30, 1999   


    Purpose   
    =======   

    DLASWP performs a series of row interchanges on the matrix A.   
    One row interchange is initiated for each of rows K1 through K2 of A.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The number of columns of the matrix A.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the matrix of column dimension N to which the row   
            interchanges will be applied.   
            On exit, the permuted matrix.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.   

    K1      (input) INTEGER   
            The first element of IPIV for which a row interchange will   
            be done.   

    K2      (input) INTEGER   
            The last element of IPIV for which a row interchange will   
            be done.   

    IPIV    (input) INTEGER array, dimension (M*abs(INCX))   
            The vector of pivot indices.  Only the elements in positions   
            K1 through K2 of IPIV are accessed.   
            IPIV(K) = L implies rows K and L are to be interchanged.   

    INCX    (input) INTEGER   
            The increment between successive values of IPIV.  If IPIV   
            is negative, the pivots are applied in reverse order.   

    Further Details   
    ===============   

    Modified by   
     R. C. Whaley, Computer Science Dept., Univ. of Tenn., Knoxville, USA   

   =====================================================================   


       Interchange row I with row IPIV(I) for each of rows K1 through K2.   

       Parameter adjustments */
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    /* Local variables */
    static doublereal temp;
    static integer i__, j, k, i1, i2, n32, ip, ix, ix0, inc;
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]

    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --ipiv;

    /* Function Body */
    if (*incx > 0) {
	ix0 = *k1;
	i1 = *k1;
	i2 = *k2;
	inc = 1;
    } else if (*incx < 0) {
	ix0 = (1 - *k2) * *incx + 1;
	i1 = *k2;
	i2 = *k1;
	inc = -1;
    } else {
	return 0;
    }

    n32 = *n / 32 << 5;
    if (n32 != 0) {
	i__1 = n32;
	for (j = 1; j <= i__1; j += 32) {
	    ix = ix0;
	    i__2 = i2;
	    i__3 = inc;
	    for (i__ = i1; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__3) 
		    {
		ip = ipiv[ix];
		if (ip != i__) {
		    i__4 = j + 31;
		    for (k = j; k <= i__4; ++k) {
			temp = a_ref(i__, k);
			a_ref(i__, k) = a_ref(ip, k);
			a_ref(ip, k) = temp;
/* L10: */
		    }
		}
		ix += *incx;
/* L20: */
	    }
/* L30: */
	}
    }
    if (n32 != *n) {
	++n32;
	ix = ix0;
	i__1 = i2;
	i__3 = inc;
	for (i__ = i1; i__3 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__3) {
	    ip = ipiv[ix];
	    if (ip != i__) {
		i__2 = *n;
		for (k = n32; k <= i__2; ++k) {
		    temp = a_ref(i__, k);
		    a_ref(i__, k) = a_ref(ip, k);
		    a_ref(ip, k) = temp;
/* L40: */
		}
	    }
	    ix += *incx;
/* L50: */
	}
    }

    return 0;

/*     End of DLASWP */

} /* dlaswp_ */

#undef a_ref


#ifdef _cpluscplus
}
#endif
