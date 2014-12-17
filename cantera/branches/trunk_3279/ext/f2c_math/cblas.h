// -*- C++ -*-

// ============================================= //
//  die double-Versionen der BLAS Level 1 und 2  //
// ============================================= //

#ifndef CBLAS1_H
// ============================================================================

// generate a plane rotation
void drotg( double *a, double *b, double *c, double *s );


#if 0
// generate a modified plane rotation
void drotmg( double *d1, double *d2, double *a, double b, double *param );
#endif

// apply a plane rotation
void drot( int n, double *x, int incx, double *y, int incy, double c,
	    double s );


#if 0
// apply a modified plane rotation
void drotm( int n, double *x, int incx, double *y, int incy, double *param );
#endif


// x <=> y
void dswap( int n, double *x, int incx, double *y, int incy );


// x <= a*x
void dscal( int n, double alpha, double *x, int incx );


// y <= x
void dcopy( int n, const double *x, int incx, double *y, int incy );


// y <= a*x+y
void daxpy( int n, double alpha, const double *x, int incx, double *y,
	    int incy );


// dot <= x^T*y
double ddot( int n, const double *x, int incx, const double *y, int incy );


// dnrm2 <= |x|_2
double dnrm2( int n, const double *x, int incx );


// asum <= |x|_1
double dasum( int n, const double *x, int incx );


// idamax <= first k such that |x_k| = max|x_i|
int idamax( int n, const double *x, int incx );

// ============================================================================
#endif // CBLAS1_H


#ifndef CBLAS2_H
// ============================================================================


enum MatrixTranspose { NoTranspose=0, Transpose=1, ConjugateTranspose=2 };
enum MatrixTriangle { UpperTriangle=0, LowerTriangle=1 };
enum MatrixUnitTriangular { UnitTriangular=0, NotUnitTriangular=1 };


// ============================================================================


// y <= alpha*A*x + beta*y,  y <= alpha*A^T*x + beta*y,   A-(m,n)
void dgemv( MatrixTranspose trans, int m, int n, double alpha,
	    const double *A, int ldA, const double *x, int incx,
	    double beta, double *y, int incy );


// y <= alpha*A*x + beta*y,  y <= alpha*A^T*x + beta*y,   A-(m,n)
void dgbmv( MatrixTranspose trans, int m, int n, int kl, int ku, double alpha,
	    const double *A, int ldA, const double *x, int incx, double *beta,
	    double *y, int incy );


// y <= alpha*A*x + beta*y
void dsymv( MatrixTriangle uplo, int n, double alpha, const double *A, int ldA,
	    const double *x, int incx, double beta, double *y, int incy );


// y <= alpha*A*x + beta*y
void dsbmv( MatrixTriangle uplo, int n, int k, double alpha, double *A,
	    int ldA, const double *x, int incx, double beta, double *y,
	    int *incy );


// y <= alpha*A*x + beta*y
void dspmv( MatrixTriangle uplo, int n, double alpha, const double *AP,
	    const double *x, int incx, double beta, double *y, int incy );


// x <= A*x,  x <= A^T*x
void dtrmv( MatrixTriangle uplo, MatrixTranspose trans,
	    MatrixUnitTriangular diag, int n, const double *A, int ldA,
	    double *x, int incx );


// x <= A*x,  x <= A^T*x
void dtbmv( MatrixTriangle uplo, MatrixTranspose trans,
	    MatrixUnitTriangular diag, int n, int k, const double *A, int ldA,
	    double *x, int incx );


// x <= A*x,  x <= A^T*x
void dtpmv( MatrixTriangle uplo, MatrixTranspose trans,
	    MatrixUnitTriangular diag, int n, int k, const double *AP,
	    double *x, int incx );


// x <= A^{-1}*x,  x <= A^{-T}*x
void dtrsv( MatrixTriangle uplo, MatrixTranspose trans,
	    MatrixUnitTriangular diag, int n, const double *A, int ldA,
	    double *x, int incx );


// x <= A^{-1}*x,  x <= A^{-T}*x
void dtbsv( MatrixTriangle uplo, MatrixTranspose trans,
	    MatrixUnitTriangular diag, int n, int k, const double *A, int ldA,
	    double *x, int incx );


// x <= A^{-1}*x,  x <= A^{-T}*x
void dtpsv( MatrixTriangle uplo, MatrixTranspose trans,
	    MatrixUnitTriangular diag, int n, int k, const double *AP,
	    double *x, int incx );


// A <= alpha*x*y^T + A,   A-(m,n)
void dger( int m, int n, double alpha, const double *x, int incx,
	   const double *y, int incy, double *A, int ldA );


// A <= alpha*x*x^T + A
void dsyr( MatrixTriangle uplo, int n, double alpha, const double *x,
	   int incx, double *A, int ldA );


// A <= alpha*x*x^T + A
void dspr( MatrixTriangle uplo, int n, double alpha, const double *x,
	   int incx, double *AP );


// A <= alpha*x*y^T + alpha*y*x^T + A
void dsyr2( MatrixTriangle uplo, int n, double alpha, const double *x,
	    int incx, const double *y, int incy, double *A, int ldA );


// A <= alpha*x*y^T + alpha*y*x^T + A
void dspr2( MatrixTriangle uplo, int n, double alpha, const double *x,
	    int incx, const double *y, int incy, double *AP );

// ============================================================================
#endif // CBLAS2_H


#ifndef BLAS1_H
#define BLAS1_H
// ============================================================================

// generate a plane rotation
extern "C"
void drotg_( double *a, double *b, double *c, double *s );


#if 0
// generate a modified plane rotation
extern "C"
void drotmg_( double *d1, double *d2, double *a, double *b, double *param );
#endif


// apply a plane rotation
extern "C"
void drot_( int *n, double *x, int *incx, double *y, int *incy,
	    double *c, double *s );


#if 0
// apply a modified plane rotation
extern "C"
void drotm_( int *n, double *x, int *incx, double *y, int *incy,
	     double *param );
#endif


// x <=> y
extern "C"
void dswap_( const int *n, double *x, const int *incx, double *y,
	     const int *incy );

// x <= a*x
extern "C"
void dscal_( const int *n, const double *alpha, double *x, const int *incx );


// y <= x
extern "C"
void dcopy_( const int *n, const double *x, const int *incx, double *y,
	     const int *incy );


// y <= a*x+y
extern "C"
void daxpy_( const int *n, const double *alpha, const double *x,
	     const int *incx, double *y, const int *incy );


// dot <= x^T*y
extern "C"
double ddot_( const int *n, const double *x, const int *incx, const double *y,
	      const int *incy );


// dnrm2 <= |x|_2
extern "C"
double dnrm2_( const int *n, const double *x, const int *incx );


// asum <= |x|_1
extern "C"
double dasum_( const int *n, const double *x, const int *incx );


// idamax <= first k such that |x_k| = max|x_i|
extern "C"
int idamax_( const int *n, const double *x, const int *incx );


// ============================================================================
#endif // BLAS1_H



#ifndef CBLAS1_H
#define CBLAS1_H
// ============================================================================


#ifdef __linux__ // muss dnorm2 f"ur linux neu implementieren
#  include <math.h>
#endif

inline
void drotg( double *a, double *b, double *c, double *s ) {
  drotg_(a,b,c,s);
}

#if 0
inline
void drotmg( double *d1, double *d2, double *a, double b, double *param ) {
  drotmg_(d1,d2,a,&b,param);
}
#endif

inline
void drot( int n, double *x, int incx, double *y, int incy, double c,
	   double s ) {
  drot_(&n,x,&incx,y,&incy,&c,&s);
}

#if 0
inline
void drotm( int n, double *x, int incx, double *y, int incy, double *param ) {
  drotm_(&n,x,&incx,y,&incy,param);
}
#endif

inline
void dswap( int n, double *x, int incx, double *y, int incy ) {
  dswap_(&n,x,&incx,y,&incy);
}

inline
void dscal( int n, double alpha, double *x, int incx ) {
  int nn = n;
  int incxx = incx;
  double aa = alpha;
  dscal_(&nn,&aa,x,&incxx);
}

inline
void dcopy( int n, const double *x, int incx, double *y, int incy ) {
  int nn = n;
  int incxx = incx;
  int incyy = incy;
  dcopy_(&nn,x,&incxx,y,&incyy);
}

inline
void daxpy( int n, double alpha, const double *x, int incx, double *y,
	    int incy ) {
  double aa = alpha;
  int incxx = incx;
  int incyy = incy;
  daxpy_(&n,&aa,x,&incxx,y,&incyy);
}

inline
double ddot( int n, const double *x, int incx, const double *y, int incy ) {
  int nn = n;
  int incxx = incx;
  int incyy = incy;
  return ddot_(&nn,x,&incxx,y,&incyy);
}

inline
double dnrm2( int n, const double *x, int incx ) {
  int nn = n;
  int incxx = incx;
#ifdef __linux__ //  fehlerhafte Berechnung
  double d=0.;
  while ( nn-- )
    d+=(*x)*(*x), x+=incxx;
  return sqrt(d);
#else            // unter nicht-Linux korrekt
  return dnrm2_(&nn,x,&incxx);
#endif
}

inline
double dasum( int n, const double *x, int incx ) {
  return dasum_(&n,x,&incx);
}

inline
int idamax( int n, const double *x, int incx ) {
  return idamax_(&n,x,&incx);
}


// ============================================================================
#endif // CBLAS1_H


#ifndef BLAS2_H
#define BLAS2_H
// ============================================================================

// y <= alpha*A*x + beta*y,  y <= alpha*A^T*x + beta*y,   A-(m,n)
//extern "C"
//void dgemv_( const char *trans, const int *m, const int *n,
//	     const double *alpha, const double *A, const int *ldA,
//	     const double *x, const int *incx,
//	     const double *beta, double *y, const int *incy );


// y <= alpha*A*x + beta*y,  y <= alpha*A^T*x + beta*y,   A-(m,n)
extern "C"
void dgbmv_( const char *trans, const int *m, const int *n, const int *kl,
	     const int *ku, const double *alpha, const double *A,
	     const int *ldA, const double *x, const int *incx,
	     const double *beta, double *y, const int *incy );


// y <= alpha*A*x + beta*y
extern "C"
void dsymv_( const char *uplo, const int *n, const double *alpha,
	     const double *A, const int *ldA, const double *x, const int *incx,
	     const double *beta, double *y, const int *incy );


// y <= alpha*A*x + beta*y
extern "C"
void dsbmv_( const char *uplo, const int *n, const int *k, const double *alpha,
	     const double *A, const int *ldA, const double *x, const int *incx,
	     const double *beta, double *y, const int *incy );


// y <= alpha*A*x + beta*y
extern "C"
void dspmv_( const char *uplo, const int *n, const double *alpha,
	     const double *AP, const double *x, const int *incx,
	     const double *beta, double *y, const int *incy );


// x <= A*x,  x <= A^T*x
extern "C"
void dtrmv_( const char *uplo, const char *trans, const char *diag,
	     const int *n, const double *A, const int *ldA,
	     double *x, const int *incx );


// x <= A*x,  x <= A^T*x
extern "C"
void dtbmv_( const char *uplo, const char *trans, const char *diag,
	     const int *n, const int *k, const double *A, const int *ldA,
	     double *x, const int *incx );


// x <= A*x,  x <= A^T*x
extern "C"
void dtpmv_( const char *uplo, const char *trans, const char *diag,
	    const int *n, const double *AP, double *x, const int *incx );


// x <= A^{-1}*x,  x <= A^{-T}*x
extern "C"
void dtrsv_( const char *uplo, const char *trans, const char *diag,
	     const int *n, const double *A, const int *ldA,
	     double *x, const int *incx );


// x <= A^{-1}*x,  x <= A^{-T}*x
extern "C"
void dtbsv_( const char *uplo, const char *trans, const char *diag,
	     const int *n, const int *k, const double *A, const int *ldA,
	     double *x, const int *incx );


// x <= A^{-1}*x,  x <= A^{-T}*x
extern "C"
void dtpsv_( const char *uplo, const char *trans, const char *diag,
	     const int *n, const double *AP, double *x, const int *incx );


// A <= alpha*x*y^T + A,   A-(m,n)
extern "C"
void dger_( const int *m, const int *n, const double *alpha, const double *x,
	    const int *incx, const double *y, const int *incy, double *A,
	    const int *ldA );


// A <= alpha*x*x^T + A
extern "C"
void dsyr_( const char *uplo, const int *n, const double *alpha,
	    const double *x, const int *incx, double *A, const int *ldA );


// A <= alpha*x*x^T + A
extern "C"
void dspr_( const char *uplo, const int *n, const double *alpha,
	    const double *x, const int *incx, double *AP );


// A <= alpha*x*y^T + alpha*y*x^T + A
extern "C"
void dsyr2_( const char *uplo, const int *n, const double *alpha,
	     const double *x, const int *incx, const double *y,
	     const int *incy, double *A, const int *ldA );


// A <= alpha*x*y^T + alpha*y*x^T + A
extern "C"
void dspr2_( const char *uplo, const int *n, const double *alpha,
	     const double *x, const int *incx, const double *y,
	     const int *incy, double *AP );

// ============================================================================
#endif // BLAS2_H


#ifndef CBLAS2_H
#define CBLAS2_H
// ============================================================================


// y <= alpha*A*x + beta*y,  y <= alpha*A^T*x + beta*y,   A-(m,n)
inline
void dgemv( MatrixTranspose trans, int m, int n, double alpha,
	    const double *A, int ldA, const double *x, int incx,
	    double beta, double *y, int incy ) {
  const char *T[3] = { "N", "T", 0 };
  int mm = m;
  int nn = n;
  double aa = alpha;
  double bb = beta;
  int ldaa = ldA;
  int incxx = incx;
  int incyy = incy;
  dgemv_(T[(int)trans],&mm,&nn,&aa,A,&ldaa,x,&incxx,&bb,y,&incyy,1);
}


// y <= alpha*A*x + beta*y,  y <= alpha*A^T*x + beta*y,   A-(m,n)
inline
void dgbmv( MatrixTranspose trans, int m, int n, int kl, int ku, double alpha,
	    const double *A, int ldA, const double *x, int incx, double beta,
	    double *y, int incy ) {
  const char *T[3] = { "N", "T" };
  dgbmv_(T[(int)trans],&m,&n,&kl,&ku,&alpha,A,&ldA,x,&incx,&beta,y,&incy);
}

// y <= alpha*A*x + beta*y
inline
void dsymv( MatrixTriangle uplo, int n, double alpha, const double *A, int ldA,
	    const double *x, int incx, double beta, double *y, int incy ) {
  const char *UL[2] = { "U", "L" };
  dsymv_(UL[(int)uplo],&n,&alpha,A,&ldA,x,&incx,&beta,y,&incy);
}


// y <= alpha*A*x + beta*y
inline
void dsbmv( MatrixTriangle uplo, int n, int k, double alpha, double *A,
	    int ldA, const double *x, int incx, double beta, double *y,
	    int incy ) {
  const char *UL[2] = { "U", "L" };
  dsbmv_(UL[(int)uplo],&n,&k,&alpha,A,&ldA,x,&incx,&beta,y,&incy);
}


// y <= alpha*A*x + beta*y
inline
void dspmv( MatrixTriangle uplo, int n, double alpha, const double *AP,
	    const double *x, int incx, double beta, double *y, int incy ) {
  const char *UL[2] = { "U", "L" };
  dspmv_(UL[(int)uplo],&n,&alpha,AP,x,&incx,&beta,y,&incy);
}


// x <= A*x,  x <= A^T*x
inline
void dtrmv( MatrixTriangle uplo, MatrixTranspose trans,
	    MatrixUnitTriangular diag, int n, const double *A, int ldA,
	    double *x, int incx ) {
  const char *UL[2] = { "U", "L" };
  const char *T[3]  = { "N", "T", 0 };
  const char *D[2]  = { "U", "N" };
  dtrmv_(UL[(int)uplo],T[(int)trans],D[(int)diag],&n,A,&ldA,x,&incx);
}


// x <= A*x,  x <= A^T*x
inline
void dtbmv( MatrixTriangle uplo, MatrixTranspose trans,
	    MatrixUnitTriangular diag, int n, int k, const double *A, int ldA,
	    double *x, int incx ) {
  const char *UL[2] = { "U", "L" };
  const char *T[3]  = { "N", "T", 0 };
  const char *D[2]  = { "U", "N" };
  dtbmv_(UL[(int)uplo],T[(int)trans],D[(int)diag],&n,&k,A,&ldA,x,&incx);
}


// x <= A*x,  x <= A^T*x
inline
void dtpmv( MatrixTriangle uplo, MatrixTranspose trans,
	    MatrixUnitTriangular diag, int n, const double *AP,
	    double *x, int incx ) {
  const char *UL[2] = { "U", "L" };
  const char *T[3]  = { "N", "T", 0 };
  const char *D[2]  = { "U", "N" };
  dtpmv_(UL[(int)uplo],T[(int)trans],D[(int)diag],&n,AP,x,&incx);
}


// x <= A^{-1}*x,  x <= A^{-T}*x
inline
void dtrsv( MatrixTriangle uplo, MatrixTranspose trans,
	    MatrixUnitTriangular diag, int n, const double *A, int ldA,
	    double *x, int incx ) {
  const char *UL[2] = { "U", "L" };
  const char *T[3]  = { "N", "T", 0 };
  const char *D[2]  = { "U", "N" };
  dtrsv_(UL[(int)uplo],T[(int)trans],D[(int)diag],&n,A,&ldA,x,&incx);
}


// x <= A^{-1}*x,  x <= A^{-T}*x
inline
void dtbsv( MatrixTriangle uplo, MatrixTranspose trans,
	    MatrixUnitTriangular diag, int n, int k, const double *A, int ldA,
	    double *x, int incx ) {
  const char *UL[2] = { "U", "L" };
  const char *T[3]  = { "N", "T", 0 };
  const char *D[2]  = { "U", "N" };
  dtbsv_(UL[(int)uplo],T[(int)trans],D[(int)diag],&n,&k,A,&ldA,x,&incx);
}


// x <= A^{-1}*x,  x <= A^{-T}*x
inline
void dtpsv( MatrixTriangle uplo, MatrixTranspose trans,
	    MatrixUnitTriangular diag, int n, const double *AP,
	    double *x, int incx ) {
  const char *UL[2] = { "U", "L" };
  const char *T[3]  = { "N", "T", 0 };
  const char *D[2]  = { "U", "N" };
  int nn = n;
  int incxx = incx;
  dtpsv_(UL[(int)uplo],T[(int)trans],D[(int)diag],&nn,AP,x,&incxx);
}


// A <= alpha*x*y^T + A,   A-(m,n)
inline
void dger( int m, int n, double alpha, const double *x, int incx,
	   const double *y, int incy, double *A, int ldA ) {
  dger_(&m,&n,&alpha,x,&incx,y,&incy,A,&ldA);
}


// A <= alpha*x*x^T + A
inline
void dsyr( MatrixTriangle uplo, int n, double alpha, const double *x,
	   int incx, double *A, int ldA ) {
  const char *UL[2] = { "U", "L" };
  dsyr_(UL[(int)uplo],&n,&alpha,x,&incx,A,&ldA);
}


// A <= alpha*x*x^T + A
inline
void dspr( MatrixTriangle uplo, int n, double alpha, const double *x,
	   int incx, double *AP ) {
  const char *UL[2] = { "U", "L" };
  dspr_(UL[(int)uplo],&n,&alpha,x,&incx,AP);
}


// A <= alpha*x*y^T + alpha*y*x^T + A
inline
void dsyr2( MatrixTriangle uplo, int n, double alpha, const double *x,
	    int incx, const double *y, int incy, double *A, int ldA ) {
  const char *UL[2] = { "U", "L" };
  dsyr2_(UL[(int)uplo],&n,&alpha,x,&incx,y,&incy,A,&ldA);
}


// A <= alpha*x*y^T + alpha*y*x^T + A
inline
void dspr2( MatrixTriangle uplo, int n, double alpha, const double *x,
	    int incx, const double *y, int incy, double *AP ) {
  const char *UL[2] = { "U", "L" };
  dspr2_(UL[(int)uplo],&n,&alpha,x,&incx,y,&incy,AP);
}


// ============================================================================


#endif // CBLAS2_H
