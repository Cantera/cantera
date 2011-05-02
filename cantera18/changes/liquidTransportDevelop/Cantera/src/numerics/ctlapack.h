/**
 * @file ctlapack.h
 */

/* $Author$
 * $Revision$
 * $Date$
 */

// Copyright 2001 California Institute of Technology.  

#ifndef CT_CTLAPACK_H
#define CT_CTLAPACK_H

#ifdef DARWIN
#undef USE_CBLAS
#undef NO_FTN_STRING_LEN_AT_END
#endif

#include "ct_defs.h"

//#include <vecLib/cblas.h>

// map BLAS names to names with or without a trailing underscore.
#ifndef LAPACK_FTN_TRAILING_UNDERSCORE

#define _DGEMV_   dgemv
#define _DGETRF_  dgetrf
#define _DGETRS_  dgetrs
#define _DGETRI_  dgetri
#define _DGELSS_  dgelss
#define _DGBSV_   dgbsv
#define _DGBTRF_  dgbtrf
#define _DGBTRS_  dgbtrs

#define _DSCAL_   dscal

#define _DGEQRF_  dgeqrf
#define _DORMQR_  dormqr
#define _DTRTRS_  dtrtrs
#define _DTRCON_  dtrcon
#define _DPOTRF_  dpotrf
#define _DPOTRS_  dpotrs

#else

#define _DGEMV_   dgemv_
#define _DGETRF_  dgetrf_
#define _DGETRS_  dgetrs_
#define _DGETRI_  dgetri_
#define _DGELSS_  dgelss_
#define _DGBSV_   dgbsv_
#define _DGBTRF_  dgbtrf_
#define _DGBTRS_  dgbtrs_

#define _DSCAL_   dscal_

#define _DGEQRF_  dgeqrf_
#define _DORMQR_  dormqr_
#define _DTRTRS_  dtrtrs_
#define _DTRCON_  dtrcon_

#define _DPOTRF_  dpotrf_
#define _DPOTRS_  dpotrs_

#endif


namespace ctlapack {
    typedef enum {Transpose = 1, NoTranspose = 0} transpose_t;
    typedef enum {ColMajor = 1, RowMajor = 0} storage_t;
    typedef enum {UpperTriangular = 0, LowerTriangular = 1} upperlower_t;
    typedef enum {Left = 0, Right = 1} side_t;
}
const char no_yes[2] = {'N', 'T'};
const char upper_lower[2] = {'U', 'L'};
const char left_right[2] = {'L', 'R'};

#ifdef USE_CBLAS
#include <Accelerate.h>
const CBLAS_ORDER cblasOrder[2] = { CblasRowMajor, CblasColMajor };
const CBLAS_TRANSPOSE cblasTrans[2] = { CblasNoTrans, CblasTrans };
#endif


// C interfaces for Fortran Lapack routines 
extern "C" {

#ifdef LAPACK_FTN_STRING_LEN_AT_END
    int _DGEMV_(const char* transpose,
        const integer* m, const integer* n, const doublereal* alpha,
        const doublereal* a, const integer* lda, const doublereal* x,
        const integer* incX, const doublereal* beta, doublereal* y,
        const integer* incY, ftnlen trsize);
#else

    int _DGEMV_(const char* transpose, ftnlen trsize, 
        const integer* m, const integer* n, const doublereal* alpha,
        const doublereal* a, const integer* lda, const doublereal* x,
        const integer* incX, const doublereal* beta, doublereal* y,
        const integer* incY);
#endif

    int _DGETRF_(const integer* m, const integer* n, 
        doublereal* a, integer* lda, integer* ipiv, 
        integer* info);

#ifdef LAPACK_FTN_STRING_LEN_AT_END

    int _DGETRS_(const char* transpose, const integer* n, 
        const integer* nrhs, doublereal* a, const integer* lda, 
        integer* ipiv, doublereal* b, const integer* ldb, 
        integer* info, ftnlen trsize);

#else

    int _DGETRS_(const char* transpose, ftnlen trsize, const integer* n, 
        const integer* nrhs, const doublereal* a, const integer* lda, 
        integer* ipiv, doublereal* b, const integer* ldb, integer* info);

#endif

    int _DGETRI_(const integer* n, doublereal* a, const integer* lda,
        integer* ipiv, doublereal* work, integer* lwork, integer* info);

    int _DGELSS_(integer *m, integer *n, integer *nrhs, 
        doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
        s, doublereal *rcond, integer *rank, doublereal *work, integer *lwork,
        integer *info);


    int _DGBSV_(integer *n, integer *kl, integer *ku, integer *nrhs,
        doublereal *a, integer *lda, integer *ipiv, doublereal *b,
        integer *ldb, integer *info);

    int _DGBTRF_(integer* m, integer *n, integer *kl, integer *ku,
        doublereal *a, integer *lda, integer *ipiv, integer *info);

#ifdef LAPACK_FTN_STRING_LEN_AT_END
    int _DGBTRS_(const char* trans, integer *n, integer *kl, integer *ku, 
        integer *nrhs, doublereal *a, integer *lda, integer *ipiv, 
        doublereal *b, integer *ldb, integer *info, ftnlen trsize);
#else
    int _DGBTRS_(const char* trans, ftnlen trsize, 
        integer *n, integer *kl, integer *ku, 
        integer *nrhs, doublereal *a, integer *lda, integer *ipiv, 
        doublereal *b, integer *ldb, integer *info);
#endif

    int _DSCAL_(integer *n, doublereal *da, doublereal *dx, integer *incx);
  void cblas_dscal(const int N, const double alpha, double *X, const int incX);


  int _DGEQRF_(const integer* m, const integer* n, doublereal* a, const integer* lda, 
	       doublereal* tau, doublereal* work,  const integer *lwork, integer *info);

#ifdef LAPACK_FTN_STRING_LEN_AT_END
  int _DORMQR_(const char* side, const char* trans, const integer* m,  const integer* n, 
	       const integer * k, doublereal* a, const integer* lda, 
	       doublereal* tau, doublereal* c, const integer* ldc,
	       doublereal* work,  const integer *lwork, integer *info, ftnlen sisize, ftnlen trsize);
#else
  int _DORMQR_(const char* side, ftnlen sisize, const char* trans, ftnlen trsize, const integer* m, 
	       const integer* n, const integer * k, doublereal* a, const integer* lda, 
	       doublereal* tau,doublereal* c, const integer* ldc,
	       doublereal* work,  const integer *lwork, integer *info);
#endif

#ifdef LAPACK_FTN_STRING_LEN_AT_END
  int _DTRTRS_(const char* uplo, const char* trans, const char *diag,  const integer* n, 
	       const integer * nrhs, doublereal* a, const integer* lda, 
	       doublereal* b, const integer* ldb, integer *info, 
	       ftnlen upsize, ftnlen trsize, ftnlen disize);
#else
  int _DTRTRS_(const char* uplo, ftnlen upsize, const char* trans, ftnlen trsize, const char *diag,
	       ftnlen disize, const integer* n, const integer * nrhs, doublereal* a, const integer* lda, 
	       doublereal* b, const integer* ldb, integer *info); 
#endif


#ifdef LAPACK_FTN_STRING_LEN_AT_END
  int _DTRCON_(const char* norm, const char* uplo, const char *diag,  const integer* n, 
	       doublereal* a, const integer* lda, const doublereal *rcond,
	       doublereal* work, const integer* iwork, integer *info, ftnlen nosize, 
	       ftnlen upsize, ftnlen disize);
#else
  int _DTRCON_(const char* norm, ftnlen nosize, const char* uplo, ftnlen upsize, const char *diag,
	       ftnlen disize, const integer* n, doublereal* a, const integer* lda, const doublereal *rcond,
	       doublereal* work, const integer* iwork, integer *info); 
#endif


#ifdef LAPACK_FTN_STRING_LEN_AT_END
  int _DPOTRF_(const char* uplo,  const integer* n,  doublereal* a, const integer* lda, integer *info, 
	       ftnlen upsize);
#else
  int _DPOTRF_(const char* uplo, ftnlen upsize, const integer* n,  doublereal* a, const integer* lda,
	       integer *info );
#endif

#ifdef LAPACK_FTN_STRING_LEN_AT_END
  int _DPOTRS_(const char* uplo,  const integer* n,  const integer* nrhs, doublereal* a, const integer* lda, 
	       doublereal* b, const integer* ldb, integer *info, ftnlen upsize);
#else
  int _DPOTRS_(const char* uplo, ftnlen upsize, const integer* n,  const integer* nrhs, doublereal* a, const integer* lda,
	       doublereal* b, const integer* ldb, integer *info);
#endif


}
//#endif

namespace Cantera {
  //====================================================================================================================
    inline void ct_dgemv(ctlapack::storage_t storage, 
        ctlapack::transpose_t trans, 
        int m, int n, doublereal alpha, const doublereal* a, int lda, 
        const doublereal* x, int incX, doublereal beta, 
        doublereal* y, int incY) 
    {
#ifdef USE_CBLAS
        cblas_dgemv(cblasOrder[storage], cblasTrans[trans], m, n, alpha, 
            a, lda, x, incX, beta, y, incY);
#else
        integer f_m = m, f_n = n, f_lda = lda, f_incX = incX, f_incY = incY;
        doublereal f_alpha = alpha, f_beta = beta;
        ftnlen trsize = 1;
#ifdef NO_FTN_STRING_LEN_AT_END
        _DGEMV_(&no_yes[trans], &f_m, &f_n, &f_alpha, a,
            &f_lda, x, &f_incX, &f_beta, y, &f_incY);
#else
#ifdef LAPACK_FTN_STRING_LEN_AT_END
        _DGEMV_(&no_yes[trans], &f_m, &f_n, &f_alpha, a,
            &f_lda, x, &f_incX, &f_beta, y, &f_incY, trsize);
#else
        _DGEMV_(&no_yes[trans], trsize, &f_m, &f_n, &f_alpha, a,
            &f_lda, x, &f_incX, &f_beta, y, &f_incY);
#endif 
#endif
#endif

    }
        
  //====================================================================================================================
    inline void ct_dgbsv(int n, int kl, int ku, int nrhs,  
        doublereal* a, int lda, integer* ipiv, doublereal* b, int ldb, 
        int& info) {
        integer f_n = n, f_kl = kl, f_ku = ku, f_nrhs = nrhs, f_lda = lda,
               f_ldb = ldb, f_info = info;
        _DGBSV_(&f_n, &f_kl, &f_ku, &f_nrhs, a, &f_lda, ipiv, 
            b, &f_ldb, &f_info);
        info = f_info;
    }
 //====================================================================================================================
    inline void ct_dgbtrf(int m, int n, int kl, int ku, 
        doublereal* a, int lda, integer* ipiv, int& info) {
        integer f_m = m, f_n = n, f_kl = kl, f_ku = ku, 
             f_lda = lda, f_info = info;
        _DGBTRF_(&f_m, &f_n, &f_kl, &f_ku, a, &f_lda, ipiv, &f_info);
        info = f_info;
    }
  //====================================================================================================================
    inline void ct_dgbtrs(ctlapack::transpose_t trans, int n, 
        int kl, int ku, int nrhs, doublereal* a, int lda, 
        integer* ipiv, doublereal* b, int ldb, int& info) {
        integer f_n = n, f_kl = kl, f_ku = ku, f_nrhs = nrhs, f_lda = lda,
               f_ldb = ldb, f_info = info;
        char tr = no_yes[trans];
#ifdef NO_FTN_STRING_LEN_AT_END
        _DGBTRS_(&tr, &f_n, &f_kl, &f_ku, &f_nrhs, a, &f_lda, ipiv, 
            b, &f_ldb, &f_info);
#else
        ftnlen trsize = 1;
#ifdef LAPACK_FTN_STRING_LEN_AT_END
        _DGBTRS_(&tr, &f_n, &f_kl, &f_ku, &f_nrhs, a, &f_lda, ipiv, 
            b, &f_ldb, &f_info, trsize);
#else
        _DGBTRS_(&tr, trsize, &f_n, &f_kl, &f_ku, &f_nrhs, a, &f_lda, ipiv, 
            b, &f_ldb, &f_info);
#endif
#endif
        info = f_info;
    }
 //====================================================================================================================
    inline void ct_dgetrf(int m, int n, 
        doublereal* a, int lda, integer* ipiv, int& info) {
        integer mm = m;
        integer nn = n;
        integer ldaa = lda;
        integer infoo = info;
        _DGETRF_(&mm, &nn, a, &ldaa, ipiv, &infoo);
        info = infoo;
    }
  //====================================================================================================================
    inline void ct_dgetrs(ctlapack::transpose_t trans, int n, 
        int nrhs, doublereal* a, int lda,
        integer* ipiv, doublereal* b, int ldb, int& info) 
    {
        integer f_n = n, f_lda = lda, f_nrhs = nrhs, f_ldb = ldb, 
               f_info = info;
        char tr = no_yes[trans];

#ifdef NO_FTN_STRING_LEN_AT_END
        _DGETRS_(&tr, &f_n, &f_nrhs, a, &f_lda, ipiv, b, &f_ldb, 
            &f_info);
#else
        ftnlen trsize = 1;
#ifdef LAPACK_FTN_STRING_LEN_AT_END
        _DGETRS_(&tr, &f_n, &f_nrhs, a, &f_lda, ipiv, b, &f_ldb, &f_info, trsize);
#else
        _DGETRS_(&tr, trsize, &f_n, &f_nrhs, a, &f_lda, ipiv, b, &f_ldb, &f_info);
#endif
#endif
        info = f_info;
    }
  //====================================================================================================================
    inline void ct_dgetri(int n, doublereal* a, int lda, integer* ipiv,
        doublereal* work, int lwork, int& info) {
        integer f_n = n, f_lda = lda, f_lwork = lwork, f_info = info;
        _DGETRI_(&f_n, a, &f_lda, ipiv, work, &f_lwork, &f_info);
    }
  //====================================================================================================================
    inline void ct_dgelss(int m, int n, int nrhs, doublereal* a, 
        int lda, doublereal* b, int ldb, doublereal* s,
        doublereal rcond, int& rank, doublereal* work, int lwork,
        int& info) {
        doublereal f_rcond = rcond;
        integer f_m = m, f_n = n, f_nrhs = nrhs, f_lda = lda, f_ldb = ldb,
             f_rank = rank, f_info = info, f_lwork = lwork;
        //f_lwork = 2*(3*min(m,n) + max(2*min(m,n), max(m,n)));
        _DGELSS_(&f_m, &f_n, &f_nrhs, a, &f_lda, b, &f_ldb, s, &f_rcond,
            &f_rank, work, &f_lwork, &f_info);
        info = f_info;
        rank = f_rank;
    }

  inline void ct_dscal(int n, doublereal da, doublereal* dx, int incx) {
        //integer f_n = n, f_incx = incx;
        //_DSCAL_(&f_n, &da, dx, &f_incx);
    cblas_dscal(n, da, dx, incx);
  }
  //====================================================================================================================
  inline void ct_dgeqrf(int m, int n, doublereal* a, int lda, doublereal *tau,
			doublereal* work, int lwork, int &info) {
    integer f_m = m;
    integer f_n = n;
    integer f_lda = lda;
    integer f_lwork = lwork;
    integer f_info = info;
    _DGEQRF_(&f_m, &f_n, a, &f_lda, tau, work, &f_lwork, &f_info);
    info = f_info;
  }
  //====================================================================================================================
  inline void ct_dormqr(ctlapack::side_t rlside, ctlapack::transpose_t trans, int m,
                        int n, int k, doublereal* a, int lda, doublereal *tau, doublereal *c, int ldc, 
                        doublereal *work, int lwork, int &info) {
    char side = left_right[rlside];
    char tr = no_yes[trans];
    integer f_m = m;
    integer f_n = n;
    integer f_k = k;
    integer f_lwork = lwork;
    integer f_lda = lda;
    integer f_ldc = ldc;
    integer f_info = info;
#ifdef NO_FTN_STRING_LEN_AT_END
    _DORMQR_(&side, &tr, &f_m, &f_n,  &f_k, a, &f_lda, tau, c, &f_ldc, work, &f_lwork, &f_info);
#else
    ftnlen trsize = 1;
#ifdef LAPACK_FTN_STRING_LEN_AT_END
    _DORMQR_(&side, &tr, &f_m, &f_n,  &f_k, a, &f_lda, tau, c, &f_ldc, work, &f_lwork, &f_info, trsize, trsize);
#else
    _DORMQR_(&side, trsize, &tr, trsize, &f_m, &f_n,  &f_k, a, &f_lda, tau, c, &f_ldc, work, &f_lwork, &f_info);
#endif
#endif
    info = f_info;
  }
  //====================================================================================================================
  inline void ct_dtrtrs(ctlapack::upperlower_t uplot, ctlapack::transpose_t trans, const char *diag,
                        int n, int nrhs, doublereal* a, int lda, doublereal *b, int ldb, int &info) {
    char uplo = upper_lower[uplot];
    char tr = no_yes[trans];
    char dd = 'N';
    if (diag) {
      dd = diag[0];
    }
    integer f_n = n;
    integer f_nrhs = nrhs;
    integer f_lda = lda;
    integer f_ldb = ldb;
    integer f_info = info;
#ifdef NO_FTN_STRING_LEN_AT_END
    _DTRTRS_(&uplo, &tr, &dd, &f_n, &f_nrhs, a, &f_lda, b, &f_ldb, &f_info);
#else
    ftnlen trsize = 1;
#ifdef LAPACK_FTN_STRING_LEN_AT_END
    _DTRTRS_(&uplo, &tr, &dd, &f_n, &f_nrhs, a, &f_lda, b, &f_ldb, &f_info, trsize, trsize, trsize);
#else
    _DTRTRS_(&uplo, trsize, &tr, trsize, &dd, trsize, &f_n, &f_nrhs, a, &f_lda, b, &f_ldb, &f_info);
#endif
#endif
    info = f_info;
  }
  //====================================================================================================================
  //!
  /*!
   *  @param work   Must be dimensioned equal to greater than 3N
   *  @param iwork  Must be dimensioned equal to or greater than N
   */
  inline doublereal ct_dtrcon(const char *norm, ctlapack::upperlower_t uplot,  const char *diag,
			      int n, doublereal* a, int lda, doublereal *work, int *iwork, int &info) {
    char uplo = upper_lower[uplot];
    char dd = 'N';
    if (diag) {
      dd = diag[0];
    }
    char nn = '1';
    if (norm) {
      nn = norm[0];
    }
    integer f_n = n;
    integer f_lda = lda;
    integer f_info = info;
    doublereal rcond;
#ifdef NO_FTN_STRING_LEN_AT_END
    _DTRCON_(&nn, &uplo, &dd, &f_n, a, &f_lda, &rcond, work, iwork, &f_info);
#else
    ftnlen trsize = 1;
#ifdef LAPACK_FTN_STRING_LEN_AT_END
    _DTRCON_(&nn, &uplo, &dd, &f_n, a, &f_lda, &rcond, work, iwork, &f_info, trsize, trsize, trsize);
#else
    _DTRCON_(&nn, trsize, &uplo, trsize, &dd, trsize, &f_n, a, &f_lda, &rcond, work, iwork, &f_info);
#endif
#endif
    info = f_info;
    return rcond;
  }
  //====================================================================================================================
  //!
  /*!
   *  @param work   Must be dimensioned equal to greater than 3N
   *  @param iwork  Must be dimensioned equal to or greater than N
   */
  inline void ct_dpotrf(ctlapack::upperlower_t uplot, int n, doublereal* a, int lda, int &info) {
    char uplo = upper_lower[uplot];
    integer f_n = n;
    integer f_lda = lda;
    integer f_info = info;

#ifdef NO_FTN_STRING_LEN_AT_END
    _DPOTRF_(&uplo, &f_n, a, &f_lda, &f_info);
#else
    ftnlen trsize = 1;
#ifdef LAPACK_FTN_STRING_LEN_AT_END
    _DPOTRF_(&uplo, &f_n, a, &f_lda, &f_info, trsize);
#else
    _DPOTRF_(&uplo, trsize, &f_n, a, &f_lda, &f_info);
#endif
#endif
    info = f_info;
    return;
  }
  //====================================================================================================================
  //!
  /*!
   */
  inline void ct_dpotrs(ctlapack::upperlower_t uplot, int n, int nrhs, doublereal* a, int lda, 
			doublereal* b, int ldb, int &info) {
    char uplo = upper_lower[uplot];
    integer f_n = n;
    integer f_nrhs = nrhs;
    integer f_lda = lda;
    integer f_ldb = ldb;
    integer f_info = info;

#ifdef NO_FTN_STRING_LEN_AT_END
    _DPOTRS_(&uplo, &f_n, &f_nrhs, a, &f_lda, b, &f_ldb, &f_info);
#else
    ftnlen trsize = 1;
#ifdef LAPACK_FTN_STRING_LEN_AT_END
    _DPOTRS_(&uplo, &f_n, &f_nrhs, a, &f_lda, b, &f_ldb, &f_info, trsize);
#else
    _DPOTRS_(&uplo, trsize, &f_n, &f_nrhs, a, &f_lda, b, &f_ldb, &f_info);
#endif
#endif
    info = f_info;
    return;
  }
  //====================================================================================================================


}

#endif
