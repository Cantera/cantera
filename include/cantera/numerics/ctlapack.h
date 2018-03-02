/**
 * @file ctlapack.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_CTLAPACK_H
#define CT_CTLAPACK_H

#include "cantera/base/ct_defs.h"

// map BLAS names to names with or without a trailing underscore.
#ifndef LAPACK_FTN_TRAILING_UNDERSCORE

#define _DGEMV_   dgemv
#define _DGETRF_  dgetrf
#define _DGETRS_  dgetrs
#define _DGETRI_  dgetri
#define _DGELSS_  dgelss
#define _DGBCON_  dgbcon
#define _DGBSV_   dgbsv
#define _DGBTRF_  dgbtrf
#define _DGBTRS_  dgbtrs
#define _DGECON_  dgecon
#define _DLANGE_  dlange
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
#define _DGBCON_  dgbcon_
#define _DGBSV_   dgbsv_
#define _DGBTRF_  dgbtrf_
#define _DGBTRS_  dgbtrs_
#define _DGECON_  dgecon_
#define _DLANGE_  dlange_
#define _DSCAL_   dscal_
#define _DGEQRF_  dgeqrf_
#define _DORMQR_  dormqr_
#define _DTRTRS_  dtrtrs_
#define _DTRCON_  dtrcon_
#define _DPOTRF_  dpotrf_
#define _DPOTRS_  dpotrs_

#endif


namespace ctlapack
{
typedef enum {Transpose = 1, NoTranspose = 0} transpose_t;
typedef enum {ColMajor = 1, RowMajor = 0} storage_t;
typedef enum {UpperTriangular = 0, LowerTriangular = 1} upperlower_t;
typedef enum {Left = 0, Right = 1} side_t;
}
const char no_yes[2] = {'N', 'T'};
const char upper_lower[2] = {'U', 'L'};
const char left_right[2] = {'L', 'R'};

// C interfaces for Fortran Lapack routines
extern "C" {

#ifdef LAPACK_FTN_STRING_LEN_AT_END
    int _DGEMV_(const char* transpose,
                const integer* m, const integer* n, const double* alpha,
                const double* a, const integer* lda, const double* x,
                const integer* incX, const double* beta, double* y,
                const integer* incY, ftnlen trsize);
#else
    int _DGEMV_(const char* transpose, ftnlen trsize,
                const integer* m, const integer* n, const double* alpha,
                const double* a, const integer* lda, const double* x,
                const integer* incX, const double* beta, double* y,
                const integer* incY);
#endif

    int _DGETRF_(const integer* m, const integer* n,
                 double* a, integer* lda, integer* ipiv,
                 integer* info);

#ifdef LAPACK_FTN_STRING_LEN_AT_END
    int _DGETRS_(const char* transpose, const integer* n,
                 const integer* nrhs, double* a, const integer* lda,
                 integer* ipiv, double* b, const integer* ldb,
                 integer* info, ftnlen trsize);
#else
    int _DGETRS_(const char* transpose, ftnlen trsize, const integer* n,
                 const integer* nrhs, const double* a, const integer* lda,
                 integer* ipiv, double* b, const integer* ldb, integer* info);
#endif

    int _DGETRI_(const integer* n, double* a, const integer* lda,
                 integer* ipiv, double* work, integer* lwork, integer* info);

    int _DGELSS_(const integer* m, const integer* n, const integer* nrhs,
                 double* a, const integer* lda, double* b,
                 const integer* ldb, double* s, const double* rcond,
                 integer* rank, double* work, integer* lwork, integer* info);

    int _DGBSV_(integer* n, integer* kl, integer* ku, integer* nrhs,
                double* a, integer* lda, integer* ipiv, double* b,
                integer* ldb, integer* info);

    int _DGBTRF_(integer* m, integer* n, integer* kl, integer* ku,
                 double* a, integer* lda, integer* ipiv, integer* info);

#ifdef LAPACK_FTN_STRING_LEN_AT_END
    int _DGBTRS_(const char* trans, integer* n, integer* kl, integer* ku,
                 integer* nrhs, double* a, integer* lda, integer* ipiv,
                 double* b, integer* ldb, integer* info, ftnlen trsize);
#else
    int _DGBTRS_(const char* trans, ftnlen trsize,
                 integer* n, integer* kl, integer* ku,
                 integer* nrhs, double* a, integer* lda, integer* ipiv,
                 double* b, integer* ldb, integer* info);
#endif

    int _DSCAL_(integer* n, double* da, double* dx, integer* incx);

    int _DGEQRF_(const integer* m, const integer* n, double* a, const integer* lda,
                 double* tau, double* work, const integer* lwork, integer* info);

#ifdef LAPACK_FTN_STRING_LEN_AT_END
    int _DORMQR_(const char* side, const char* trans, const integer* m, const integer* n,
                 const integer* k, double* a, const integer* lda,
                 double* tau, double* c, const integer* ldc,
                 double* work, const integer* lwork, integer* info, ftnlen sisize, ftnlen trsize);
#else
    int _DORMQR_(const char* side, ftnlen sisize, const char* trans, ftnlen trsize, const integer* m,
                 const integer* n, const integer* k, double* a, const integer* lda,
                 double* tau,double* c, const integer* ldc,
                 double* work, const integer* lwork, integer* info);
#endif

#ifdef LAPACK_FTN_STRING_LEN_AT_END
    int _DTRTRS_(const char* uplo, const char* trans, const char* diag, const integer* n,
                 const integer* nrhs, double* a, const integer* lda,
                 double* b, const integer* ldb, integer* info,
                 ftnlen upsize, ftnlen trsize, ftnlen disize);
#else
    int _DTRTRS_(const char* uplo, ftnlen upsize, const char* trans, ftnlen trsize, const char* diag,
                 ftnlen disize, const integer* n, const integer* nrhs, double* a, const integer* lda,
                 double* b, const integer* ldb, integer* info);
#endif


#ifdef LAPACK_FTN_STRING_LEN_AT_END
    int _DTRCON_(const char* norm, const char* uplo, const char* diag, const integer* n,
                 double* a, const integer* lda, const double* rcond,
                 double* work, const integer* iwork, integer* info, ftnlen nosize,
                 ftnlen upsize, ftnlen disize);
#else
    int _DTRCON_(const char* norm, ftnlen nosize, const char* uplo, ftnlen upsize, const char* diag,
                 ftnlen disize, const integer* n, double* a, const integer* lda, const double* rcond,
                 double* work, const integer* iwork, integer* info);
#endif


#ifdef LAPACK_FTN_STRING_LEN_AT_END
    int _DPOTRF_(const char* uplo, const integer* n, double* a, const integer* lda, integer* info,
                 ftnlen upsize);
#else
    int _DPOTRF_(const char* uplo, ftnlen upsize, const integer* n, double* a, const integer* lda,
                 integer* info);
#endif

#ifdef LAPACK_FTN_STRING_LEN_AT_END
    int _DPOTRS_(const char* uplo, const integer* n, const integer* nrhs, double* a, const integer* lda,
                 double* b, const integer* ldb, integer* info, ftnlen upsize);
#else
    int _DPOTRS_(const char* uplo, ftnlen upsize, const integer* n, const integer* nrhs, double* a, const integer* lda,
                 double* b, const integer* ldb, integer* info);
#endif

#ifdef LAPACK_FTN_STRING_LEN_AT_END
    int _DGECON_(const char* norm, const integer* n, double* a, const integer* lda,
                 const double* rnorm, const double* rcond,
                 double* work, const integer* iwork, integer* info, ftnlen nosize);
#else
    int _DGECON_(const char* norm, ftnlen nosize, const integer* n, double* a, const integer* lda,
                 const double* rnorm, const double* rcond,
                 double* work, const integer* iwork, integer* info);
#endif

#ifdef LAPACK_FTN_STRING_LEN_AT_END
    int _DGBCON_(const char* norm, const integer* n, integer* kl, integer* ku, double* ab, const integer* ldab,
                 const integer* ipiv, const double* anorm, const double* rcond,
                 double* work, const integer* iwork, integer* info, ftnlen nosize);
#else
    int _DGBCON_(const char* norm, ftnlen nosize, const integer* n, integer* kl, integer* ku, double* ab, const integer* ldab,
                 const integer* ipiv, const double* anorm, const double* rcond,
                 double* work, const integer* iwork, integer* info);
#endif

#ifdef LAPACK_FTN_STRING_LEN_AT_END
    double _DLANGE_(const char* norm, const integer* m, const integer* n, double* a, const integer* lda,
                        double* work, ftnlen nosize);
#else
    double _DLANGE_(const char* norm, ftnlen nosize, const integer* m, const integer* n, double* a, const integer* lda,
                        double* work);
#endif
}

namespace Cantera
{

inline void ct_dgemv(ctlapack::storage_t storage,
                     ctlapack::transpose_t trans,
                     int m, int n, double alpha, const double* a, int lda,
                     const double* x, int incX, double beta,
                     double* y, int incY)
{
    integer f_m = m, f_n = n, f_lda = lda, f_incX = incX, f_incY = incY;
    double f_alpha = alpha, f_beta = beta;
    ftnlen trsize = 1;
#ifdef LAPACK_FTN_STRING_LEN_AT_END
    _DGEMV_(&no_yes[trans], &f_m, &f_n, &f_alpha, a,
            &f_lda, x, &f_incX, &f_beta, y, &f_incY, trsize);
#else
    _DGEMV_(&no_yes[trans], trsize, &f_m, &f_n, &f_alpha, a,
            &f_lda, x, &f_incX, &f_beta, y, &f_incY);
#endif
}

inline void ct_dgbsv(int n, int kl, int ku, int nrhs,
                     double* a, int lda, integer* ipiv, double* b, int ldb,
                     int& info)
{
    integer f_n = n, f_kl = kl, f_ku = ku, f_nrhs = nrhs, f_lda = lda,
            f_ldb = ldb, f_info = 0;
    _DGBSV_(&f_n, &f_kl, &f_ku, &f_nrhs, a, &f_lda, ipiv,
            b, &f_ldb, &f_info);
    info = f_info;
}

inline void ct_dgelss(size_t m, size_t n, size_t nrhs, double* a,
                      size_t lda, double* b, size_t ldb, double* s,
                      double rcond, size_t& rank, double* work,
                      int& lwork, int& info)
{
    integer f_m = static_cast<integer>(m);
    integer f_n = static_cast<integer>(n);
    integer f_nrhs = static_cast<integer>(nrhs);
    integer f_lda = static_cast<integer>(lda);
    integer f_ldb = static_cast<integer>(ldb);
    integer f_lwork = static_cast<integer>(lwork);
    integer f_rank, f_info;

    _DGELSS_(&f_m, &f_n, &f_nrhs, a, &f_lda, b, &f_ldb, s, &rcond, &f_rank,
            work, &f_lwork, &f_info);

    info = static_cast<int>(f_info);
    rank = static_cast<size_t>(f_rank);
    lwork = static_cast<int>(f_lwork);
}

inline void ct_dgbtrf(size_t m, size_t n, size_t kl, size_t ku,
                      double* a, size_t lda, integer* ipiv, int& info)
{
    integer f_m = (int) m;
    integer f_n = (int) n;
    integer f_kl = (int) kl;
    integer f_ku = (int) ku;
    integer f_lda = (int) lda;
    integer f_info = 0;
    _DGBTRF_(&f_m, &f_n, &f_kl, &f_ku, a, &f_lda, ipiv, &f_info);
    info = f_info;
}

inline void ct_dgbtrs(ctlapack::transpose_t trans, size_t n,
                      size_t kl, size_t ku, size_t nrhs, double* a, size_t lda,
                      integer* ipiv, double* b, size_t ldb, int& info)
{
    integer f_n = (int) n;
    integer f_kl = (int) kl;
    integer f_ku = (int) ku;
    integer f_nrhs = (int) nrhs;
    integer f_lda = (int) lda;
    integer f_ldb = (int) ldb;
    integer f_info = 0;
    char tr = no_yes[trans];
    ftnlen trsize = 1;
#ifdef LAPACK_FTN_STRING_LEN_AT_END
    _DGBTRS_(&tr, &f_n, &f_kl, &f_ku, &f_nrhs, a, &f_lda, ipiv,
             b, &f_ldb, &f_info, trsize);
#else
    _DGBTRS_(&tr, trsize, &f_n, &f_kl, &f_ku, &f_nrhs, a, &f_lda, ipiv,
             b, &f_ldb, &f_info);
#endif
    info = f_info;
}

inline void ct_dgetrf(size_t m, size_t n,
                      double* a, size_t lda, integer* ipiv, int& info)
{
    integer mm = (int) m;
    integer nn = (int) n;
    integer ldaa = (int) lda;
    integer infoo = 0;
    _DGETRF_(&mm, &nn, a, &ldaa, ipiv, &infoo);
    info = infoo;
}

inline void ct_dgetrs(ctlapack::transpose_t trans, size_t n,
                      size_t nrhs, double* a, size_t lda,
                      integer* ipiv, double* b, size_t ldb, int& info)
{
    integer f_n = (int) n;
    integer f_lda = (int) lda;
    integer f_nrhs = (int) nrhs;
    integer f_ldb = (int) ldb;
    integer f_info = 0;
    char tr = no_yes[trans];
    ftnlen trsize = 1;
#ifdef LAPACK_FTN_STRING_LEN_AT_END
    _DGETRS_(&tr, &f_n, &f_nrhs, a, &f_lda, ipiv, b, &f_ldb, &f_info, trsize);
#else
    _DGETRS_(&tr, trsize, &f_n, &f_nrhs, a, &f_lda, ipiv, b, &f_ldb, &f_info);
#endif
    info = f_info;
}

inline void ct_dgetri(int n, double* a, int lda, integer* ipiv,
                      double* work, int lwork, int& info)
{
    integer f_n = n, f_lda = lda, f_lwork = lwork, f_info = 0;
    _DGETRI_(&f_n, a, &f_lda, ipiv, work, &f_lwork, &f_info);
}

inline void ct_dscal(int n, double da, double* dx, int incx)
{
    integer f_n = n, f_incx = incx;
    _DSCAL_(&f_n, &da, dx, &f_incx);
}

inline void ct_dgeqrf(size_t m, size_t n, double* a, size_t lda, double* tau,
                      double* work, size_t lwork, int& info)
{
    integer f_m = static_cast<integer>(m);
    integer f_n = static_cast<integer>(n);
    integer f_lda = static_cast<integer>(lda);
    integer f_lwork = static_cast<integer>(lwork);
    integer f_info = 0;
    _DGEQRF_(&f_m, &f_n, a, &f_lda, tau, work, &f_lwork, &f_info);
    info = f_info;
}

inline void ct_dormqr(ctlapack::side_t rlside, ctlapack::transpose_t trans, size_t m,
                      size_t n, size_t k, double* a, size_t lda, double* tau, double* c, size_t ldc,
                      double* work, size_t lwork, int& info)
{
    char side = left_right[rlside];
    char tr = no_yes[trans];
    integer f_m = static_cast<integer>(m);
    integer f_n = static_cast<integer>(n);
    integer f_k = static_cast<integer>(k);
    integer f_lwork = static_cast<integer>(lwork);
    integer f_lda = static_cast<integer>(lda);
    integer f_ldc = static_cast<integer>(ldc);
    integer f_info = 0;
    ftnlen trsize = 1;
#ifdef LAPACK_FTN_STRING_LEN_AT_END
    _DORMQR_(&side, &tr, &f_m, &f_n, &f_k, a, &f_lda, tau, c, &f_ldc, work, &f_lwork, &f_info, trsize, trsize);
#else
    _DORMQR_(&side, trsize, &tr, trsize, &f_m, &f_n, &f_k, a, &f_lda, tau, c, &f_ldc, work, &f_lwork, &f_info);
#endif
    info = f_info;
}

inline void ct_dtrtrs(ctlapack::upperlower_t uplot, ctlapack::transpose_t trans, const char* diag,
                      size_t n, size_t nrhs, double* a, size_t lda, double* b, size_t ldb, int& info)
{
    char uplo = upper_lower[uplot];
    char tr = no_yes[trans];
    char dd = 'N';
    if (diag) {
        dd = diag[0];
    }
    integer f_n = static_cast<integer>(n);
    integer f_nrhs = static_cast<integer>(nrhs);
    integer f_lda = static_cast<integer>(lda);
    integer f_ldb = static_cast<integer>(ldb);
    integer f_info = 0;
    ftnlen trsize = 1;
#ifdef LAPACK_FTN_STRING_LEN_AT_END
    _DTRTRS_(&uplo, &tr, &dd, &f_n, &f_nrhs, a, &f_lda, b, &f_ldb, &f_info, trsize, trsize, trsize);
#else
    _DTRTRS_(&uplo, trsize, &tr, trsize, &dd, trsize, &f_n, &f_nrhs, a, &f_lda, b, &f_ldb, &f_info);
#endif
    info = f_info;
}

/*!
 *  @param work   Must be dimensioned equal to greater than 3N
 *  @param iwork  Must be dimensioned equal to or greater than N
 */
inline double ct_dtrcon(const char* norm, ctlapack::upperlower_t uplot, const char* diag,
                            size_t n, double* a, size_t lda, double* work, int* iwork, int& info)
{
    char uplo = upper_lower[uplot];
    char dd = 'N';
    if (diag) {
        dd = diag[0];
    }
    char nn = '1';
    if (norm) {
        nn = norm[0];
    }
    integer f_n = static_cast<integer>(n);
    integer f_lda = static_cast<integer>(lda);
    integer f_info = 0;
    double rcond = 0.0;
    ftnlen trsize = 1;
#ifdef LAPACK_FTN_STRING_LEN_AT_END
    _DTRCON_(&nn, &uplo, &dd, &f_n, a, &f_lda, &rcond, work, iwork, &f_info, trsize, trsize, trsize);
#else
    _DTRCON_(&nn, trsize, &uplo, trsize, &dd, trsize, &f_n, a, &f_lda, &rcond, work, iwork, &f_info);
#endif
    info = f_info;
    return rcond;
}

inline void ct_dpotrf(ctlapack::upperlower_t uplot, size_t n, double* a, size_t lda, int& info)
{
    char uplo = upper_lower[uplot];
    integer f_n = static_cast<integer>(n);
    integer f_lda = static_cast<integer>(lda);
    integer f_info = 0;
    ftnlen trsize = 1;
#ifdef LAPACK_FTN_STRING_LEN_AT_END
    _DPOTRF_(&uplo, &f_n, a, &f_lda, &f_info, trsize);
#else
    _DPOTRF_(&uplo, trsize, &f_n, a, &f_lda, &f_info);
#endif
    info = f_info;
    return;
}

inline void ct_dpotrs(ctlapack::upperlower_t uplot, size_t n, size_t nrhs, double* a, size_t lda,
                      double* b, size_t ldb, int& info)
{
    char uplo = upper_lower[uplot];
    integer f_n = static_cast<integer>(n);
    integer f_nrhs = static_cast<integer>(nrhs);
    integer f_lda = static_cast<integer>(lda);
    integer f_ldb = static_cast<integer>(ldb);
    integer f_info = 0;
    ftnlen trsize = 1;
#ifdef LAPACK_FTN_STRING_LEN_AT_END
    _DPOTRS_(&uplo, &f_n, &f_nrhs, a, &f_lda, b, &f_ldb, &f_info, trsize);
#else
    _DPOTRS_(&uplo, trsize, &f_n, &f_nrhs, a, &f_lda, b, &f_ldb, &f_info);
#endif
    info = f_info;
    return;
}

inline double ct_dgecon(const char norm, size_t n, double* a, size_t lda, double anorm,
                            double* work, int* iwork, int& info)
{
    char cnorm = '1';
    if (norm) {
        cnorm = norm;
    }
    integer f_n = static_cast<integer>(n);
    integer f_lda = static_cast<integer>(lda);
    integer f_info = 0;
    double rcond = 0.0;
    ftnlen trsize = 1;
#ifdef LAPACK_FTN_STRING_LEN_AT_END
    _DGECON_(&cnorm, &f_n, a, &f_lda, &anorm, &rcond, work, iwork, &f_info, trsize);
#else
    _DGECON_(&cnorm, trsize, &f_n, a, &f_lda, &anorm, &rcond, work, iwork, &f_info);
#endif
    info = f_info;
    return rcond;
}

inline double ct_dgbcon(const char norm, size_t n, size_t kl, size_t ku,
                            double* a, size_t ldab, int* ipiv, double anorm,
                            double* work, int* iwork, int& info)
{
    char cnorm = '1';
    if (norm) {
        cnorm = norm;
    }
    integer f_n = static_cast<integer>(n);
    integer f_kl = static_cast<integer>(kl);
    integer f_ku = static_cast<integer>(ku);
    integer f_ldab = static_cast<integer>(ldab);
    integer f_info = 0;
    double rcond = 0.0;
    ftnlen trsize = 1;
#ifdef LAPACK_FTN_STRING_LEN_AT_END
    _DGBCON_(&cnorm, &f_n, &f_kl, &f_ku, a, &f_ldab, ipiv, &anorm, &rcond, work, iwork, &f_info, trsize);
#else
    _DGBCON_(&cnorm, trsize, &f_n, &f_kl, &f_ku, a, &f_ldab, ipiv, &anorm, &rcond, work, iwork, &f_info);
#endif
    info = f_info;
    return rcond;
}

inline double ct_dlange(const char norm, size_t m, size_t n, double* a, size_t lda,
                            double* work)
{
    char cnorm = '1';
    if (norm) {
        cnorm = norm;
    }
    integer f_m = static_cast<integer>(m);
    integer f_n = static_cast<integer>(n);
    integer f_lda = static_cast<integer>(lda);
    double anorm;
    ftnlen trsize = 1;
#ifdef LAPACK_FTN_STRING_LEN_AT_END
    anorm = _DLANGE_(&cnorm, &f_m, &f_n, a, &f_lda, work, trsize);
#else
    anorm = _DLANGE_(&cnorm, trsize, &f_m, &f_n, a, &f_lda, work);
#endif
    return anorm;
}

}

#endif
