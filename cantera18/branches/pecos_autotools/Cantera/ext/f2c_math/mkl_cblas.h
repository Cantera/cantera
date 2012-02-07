/*
//               INTEL CORPORATION PROPRIETARY INFORMATION
//  This software is supplied under the terms of a license agreement or
//  nondisclosure agreement with Intel Corporation and may not be copied
//  or disclosed except in accordance with the terms of that agreement.
//    Copyright 1999, 2000 Intel Corporation. All Rights Reserved.
//
// File    : mkl_cblas.h
// Purpose : MKL CBLAS interface
// Author  : Shemyakin Andrey
*/

#ifndef __MKL_CBLAS_H__
#define __MKL_CBLAS_H__
#include <stddef.h>

#ifdef __cplusplus
extern "C" {            /* Assume C declarations for C++ */
#endif /* __cplusplus */

/*
 * Enumerated and derived types
 */
#define CBLAS_INDEX size_t  /* this may vary between platforms */

typedef enum {CblasRowMajor=101, CblasColMajor=102} CBLAS_ORDER;
typedef enum {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113} CBLAS_TRANSPOSE;
typedef enum {CblasUpper=121, CblasLower=122} CBLAS_UPLO;
typedef enum {CblasNonUnit=131, CblasUnit=132} CBLAS_DIAG;
typedef enum {CblasLeft=141, CblasRight=142} CBLAS_SIDE;

/*
 * ===========================================================================
 * Prototypes for level 1 BLAS functions (complex are recast as routines)
 * ===========================================================================
 */

float  cblas_sdot(const int N, const float  *X, const int incX,
                  const float  *Y, const int incY);
float  cblas_sdoti(const int N, const float *X, const int *indx,
                   const float *Y);
double cblas_ddot(const int N, const double *X, const int incX,
                  const double *Y, const int incY);
double cblas_ddoti(const int N, const double *X, const int *indx,
                   const double *Y);

/*
 * Functions having prefixes Z and C only
 */
void   cblas_cdotu_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotu);
void   cblas_cdotui_sub(const int N, const void *X, const int *indx,
                        const void *Y, void *dotui);
void   cblas_cdotc_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotc);
void   cblas_cdotci_sub(const int N, const void *X, const int *indx,
                        const void *Y, void *dotui);

void   cblas_zdotu_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotu);
void   cblas_zdotui_sub(const int N, const void *X, const int *indx,
                        const void *Y, void *dotui);
void   cblas_zdotc_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotc);
void   cblas_zdotci_sub(const int N, const void *X, const int *indx,
                        const void *Y, void *dotui);


/*
 * Functions having prefixes S D SC DZ
 */
float  cblas_snrm2(const int N, const float *X, const int incX);
float  cblas_sasum(const int N, const float *X, const int incX);

double cblas_dnrm2(const int N, const double *X, const int incX);
double cblas_dasum(const int N, const double *X, const int incX);

float  cblas_scnrm2(const int N, const void *X, const int incX);
float  cblas_scasum(const int N, const void *X, const int incX);

double cblas_dznrm2(const int N, const void *X, const int incX);
double cblas_dzasum(const int N, const void *X, const int incX);


/*
 * Functions having standard 4 prefixes (S D C Z)
 */
CBLAS_INDEX cblas_isamax(const int N, const float  *X, const int incX);
CBLAS_INDEX cblas_idamax(const int N, const double *X, const int incX);
CBLAS_INDEX cblas_icamax(const int N, const void   *X, const int incX);
CBLAS_INDEX cblas_izamax(const int N, const void   *X, const int incX);
CBLAS_INDEX cblas_isamin(const int N, const float  *X, const int incX);
CBLAS_INDEX cblas_idamin(const int N, const double *X, const int incX);
CBLAS_INDEX cblas_icamin(const int N, const void   *X, const int incX);
CBLAS_INDEX cblas_izamin(const int N, const void   *X, const int incX);

/*
 * ===========================================================================
 * Prototypes for level 1 BLAS routines
 * ===========================================================================
 */

/*
 * Routines with standard 4 prefixes (s, d, c, z)
 */
void cblas_sswap(const int N, float *X, const int incX,
                 float *Y, const int incY);
void cblas_scopy(const int N, const float *X, const int incX,
                 float *Y, const int incY);
void cblas_saxpy(const int N, const float alpha, const float *X,
                 const int incX, float *Y, const int incY);
void cblas_saxpyi(const int N, const float alpha, const float *X,
                 const int *indx, float *Y);
void cblas_sgthr(const int N, const float *Y, float *X, 
                                 const int *indx);
void cblas_sgthrz(const int N, float *Y, float *X, 
                                 const int *indx);
void cblas_ssctr(const int N, const float *X, const int *indx, 
                                 float *Y);
void cblas_srotg(float *a, float *b, float *c, float *s);

void cblas_dswap(const int N, double *X, const int incX,
                 double *Y, const int incY);
void cblas_dcopy(const int N, const double *X, const int incX,
                 double *Y, const int incY);
void cblas_daxpy(const int N, const double alpha, const double *X,
                 const int incX, double *Y, const int incY);
void cblas_daxpyi(const int N, const double alpha, const double *X,
                 const int *indx, double *Y);
void cblas_dgthr(const int N, const double *Y, double *X, 
                                 const int *indx);
void cblas_dgthrz(const int N, double *Y, double *X, 
                                 const int *indx);
void cblas_dsctr(const int N, const double *X, const int *indx, 
                                 double *Y);
void cblas_drotg(double *a, double *b, double *c, double *s);

void cblas_cswap(const int N, void *X, const int incX,
                 void *Y, const int incY);
void cblas_ccopy(const int N, const void *X, const int incX,
                 void *Y, const int incY);
void cblas_caxpy(const int N, const void *alpha, const void *X,
                 const int incX, void *Y, const int incY);
void cblas_caxpyi(const int N, const void *alpha, const void *X,
                 const int *indx, void *Y);
void cblas_cgthr(const int N, const void *Y, void *X, 
                                 const int *indx);
void cblas_cgthrz(const int N, void *Y, void *X, 
                                 const int *indx);
void cblas_csctr(const int N, const void *X, const int *indx, 
                                 void *Y);
void cblas_crotg(void *a, const void *b, float *c, void *s);

void cblas_zswap(const int N, void *X, const int incX,
                 void *Y, const int incY);
void cblas_zcopy(const int N, const void *X, const int incX,
                 void *Y, const int incY);
void cblas_zaxpy(const int N, const void *alpha, const void *X,
                 const int incX, void *Y, const int incY);
void cblas_zaxpyi(const int N, const void *alpha, const void *X,
                 const int *indx, void *Y);
void cblas_zgthr(const int N, const void *Y, void *X, 
                                 const int *indx);
void cblas_zgthrz(const int N, void *Y, void *X, 
                                 const int *indx);
void cblas_zsctr(const int N, const void *X, const int *indx, 
                                 void *Y);
void cblas_zrotg(void *a, const void *b, double *c, void *s);

/*
 * Routines with S and D prefix only
 */
void cblas_srotmg(float *d1, float *d2, float *b1, const float *b2, float *P);
void cblas_srot(const int N, float *X, const int incX,
                float *Y, const int incY, const float c, const float s);
void cblas_sroti(const int N, float *X, const int *indx,
                float *Y, const float c, const float s);
void cblas_srotm(const int N, float *X, const int incX,
                float *Y, const int incY, const float *P);

void cblas_drotmg(double *d1, double *d2, double *b1, const double *b2, double *P);
void cblas_drot(const int N, double *X, const int incX,
                double *Y, const int incY, const double c, const double  s);
void cblas_drotm(const int N, double *X, const int incX,
                double *Y, const int incY, const double *P);
void cblas_droti(const int N, double *X, const int *indx,
                double *Y, const double c, const double s);

/*
 * Routines with CS and ZD prefix only
 */
void cblas_csrot(const int N, void *X, const int incX,
                 void *Y, const int incY, const float c, const float s);
void cblas_zdrot(const int N, void *X, const int incX,
                 void *Y, const int incY, const double c, const double  s);

/*
 * Routines with S D C Z CS and ZD prefixes
 */
void cblas_sscal(const int N, const float alpha, float *X, const int incX);
void cblas_dscal(const int N, const double alpha, double *X, const int incX);
void cblas_cscal(const int N, const void *alpha, void *X, const int incX);
void cblas_zscal(const int N, const void *alpha, void *X, const int incX);
void cblas_csscal(const int N, const float alpha, void *X, const int incX);
void cblas_zdscal(const int N, const double alpha, void *X, const int incX);

/*
 * ===========================================================================
 * Prototypes for level 2 BLAS
 * ===========================================================================
 */

/*
 * Routines with standard 4 prefixes (S, D, C, Z)
 */
void cblas_sgemv(const  CBLAS_ORDER order,
                 const  CBLAS_TRANSPOSE TransA, const int M, const int N,
                 const float alpha, const float *A, const int lda,
                 const float *X, const int incX, const float beta,
                 float *Y, const int incY);
void cblas_sgbmv(const  CBLAS_ORDER order,
                 const  CBLAS_TRANSPOSE TransA, const int M, const int N,
                 const int KL, const int KU, const float alpha,
                 const float *A, const int lda, const float *X,
                 const int incX, const float beta, float *Y, const int incY);
void cblas_strmv(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo,
                 const  CBLAS_TRANSPOSE TransA, const  CBLAS_DIAG Diag,
                 const int N, const float *A, const int lda,
                 float *X, const int incX);
void cblas_stbmv(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo,
                 const  CBLAS_TRANSPOSE TransA, const  CBLAS_DIAG Diag,
                 const int N, const int K, const float *A, const int lda,
                 float *X, const int incX);
void cblas_stpmv(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo,
                 const  CBLAS_TRANSPOSE TransA, const  CBLAS_DIAG Diag,
                 const int N, const float *Ap, float *X, const int incX);
void cblas_strsv(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo,
                 const  CBLAS_TRANSPOSE TransA, const  CBLAS_DIAG Diag,
                 const int N, const float *A, const int lda, float *X,
                 const int incX);
void cblas_stbsv(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo,
                 const  CBLAS_TRANSPOSE TransA, const  CBLAS_DIAG Diag,
                 const int N, const int K, const float *A, const int lda,
                 float *X, const int incX);
void cblas_stpsv(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo,
                 const  CBLAS_TRANSPOSE TransA, const  CBLAS_DIAG Diag,
                 const int N, const float *Ap, float *X, const int incX);

void cblas_dgemv(const  CBLAS_ORDER order,
                 const  CBLAS_TRANSPOSE TransA, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 const double *X, const int incX, const double beta,
                 double *Y, const int incY);
void cblas_dgbmv(const  CBLAS_ORDER order,
                 const  CBLAS_TRANSPOSE TransA, const int M, const int N,
                 const int KL, const int KU, const double alpha,
                 const double *A, const int lda, const double *X,
                 const int incX, const double beta, double *Y, const int incY);
void cblas_dtrmv(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo,
                 const  CBLAS_TRANSPOSE TransA, const  CBLAS_DIAG Diag,
                 const int N, const double *A, const int lda,
                 double *X, const int incX);
void cblas_dtbmv(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo,
                 const  CBLAS_TRANSPOSE TransA, const  CBLAS_DIAG Diag,
                 const int N, const int K, const double *A, const int lda,
                 double *X, const int incX);
void cblas_dtpmv(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo,
                 const  CBLAS_TRANSPOSE TransA, const  CBLAS_DIAG Diag,
                 const int N, const double *Ap, double *X, const int incX);
void cblas_dtrsv(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo,
                 const  CBLAS_TRANSPOSE TransA, const  CBLAS_DIAG Diag,
                 const int N, const double *A, const int lda, double *X,
                 const int incX);
void cblas_dtbsv(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo,
                 const  CBLAS_TRANSPOSE TransA, const  CBLAS_DIAG Diag,
                 const int N, const int K, const double *A, const int lda,
                 double *X, const int incX);
void cblas_dtpsv(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo,
                 const  CBLAS_TRANSPOSE TransA, const  CBLAS_DIAG Diag,
                 const int N, const double *Ap, double *X, const int incX);

void cblas_cgemv(const  CBLAS_ORDER order,
                 const  CBLAS_TRANSPOSE TransA, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 const void *X, const int incX, const void *beta,
                 void *Y, const int incY);
void cblas_cgbmv(const  CBLAS_ORDER order,
                 const  CBLAS_TRANSPOSE TransA, const int M, const int N,
                 const int KL, const int KU, const void *alpha,
                 const void *A, const int lda, const void *X,
                 const int incX, const void *beta, void *Y, const int incY);
void cblas_ctrmv(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo,
                 const  CBLAS_TRANSPOSE TransA, const  CBLAS_DIAG Diag,
                 const int N, const void *A, const int lda,
                 void *X, const int incX);
void cblas_ctbmv(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo,
                 const  CBLAS_TRANSPOSE TransA, const  CBLAS_DIAG Diag,
                 const int N, const int K, const void *A, const int lda,
                 void *X, const int incX);
void cblas_ctpmv(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo,
                 const  CBLAS_TRANSPOSE TransA, const  CBLAS_DIAG Diag,
                 const int N, const void *Ap, void *X, const int incX);
void cblas_ctrsv(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo,
                 const  CBLAS_TRANSPOSE TransA, const  CBLAS_DIAG Diag,
                 const int N, const void *A, const int lda, void *X,
                 const int incX);
void cblas_ctbsv(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo,
                 const  CBLAS_TRANSPOSE TransA, const  CBLAS_DIAG Diag,
                 const int N, const int K, const void *A, const int lda,
                 void *X, const int incX);
void cblas_ctpsv(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo,
                 const  CBLAS_TRANSPOSE TransA, const  CBLAS_DIAG Diag,
                 const int N, const void *Ap, void *X, const int incX);

void cblas_zgemv(const  CBLAS_ORDER order,
                 const  CBLAS_TRANSPOSE TransA, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 const void *X, const int incX, const void *beta,
                 void *Y, const int incY);
void cblas_zgbmv(const  CBLAS_ORDER order,
                 const  CBLAS_TRANSPOSE TransA, const int M, const int N,
                 const int KL, const int KU, const void *alpha,
                 const void *A, const int lda, const void *X,
                 const int incX, const void *beta, void *Y, const int incY);
void cblas_ztrmv(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo,
                 const  CBLAS_TRANSPOSE TransA, const  CBLAS_DIAG Diag,
                 const int N, const void *A, const int lda,
                 void *X, const int incX);
void cblas_ztbmv(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo,
                 const  CBLAS_TRANSPOSE TransA, const  CBLAS_DIAG Diag,
                 const int N, const int K, const void *A, const int lda,
                 void *X, const int incX);
void cblas_ztpmv(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo,
                 const  CBLAS_TRANSPOSE TransA, const  CBLAS_DIAG Diag,
                 const int N, const void *Ap, void *X, const int incX);
void cblas_ztrsv(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo,
                 const  CBLAS_TRANSPOSE TransA, const  CBLAS_DIAG Diag,
                 const int N, const void *A, const int lda, void *X,
                 const int incX);
void cblas_ztbsv(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo,
                 const  CBLAS_TRANSPOSE TransA, const  CBLAS_DIAG Diag,
                 const int N, const int K, const void *A, const int lda,
                 void *X, const int incX);
void cblas_ztpsv(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo,
                 const  CBLAS_TRANSPOSE TransA, const  CBLAS_DIAG Diag,
                 const int N, const void *Ap, void *X, const int incX);


/*
 * Routines with S and D prefixes only
 */
void cblas_ssymv(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo,
                 const int N, const float alpha, const float *A,
                 const int lda, const float *X, const int incX,
                 const float beta, float *Y, const int incY);
void cblas_ssbmv(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo,
                 const int N, const int K, const float alpha, const float *A,
                 const int lda, const float *X, const int incX,
                 const float beta, float *Y, const int incY);
void cblas_sspmv(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo,
                 const int N, const float alpha, const float *Ap,
                 const float *X, const int incX,
                 const float beta, float *Y, const int incY);
void cblas_sger(const  CBLAS_ORDER order, const int M, const int N,
                const float alpha, const float *X, const int incX,
                const float *Y, const int incY, float *A, const int lda);
void cblas_ssyr(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo,
                const int N, const float alpha, const float *X,
                const int incX, float *A, const int lda);
void cblas_sspr(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo,
                const int N, const float alpha, const float *X,
                const int incX, float *Ap);
void cblas_ssyr2(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo,
                const int N, const float alpha, const float *X,
                const int incX, const float *Y, const int incY, float *A,
                const int lda);
void cblas_sspr2(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo,
                const int N, const float alpha, const float *X,
                const int incX, const float *Y, const int incY, float *A);

void cblas_dsymv(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo,
                 const int N, const double alpha, const double *A,
                 const int lda, const double *X, const int incX,
                 const double beta, double *Y, const int incY);
void cblas_dsbmv(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo,
                 const int N, const int K, const double alpha, const double *A,
                 const int lda, const double *X, const int incX,
                 const double beta, double *Y, const int incY);
void cblas_dspmv(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo,
                 const int N, const double alpha, const double *Ap,
                 const double *X, const int incX,
                 const double beta, double *Y, const int incY);
void cblas_dger(const  CBLAS_ORDER order, const int M, const int N,
                const double alpha, const double *X, const int incX,
                const double *Y, const int incY, double *A, const int lda);
void cblas_dsyr(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo,
                const int N, const double alpha, const double *X,
                const int incX, double *A, const int lda);
void cblas_dspr(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo,
                const int N, const double alpha, const double *X,
                const int incX, double *Ap);
void cblas_dsyr2(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo,
                const int N, const double alpha, const double *X,
                const int incX, const double *Y, const int incY, double *A,
                const int lda);
void cblas_dspr2(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo,
                const int N, const double alpha, const double *X,
                const int incX, const double *Y, const int incY, double *A);


/*
 * Routines with C and Z prefixes only
 */
void cblas_chemv(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo,
                 const int N, const void *alpha, const void *A,
                 const int lda, const void *X, const int incX,
                 const void *beta, void *Y, const int incY);
void cblas_chbmv(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo,
                 const int N, const int K, const void *alpha, const void *A,
                 const int lda, const void *X, const int incX,
                 const void *beta, void *Y, const int incY);
void cblas_chpmv(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo,
                 const int N, const void *alpha, const void *Ap,
                 const void *X, const int incX,
                 const void *beta, void *Y, const int incY);
void cblas_cgeru(const  CBLAS_ORDER order, const int M, const int N,
                 const void *alpha, const void *X, const int incX,
                 const void *Y, const int incY, void *A, const int lda);
void cblas_cgerc(const  CBLAS_ORDER order, const int M, const int N,
                 const void *alpha, const void *X, const int incX,
                 const void *Y, const int incY, void *A, const int lda);
void cblas_cher(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo,
                const int N, const float alpha, const void *X, const int incX,
                void *A, const int lda);
void cblas_chpr(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo,
                const int N, const float alpha, const void *X,
                const int incX, void *A);
void cblas_cher2(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo, const int N,
                const void *alpha, const void *X, const int incX,
                const void *Y, const int incY, void *A, const int lda);
void cblas_chpr2(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo, const int N,
                const void *alpha, const void *X, const int incX,
                const void *Y, const int incY, void *Ap);

void cblas_zhemv(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo,
                 const int N, const void *alpha, const void *A,
                 const int lda, const void *X, const int incX,
                 const void *beta, void *Y, const int incY);
void cblas_zhbmv(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo,
                 const int N, const int K, const void *alpha, const void *A,
                 const int lda, const void *X, const int incX,
                 const void *beta, void *Y, const int incY);
void cblas_zhpmv(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo,
                 const int N, const void *alpha, const void *Ap,
                 const void *X, const int incX,
                 const void *beta, void *Y, const int incY);
void cblas_zgeru(const  CBLAS_ORDER order, const int M, const int N,
                 const void *alpha, const void *X, const int incX,
                 const void *Y, const int incY, void *A, const int lda);
void cblas_zgerc(const  CBLAS_ORDER order, const int M, const int N,
                 const void *alpha, const void *X, const int incX,
                 const void *Y, const int incY, void *A, const int lda);
void cblas_zher(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo,
                const int N, const double alpha, const void *X, const int incX,
                void *A, const int lda);
void cblas_zhpr(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo,
                const int N, const double alpha, const void *X,
                const int incX, void *A);
void cblas_zher2(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo, const int N,
                const void *alpha, const void *X, const int incX,
                const void *Y, const int incY, void *A, const int lda);
void cblas_zhpr2(const  CBLAS_ORDER order, const  CBLAS_UPLO Uplo, const int N,
                const void *alpha, const void *X, const int incX,
                const void *Y, const int incY, void *Ap);

/*
 * ===========================================================================
 * Prototypes for level 3 BLAS
 * ===========================================================================
 */

/*
 * Routines with standard 4 prefixes (S, D, C, Z)
 */
void cblas_sgemm(const  CBLAS_ORDER Order, const  CBLAS_TRANSPOSE TransA,
                 const  CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const float alpha, const float *A,
                 const int lda, const float *B, const int ldb,
                 const float beta, float *C, const int ldc);
void cblas_ssymm(const  CBLAS_ORDER Order, const  CBLAS_SIDE Side,
                 const  CBLAS_UPLO Uplo, const int M, const int N,
                 const float alpha, const float *A, const int lda,
                 const float *B, const int ldb, const float beta,
                 float *C, const int ldc);
void cblas_ssyrk(const  CBLAS_ORDER Order, const  CBLAS_UPLO Uplo,
                 const  CBLAS_TRANSPOSE Trans, const int N, const int K,
                 const float alpha, const float *A, const int lda,
                 const float beta, float *C, const int ldc);
void cblas_ssyr2k(const  CBLAS_ORDER Order, const  CBLAS_UPLO Uplo,
                  const  CBLAS_TRANSPOSE Trans, const int N, const int K,
                  const float alpha, const float *A, const int lda,
                  const float *B, const int ldb, const float beta,
                  float *C, const int ldc);
void cblas_strmm(const  CBLAS_ORDER Order, const  CBLAS_SIDE Side,
                 const  CBLAS_UPLO Uplo, const  CBLAS_TRANSPOSE TransA,
                 const  CBLAS_DIAG Diag, const int M, const int N,
                 const float alpha, const float *A, const int lda,
                 float *B, const int ldb);
void cblas_strsm(const  CBLAS_ORDER Order, const  CBLAS_SIDE Side,
                 const  CBLAS_UPLO Uplo, const  CBLAS_TRANSPOSE TransA,
                 const  CBLAS_DIAG Diag, const int M, const int N,
                 const float alpha, const float *A, const int lda,
                 float *B, const int ldb);

void cblas_dgemm(const  CBLAS_ORDER Order, const  CBLAS_TRANSPOSE TransA,
                 const  CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const double alpha, const double *A,
                 const int lda, const double *B, const int ldb,
                 const double beta, double *C, const int ldc);
void cblas_dsymm(const  CBLAS_ORDER Order, const  CBLAS_SIDE Side,
                 const  CBLAS_UPLO Uplo, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 const double *B, const int ldb, const double beta,
                 double *C, const int ldc);
void cblas_dsyrk(const  CBLAS_ORDER Order, const  CBLAS_UPLO Uplo,
                 const  CBLAS_TRANSPOSE Trans, const int N, const int K,
                 const double alpha, const double *A, const int lda,
                 const double beta, double *C, const int ldc);
void cblas_dsyr2k(const  CBLAS_ORDER Order, const  CBLAS_UPLO Uplo,
                  const  CBLAS_TRANSPOSE Trans, const int N, const int K,
                  const double alpha, const double *A, const int lda,
                  const double *B, const int ldb, const double beta,
                  double *C, const int ldc);
void cblas_dtrmm(const  CBLAS_ORDER Order, const  CBLAS_SIDE Side,
                 const  CBLAS_UPLO Uplo, const  CBLAS_TRANSPOSE TransA,
                 const  CBLAS_DIAG Diag, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 double *B, const int ldb);
void cblas_dtrsm(const  CBLAS_ORDER Order, const  CBLAS_SIDE Side,
                 const  CBLAS_UPLO Uplo, const  CBLAS_TRANSPOSE TransA,
                 const  CBLAS_DIAG Diag, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 double *B, const int ldb);

void cblas_cgemm(const  CBLAS_ORDER Order, const  CBLAS_TRANSPOSE TransA,
                 const  CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const void *alpha, const void *A,
                 const int lda, const void *B, const int ldb,
                 const void *beta, void *C, const int ldc);
void cblas_csymm(const  CBLAS_ORDER Order, const  CBLAS_SIDE Side,
                 const  CBLAS_UPLO Uplo, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 const void *B, const int ldb, const void *beta,
                 void *C, const int ldc);
void cblas_csyrk(const  CBLAS_ORDER Order, const  CBLAS_UPLO Uplo,
                 const  CBLAS_TRANSPOSE Trans, const int N, const int K,
                 const void *alpha, const void *A, const int lda,
                 const void *beta, void *C, const int ldc);
void cblas_csyr2k(const  CBLAS_ORDER Order, const  CBLAS_UPLO Uplo,
                  const  CBLAS_TRANSPOSE Trans, const int N, const int K,
                  const void *alpha, const void *A, const int lda,
                  const void *B, const int ldb, const void *beta,
                  void *C, const int ldc);
void cblas_ctrmm(const  CBLAS_ORDER Order, const  CBLAS_SIDE Side,
                 const  CBLAS_UPLO Uplo, const  CBLAS_TRANSPOSE TransA,
                 const  CBLAS_DIAG Diag, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 void *B, const int ldb);
void cblas_ctrsm(const  CBLAS_ORDER Order, const  CBLAS_SIDE Side,
                 const  CBLAS_UPLO Uplo, const  CBLAS_TRANSPOSE TransA,
                 const  CBLAS_DIAG Diag, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 void *B, const int ldb);

void cblas_zgemm(const  CBLAS_ORDER Order, const  CBLAS_TRANSPOSE TransA,
                 const  CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const void *alpha, const void *A,
                 const int lda, const void *B, const int ldb,
                 const void *beta, void *C, const int ldc);
void cblas_zsymm(const  CBLAS_ORDER Order, const  CBLAS_SIDE Side,
                 const  CBLAS_UPLO Uplo, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 const void *B, const int ldb, const void *beta,
                 void *C, const int ldc);
void cblas_zsyrk(const  CBLAS_ORDER Order, const  CBLAS_UPLO Uplo,
                 const  CBLAS_TRANSPOSE Trans, const int N, const int K,
                 const void *alpha, const void *A, const int lda,
                 const void *beta, void *C, const int ldc);
void cblas_zsyr2k(const  CBLAS_ORDER Order, const  CBLAS_UPLO Uplo,
                  const  CBLAS_TRANSPOSE Trans, const int N, const int K,
                  const void *alpha, const void *A, const int lda,
                  const void *B, const int ldb, const void *beta,
                  void *C, const int ldc);
void cblas_ztrmm(const  CBLAS_ORDER Order, const  CBLAS_SIDE Side,
                 const  CBLAS_UPLO Uplo, const  CBLAS_TRANSPOSE TransA,
                 const  CBLAS_DIAG Diag, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 void *B, const int ldb);
void cblas_ztrsm(const  CBLAS_ORDER Order, const  CBLAS_SIDE Side,
                 const  CBLAS_UPLO Uplo, const  CBLAS_TRANSPOSE TransA,
                 const  CBLAS_DIAG Diag, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 void *B, const int ldb);


/*
 * Routines with prefixes C and Z only
 */
void cblas_chemm(const  CBLAS_ORDER Order, const  CBLAS_SIDE Side,
                 const  CBLAS_UPLO Uplo, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 const void *B, const int ldb, const void *beta,
                 void *C, const int ldc);
void cblas_cherk(const  CBLAS_ORDER Order, const  CBLAS_UPLO Uplo,
                 const  CBLAS_TRANSPOSE Trans, const int N, const int K,
                 const float alpha, const void *A, const int lda,
                 const float beta, void *C, const int ldc);
void cblas_cher2k(const  CBLAS_ORDER Order, const  CBLAS_UPLO Uplo,
                  const  CBLAS_TRANSPOSE Trans, const int N, const int K,
                  const void *alpha, const void *A, const int lda,
                  const void *B, const int ldb, const float beta,
                  void *C, const int ldc);

void cblas_zhemm(const  CBLAS_ORDER Order, const  CBLAS_SIDE Side,
                 const  CBLAS_UPLO Uplo, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 const void *B, const int ldb, const void *beta,
                 void *C, const int ldc);
void cblas_zherk(const  CBLAS_ORDER Order, const  CBLAS_UPLO Uplo,
                 const  CBLAS_TRANSPOSE Trans, const int N, const int K,
                 const double alpha, const void *A, const int lda,
                 const double beta, void *C, const int ldc);
void cblas_zher2k(const  CBLAS_ORDER Order, const  CBLAS_UPLO Uplo,
                  const  CBLAS_TRANSPOSE Trans, const int N, const int K,
                  const void *alpha, const void *A, const int lda,
                  const void *B, const int ldb, const double beta,
                  void *C, const int ldc);


#ifdef __cplusplus
}
#endif    /* __cplusplus */

#endif    /* __MKL_CBLAS_H__ */
