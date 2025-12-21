#ifndef _SUNDIALS_LAPACK_H
#define _SUNDIALS_LAPACK_H

#include <sundials/sundials_types.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

#define dgbtrf_f77 dgbtrf_
#define dgbtrs_f77 dgbtrs_
#define dgetrf_f77 dgetrf_
#define dgetrs_f77 dgetrs_
#define dgeev_f77  dgeev_

extern void dgbtrf_f77(const sunindextype* m, const sunindextype* n,
                       const sunindextype* kl, const sunindextype* ku,
                       double* ab, sunindextype* ldab, sunindextype* ipiv,
                       sunindextype* info);

extern void dgbtrs_f77(const char* trans, const sunindextype* n,
                       const sunindextype* kl, const sunindextype* ku,
                       const sunindextype* nrhs, double* ab,
                       const sunindextype* ldab, sunindextype* ipiv, double* b,
                       const sunindextype* ldb, sunindextype* info);

extern void dgetrf_f77(const sunindextype* m, const sunindextype* n, double* a,
                       sunindextype* lda, sunindextype* ipiv, sunindextype* info);

extern void dgetrs_f77(const char* trans, const sunindextype* n,
                       const sunindextype* nrhs, double* a,
                       const sunindextype* lda, sunindextype* ipiv, double* b,
                       const sunindextype* ldb, sunindextype* info);

extern void dgeev_f77(const char* jobvl, const char* jobvr, const sunindextype* n,
                      double* a, const sunindextype* lda, double* wr, double* wi,
                      double* vl, const sunindextype* ldvl, double* vr,
                      const sunindextype* ldvr, double* work,
                      const sunindextype* lwork, sunindextype* info);

#ifdef __cplusplus
}
#endif

#endif
