
// Copyright 2001 California Institute of Technology.
// All rights reserved.

#include "cantera/base/ct_defs.h"

#ifndef LAPACK_H
#define LAPACK_H

#if defined(NEEDS_F77_TRANSLATION)

#if defined(F77EXTERNS_UPPERCASE_NOTRAILINGBAR)
#define dgetrs_ DGETRS
#define dgetrf_ DGETRF
#define dgetri_ DGETRI
#define dvode_ DVODE
#define ddassl_ DDASSL
#define simplx_ SIMPLX
#define splin2_ SPLIN2
#define splie2_ SPLIE2
#endif

#endif

extern "C"
{

    void dcopy_(const integer* n, const doublereal* dx, const integer* incx, doublereal* dy, const integer* incy);

    doublereal ddot_(integer* n, doublereal* dx, integer* incx, doublereal* dy, integer* incy);

    void daxpy_(integer* n, doublereal* a, doublereal* x, integer* incx,
                doublereal* y, integer* incy);

    void dscal_(integer* n, doublereal* da, doublereal* dx, integer* incx);

    integer idamax_(integer* n, doublereal* a, integer* incx);
}


#endif
