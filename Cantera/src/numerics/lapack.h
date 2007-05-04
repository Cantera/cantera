
// Copyright 2001 California Institute of Technology.  
// All rights reserved.

#include "ct_defs.h"

#ifndef LAPACK_H
#define LAPACK_H

#if defined(NEEDS_F77_TRANSLATION)

#if defined(F77EXTERNS_UPPERCASE_NOTRAILINGBAR)
#define dgelss_ DGELSS
#define dgetrs_ DGETRS
#define dgetrf_ DGETRF
#define dgetri_ DGETRI
#define dvode_ DVODE
#define ddassl_ DDASSL
#define simplx_ SIMPLX
#define splin2_ SPLIN2
#define splie2_ SPLIE2
#define dgelss_ DGELSS
#endif

#endif

extern "C" 
{

  //  /* Subroutine */ int dgelss_(integer *m, integer *n, integer *nrhs, 
  //	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
  //	s, doublereal *rcond, integer *rank, doublereal *work, integer *lwork,
  //	 integer *info);

  ///* Subroutine */ int dgetrs_(char *trans, integer *n, integer *nrhs, 
//	doublereal *a, integer *lda, integer *ipiv, doublereal *b, integer *
//	ldb, integer *info, int len);

///* Subroutine */ int dgetrf_(integer *m, integer *n, doublereal *a, integer *
//	lda, integer *ipiv, integer *info);

///* Subroutine */ int dgetri_(integer *n, doublereal *a, integer *
//	lda, integer *ipiv, doublereal *work, integer *lwork, integer *info);

//  /* Subroutine */ int dgemv_(char *trans, integer *m, integer *n, doublereal *
//  	alpha, doublereal *a, integer *lda, doublereal *x, integer *incx, 
//  	doublereal *beta, doublereal *y, integer *incy, int len);

void dcopy_(const integer *n, const doublereal *dx, const integer *incx, doublereal *dy, const integer *incy); 
  //doublereal ddot_(integer *n, doublereal *dx, integer *incx, doublereal *dy, integer *incy); 

doublereal ddot_(integer *n, doublereal *dx, integer *incx, doublereal *dy, integer *incy); 

void daxpy_(integer* n, doublereal* a, doublereal* x, integer* incx, 
					 doublereal* y, integer* incy);
void dscal_(integer *n, doublereal *da, doublereal *dx, integer *incx);
integer idamax_(integer* n, doublereal* a, integer* incx);
  
}


#endif
