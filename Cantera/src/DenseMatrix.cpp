/**
 *  @file DenseMatrix.cpp
 */

// Copyright 2001  California Institute of Technology

#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include <iostream>
#include <vector>
using namespace std;

#include "ct_defs.h"
#include "ctlapack.h"
#include "utilities.h"
#include "DenseMatrix.h"

namespace Cantera {

        /// assignment.
        DenseMatrix& DenseMatrix::operator=(const DenseMatrix& y) {
            if (&y == this) return *this;
            Array2D::operator=(y);
            m_ipiv = y.ipiv();
	    return *this;
        }

        void DenseMatrix::resize(int n, int m, doublereal v) {
            Array2D::resize(n,m,v);
            m_ipiv.resize( max(n,m) );
        }

        void DenseMatrix::mult(const double* b, double* prod) const {
            ct_dgemv(ctlapack::ColMajor, ctlapack::NoTranspose, nRows(), 
                nRows(), 1.0, begin(), nRows(), b, 1, 0.0, prod, 1);
        }

        void DenseMatrix::leftMult(const double* b, double* prod) const {
            int nc = nColumns();
            int nr = nRows();
            int n, i;
            double sum = 0.0;
            for (n = 0; n < nc; n++) {
                sum = 0.0;
                for (i = 0; i < nr; i++) {
                    sum += value(i,n)*b[i];
                }
                prod[n] = sum;
            }
        }

    /**
     * Solve Ax = b. Array b is overwritten on exit with x.
     */
    int solve(DenseMatrix& A, double* b) {
        int info=0;
        ct_dgetrf(A.nRows(), A.nColumns(), A.begin(), A.nRows(), 
            A.ipiv().begin(), info);
        if (info != 0) 
            throw CanteraError("DenseMatrix::solve",
                "DGETRF returned INFO = "+int2str(info));
        ct_dgetrs(ctlapack::NoTranspose, A.nRows(), 1, A.begin(), A.nRows(), 
            A.ipiv().begin(), b, A.nColumns(), info);
        if (info != 0) 
            throw CanteraError("DenseMatrix::solve",
                "DGETRS returned INFO = "+int2str(info));
        return 0;
    }

    /**  Solve Ax = b for multiple right-hand-side vectors. */
    int solve(DenseMatrix& A, DenseMatrix& b) {
        int info=0;
        ct_dgetrf(A.nRows(), A.nColumns(), A.begin(), A.nRows(), 
            A.ipiv().begin(), info);
        if (info != 0) 
            throw CanteraError("DenseMatrix::solve",
                "DGETRF returned INFO = "+int2str(info)); 
        ct_dgetrs(ctlapack::NoTranspose, A.nRows(), b.nColumns(), 
            A.begin(), A.nRows(), 
            A.ipiv().begin(), b.begin(), b.nRows(), info);
        if (info != 0)
            throw CanteraError("DenseMatrix::solve",
                "DGETRS returned INFO = "+int2str(info));
        return 0;
    }


    /** @todo fix lwork */
    int leastSquares(DenseMatrix& A, double* b) {
        int info = 0;
        int rank = 0;
        double rcond = -1.0;
        // fix this!
        int lwork = 6000; // 2*(3*min(m,n) + max(2*min(m,n), max(m,n)));
        vector_fp work(lwork);
        vector_fp s(min(A.nRows(),A.nColumns()));
        ct_dgelss(A.nRows(), A.nColumns(), 1, A.begin(), 
            A.nRows(), b, A.nColumns(), s.begin(),
            rcond, rank, work.begin(), work.size(), info);
        if (info != 0) 
            throw CanteraError("DenseMatrix::leaseSquares",
                "DGELSS returned INFO = "+int2str(info));
        return 0;
    }

    /**
     * Multiply \c A*b and return the result in \c prod. Uses BLAS
     * routine DGEMV.
     */
    void multiply(const DenseMatrix& A, const double* b, double* prod) {
        ct_dgemv(ctlapack::ColMajor, ctlapack::NoTranspose, 
            A.nRows(), A.nRows(), 1.0, 
            A.begin(), A.nRows(), b, 1, 0.0, prod, 1);
    }

    void increment(const DenseMatrix& A, 
        const double* b, double* prod) {
        ct_dgemv(ctlapack::ColMajor, ctlapack::NoTranspose, 
            A.nRows(), A.nRows(), 1.0, 
            A.begin(), A.nRows(), b, 1, 1.0, prod, 1);
    }


    /**
     * invert A. A is overwritten with A^-1.
     */
    int invert(DenseMatrix& A, int nn) {
        integer n = (nn > 0 ? nn : A.nRows());
        int info=0;
        ct_dgetrf(n, n, A.begin(), A.nRows(), A.ipiv().begin(), info);
        if (info != 0) 
            throw CanteraError("invert",
                "DGETRF returned INFO="+int2str(info));

        vector_fp work(n);
        integer lwork = n; 
        ct_dgetri(n, A.begin(), A.nRows(), A.ipiv().begin(), 
            work.begin(), lwork, info);
        if (info != 0) 
            throw CanteraError("invert",
                "DGETRI returned INFO="+int2str(info));
        return 0;
    }
}

