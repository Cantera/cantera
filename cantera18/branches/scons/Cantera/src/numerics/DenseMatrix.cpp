/**
 *  @file DenseMatrix.cpp
 */

// Copyright 2001  California Institute of Technology

#include "ct_defs.h"
#include "ctlapack.h"
#include "utilities.h"
#include "DenseMatrix.h"
#include "stringUtils.h"

namespace Cantera {

        /// assignment.
        DenseMatrix& DenseMatrix::operator=(const DenseMatrix& y) {
            if (&y == this) return *this;
            Array2D::operator=(y);
            m_ipiv = y.ipiv();
	    return *this;
        }

        void DenseMatrix::resize(size_t n, size_t m, doublereal v) {
            Array2D::resize(n,m,v);
            m_ipiv.resize( max(n,m) );
        }

        void DenseMatrix::mult(const double* b, double* prod) const {
            ct_dgemv(ctlapack::ColMajor, ctlapack::NoTranspose, 
                int(nRows()), int(nRows()), 1.0, ptrColumn(0),
                int(nRows()), b, 1, 0.0, prod, 1);
        }

        void DenseMatrix::leftMult(const double* b, double* prod) const {
            size_t nc = nColumns();
            size_t nr = nRows();
            size_t n, i;
            double sum = 0.0;
            for (n = 0; n < nc; n++) {
                sum = 0.0;
                for (i = 0; i < nr; i++) {
                    sum += value(i,n)*b[i];
                }
                prod[n] = sum;
            }
        }

    int solve(DenseMatrix& A, double* b) {
        int info=0;
        ct_dgetrf(int(A.nRows()), int(A.nColumns()), A.ptrColumn(0),
		    int(A.nRows()), &A.ipiv()[0], info);
        if (info != 0) 
            throw CanteraError("DenseMatrix::solve",
                "DGETRF returned INFO = "+int2str(info));
        ct_dgetrs(ctlapack::NoTranspose, int(A.nRows()), 1, A.ptrColumn(0),
            int(A.nRows()), &A.ipiv()[0], b, int(A.nColumns()), info);
        if (info != 0) 
            throw CanteraError("DenseMatrix::solve",
                "DGETRS returned INFO = "+int2str(info));
        return 0;
    }

    int solve(DenseMatrix& A, DenseMatrix& b) {
        int info=0;
        ct_dgetrf(int(A.nRows()), int(A.nColumns()), A.ptrColumn(0),
		    int(A.nRows()), &A.ipiv()[0], info);
        if (info != 0)
            throw CanteraError("DenseMatrix::solve",
                "DGETRF returned INFO = "+int2str(info)); 
        ct_dgetrs(ctlapack::NoTranspose, int(A.nRows()), int(b.nColumns()),
            A.ptrColumn(0), int(A.nRows()), &A.ipiv()[0], b.ptrColumn(0),
            int(b.nRows()), info);
        if (info != 0)
            throw CanteraError("DenseMatrix::solve",
                "DGETRS returned INFO = "+int2str(info));
        return 0;
    }

#ifdef INCL_LEAST_SQUARES
    /** @todo fix lwork */
    int leastSquares(DenseMatrix& A, double* b) {
        int info = 0;
        int rank = 0;
        double rcond = -1.0;
        // fix this!
        int lwork = 6000; // 2*(3*min(m,n) + max(2*min(m,n), max(m,n)));
        vector_fp work(lwork);
        vector_fp s(min(A.nRows(), A.nColumns()));
        ct_dgelss(int(A.nRows()), int(A.nColumns()), 1, A.ptrColumn(0),
		    int(A.nRows()), b, int(A.nColumns()), &s[0], rcond, rank,
		    &work[0], work.size(), info);
        if (info != 0) 
            throw CanteraError("DenseMatrix::leaseSquares",
                "DGELSS returned INFO = "+int2str(info));
        return 0; 
    }
#endif

    void multiply(const DenseMatrix& A, const double* b, double* prod) {
        ct_dgemv(ctlapack::ColMajor, ctlapack::NoTranspose, int(A.nRows()),
            int(A.nColumns()), 1.0, A.ptrColumn(0), int(A.nRows()), b, 1,
            0.0, prod, 1);
    }

    void increment(const DenseMatrix& A, const double* b, double* prod) {
        ct_dgemv(ctlapack::ColMajor, ctlapack::NoTranspose, int(A.nRows()),
		    int(A.nRows()), 1.0, A.ptrColumn(0), int(A.nRows()), b, 1, 1.0,
		    prod, 1);
    }

    int invert(DenseMatrix& A, size_t nn) {
        int n = int(nn > 0 ? nn : A.nRows());
        int info=0;
        ct_dgetrf(n, n, A.ptrColumn(0), int(A.nRows()), &A.ipiv()[0], info);
        if (info != 0)
            throw CanteraError("invert",
                "DGETRF returned INFO="+int2str(info));

        vector_fp work(n);
        ct_dgetri(n, A.ptrColumn(0), int(A.nRows()), &A.ipiv()[0], &work[0],
		    int(work.size()), info);
        if (info != 0)
            throw CanteraError("invert",
                "DGETRI returned INFO="+int2str(info));
        return 0;
    }
}
