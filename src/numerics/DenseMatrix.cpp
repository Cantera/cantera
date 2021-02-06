/**
 *  @file DenseMatrix.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/numerics/DenseMatrix.h"
#include "cantera/base/stringUtils.h"
#if CT_USE_LAPACK
    #include "cantera/numerics/ctlapack.h"
#else
    #include "cantera/numerics/eigen_dense.h"
#endif

namespace Cantera
{

DenseMatrix::DenseMatrix() :
    m_useReturnErrorCode(0),
    m_printLevel(0)
{
}

DenseMatrix::DenseMatrix(size_t n, size_t m, doublereal v) :
    Array2D(n, m, v),
    m_useReturnErrorCode(0),
    m_printLevel(0)
{
    m_ipiv.resize(std::max(n, m));
    m_colPts.resize(m);
    if (!m_data.empty()) {
        for (size_t j = 0; j < m; j++) {
            m_colPts[j] = &m_data[m_nrows*j];
        }
    }
}

DenseMatrix::DenseMatrix(const DenseMatrix& y) :
    Array2D(y),
    m_useReturnErrorCode(0),
    m_printLevel(0)
{
    m_ipiv = y.ipiv();
    m_colPts.resize(m_ncols);
    if (!m_data.empty()) {
        for (size_t j = 0; j < m_ncols; j++) {
            m_colPts[j] = &m_data[m_nrows*j];
        }
    }
}

DenseMatrix& DenseMatrix::operator=(const DenseMatrix& y)
{
    if (&y == this) {
        return *this;
    }
    Array2D::operator=(y);
    m_ipiv = y.ipiv();
    m_colPts.resize(m_ncols);
    for (size_t j = 0; j < m_ncols; j++) {
        m_colPts[j] = &m_data[m_nrows*j];
    }
    m_useReturnErrorCode = y.m_useReturnErrorCode;
    m_printLevel = y.m_printLevel;
    return *this;
}

void DenseMatrix::resize(size_t n, size_t m, doublereal v)
{
    Array2D::resize(n,m,v);
    m_ipiv.resize(std::max(n,m));
    m_colPts.resize(m_ncols);
    if (!m_data.empty()) {
        for (size_t j = 0; j < m_ncols; j++) {
            m_colPts[j] = &m_data[m_nrows*j];
        }
    }
}

doublereal* const* DenseMatrix::colPts()
{
    return &m_colPts[0];
}

const doublereal* const* DenseMatrix::const_colPts() const
{
    return &m_colPts[0];
}

void DenseMatrix::mult(const double* b, double* prod) const
{
#if CT_USE_LAPACK
    ct_dgemv(ctlapack::ColMajor, ctlapack::NoTranspose,
             static_cast<int>(nRows()),
             static_cast<int>(nColumns()), 1.0, ptrColumn(0),
             static_cast<int>(nRows()), b, 1, 0.0, prod, 1);
#else
    MappedMatrix mat(const_cast<double*>(m_data.data()), nRows(), nColumns());
    MappedVector bm(const_cast<double*>(b), nColumns());
    MappedVector pm(prod, nRows());
    pm = mat * bm;
#endif
}

void DenseMatrix::mult(const DenseMatrix& B, DenseMatrix& prod) const
{
    if (nColumns() != B.nRows()) {
        throw CanteraError("DenseMatrix::mult",
            "Inner matrix dimensions do not agree: {} != {}", nColumns(), B.nRows());
    }
    if (nRows() != prod.nRows() || B.nColumns() != prod.nColumns()) {
        throw CanteraError("DenseMatrix::mult",
            "Output matrix has wrong dimensions: {}x{} != {}x{}",
            prod.nRows(), prod.nColumns(), nRows(), B.nColumns());
    }
    const doublereal* const* bcols = B.const_colPts();
    doublereal* const* prodcols = prod.colPts();
    for (size_t col=0; col < B.nColumns(); ++col) {
        // Loop over ncols multiplying A*column of B and storing in
        // corresponding prod column
        mult(bcols[col], prodcols[col]);
    }
}

void DenseMatrix::leftMult(const double* const b, double* const prod) const
{
    for (size_t n = 0; n < nColumns(); n++) {
        double sum = 0.0;
        for (size_t i = 0; i < nRows(); i++) {
            sum += value(i,n)*b[i];
        }
        prod[n] = sum;
    }
}

vector_int& DenseMatrix::ipiv()
{
    return m_ipiv;
}

int solve(DenseMatrix& A, double* b, size_t nrhs, size_t ldb)
{
    if (A.nColumns() != A.nRows()) {
        if (A.m_printLevel) {
            writelogf("solve(DenseMatrix& A, double* b): Can only solve a square matrix\n");
        }
        throw CanteraError("solve(DenseMatrix& A, double* b)", "Can only solve a square matrix");
    }

    int info = 0;
    if (ldb == 0) {
        ldb = A.nColumns();
    }
    #if CT_USE_LAPACK
        ct_dgetrf(A.nRows(), A.nColumns(), A.ptrColumn(0),
                  A.nRows(), &A.ipiv()[0], info);
        if (info > 0) {
            if (A.m_printLevel) {
                writelogf("solve(DenseMatrix& A, double* b): DGETRF returned INFO = %d   U(i,i) is exactly zero. The factorization has"
                          " been completed, but the factor U is exactly singular, and division by zero will occur if "
                          "it is used to solve a system of equations.\n", info);
            }
            if (!A.m_useReturnErrorCode) {
                throw CanteraError("solve(DenseMatrix& A, double* b)",
                                    "DGETRF returned INFO = {}. U(i,i) is exactly zero. The factorization has"
                                    " been completed, but the factor U is exactly singular, and division by zero will occur if "
                                    "it is used to solve a system of equations.", info);
            }
            return info;
        } else if (info < 0) {
            if (A.m_printLevel) {
                writelogf("solve(DenseMatrix& A, double* b): DGETRF returned INFO = %d. The argument i has an illegal value\n", info);
            }

            throw CanteraError("solve(DenseMatrix& A, double* b)",
                               "DGETRF returned INFO = {}. The argument i has an illegal value", info);
        }

        ct_dgetrs(ctlapack::NoTranspose, A.nRows(), nrhs, A.ptrColumn(0),
                  A.nRows(), &A.ipiv()[0], b, ldb, info);
        if (info != 0) {
            if (A.m_printLevel) {
                writelog("solve(DenseMatrix& A, double* b): DGETRS returned INFO = {}\n", info);
            }
            if (info < 0 || !A.m_useReturnErrorCode) {
                throw CanteraError("solve(DenseMatrix& A, double* b)",
                                   "DGETRS returned INFO = {}", info);
            }
        }
    #else
        MappedMatrix Am(&A(0,0), A.nRows(), A.nColumns());
        #ifdef NDEBUG
            auto lu = Am.partialPivLu();
        #else
            auto lu = Am.fullPivLu();
            if (lu.nonzeroPivots() < static_cast<long int>(A.nColumns())) {
                throw CanteraError("solve(DenseMatrix& A, double* b)",
                    "Matrix appears to be rank-deficient: non-zero pivots = {}; columns = {}",
                    lu.nonzeroPivots(), A.nColumns());
            }
        #endif
        for (size_t i = 0; i < nrhs; i++) {
            MappedVector bm(b + ldb*i, A.nColumns());
            bm = lu.solve(bm);
        }
    #endif
    return info;
}

int solve(DenseMatrix& A, DenseMatrix& b)
{
    return solve(A, b.ptrColumn(0), b.nColumns(), b.nRows());
}

void multiply(const DenseMatrix& A, const double* const b, double* const prod)
{
    A.mult(b, prod);
}

void increment(const DenseMatrix& A, const double* b, double* prod)
{
    #if CT_USE_LAPACK
        ct_dgemv(ctlapack::ColMajor, ctlapack::NoTranspose,
                 static_cast<int>(A.nRows()), static_cast<int>(A.nColumns()), 1.0,
                 A.ptrColumn(0), static_cast<int>(A.nRows()), b, 1, 1.0, prod, 1);
    #else
        MappedMatrix Am(&const_cast<DenseMatrix&>(A)(0,0), A.nRows(), A.nColumns());
        MappedVector bm(const_cast<double*>(b), A.nColumns());
        MappedVector pm(prod, A.nRows());
        pm += Am * bm;
    #endif
}

int invert(DenseMatrix& A, size_t nn)
{
    int info=0;
    #if CT_USE_LAPACK
        integer n = static_cast<int>(nn != npos ? nn : A.nRows());
        ct_dgetrf(n, n, A.ptrColumn(0), static_cast<int>(A.nRows()),
                  &A.ipiv()[0], info);
        if (info != 0) {
            if (A.m_printLevel) {
                writelogf("invert(DenseMatrix& A, size_t nn): DGETRS returned INFO = %d\n", info);
            }
            if (! A.m_useReturnErrorCode) {
                throw CanteraError("invert(DenseMatrix& A, size_t nn)", "DGETRS returned INFO = {}", info);
            }
            return info;
        }

        vector_fp work(n);
        integer lwork = static_cast<int>(work.size());
        ct_dgetri(n, A.ptrColumn(0), static_cast<int>(A.nRows()),
                  &A.ipiv()[0], &work[0], lwork, info);
        if (info != 0) {
            if (A.m_printLevel) {
                writelogf("invert(DenseMatrix& A, size_t nn): DGETRS returned INFO = %d\n", info);
            }
            if (! A.m_useReturnErrorCode) {
                throw CanteraError("invert(DenseMatrix& A, size_t nn)", "DGETRI returned INFO={}", info);
            }
        }
    #else
        MappedMatrix Am(&A(0,0), A.nRows(), A.nColumns());
        if (nn == npos) {
            Am = Am.inverse();
        } else {
            Am.topLeftCorner(nn, nn) = Am.topLeftCorner(nn, nn).inverse();
        }
    #endif
    return info;
}

}
