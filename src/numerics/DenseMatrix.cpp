/**
 *  @file DenseMatrix.cpp
 */

// Copyright 2001  California Institute of Technology

#include "cantera/numerics/ctlapack.h"
#include "cantera/numerics/DenseMatrix.h"
#include "cantera/base/stringUtils.h"

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
            m_colPts[j] = &(m_data[m_nrows*j]);
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
            m_colPts[j] = &(m_data[m_nrows*j]);
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
        m_colPts[j] = &(m_data[m_nrows*j]);
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
            m_colPts[j] = &(m_data[m_nrows*j]);
        }
    }
}

doublereal*   const* DenseMatrix::colPts()
{
    return &(m_colPts[0]);
}

const doublereal*   const* DenseMatrix::const_colPts() const
{
    return &(m_colPts[0]);
}

void DenseMatrix::mult(const double* b, double* prod) const
{
    ct_dgemv(ctlapack::ColMajor, ctlapack::NoTranspose,
             static_cast<int>(nRows()),
             static_cast<int>(nRows()), 1.0, ptrColumn(0),
             static_cast<int>(nRows()), b, 1, 0.0, prod, 1);
}

void DenseMatrix::mult(const DenseMatrix& B, DenseMatrix& prod) const
{
    if (m_ncols != B.nColumns() || m_nrows != B.nRows() || m_ncols != m_nrows || m_ncols != prod.nColumns() || m_nrows != prod.nColumns()) {
        throw CanteraError("mult(const DenseMatrix &B, DenseMatrix &prod)",
                           "Cannot multiply matrices that are not square and/or not the same size.");
    }
    const doublereal* const* bcols = B.const_colPts();
    doublereal* const* prodcols = prod.colPts();
    for (size_t col=0; col < m_ncols; ++col) {
        // Loop over ncols multiplying A*column of B and storing in corresponding prod column
        mult(bcols[col], prodcols[col]);
    }
}

void DenseMatrix::leftMult(const double* const b, double* const prod) const
{
    size_t nc = nColumns();
    size_t nr = nRows();
    double sum = 0.0;
    for (size_t n = 0; n < nc; n++) {
        sum = 0.0;
        for (size_t i = 0; i < nr; i++) {
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
    int info = 0;
    if (A.nColumns() != A.nRows()) {
        if (A.m_printLevel) {
            writelogf("solve(DenseMatrix& A, double* b): Can only solve a square matrix\n");
        }
        throw CELapackError("solve(DenseMatrix& A, double* b)", "Can only solve a square matrix");
    }
    ct_dgetrf(A.nRows(), A.nColumns(), A.ptrColumn(0),
              A.nRows(), &A.ipiv()[0], info);
    if (info > 0) {
        if (A.m_printLevel) {
            writelogf("solve(DenseMatrix& A, double* b): DGETRF returned INFO = %d   U(i,i) is exactly zero. The factorization has"
                      " been completed, but the factor U is exactly singular, and division by zero will occur if "
                      "it is used to solve a system of equations.\n", info);
        }
        if (!A.m_useReturnErrorCode) {
            throw CELapackError("solve(DenseMatrix& A, double* b)",
                                "DGETRF returned INFO = "+int2str(info) + ".   U(i,i) is exactly zero. The factorization has"
                                " been completed, but the factor U is exactly singular, and division by zero will occur if "
                                "it is used to solve a system of equations.");
        }
        return info;
    } else if (info < 0) {
        if (A.m_printLevel) {
            writelogf("solve(DenseMatrix& A, double* b): DGETRF returned INFO = %d. The argument i has an illegal value\n", info);
        }

        throw CELapackError("solve(DenseMatrix& A, double* b)",
                            "DGETRF returned INFO = "+int2str(info) + ". The argument i has an illegal value");
    }

    if (ldb == 0) {
        ldb = A.nColumns();
    }
    ct_dgetrs(ctlapack::NoTranspose, A.nRows(), nrhs, A.ptrColumn(0),
              A.nRows(), &A.ipiv()[0], b, ldb, info);
    if (info != 0) {
        if (A.m_printLevel) {
            writelogf("solve(DenseMatrix& A, double* b): DGETRS returned INFO = %d\n", info);
        }
        if (info < 0 || !A.m_useReturnErrorCode) {
            throw CELapackError("solve(DenseMatrix& A, double* b)", "DGETRS returned INFO = "+int2str(info));
        }
    }
    return info;
}

int solve(DenseMatrix& A, DenseMatrix& b)
{
    return solve(A, b.ptrColumn(0), b.nColumns(), b.nRows());
}

void multiply(const DenseMatrix& A, const double* const b, double* const prod)
{
    ct_dgemv(ctlapack::ColMajor, ctlapack::NoTranspose,
             static_cast<int>(A.nRows()), static_cast<int>(A.nColumns()), 1.0,
             A.ptrColumn(0), static_cast<int>(A.nRows()), b, 1, 0.0, prod, 1);
}

void increment(const DenseMatrix& A, const double* b, double* prod)
{
    ct_dgemv(ctlapack::ColMajor, ctlapack::NoTranspose,
             static_cast<int>(A.nRows()), static_cast<int>(A.nRows()), 1.0,
             A.ptrColumn(0), static_cast<int>(A.nRows()), b, 1, 1.0, prod, 1);
}

int invert(DenseMatrix& A, size_t nn)
{
    integer n = static_cast<int>(nn != npos ? nn : A.nRows());
    int info=0;
    ct_dgetrf(n, n, A.ptrColumn(0), static_cast<int>(A.nRows()),
              &A.ipiv()[0], info);
    if (info != 0) {
        if (A.m_printLevel) {
            writelogf("invert(DenseMatrix& A, int nn): DGETRS returned INFO = %d\n", info);
        }
        if (! A.m_useReturnErrorCode) {
            throw CELapackError("invert(DenseMatrix& A, int nn)", "DGETRS returned INFO = "+int2str(info));
        }
        return info;
    }

    vector_fp work(n);
    integer lwork = static_cast<int>(work.size());
    ct_dgetri(n, A.ptrColumn(0), static_cast<int>(A.nRows()),
              &A.ipiv()[0],  &work[0], lwork, info);
    if (info != 0) {
        if (A.m_printLevel) {
            writelogf("invert(DenseMatrix& A, int nn): DGETRS returned INFO = %d\n", info);
        }
        if (! A.m_useReturnErrorCode) {
            throw CELapackError("invert(DenseMatrix& A, int nn)", "DGETRI returned INFO="+int2str(info));
        }
    }
    return info;
}

}
