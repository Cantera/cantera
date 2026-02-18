/**
 *  @file DenseMatrix.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/numerics/DenseMatrix.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/global.h"
#if CT_USE_LAPACK
    #include "cantera/numerics/ctlapack.h"
#else
    #include "cantera/numerics/eigen_dense.h"
#endif

namespace Cantera
{

DenseMatrix::DenseMatrix(size_t n, size_t m, double v) :
    Array2D(n, m, v)
{
    m_ipiv.resize(std::max(n, m));
}

DenseMatrix::DenseMatrix(const DenseMatrix& y) :
    Array2D(y)
{
    m_ipiv = y.ipiv();
}

DenseMatrix& DenseMatrix::operator=(const DenseMatrix& y)
{
    if (&y == this) {
        return *this;
    }
    Array2D::operator=(y);
    m_ipiv = y.ipiv();
    return *this;
}

void DenseMatrix::resize(size_t n, size_t m, double v)
{
    Array2D::resize(n,m,v);
    m_ipiv.resize(std::max(n,m));
}

void DenseMatrix::mult(span<const double> b, span<double> prod) const
{
    checkArraySize("DenseMatrix::mult", b.size(), nColumns());
    checkArraySize("DenseMatrix::mult", prod.size(), nRows());
#if CT_USE_LAPACK
    ct_dgemv(ctlapack::ColMajor, ctlapack::NoTranspose,
             static_cast<int>(nRows()),
             static_cast<int>(nColumns()), 1.0, col(0).data(),
             static_cast<int>(nRows()), b.data(), 1, 0.0, prod.data(), 1);
#else
    MappedMatrix mat(const_cast<double*>(m_data.data()), nRows(), nColumns());
    MappedVector bm(const_cast<double*>(b.data()), nColumns());
    MappedVector pm(prod.data(), nRows());
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
    for (size_t col=0; col < B.nColumns(); ++col) {
        // Loop over ncols multiplying A*column of B and storing in
        // corresponding prod column
        mult(B.col(col), prod.col(col));
    }
}

void DenseMatrix::leftMult(span<const double> b, span<double> prod) const
{
    checkArraySize("DenseMatrix::leftMult", b.size(), nRows());
    checkArraySize("DenseMatrix::leftMult", prod.size(), nColumns());
    for (size_t n = 0; n < nColumns(); n++) {
        double sum = 0.0;
        for (size_t i = 0; i < nRows(); i++) {
            sum += value(i,n)*b[i];
        }
        prod[n] = sum;
    }
}

vector<int>& DenseMatrix::ipiv()
{
    return m_ipiv;
}

void solve(DenseMatrix& A, span<double> b, size_t nrhs, size_t ldb)
{
    if (A.nColumns() != A.nRows()) {
        throw CanteraError("solve(DenseMatrix& A, span<double> b)",
            "Can only solve a square matrix");
    }

    int info = 0;
    if (ldb == 0) {
        ldb = A.nColumns();
    }
    checkArraySize("solve(DenseMatrix& A, span<double> b)", b.size(), ldb * nrhs);
    #if CT_USE_LAPACK
        ct_dgetrf(A.nRows(), A.nColumns(), A.data().data(),
                  A.nRows(), &A.ipiv()[0], info);
        if (info > 0) {
            throw CanteraError("solve(DenseMatrix& A, span<double> b)",
                                "DGETRF returned INFO = {}. U(i,i) is exactly zero. The factorization has"
                                " been completed, but the factor U is exactly singular, and division by zero will occur if "
                                "it is used to solve a system of equations.", info);
        } else if (info < 0) {
            throw CanteraError("solve(DenseMatrix& A, span<double> b)",
                               "DGETRF returned INFO = {}. The argument i has an illegal value", info);
        }

        ct_dgetrs(ctlapack::NoTranspose, A.nRows(), nrhs, A.data().data(),
                  A.nRows(), &A.ipiv()[0], b.data(), ldb, info);
        if (info != 0) {
            throw CanteraError("solve(DenseMatrix& A, span<double> b)",
                                "DGETRS returned INFO = {}", info);
        }
    #else
        MappedMatrix Am(&A(0,0), A.nRows(), A.nColumns());
        #ifdef NDEBUG
            auto lu = Am.partialPivLu();
        #else
            auto lu = Am.fullPivLu();
            if (lu.nonzeroPivots() < static_cast<long int>(A.nColumns())) {
                throw CanteraError("solve(DenseMatrix& A, span<double> b)",
                    "Matrix appears to be rank-deficient: non-zero pivots = {}; columns = {}",
                    lu.nonzeroPivots(), A.nColumns());
            }
        #endif
        for (size_t i = 0; i < nrhs; i++) {
            MappedVector bm(b.data() + ldb*i, A.nColumns());
            bm = lu.solve(bm);
        }
    #endif
}

void solve(DenseMatrix& A, DenseMatrix& b)
{
    solve(A, b.data(), b.nColumns(), b.nRows());
}

void multiply(const DenseMatrix& A, span<const double> b, span<double> prod)
{
    A.mult(b, prod);
}

void increment(const DenseMatrix& A, span<const double> b, span<double> prod)
{
    checkArraySize("increment", b.size(), A.nColumns());
    checkArraySize("increment", prod.size(), A.nRows());
    #if CT_USE_LAPACK
        ct_dgemv(ctlapack::ColMajor, ctlapack::NoTranspose,
                 static_cast<int>(A.nRows()), static_cast<int>(A.nColumns()), 1.0,
                 A.data().data(), static_cast<int>(A.nRows()), b.data(), 1, 1.0, prod.data(), 1);
    #else
        MappedMatrix Am(&const_cast<DenseMatrix&>(A)(0,0), A.nRows(), A.nColumns());
        MappedVector bm(const_cast<double*>(b.data()), A.nColumns());
        MappedVector pm(prod.data(), A.nRows());
        pm += Am * bm;
    #endif
}

void invert(DenseMatrix& A, size_t nn)
{
    int info=0;
    #if CT_USE_LAPACK
        integer n = static_cast<int>(nn != npos ? nn : A.nRows());
        ct_dgetrf(n, n, A.data().data(), static_cast<int>(A.nRows()),
                  &A.ipiv()[0], info);
        if (info != 0) {
            throw CanteraError("invert(DenseMatrix& A, size_t nn)", "DGETRF returned INFO = {}", info);
        }

        vector<double> work(n);
        integer lwork = static_cast<int>(work.size());
        ct_dgetri(n, A.data().data(), static_cast<int>(A.nRows()),
                  &A.ipiv()[0], &work[0], lwork, info);
        if (info != 0) {
            throw CanteraError("invert(DenseMatrix& A, size_t nn)", "DGETRI returned INFO={}", info);
        }
    #else
        MappedMatrix Am(&A(0,0), A.nRows(), A.nColumns());
        if (nn == npos) {
            Am = Am.inverse();
        } else {
            Am.topLeftCorner(nn, nn) = Am.topLeftCorner(nn, nn).inverse();
        }
    #endif
}

}
