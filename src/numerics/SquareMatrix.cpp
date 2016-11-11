//! @file SquareMatrix.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/base/stringUtils.h"
#include "cantera/numerics/SquareMatrix.h"
#if CT_USE_LAPACK
    #include "cantera/numerics/ctlapack.h"
#endif

using namespace std;

namespace Cantera
{

SquareMatrix::SquareMatrix() :
    a1norm_(0.0),
    useQR_(0)
{
    warn_deprecated("class SquareMatrix",
        "Use class DenseMatrix instead. To be removed after Cantera 2.3.");
}

SquareMatrix::SquareMatrix(size_t n, doublereal v)  :
    DenseMatrix(n, n, v),
    a1norm_(0.0),
    useQR_(0)
{
    warn_deprecated("class SquareMatrix",
        "Use class DenseMatrix instead. To be removed after Cantera 2.3.");
}

SquareMatrix::SquareMatrix(const SquareMatrix& y) :
    DenseMatrix(y),
    GeneralMatrix(y),
    a1norm_(y.a1norm_),
    useQR_(y.useQR_)
{
}

SquareMatrix& SquareMatrix::operator=(const SquareMatrix& y)
{
    if (&y == this) {
        return *this;
    }
    DenseMatrix::operator=(y);
    GeneralMatrix::operator=(y);
    a1norm_ = y.a1norm_;
    useQR_ = y.useQR_;
    return *this;
}

int SquareMatrix::solve(doublereal* b, size_t nrhs, size_t ldb)
{
#if CT_USE_LAPACK
    if (useQR_) {
        return solveQR(b);
    }
    int info=0;

    // Check to see whether the matrix has been factored.
    if (!m_factored) {
        int retn = factor();
        if (retn) {
            return retn;
        }
    }
    if (ldb == 0) {
        ldb = nColumns();
    }

    // Solve the factored system
    ct_dgetrs(ctlapack::NoTranspose, static_cast<int>(nRows()),
              nrhs, &*begin(), static_cast<int>(nRows()),
              ipiv().data(), b, ldb, info);
    if (info != 0) {
        if (m_printLevel) {
            writelogf("SquareMatrix::solve(): DGETRS returned INFO = %d\n", info);
        }
        if (! m_useReturnErrorCode) {
            throw CanteraError("SquareMatrix::solve()", "DGETRS returned INFO = {}", info);
        }
    }
    return info;
#else
    throw CanteraError("SquareMatrix::solve",
                       "Not Implemented when LAPACK is not available");
#endif
}

void SquareMatrix::zero()
{
    m_data.assign(m_data.size(), 0.0);
}

void SquareMatrix::resize(size_t n, size_t m, doublereal v)
{
    DenseMatrix::resize(n, m, v);
}

void SquareMatrix::mult(const doublereal* b, doublereal* prod) const
{
    DenseMatrix::mult(b, prod);
}

void SquareMatrix::mult(const DenseMatrix& b, DenseMatrix& prod) const
{
    DenseMatrix::mult(b, prod);
}

void SquareMatrix::leftMult(const doublereal* const b, doublereal* const prod) const
{
    DenseMatrix::leftMult(b, prod);
}

int SquareMatrix::factor()
{
#if CT_USE_LAPACK
    if (useQR_) {
        return factorQR();
    }
    a1norm_ = ct_dlange('1', m_nrows, m_nrows, &*begin(), m_nrows, 0);
    integer n = static_cast<int>(nRows());
    int info=0;
    m_factored = 1;
    ct_dgetrf(n, n, &*begin(), static_cast<int>(nRows()), ipiv().data(), info);
    if (info != 0) {
        if (m_printLevel) {
            writelogf("SquareMatrix::factor(): DGETRS returned INFO = %d\n", info);
        }
        if (! m_useReturnErrorCode) {
            throw CanteraError("SquareMatrix::factor()", "DGETRS returned INFO = {}", info);
        }
    }
    return info;
#else
    throw CanteraError("SquareMatrix::factor",
                   "Not Implemented when LAPACK is not available");
#endif
}

void SquareMatrix::setFactorFlag()
{
    m_factored = 1;
}

int SquareMatrix::factorQR()
{
#if CT_USE_LAPACK
    if (tau.size() < m_nrows) {
        tau.resize(m_nrows, 0.0);
        work.resize(8 * m_nrows, 0.0);
    }
    a1norm_ = ct_dlange('1', m_nrows, m_nrows, &*begin(), m_nrows, work.data());
    int info = 0;
    m_factored = 2;
    size_t lwork = work.size();
    ct_dgeqrf(m_nrows, m_nrows, &*begin(), m_nrows, tau.data(), work.data(), lwork, info);
    if (info != 0) {
        if (m_printLevel) {
            writelogf("SquareMatrix::factorQR(): DGEQRF returned INFO = %d\n", info);
        }
        if (! m_useReturnErrorCode) {
            throw CanteraError("SquareMatrix::factorQR()", "DGEQRF returned INFO = {}", info);
        }
    }
    size_t lworkOpt = static_cast<size_t>(work[0]);
    if (lworkOpt > lwork) {
        work.resize(lworkOpt);
    }
    return info;
#else
    throw CanteraError("SquareMatrix::factorQR",
                   "Not Implemented when LAPACK is not available");
#endif
}

int SquareMatrix::solveQR(doublereal* b)
{
#if CT_USE_LAPACK
    int info=0;

    // Check to see whether the matrix has been factored.
    if (!m_factored) {
        int retn = factorQR();
        if (retn) {
            return retn;
        }
    }

    size_t lwork = work.size();
    if (lwork < m_nrows) {
        work.resize(8 * m_nrows, 0.0);
        lwork = 8 * m_nrows;
    }

    // Solve the factored system
    ct_dormqr(ctlapack::Left, ctlapack::Transpose, m_nrows, 1, m_nrows, &*begin(), m_nrows, tau.data(), b, m_nrows,
              work.data(), lwork, info);
    if (info != 0) {
        if (m_printLevel) {
            writelogf("SquareMatrix::solveQR(): DORMQR returned INFO = %d\n", info);
        }
        if (! m_useReturnErrorCode) {
            throw CanteraError("SquareMatrix::solveQR()", "DORMQR returned INFO = {}", info);
        }
    }
    size_t lworkOpt = static_cast<size_t>(work[0]);
    if (lworkOpt > lwork) {
        work.resize(lworkOpt);
    }

    char dd = 'N';
    ct_dtrtrs(ctlapack::UpperTriangular, ctlapack::NoTranspose, &dd, m_nrows, 1, &*begin(), m_nrows, b,
              m_nrows, info);
    if (info != 0) {
        if (m_printLevel) {
            writelogf("SquareMatrix::solveQR(): DTRTRS returned INFO = %d\n", info);
        }
        if (! m_useReturnErrorCode) {
            throw CanteraError("SquareMatrix::solveQR()", "DTRTRS returned INFO = {}", info);
        }
    }
    return info;
#else
    throw CanteraError("SquareMatrix::solveQR",
                   "Not Implemented when LAPACK is not available");
#endif
}

doublereal SquareMatrix::rcond(doublereal anorm)
{
#if CT_USE_LAPACK
    if (iwork_.size() < m_nrows) {
        iwork_.resize(m_nrows);
    }
    if (work.size() <4 * m_nrows) {
        work.resize(4 * m_nrows);
    }
    doublereal rcond = 0.0;
    if (m_factored != 1) {
        throw CanteraError("SquareMatrix::rcond()", "matrix isn't factored correctly");
    }

    int rinfo = 0;
    rcond = ct_dgecon('1', m_nrows, &*begin(), m_nrows, anorm, work.data(),
                      iwork_.data(), rinfo);
    if (rinfo != 0) {
        if (m_printLevel) {
            writelogf("SquareMatrix::rcond(): DGECON returned INFO = %d\n", rinfo);
        }
        if (! m_useReturnErrorCode) {
            throw CanteraError("SquareMatrix::rcond()", "DGECON returned INFO = {}", rinfo);
        }
    }
    return rcond;
#else
    throw CanteraError("SquareMatrix::rcond",
                   "Not Implemented when LAPACK is not available");
#endif
}

doublereal SquareMatrix::oneNorm() const
{
    return a1norm_;
}

doublereal SquareMatrix::rcondQR()
{
#if CT_USE_LAPACK
    if (iwork_.size() < m_nrows) {
        iwork_.resize(m_nrows);
    }
    if (work.size() <3 * m_nrows) {
        work.resize(3 * m_nrows);
    }
    doublereal rcond = 0.0;
    if (m_factored != 2) {
        throw CanteraError("SquareMatrix::rcondQR()", "matrix isn't factored correctly");
    }

    int rinfo = 0;
    rcond = ct_dtrcon(0, ctlapack::UpperTriangular, 0, m_nrows, &*begin(), m_nrows, work.data(),
                       iwork_.data(), rinfo);
    if (rinfo != 0) {
        if (m_printLevel) {
            writelogf("SquareMatrix::rcondQR(): DTRCON returned INFO = %d\n", rinfo);
        }
        if (! m_useReturnErrorCode) {
            throw CanteraError("SquareMatrix::rcondQR()", "DTRCON returned INFO = {}", rinfo);
        }
    }
    return rcond;
#else
    throw CanteraError("SquareMatrix::rcondQR",
                       "Not Implemented when LAPACK is not available");
#endif
}

void SquareMatrix::useFactorAlgorithm(int fAlgorithm)
{
    useQR_ = fAlgorithm;
}

int SquareMatrix::factorAlgorithm() const
{
    return (int) useQR_;
}

doublereal* SquareMatrix::ptrColumn(size_t j)
{
    return Array2D::ptrColumn(j);
}

size_t SquareMatrix::nRows() const
{
    return m_nrows;
}

size_t SquareMatrix::nRowsAndStruct(size_t* const iStruct) const
{
    warn_deprecated("SquareMatrix::nRowsAndStruct",
                    "To be removed after Cantera 2.3.");
    return m_nrows;
}

GeneralMatrix* SquareMatrix::duplMyselfAsGeneralMatrix() const
{
    return new SquareMatrix(*this);
}

vector_fp::iterator SquareMatrix::begin()
{
    return m_data.begin();
}

vector_fp::const_iterator SquareMatrix::begin() const
{
    return m_data.begin();
}

doublereal* const* SquareMatrix::colPts()
{
    return DenseMatrix::colPts();
}

size_t SquareMatrix::checkRows(doublereal& valueSmall) const
{
    valueSmall = 1.0E300;
    size_t iSmall = npos;
    for (size_t i = 0; i < m_nrows; i++) {
        double valueS = 0.0;
        for (size_t j = 0; j < m_nrows; j++) {
            valueS = std::max(fabs(value(i,j)), valueS);
        }
        if (valueS < valueSmall) {
            iSmall = i;
            valueSmall = valueS;
        }
    }
    return iSmall;
}

size_t SquareMatrix::checkColumns(doublereal& valueSmall) const
{
    valueSmall = 1.0E300;
    size_t jSmall = npos;
    for (size_t j = 0; j < m_nrows; j++) {
        double valueS = 0.0;
        for (size_t i = 0; i < m_nrows; i++) {
            valueS = std::max(fabs(value(i,j)), valueS);
        }
        if (valueS < valueSmall) {
            jSmall = j;
            valueSmall = valueS;
        }
    }
    return jSmall;
}

}
