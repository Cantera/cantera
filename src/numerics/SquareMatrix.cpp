/**
 *  @file SquareMatrix.cpp
 */

/*
 * Copyright 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "cantera/base/stringUtils.h"
#include "cantera/numerics/ctlapack.h"
#include "cantera/numerics/SquareMatrix.h"

using namespace std;

namespace Cantera
{

SquareMatrix::SquareMatrix() :
    GeneralMatrix(0),
    a1norm_(0.0),
    useQR_(0)
{
}

SquareMatrix::SquareMatrix(size_t n, doublereal v)  :
    DenseMatrix(n, n, v),
    GeneralMatrix(0),
    a1norm_(0.0),
    useQR_(0)

{
}

SquareMatrix::SquareMatrix(const SquareMatrix& y) :
    DenseMatrix(y),
    GeneralMatrix(0),
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
    if (useQR_) {
        return solveQR(b);
    }
    int info=0;
    /*
     * Check to see whether the matrix has been factored.
     */
    if (!m_factored) {
        int retn = factor();
        if (retn) {
            return retn;
        }
    }
    if (ldb == 0) {
        ldb = nColumns();
    }
    /*
     * Solve the factored system
     */
    ct_dgetrs(ctlapack::NoTranspose, static_cast<int>(nRows()),
              nrhs, &(*(begin())), static_cast<int>(nRows()),
              DATA_PTR(ipiv()), b, ldb, info);
    if (info != 0) {
        if (m_printLevel) {
            writelogf("SquareMatrix::solve(): DGETRS returned INFO = %d\n", info);
        }
        if (! m_useReturnErrorCode) {
            throw CELapackError("SquareMatrix::solve()", "DGETRS returned INFO = " + int2str(info));
        }
    }
    return info;
}

void SquareMatrix::zero()
{
    size_t n = nRows();
    if (n > 0) {
        size_t nn = n * n;
        double* sm = &m_data[0];
        /*
         * Using memset is the fastest way to zero a contiguous
         * section of memory.
         */
        (void) memset((void*) sm, 0, nn * sizeof(double));
    }
}

void SquareMatrix::resize(size_t n, size_t m, doublereal v)
{
    DenseMatrix::resize(n, m, v);
}

void  SquareMatrix::mult(const doublereal* b, doublereal* prod) const
{
    DenseMatrix::mult(b, prod);
}

void SquareMatrix::mult(const DenseMatrix& b, DenseMatrix& prod) const
{
    DenseMatrix::mult(b, prod);
}

void  SquareMatrix::leftMult(const doublereal* const b, doublereal* const prod) const
{
    DenseMatrix::leftMult(b, prod);
}

int SquareMatrix::factor()
{
    if (useQR_) {
        return factorQR();
    }
    a1norm_ = ct_dlange('1', m_nrows, m_nrows, &(*(begin())), m_nrows, 0);
    integer n = static_cast<int>(nRows());
    int info=0;
    m_factored = 1;
    ct_dgetrf(n, n, &(*(begin())), static_cast<int>(nRows()), DATA_PTR(ipiv()), info);
    if (info != 0) {
        if (m_printLevel) {
            writelogf("SquareMatrix::factor(): DGETRS returned INFO = %d\n", info);
        }
        if (! m_useReturnErrorCode) {
            throw CELapackError("SquareMatrix::factor()", "DGETRS returned INFO = "+int2str(info));
        }
    }
    return info;
}

void SquareMatrix::setFactorFlag()
{
    m_factored = 1;
}

int SquareMatrix::factorQR()
{
    if (tau.size() < m_nrows)  {
        tau.resize(m_nrows, 0.0);
        work.resize(8 * m_nrows, 0.0);
    }
    a1norm_ = ct_dlange('1', m_nrows, m_nrows, &(*(begin())), m_nrows, DATA_PTR(work));
    int info = 0;
    m_factored = 2;
    size_t lwork = work.size();
    ct_dgeqrf(m_nrows, m_nrows, &(*(begin())), m_nrows, DATA_PTR(tau), DATA_PTR(work), lwork, info);
    if (info != 0) {
        if (m_printLevel) {
            writelogf("SquareMatrix::factorQR(): DGEQRF returned INFO = %d\n", info);
        }
        if (! m_useReturnErrorCode) {
            throw CELapackError("SquareMatrix::factorQR()", "DGEQRF returned INFO = " + int2str(info));
        }
    }
    size_t lworkOpt = static_cast<size_t>(work[0]);
    if (lworkOpt > lwork) {
        work.resize(lworkOpt);
    }


    return info;
}

int SquareMatrix::solveQR(doublereal* b)
{
    int info=0;
    /*
     * Check to see whether the matrix has been factored.
     */
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

    /*
     * Solve the factored system
     */
    ct_dormqr(ctlapack::Left, ctlapack::Transpose, m_nrows, 1, m_nrows, &(*(begin())), m_nrows, DATA_PTR(tau), b, m_nrows,
              DATA_PTR(work), lwork, info);
    if (info != 0) {
        if (m_printLevel) {
            writelogf("SquareMatrix::solveQR(): DORMQR returned INFO = %d\n", info);
        }
        if (! m_useReturnErrorCode) {
            throw CELapackError("SquareMatrix::solveQR()", "DORMQR returned INFO = " + int2str(info));
        }
    }
    size_t lworkOpt = static_cast<size_t>(work[0]);
    if (lworkOpt > lwork) {
        work.resize(lworkOpt);
    }

    char dd = 'N';

    ct_dtrtrs(ctlapack::UpperTriangular, ctlapack::NoTranspose, &dd, m_nrows, 1,  &(*(begin())), m_nrows, b,
              m_nrows, info);
    if (info != 0) {
        if (m_printLevel) {
            writelogf("SquareMatrix::solveQR(): DTRTRS returned INFO = %d\n", info);
        }
        if (! m_useReturnErrorCode) {
            throw CELapackError("SquareMatrix::solveQR()", "DTRTRS returned INFO = " + int2str(info));
        }
    }

    return info;
}

doublereal SquareMatrix::rcond(doublereal anorm)
{

    if (iwork_.size() < m_nrows) {
        iwork_.resize(m_nrows);
    }
    if (work.size() <4 * m_nrows) {
        work.resize(4 * m_nrows);
    }
    doublereal rcond = 0.0;
    if (m_factored != 1) {
        throw CELapackError("SquareMatrix::rcond()", "matrix isn't factored correctly");
    }

    int rinfo = 0;
    rcond = ct_dgecon('1', m_nrows, &(*(begin())), m_nrows, anorm, DATA_PTR(work),
                      DATA_PTR(iwork_), rinfo);
    if (rinfo != 0) {
        if (m_printLevel) {
            writelogf("SquareMatrix::rcond(): DGECON returned INFO = %d\n", rinfo);
        }
        if (! m_useReturnErrorCode) {
            throw CELapackError("SquareMatrix::rcond()", "DGECON returned INFO = " + int2str(rinfo));
        }
    }
    return rcond;
}

doublereal SquareMatrix::oneNorm() const
{
    return a1norm_;
}

doublereal SquareMatrix::rcondQR()
{

    if (iwork_.size() < m_nrows) {
        iwork_.resize(m_nrows);
    }
    if (work.size() <3 * m_nrows) {
        work.resize(3 * m_nrows);
    }
    doublereal rcond = 0.0;
    if (m_factored != 2) {
        throw CELapackError("SquareMatrix::rcondQR()", "matrix isn't factored correctly");
    }

    int rinfo = 0;
    rcond =  ct_dtrcon(0, ctlapack::UpperTriangular, 0, m_nrows, &(*(begin())), m_nrows, DATA_PTR(work),
                       DATA_PTR(iwork_), rinfo);
    if (rinfo != 0) {
        if (m_printLevel) {
            writelogf("SquareMatrix::rcondQR(): DTRCON returned INFO = %d\n", rinfo);
        }
        if (! m_useReturnErrorCode) {
            throw CELapackError("SquareMatrix::rcondQR()", "DTRCON returned INFO = " + int2str(rinfo));
        }
    }
    return rcond;
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

void  SquareMatrix::copyData(const GeneralMatrix& y)
{
    const SquareMatrix* yy_ptr = dynamic_cast<const SquareMatrix*>(& y);
    Array2D::copyData(*yy_ptr);
}

size_t  SquareMatrix::nRows() const
{
    return m_nrows;
}

size_t SquareMatrix::nRowsAndStruct(size_t* const iStruct) const
{
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

doublereal*   const* SquareMatrix::colPts()
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
