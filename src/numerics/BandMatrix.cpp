//! @file BandMatrix.cpp Banded matrices.

// Copyright 2001  California Institute of Technology

#include "cantera/numerics/BandMatrix.h"
#include "cantera/numerics/ctlapack.h"
#include "cantera/base/utilities.h"
#include "cantera/base/stringUtils.h"

#include <cstring>
#include <fstream>

using namespace std;

namespace Cantera
{

BandMatrix::BandMatrix() :
    GeneralMatrix(1),
    m_n(0),
    m_kl(0),
    m_ku(0),
    m_zero(0.0)
{
}

BandMatrix::BandMatrix(size_t n, size_t kl, size_t ku, doublereal v)   :
    GeneralMatrix(1),
    m_n(n),
    m_kl(kl),
    m_ku(ku),
    m_zero(0.0)
{
    data.resize(n*(2*kl + ku + 1));
    ludata.resize(n*(2*kl + ku + 1));
    fill(data.begin(), data.end(), v);
    fill(ludata.begin(), ludata.end(), 0.0);
    m_ipiv.resize(m_n);
    m_colPtrs.resize(n);
    size_t ldab = (2*kl + ku + 1);
    for (size_t j = 0; j < n; j++) {
        m_colPtrs[j] = &data[ldab * j];
    }
}

BandMatrix::BandMatrix(const BandMatrix& y) :
    GeneralMatrix(y),
    m_n(0),
    m_kl(0),
    m_ku(0),
    m_zero(0.0)
{
    m_n = y.m_n;
    m_kl = y.m_kl;
    m_ku = y.m_ku;
    data = y.data;
    ludata = y.ludata;
    m_ipiv = y.m_ipiv;
    m_colPtrs.resize(m_n);
    size_t ldab = (2 *m_kl + m_ku + 1);
    for (size_t j = 0; j < m_n; j++) {
        m_colPtrs[j] = &data[ldab * j];
    }
}

BandMatrix& BandMatrix::operator=(const BandMatrix& y)
{
    if (&y == this) {
        return *this;
    }
    GeneralMatrix::operator=(y);
    m_n = y.m_n;
    m_kl = y.m_kl;
    m_ku = y.m_ku;
    m_ipiv = y.m_ipiv;
    data = y.data;
    ludata = y.ludata;
    m_colPtrs.resize(m_n);
    size_t ldab = (2 * m_kl + m_ku + 1);
    for (size_t j = 0; j < m_n; j++) {
        m_colPtrs[j] = &data[ldab * j];
    }
    return *this;
}

void BandMatrix::resize(size_t n, size_t kl, size_t ku, doublereal v)
{
    m_n = n;
    m_kl = kl;
    m_ku = ku;
    data.resize(n*(2*kl + ku + 1));
    ludata.resize(n*(2*kl + ku + 1));
    m_ipiv.resize(m_n);
    fill(data.begin(), data.end(), v);
    m_colPtrs.resize(m_n);
    size_t ldab = (2 * m_kl + m_ku + 1);
    for (size_t j = 0; j < n; j++) {
        m_colPtrs[j] = &data[ldab * j];
    }
    m_factored = false;
}

void BandMatrix::bfill(doublereal v)
{
    std::fill(data.begin(), data.end(), v);
    m_factored = false;
}

void BandMatrix::zero()
{
    std::fill(data.begin(), data.end(), 0.0);
    m_factored = false;
}

doublereal& BandMatrix::operator()(size_t i, size_t j)
{
    return value(i,j);
}

doublereal BandMatrix::operator()(size_t i, size_t j) const
{
    return value(i,j);
}

doublereal& BandMatrix::value(size_t i, size_t j)
{
    m_factored = false;
    if (i + m_ku < j || i > j + m_kl) {
        return m_zero;
    }
    return data[index(i,j)];
}

doublereal BandMatrix::value(size_t i, size_t j) const
{
    if (i + m_ku < j || i > j + m_kl) {
        return 0.0;
    }
    return data[index(i,j)];
}

size_t BandMatrix::index(size_t i, size_t j) const
{
    return (2*m_kl + m_ku)*j + m_kl + m_ku + i;
}

doublereal BandMatrix::_value(size_t i, size_t j) const
{
    return data[index(i,j)];
}

size_t BandMatrix::nRows() const
{
    return m_n;
}

size_t BandMatrix::nRowsAndStruct(size_t* const iStruct) const
{
    if (iStruct) {
        iStruct[0] = m_kl;
        iStruct[1] = m_ku;
    }
    return m_n;
}

size_t BandMatrix::nColumns() const
{
    return m_n;
}

size_t BandMatrix::nSubDiagonals() const
{
    return m_kl;
}

size_t BandMatrix::nSuperDiagonals() const
{
    return m_ku;
}

size_t BandMatrix::ldim() const
{
    return 2*m_kl + m_ku + 1;
}

vector_int& BandMatrix::ipiv()
{
    return m_ipiv;
}

void BandMatrix::mult(const doublereal* b, doublereal* prod) const
{
    for (size_t m = 0; m < m_n; m++) {
        double sum = 0.0;
        size_t start = (m >= m_kl) ? m - m_kl : 0;
        size_t stop = std::min(m + m_ku + 1, m_n);
        for (size_t j = start; j < stop; j++) {
            sum += _value(m,j) * b[j];
        }
        prod[m] = sum;
    }
}

void BandMatrix::leftMult(const doublereal* const b, doublereal* const prod) const
{
    for (size_t n = 0; n < m_n; n++) {
        double sum = 0.0;
        size_t start = (n >= m_ku) ? n - m_ku : 0;
        size_t stop = std::min(n + m_kl + 1, m_n);
        for (size_t i = start; i < stop; i++) {
            sum += _value(i,n) * b[i];
        }
        prod[n] = sum;
    }
}

int BandMatrix::factor()
{
    int info=0;
    ludata = data;
    ct_dgbtrf(nRows(), nColumns(), nSubDiagonals(), nSuperDiagonals(),
              ludata.data(), ldim(), ipiv().data(), info);

    // if info = 0, LU decomp succeeded.
    if (info == 0) {
        m_factored = true;
    } else {
        m_factored = false;
        ofstream fout("bandmatrix.csv");
        fout << *this << endl;
    }
    return info;
}

int BandMatrix::solve(const doublereal* const b, doublereal* const x)
{
    copy(b, b + m_n, x);
    return solve(x);
}

int BandMatrix::solve(doublereal* b, size_t nrhs, size_t ldb)
{
    int info = 0;
    if (!m_factored) {
        info = factor();
    }
    if (ldb == 0) {
        ldb = nColumns();
    }
    if (info == 0) {
        ct_dgbtrs(ctlapack::NoTranspose, nColumns(), nSubDiagonals(),
                  nSuperDiagonals(), nrhs, ludata.data(), ldim(),
                  ipiv().data(), b, ldb, info);
    }

    // error handling
    if (info != 0) {
        ofstream fout("bandmatrix.csv");
        fout << *this << endl;
    }
    return info;
}

vector_fp::iterator BandMatrix::begin()
{
    m_factored = false;
    return data.begin();
}

vector_fp::iterator BandMatrix::end()
{
    m_factored = false;
    return data.end();
}

vector_fp::const_iterator BandMatrix::begin() const
{
    return data.begin();
}

vector_fp::const_iterator BandMatrix::end() const
{
    return data.end();
}

ostream& operator<<(ostream& s, const BandMatrix& m)
{
    size_t nr = m.nRows();
    size_t nc = m.nColumns();
    for (size_t i = 0; i < nr; i++) {
        for (size_t j = 0; j < nc; j++) {
            s << m(i,j) << ", ";
        }
        s << endl;
    }
    return s;
}

doublereal BandMatrix::rcond(doublereal a1norm)
{
    int printLevel = 0;
    int useReturnErrorCode = 0;
    if (iwork_.size() < m_n) {
        iwork_.resize(m_n);
    }
    if (work_.size() < 3 * m_n) {
        work_.resize(3 * m_n);
    }
    doublereal rcond = 0.0;
    if (m_factored != 1) {
        throw CanteraError("BandMatrix::rcond()", "matrix isn't factored correctly");
    }

    size_t ldab = (2 *m_kl + m_ku + 1);
    int rinfo = 0;
    rcond = ct_dgbcon('1', m_n, m_kl, m_ku, ludata.data(), ldab, m_ipiv.data(), a1norm, work_.data(),
                      iwork_.data(), rinfo);
    if (rinfo != 0) {
        if (printLevel) {
            writelogf("BandMatrix::rcond(): DGBCON returned INFO = %d\n", rinfo);
        }
        if (! useReturnErrorCode) {
            throw CanteraError("BandMatrix::rcond()", "DGBCON returned INFO = {}", rinfo);
        }
    }
    return rcond;
}

int BandMatrix::factorAlgorithm() const
{
    return 0;
}

doublereal BandMatrix::oneNorm() const
{
    double value = 0.0;
    for (size_t j = 0; j < m_n; j++) {
        double sum = 0.0;
        size_t start = (j >= m_ku) ? j - m_ku : 0;
        size_t stop = std::min(j + m_kl + 1, m_n);
        for (size_t i = start; i < stop; i++) {
            sum += std::abs(_value(i,j));
        }
        value = std::max(sum, value);
    }
    return value;
}

size_t BandMatrix::checkRows(doublereal& valueSmall) const
{
    valueSmall = 1.0E300;
    size_t iSmall = npos;
    for (size_t i = 0; i < m_n; i++) {
        double valueS = 0.0;
        size_t start = (i > m_kl) ? i - m_kl : 0;
        size_t stop = std::min(i + m_ku + 1, m_n);
        for (size_t j = start; j < stop; j++) {
            valueS = std::max(fabs(_value(i,j)), valueS);
        }
        if (valueS < valueSmall) {
            iSmall = i;
            valueSmall = valueS;
            if (valueSmall == 0.0) {
                return iSmall;
            }
        }
    }
    return iSmall;
}

size_t BandMatrix::checkColumns(doublereal& valueSmall) const
{
    valueSmall = 1.0E300;
    size_t jSmall = npos;
    for (size_t j = 0; j < m_n; j++) {
        double valueS = 0.0;
        size_t start = (j > m_ku) ? j - m_ku : 0;
        size_t stop = std::min(j + m_kl + 1, m_n);
        for (size_t i = start; i < stop; i++) {
            valueS = std::max(fabs(_value(i,j)), valueS);
        }
        if (valueS < valueSmall) {
            jSmall = j;
            valueSmall = valueS;
            if (valueSmall == 0.0) {
                return jSmall;
            }
        }
    }
    return jSmall;
}

GeneralMatrix* BandMatrix::duplMyselfAsGeneralMatrix() const
{
    return new BandMatrix(*this);
}

doublereal* BandMatrix::ptrColumn(size_t j)
{
    return m_colPtrs[j];
}

doublereal* const* BandMatrix::colPts()
{
    return &m_colPtrs[0];
}

}
