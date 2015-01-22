/**
 *  @file BandMatrix.cpp
 *
 *  Banded matrices.
 */

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
        m_colPtrs[j] = &(data[ldab * j]);
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
        m_colPtrs[j] = &(data[ldab * j]);
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
        m_colPtrs[j] = &(data[ldab * j]);
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
        m_colPtrs[j] = &(data[ldab * j]);
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
    int jj = static_cast<int>(j);
    int ii = static_cast<int>(i);
    size_t rw = (int) m_kl + (int) m_ku + (int) ii - jj;
    return (2*m_kl + m_ku + 1)*j + rw;
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

vector_int&   BandMatrix::ipiv()
{
    return m_ipiv;
}

void BandMatrix::mult(const doublereal* b, doublereal* prod) const
{
    int kl = static_cast<int>(m_kl);
    int ku = static_cast<int>(m_ku);
    int nr = static_cast<int>(nRows());
    doublereal sum = 0.0;
    for (int m = 0; m < nr; m++) {
        sum = 0.0;
        for (int j = m - kl; j <= m + ku; j++) {
            if (j >= 0 && j < (int) m_n) {
                sum += _value(m,j) * b[j];
            }
        }
        prod[m] = sum;
    }
}

void BandMatrix::leftMult(const doublereal* const b, doublereal* const prod) const
{
    int kl = static_cast<int>(m_kl);
    int ku = static_cast<int>(m_ku);
    int nc = static_cast<int>(nColumns());
    doublereal sum = 0.0;
    for (int n = 0; n < nc; n++) {
        sum = 0.0;
        for (int i = n - ku; i <= n + kl; i++) {
            if (i >= 0 && i < (int) m_n) {
                size_t ii = i;
                sum += _value(ii,n) * b[ii];
            }
        }
        prod[n] = sum;
    }
}

int BandMatrix::factor()
{
    int info=0;
    copy(data.begin(), data.end(), ludata.begin());
    ct_dgbtrf(nRows(), nColumns(), nSubDiagonals(), nSuperDiagonals(),
              DATA_PTR(ludata), ldim(), DATA_PTR(ipiv()), info);

    // if info = 0, LU decomp succeeded.
    if (info == 0) {
        m_factored = true;
    } else {
        m_factored = false;
        ofstream fout("bandmatrix.csv");
        fout << *this << endl;
        fout.close();
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
    if (info == 0)
        ct_dgbtrs(ctlapack::NoTranspose, nColumns(), nSubDiagonals(),
                  nSuperDiagonals(), nrhs, DATA_PTR(ludata), ldim(),
                  DATA_PTR(ipiv()), b, ldb, info);

    // error handling
    if (info != 0) {
        ofstream fout("bandmatrix.csv");
        fout << *this << endl;
        fout.close();
    }
    return info;
}

vector_fp::iterator  BandMatrix::begin()
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

doublereal  BandMatrix::rcond(doublereal a1norm)
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
    rcond = ct_dgbcon('1', m_n, m_kl, m_ku, DATA_PTR(ludata), ldab, DATA_PTR(m_ipiv), a1norm, DATA_PTR(work_),
                      DATA_PTR(iwork_), rinfo);
    if (rinfo != 0) {
        if (printLevel) {
            writelogf("BandMatrix::rcond(): DGBCON returned INFO = %d\n", rinfo);
        }
        if (! useReturnErrorCode) {
            throw CanteraError("BandMatrix::rcond()", "DGBCON returned INFO = " + int2str(rinfo));
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
    int ku = static_cast<int>(m_ku);
    int kl = static_cast<int>(m_kl);
    doublereal value = 0.0;
    for (int j = 0; j < (int) m_n; j++) {
        doublereal sum = 0.0;
        doublereal* colP =  m_colPtrs[j];
        for (int i = j - ku; i <= j + kl; i++) {
            sum += fabs(colP[kl + ku + i - j]);
        }
        value = std::max(sum, value);
    }
    return value;
}

size_t BandMatrix::checkRows(doublereal& valueSmall) const
{
    valueSmall = 1.0E300;
    size_t iSmall = npos;
    for (int i = 0; i < (int) m_n; i++) {
        double valueS = 0.0;
        for (int j = i - (int) m_kl; j <= i + (int) m_ku; j++) {
            if (j >= 0 && j < (int) m_n) {
                valueS = std::max(fabs(value(i,j)), valueS);
            }
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
    for (int j = 0; j < (int) m_n; j++) {
        double valueS = 0.0;
        for (int i = j - (int) m_ku; i <= j + (int) m_kl; i++) {
            if (i >= 0 && i < (int) m_n) {
                valueS = std::max(fabs(value(i,j)), valueS);
            }
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
    return &(m_colPtrs[0]);
}

void BandMatrix::copyData(const GeneralMatrix& y)
{
    warn_deprecated("BandMatrix::copyData", "To be removed after Cantera 2.2.");
    m_factored = false;
    size_t n = sizeof(doublereal) * m_n * (2 *m_kl + m_ku + 1);
    GeneralMatrix* yyPtr = const_cast<GeneralMatrix*>(&y);
    (void) memcpy(DATA_PTR(data), yyPtr->ptrColumn(0), n);
}

void BandMatrix::useFactorAlgorithm(int fAlgorithm)
{
//    useQR_ = fAlgorithm;
}


}
