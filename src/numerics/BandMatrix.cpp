//! @file BandMatrix.cpp Banded matrices.

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/numerics/BandMatrix.h"
#include "cantera/base/utilities.h"
#include "cantera/base/stringUtils.h"

#if CT_USE_LAPACK
    #include "cantera/numerics/ctlapack.h"
#else
    #if CT_SUNDIALS_USE_LAPACK
        #if CT_SUNDIALS_VERSION >= 30
            #include "sunlinsol/sunlinsol_lapackband.h"
        #else
            #include "cvodes/cvodes_lapack.h"
        #endif
    #else
        #if CT_SUNDIALS_VERSION >= 30
            #include "sunlinsol/sunlinsol_band.h"
        #else
            #include "cvodes/cvodes_dense.h"
            #include "cvodes/cvodes_band.h"
        #endif
    #endif
#endif

#include <cstring>
#include <fstream>

using namespace std;

namespace Cantera
{

// pImpl wrapper class for vector of Sundials index types to avoid needing to
// include Sundials headers in BandMatrix.h
struct BandMatrix::PivData {
#if CT_USE_LAPACK
    vector_int data;
#elif CT_SUNDIALS_VERSION >= 30
    std::vector<sunindextype> data;
#elif CT_SUNDIALS_VERSION >= 25
    std::vector<long int> data;
#else
    std::vector<int> data;
#endif
};

BandMatrix::BandMatrix() :
    m_n(0),
    m_kl(0),
    m_ku(0),
    m_zero(0.0),
    m_ipiv{new PivData()},
    m_info(0)
{
}

BandMatrix::~BandMatrix()
{
    // Needs to be defined here so m_ipiv can be deleted
}

BandMatrix::BandMatrix(size_t n, size_t kl, size_t ku, doublereal v)   :
    m_n(n),
    m_kl(kl),
    m_ku(ku),
    m_zero(0.0),
    m_ipiv{new PivData()},
    m_info(0)
{
    data.resize(n*(2*kl + ku + 1));
    ludata.resize(n*(2*kl + ku + 1));
    fill(data.begin(), data.end(), v);
    fill(ludata.begin(), ludata.end(), 0.0);
    m_ipiv->data.resize(m_n);
    m_colPtrs.resize(n);
    m_lu_col_ptrs.resize(n);
    size_t ldab = (2*kl + ku + 1);
    for (size_t j = 0; j < n; j++) {
        m_colPtrs[j] = &data[ldab * j];
        m_lu_col_ptrs[j] = &ludata[ldab * j];
    }
}

BandMatrix::BandMatrix(const BandMatrix& y) :
    GeneralMatrix(y),
    m_n(0),
    m_kl(0),
    m_ku(0),
    m_zero(0.0),
    m_ipiv{new PivData()},
    m_info(y.m_info)
{
    m_n = y.m_n;
    m_kl = y.m_kl;
    m_ku = y.m_ku;
    data = y.data;
    ludata = y.ludata;
    m_ipiv->data = y.m_ipiv->data;
    m_colPtrs.resize(m_n);
    m_lu_col_ptrs.resize(m_n);
    size_t ldab = (2 *m_kl + m_ku + 1);
    for (size_t j = 0; j < m_n; j++) {
        m_colPtrs[j] = &data[ldab * j];
        m_lu_col_ptrs[j] = &ludata[ldab * j];
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
    m_ipiv->data = y.m_ipiv->data;
    data = y.data;
    ludata = y.ludata;
    m_colPtrs.resize(m_n);
    m_lu_col_ptrs.resize(m_n);
    size_t ldab = (2 * m_kl + m_ku + 1);
    for (size_t j = 0; j < m_n; j++) {
        m_colPtrs[j] = &data[ldab * j];
        m_lu_col_ptrs[j] = &ludata[ldab * j];
    }
    m_info = y.m_info;
    return *this;
}

void BandMatrix::resize(size_t n, size_t kl, size_t ku, doublereal v)
{
    m_n = n;
    m_kl = kl;
    m_ku = ku;
    data.resize(n*(2*kl + ku + 1));
    ludata.resize(n*(2*kl + ku + 1));
    m_ipiv->data.resize(m_n);
    fill(data.begin(), data.end(), v);
    m_colPtrs.resize(m_n);
    m_lu_col_ptrs.resize(m_n);
    size_t ldab = (2 * m_kl + m_ku + 1);
    for (size_t j = 0; j < n; j++) {
        m_colPtrs[j] = &data[ldab * j];
        m_lu_col_ptrs[j] = &ludata[ldab * j];
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
    ludata = data;
#if CT_USE_LAPACK
    ct_dgbtrf(nRows(), nColumns(), nSubDiagonals(), nSuperDiagonals(),
              ludata.data(), ldim(), m_ipiv->data.data(), m_info);
#else
    long int nu = static_cast<long int>(nSuperDiagonals());
    long int nl = static_cast<long int>(nSubDiagonals());
    long int smu = nu + nl;
    m_info = bandGBTRF(m_lu_col_ptrs.data(), static_cast<long int>(nColumns()),
                       nu, nl, smu, m_ipiv->data.data());
#endif
    if (m_info != 0) {
        throw Cantera::CanteraError("BandMatrix::factor",
            "Factorization failed with DGBTRF error code {}.", m_info);
    }
    m_factored = true;
    return m_info;
}

int BandMatrix::solve(const doublereal* const b, doublereal* const x)
{
    copy(b, b + m_n, x);
    return solve(x);
}

int BandMatrix::solve(doublereal* b, size_t nrhs, size_t ldb)
{
    if (!m_factored) {
        factor();
    }
    if (ldb == 0) {
        ldb = nColumns();
    }
#if CT_USE_LAPACK
    ct_dgbtrs(ctlapack::NoTranspose, nColumns(), nSubDiagonals(),
              nSuperDiagonals(), nrhs, ludata.data(), ldim(),
              m_ipiv->data.data(), b, ldb, m_info);
#else
    long int nu = static_cast<long int>(nSuperDiagonals());
    long int nl = static_cast<long int>(nSubDiagonals());
    long int smu = nu + nl;
    double** a = m_lu_col_ptrs.data();
    bandGBTRS(a, static_cast<long int>(nColumns()), smu, nl, m_ipiv->data.data(), b);
    m_info = 0;
#endif

    if (m_info != 0) {
        throw Cantera::CanteraError("BandMatrix::solve",
            "Linear solve failed with DGBTRS error code {}.", m_info);
    }
    return m_info;
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
    for (size_t i = 0; i < m.nRows(); i++) {
        s << m(i, 0);
        for (size_t j = 1; j < m.nColumns(); j++) {
            s << ", " << m(i,j);
        }
        s << endl;
    }
    return s;
}

doublereal BandMatrix::rcond(doublereal a1norm)
{
    iwork_.resize(m_n);
    work_.resize(3 * m_n);

    if (m_factored != 1) {
        throw CanteraError("BandMatrix::rcond", "matrix isn't factored correctly");
    }

#if CT_USE_LAPACK
    size_t ldab = (2 *m_kl + m_ku + 1);
    int rinfo = 0;
    double rcond = ct_dgbcon('1', m_n, m_kl, m_ku, ludata.data(),
        ldab, m_ipiv->data.data(), a1norm, work_.data(), iwork_.data(), rinfo);
    if (rinfo != 0) {
        throw CanteraError("BandMatrix::rcond", "DGBCON returned INFO = {}", rinfo);
    }
    return rcond;
#else
    throw CanteraError("BandMatrix::rcond", "not implemented when LAPACK is missing");
#endif
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

doublereal* BandMatrix::ptrColumn(size_t j)
{
    return m_colPtrs[j];
}

doublereal* const* BandMatrix::colPts()
{
    return &m_colPtrs[0];
}

}
