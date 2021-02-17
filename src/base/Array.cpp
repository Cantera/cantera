/**
 *  @file Array.cpp Implementation file for class Cantera::Array2D
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/Array.h"
#include "cantera/base/utilities.h"

namespace Cantera
{

Array2D::Array2D()
    : m_nrows(0)
    , m_ncols(0)
{
}

Array2D::Array2D(const size_t m, const size_t n, const double v)
    : m_nrows(m)
    , m_ncols(n)
{
    m_data.assign(n*m, v);
}

Array2D::Array2D(const size_t m, const size_t n, const double* values)
    : m_nrows(m)
    , m_ncols(n)
{
    m_data.assign(values, values + n*m);
}

Array2D::Array2D(const Array2D& y)
    : m_data(y.m_data)
    , m_nrows(y.m_nrows)
    , m_ncols(y.m_ncols)
{
}

Array2D& Array2D::operator=(const Array2D& y)
{
    if (&y == this) {
        return *this;
    }
    m_nrows = y.m_nrows;
    m_ncols = y.m_ncols;
    m_data = y.m_data;
    return *this;
}

void Array2D::resize(size_t n, size_t m, double v)
{
    m_nrows = n;
    m_ncols = m;
    m_data.resize(n*m, v);
}

void Array2D::appendColumn(const vector_fp& c)
{
    m_ncols++;
    m_data.resize(m_nrows * m_ncols);
    for (size_t m = 0; m < m_nrows; m++) {
        value(m_ncols, m) = c[m];
    }
}

void Array2D::appendColumn(const double* const c)
{
    m_ncols++;
    m_data.resize(m_nrows * m_ncols);
    for (size_t m = 0; m < m_nrows; m++) {
        value(m_ncols, m) = c[m];
    }
}

void Array2D::setRow(size_t n, const double* const rw)
{
    for (size_t j = 0; j < m_ncols; j++) {
        m_data[m_nrows*j + n] = rw[j];
    }
}

void Array2D::getRow(size_t n, double* const rw)
{
    for (size_t j = 0; j < m_ncols; j++) {
        rw[j] = m_data[m_nrows*j + n];
    }
}

void Array2D::setColumn(size_t m, double* const col)
{
    for (size_t i = 0; i < m_nrows; i++) {
        m_data[m_nrows*m + i] = col[i];
    }
}

void Array2D::getColumn(size_t m, double* const col)
{
    for (size_t i = 0; i < m_nrows; i++) {
        col[i] = m_data[m_nrows*m + i];
    }
}


std::ostream& operator<<(std::ostream& s, const Array2D& m)
{
    size_t nr = m.nRows();
    size_t nc = m.nColumns();
    for (size_t i = 0; i < nr; i++) {
        s << m(i,0);
        for (size_t j = 1; j < nc; j++) {
            s << ", " << m(i,j);
        }
        s << std::endl;
    }
    return s;
}

void operator*=(Array2D& m, double a)
{
    scale(m.begin(), m.end(), m.begin(), a);
}

}
