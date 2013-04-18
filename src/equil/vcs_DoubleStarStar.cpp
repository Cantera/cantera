/**
 *  @file vcs_DoubleStarStar.cpp
 *
 *  Implementation file for class DoubleStarStar
 */
#include "cantera/equil/vcs_DoubleStarStar.h"

namespace VCSnonideal
{

DoubleStarStar::DoubleStarStar() :
    m_nrows(0),
    m_ncols(0)
{
    m_data.clear();
    m_colAddr.clear();
}

DoubleStarStar::DoubleStarStar(size_t m, size_t n, double v) :
    m_nrows(n),
    m_ncols(m)
{
    m_data.resize(n*m);
    std::fill(m_data.begin(), m_data.end(), v);
    m_colAddr.resize(m);
    for (size_t jcol = 0; jcol < m_ncols; jcol++) {
        m_colAddr[jcol] = &(m_data[jcol*m_nrows]);
    }
}

DoubleStarStar::DoubleStarStar(const DoubleStarStar& y)
{
    m_nrows = y.m_nrows;
    m_ncols = y.m_ncols;
    m_data.resize(m_nrows*m_ncols);
    m_data = y.m_data;
    m_colAddr.resize(m_ncols);
    if (!m_data.empty()) {
        for (size_t jcol = 0; jcol < m_ncols; jcol++) {
            m_colAddr[jcol] = &(m_data[jcol*m_nrows]);
        }
    }
}

DoubleStarStar& DoubleStarStar::operator=(const DoubleStarStar& y)
{
    if (&y == this) {
        return *this;
    }
    m_nrows = y.m_nrows;
    m_ncols = y.m_ncols;
    m_data.resize(m_nrows*m_ncols);
    m_data = y.m_data;
    m_colAddr.resize(m_ncols);
    if (!m_data.empty()) {
        for (size_t jcol = 0; jcol < m_ncols; jcol++) {
            m_colAddr[jcol] = &(m_data[jcol*m_nrows]);
        }
    }
    return *this;
}

void DoubleStarStar::resize(size_t m, size_t n, double v)
{
    std::vector<double> old_data;
    bool doCopy = false;
    if (m_nrows > 0 && m_ncols > 0) {
        if (m_nrows != n) {
            doCopy = true;
            old_data = m_data;
        }
    }
    m_data.resize(n*m, v);
    if (doCopy) {
        if (n >= m_nrows && m >= m_ncols) {
            for (size_t jcol = 0; jcol < m_ncols; jcol++) {
                for (size_t irow = 0; irow < m_nrows; irow++) {
                    m_data[jcol*n + irow] = old_data[jcol*m_nrows + irow];
                }
                for (size_t irow = m_nrows; irow < n; irow++) {
                    m_data[jcol*n + irow] = v;
                }
            }
            for (size_t jcol = m_ncols; jcol < m; jcol++) {
                for (size_t irow = 0; irow < n; irow++) {
                    m_data[jcol*n + irow] = v;
                }
            }
        } else {
            std::fill(m_data.begin(), m_data.end(), v);
            for (size_t jcol = 0; jcol < m_ncols; jcol++) {
                for (size_t irow = 0; irow < m_nrows; irow++) {
                    m_data[jcol*n + irow] = old_data[jcol*m_nrows + irow];
                }
            }
        }
    }
    m_nrows = n;
    m_ncols = m;
    m_colAddr.resize(m_ncols);
    if (!m_data.empty()) {
        for (size_t jcol = 0; jcol < m_ncols; jcol++) {
            m_colAddr[jcol] = &(m_data[jcol*m_nrows]);
        }
    }
}

double* DoubleStarStar::operator[](size_t jcol)
{
    return m_colAddr[jcol];
}

const double* DoubleStarStar::operator[](size_t jcol) const
{
    return (const double*) m_colAddr[jcol];
}

double* const* DoubleStarStar::baseDataAddr()
{
    return (double* const*) &(m_colAddr[0]);
}

double const* const* DoubleStarStar::constBaseDataAddr() const
{
    return (double const* const*) &(m_colAddr[0]);
}

size_t DoubleStarStar::nRows() const
{
    return m_nrows;
}

size_t DoubleStarStar::nColumns() const
{
    return m_ncols;
}

}
