/**
 *  @file vcs_IntStarStar.cpp Implementation of class IntStarStar
 */
#include "cantera/equil/vcs_IntStarStar.h"

namespace VCSnonideal
{
IntStarStar::IntStarStar() :
    m_nrows(0),
    m_ncols(0)
{
    m_data.clear();
    m_colAddr.clear();
}

IntStarStar::IntStarStar(size_t m, size_t n, int v) :
    m_nrows(n),
    m_ncols(m)
{
    m_data.resize(n*m);
    std::fill(m_data.begin(), m_data.end(), v);
    m_colAddr.resize(m);
    if (!m_data.empty()) {
        for (size_t jcol = 0; jcol < m_ncols; jcol++) {
            m_colAddr[jcol] = &(m_data[jcol*m_nrows]);
        }
    }
}

IntStarStar::IntStarStar(const IntStarStar& y)
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

IntStarStar& IntStarStar::operator=(const IntStarStar& y)
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

void IntStarStar::resize(size_t m, size_t n, int v)
{
    std::vector<int> old_data;
    bool doCopy = false;
    if (m_nrows > 0 && m_ncols > 0) {
        if (m_ncols != m) {
            doCopy = true;
            old_data = m_data;
        }
    }
    m_data.resize(n*m, v);
    if (doCopy) {
        if (n >= m_nrows && m >= m_ncols) {
            for (size_t jcol = 0; jcol < m_ncols; jcol++) {
                for (size_t irow = 0; irow < m_nrows; irow++) {
                    m_data[jcol*m + irow] = old_data[jcol*m_ncols + irow];
                }
                for (size_t irow = m_nrows; irow < n; irow++) {
                    m_data[jcol*m + irow] = v;
                }
            }
            for (size_t jcol = m_ncols; jcol < m; jcol++) {
                for (size_t irow = 0; irow < n; irow++) {
                    m_data[jcol*m + irow] = v;
                }
            }
        } else {
            std::fill(m_data.begin(), m_data.end(), v);
            for (size_t jcol = 0; jcol < m_ncols; jcol++) {
                for (size_t irow = 0; irow < m_nrows; irow++) {
                    m_data[jcol*m + irow] = old_data[jcol*m_ncols + irow];
                }
            }
        }
    }
    m_nrows = n;
    m_ncols = m;
    m_colAddr.resize(m_ncols);
    for (size_t jcol = 0; jcol < m_ncols; jcol++) {
        m_colAddr[jcol] = &(m_data[jcol*m_nrows]);
    }
}

int* IntStarStar::operator[](size_t jcol)
{
    return m_colAddr[jcol];
}

const int* IntStarStar::operator[](size_t jcol) const
{
    return (const int*) m_colAddr[jcol];
}

int* const* IntStarStar::baseDataAddr()
{
    return (int* const*) &(m_colAddr[0]);
}

size_t IntStarStar::nRows() const
{
    return m_nrows;
}

size_t IntStarStar::nColumns() const
{
    return m_ncols;
}

}
