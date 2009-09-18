/**
 *  @file DoubleStarStar.cpp
 *
 *  Header file for class DoubleStarStar
 */

/*
 *  $Author: hkmoffa $
 *  $Revision: 1.2 $
 *  $Date: 2008/08/01 16:05:38 $
 */

#include "vcs_DoubleStarStar.h"

namespace VCSnonideal {

//!Default constructor. Create an empty array.
DoubleStarStar::DoubleStarStar() :
  m_nrows(0),
  m_ncols(0)
{
  m_data.clear();
  m_colAddr.clear();
}

/*
 *  Constructor. Create an \c m by \c n array, and initialize
 *  all elements to \c v.
 */
DoubleStarStar::DoubleStarStar(int m, int n, double v) :
  m_nrows(n),
  m_ncols(m) 
{
  m_data.resize(n*m);
  std::fill(m_data.begin(), m_data.end(), v);
  m_colAddr.resize(m);
  for (int jcol = 0; jcol < m_ncols; jcol++) {
    m_colAddr[jcol] = &(m_data[jcol*m_nrows]);
  }
}

// copy constructor
DoubleStarStar::DoubleStarStar(const DoubleStarStar& y) {
  m_nrows = y.m_nrows;
  m_ncols = y.m_ncols;
  m_data.resize(m_nrows*m_ncols);
  m_data = y.m_data;
  m_colAddr.resize(m_ncols);
  for (int jcol = 0; jcol < m_ncols; jcol++) {
    m_colAddr[jcol] = &(m_data[jcol*m_nrows]);
  }
}

// assignment operator
DoubleStarStar& DoubleStarStar::operator=(const DoubleStarStar& y) {
  if (&y == this) return *this;
  m_nrows = y.m_nrows;
  m_ncols = y.m_ncols;
  m_data.resize(m_nrows*m_ncols);
  m_data = y.m_data;
  m_colAddr.resize(m_ncols);
  for (int jcol = 0; jcol < m_ncols; jcol++) {
    m_colAddr[jcol] = &(m_data[jcol*m_nrows]);
  }
  return *this;
}


// resize the array, and fill the new entries with 'v'
/*
 * @param n  This is the number of rows
 * @param m  This is the number of columns in the new matrix
 * @param v  Default fill value -> defaults to zero.
 */
void DoubleStarStar::resize(int m, int n, double v) {
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
      for (int jcol = 0; jcol < m_ncols; jcol++) {
        for (int irow = 0; irow < m_nrows; irow++) {
          m_data[jcol*n + irow] = old_data[jcol*m_nrows + irow];
        }
        for (int irow = m_nrows; irow < n; irow++) {
          m_data[jcol*n + irow] = v;
        }
      }
      for (int jcol = m_ncols; jcol < m; jcol++) {
        for (int irow = 0; irow < n; irow++) {
          m_data[jcol*n + irow] = v;
        }
      }
    } else {
      std::fill(m_data.begin(), m_data.end(), v);
      for (int jcol = 0; jcol < m_ncols; jcol++) {
        for (int irow = 0; irow < m_nrows; irow++) {
          m_data[jcol*n + irow] = old_data[jcol*m_nrows + irow];
        }
      }
    }
  }
  m_nrows = n;
  m_ncols = m;
  m_colAddr.resize(m_ncols);
  for (int jcol = 0; jcol < m_ncols; jcol++) {
    m_colAddr[jcol] = &(m_data[jcol*m_nrows]);
  }
}

double * const DoubleStarStar::operator[](int jcol) {
  return m_colAddr[jcol];
}

const double * const DoubleStarStar::operator[](int jcol) const {
  return (const double * const) m_colAddr[jcol];
}

double * const * const DoubleStarStar::baseDataAddr() {
  return (double * const * const) &(m_colAddr[0]);
}

double const * const * const DoubleStarStar::constBaseDataAddr() const {
  return (double const * const * const) &(m_colAddr[0]);
}

// Number of rows
int DoubleStarStar::nRows() const { 
  return m_nrows; 
}

// Number of columns
int DoubleStarStar::nColumns() const { 
  return m_ncols;
}

} 
  
