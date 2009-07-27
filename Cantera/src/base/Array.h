/**
 *  @file Array.h
 *  Header file for class Array2D
 */

/*
 *  $Author: hkmoffa $
 *  $Revision: 1.5 $
 *  $Date: 2009/01/15 23:05:34 $
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_ARRAY_H
#define CT_ARRAY_H

#include "ct_defs.h"
#include "ctexceptions.h"
#include "utilities.h"

namespace Cantera { 


 
  //!  A class for 2D arrays stored in column-major
  //!  (Fortran-compatible) form.
  /*!
   *  In this form, the data entry for an n row, m col 
   *  matrix is 
   *       index = i + (n-1) * j
   *  where
   *     J(i,j) = data_start + index 
   *         i = row
   *         j = column
   */
  class Array2D {

  public:

    typedef vector_fp::iterator iterator;
    typedef vector_fp::const_iterator const_iterator;

    /**
     * Default constructor. Create an empty array.
     */
    Array2D() : m_nrows(0), m_ncols(0) { m_data.clear(); }


    //!  Constructor.
    /*!
     *    Create an \c m by \c n array, and initialize
     *    all elements to \c v.
     *
     *  @param m   Number of rows
     *  @param n   Number of columns
     *  @param v   Default fill value. The default is 0.0
     */
    Array2D(const int m, const int n, const doublereal v = 0.0) 
      : m_nrows(m), m_ncols(n) {
      m_data.resize(n*m);
      std::fill(m_data.begin(), m_data.end(), v);
    }

    //!  Copy constructor
    /*!
     *   @param y  Array2D to make the copy from
     */
    Array2D(const Array2D& y) {
      m_nrows = y.m_nrows;
      m_ncols = y.m_ncols;
      m_data.resize(m_nrows*m_ncols);
      m_data = y.m_data;
    }

    //! assignment operator
    /*!
     *  @param y Array2D to get the values from
     */
    Array2D& operator=(const Array2D& y) {
      if (&y == this) return *this;
      m_nrows = y.m_nrows;
      m_ncols = y.m_ncols;
      m_data.resize(m_nrows*m_ncols);
      m_data = y.m_data;
      return *this;
    }

    //! resize the array, and fill the new entries with 'v'
    /*!
     * @param n  This is the number of rows
     * @param m  This is the number of columns in the new matrix
     * @param v  Default fill value -> defaults to zero.
     */
    void resize(int n, int m, doublereal v = 0.0) {
      m_nrows = n;
      m_ncols = m;
      m_data.resize(n*m, v);
    }

    /// append a column
    void appendColumn(const vector_fp& c) {
      m_ncols++;
      m_data.resize(m_nrows*m_ncols);
      int m;
      for (m = 0;  m < m_nrows; m++) value(m_ncols, m) = c[m];
    }

    /// append a column
    void appendColumn(doublereal* c) {
      m_ncols++;
      m_data.resize(m_nrows*m_ncols);
      int m;
      for (m = 0;  m < m_nrows; m++) value(m_ncols, m) = c[m];
    }

    /// set the nth row to array rw
    void setRow(int n, doublereal* rw) {
      for (int j = 0; j < m_ncols; j++) {
	m_data[m_nrows*j + n] = rw[j];
      }
    }

    /// get the nth row
    void getRow(int n, doublereal* rw) {
      for (int j = 0; j < m_ncols; j++) {
	rw[j] = m_data[m_nrows*j + n];
      }
    }

    /// set the values in column m to those in array col
    void setColumn(int m, doublereal* col) {
      for (int i = 0; i < m_nrows; i++) {
	m_data[m_nrows*m + i] = col[i];
      }
    }

    /// get the values in column m
    void getColumn(int m, doublereal* col) {
      for (int i = 0; i < m_nrows; i++) {
	col[i] = m_data[m_nrows*m + i];
      }
    }
                
    /**
     * Destructor. Does nothing, since no memory allocated on the
     * heap.
     */
    virtual ~Array2D(){}


    /**
     * Evaluate a*x + y.
     */
    void axpy(doublereal a, const Array2D& x, const Array2D& y) {
      iterator b = begin();
      const_iterator xb = x.begin();
      const_iterator yb = y.begin();
      for (; b != end(); ++b, ++xb, ++yb)  *b = a*(*xb) + *yb;
    }

    
    //! Allows setting elements using the syntax A(i,j) = x.
    /*!
     *  @param  i            row index
     *  @param  j            column index.
     *
     *  @return Returns a reference to A(i,j) which may be assigned.
     */ 
    doublereal& operator()( int i, int j) { return value(i,j); }

    /**
     * Allows retrieving elements using the syntax x = A(i,j).
     */ 
    doublereal operator() ( int i, int j) const {return value(i,j);}

    //! Returns a changeable reference to position in the matrix
    /*!
     * This is a key entry. Returns a reference to the matrixes (i,j)
     * element. This may be used as an L value.
     * @param i   The row index
     * @param j   The column index
     */
    doublereal& value( int i, int j) {return m_data[m_nrows*j + i];}
    doublereal value( int i, int j) const {return m_data[m_nrows*j + i];}

    /// Number of rows
    size_t nRows() const { return m_nrows; }

    /// Number of columns
    size_t nColumns() const { return m_ncols; }

    /// Return an iterator pointing to the first element
    iterator begin() { return m_data.begin(); }

    /// Return an iterator pointing past the last element
    iterator end() { return m_data.end(); }

    /// Return a const iterator pointing to the first element
    const_iterator begin() const { return m_data.begin(); }

    /// Return a const iterator pointing to past the last element
    const_iterator end() const { return m_data.end(); }

    /// Return a reference to the data vector
    vector_fp& data() { return m_data; }

    /// Return a const reference to the data vector
    const vector_fp& data() const { return m_data; }

    /// Return a pointer to the top of column j, columns are contiguous
    /// in memory
    doublereal * ptrColumn(int j) { return &(m_data[m_nrows*j]); }
    const doublereal * ptrColumn(int j) const { return &(m_data[m_nrows*j]); }

  protected:

    //! Data storred in a single array
    vector_fp m_data;

    //! Number of rows
    int m_nrows;

    //! Number of columns
    int m_ncols;
  };

  /// output the array
  inline std::ostream& operator<<(std::ostream& s, const Array2D& m) {
    int nr = static_cast<int>(m.nRows());
    int nc = static_cast<int>(m.nColumns());
    int i,j;
    for (i = 0; i < nr; i++) {
      for (j = 0; j < nc; j++) {
	s << m(i,j) << ", ";
      }
      s << std::endl;
    }
    return s;
  }

  //! Overload the times equals operator for multiplication
  //! of a matrix and a scalar.
  /*!
   *  Scaled every element of the matrix by the scalar input
   *
   *   @param m   Matrix
   *   @param a   scalar
   */
  inline void operator*=(Array2D& m, doublereal a) {
    scale(m.begin(), m.end(), m.begin(), a);
  }

  //! Overload the plus equals operator for addition
  //! of one matrix with another
  /*!
   *   Adds each element of the second matrix into the first
   *   matrix
   *
   *   @param x   First matrix
   *   @param y   Second matrix, which is a const
   */
  inline void operator+=(Array2D& x, const Array2D& y) {
    sum_each(x.begin(), x.end(), y.begin());
  }
  
}

#endif
