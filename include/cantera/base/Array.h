/**
 *  @file Array.h Header file for class Cantera::Array2D
 */
// Copyright 2001  California Institute of Technology


#ifndef CT_ARRAY_H
#define CT_ARRAY_H

#include "ct_defs.h"
#include "ctexceptions.h"
#include "utilities.h"

#include <iostream>
#include <cstring>

namespace Cantera
{



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
class Array2D
{

public:

    //! Type definition for the iterator class that is
    //! can be used by Array2D types.
    /*!
     *  this is just equal to vector_fp iterator.
     */
    typedef vector_fp::iterator iterator;


    //! Type definition for the const_iterator class that is
    //! can be used by Array2D types.
    /*!
     *  this is just equal to vector_fp const_iterator.
     */
    typedef vector_fp::const_iterator const_iterator;

    /**
     * Default constructor. Create an empty array.
     */
    Array2D() :
        m_data(0),
        m_nrows(0),
        m_ncols(0) {
    }

    //!  Constructor.
    /*!
     *    Create an \c m by \c n array, and initialize
     *    all elements to \c v.
     *
     *  @param m   Number of rows
     *  @param n   Number of columns
     *  @param v   Default fill value. The default is 0.0
     */
    Array2D(const size_t m, const size_t n, const doublereal v = 0.0)
        :  m_data(0), m_nrows(m), m_ncols(n) {
        m_data.resize(n*m);
        std::fill(m_data.begin(), m_data.end(), v);
    }

    //!  Constructor.
    /*!
     *    Create an \c m by \c n array, initialized with the contents
     *    of the array \c values.
     *
     *  @param m   Number of rows
     *  @param n   Number of columns
     *  @param values Initial values of the array. Must be of length m*n,
     *  and stored in column-major order.
     */
    Array2D(const size_t m, const size_t n, const doublereal* values)
        :  m_data(0), m_nrows(m), m_ncols(n) {
        m_data.resize(n*m);
        std::copy(values, values + m_nrows*m_ncols, m_data.begin());
    }

    //!  Copy constructor
    /*!
     *   @param y  Array2D to make the copy from
     */
    Array2D(const Array2D& y) :
        m_data(0),
        m_nrows(0),
        m_ncols(0) {
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
        if (&y == this) {
            return *this;
        }
        m_nrows = y.m_nrows;
        m_ncols = y.m_ncols;
        m_data.resize(m_nrows*m_ncols);
        m_data = y.m_data;
        return *this;
    }

    //! Resize the array, and fill the new entries with 'v'
    /*!
     * @param n  This is the number of rows
     * @param m  This is the number of columns in the new matrix
     * @param v  Default fill value -> defaults to zero.
     */
    void resize(size_t n, size_t m, doublereal v = 0.0) {
        m_nrows = n;
        m_ncols = m;
        m_data.resize(n*m, v);
    }

    //! Copy the data from one array into another without doing any checking
    /*!
     *  This differs from the assignment operator as no resizing is done and memcpy() is used.
     *  @param y Array to be copied
     */
    void copyData(const Array2D& y) {
        size_t n = sizeof(doublereal) * m_nrows * m_ncols;
        (void) memcpy(DATA_PTR(m_data), y.ptrColumn(0), n);
    }

    //! Append a column to the existing matrix using a std vector
    /*!
     *  This operation will add a column onto the existing matrix.
     *
     *  @param c  This vector<doublereal> is the entries in the
     *            column to be added. It must have a length
     *            equal to m_nrows or greater.
     */
    void appendColumn(const vector_fp& c) {
        m_ncols++;
        m_data.resize(m_nrows*m_ncols);
        size_t m;
        for (m = 0;  m < m_nrows; m++) {
            value(m_ncols, m) = c[m];
        }
    }

    //! Append a column to the existing matrix
    /*!
     *  This operation will add a column onto the existing matrix.
     *
     *  @param c  This vector of doubles is the entries in the
     *            column to be added. It must have a length
     *            equal to m_nrows or greater.
     */
    void appendColumn(const doublereal* const c) {
        m_ncols++;
        m_data.resize(m_nrows*m_ncols);
        size_t m;
        for (m = 0;  m < m_nrows; m++) {
            value(m_ncols, m) = c[m];
        }
    }

    //! Set the nth row to array rw
    /*!
     *  @param  n  Index of the row to be changed
     *  @param  rw  Vector for the row. Must have a length of m_ncols.
     */
    void setRow(size_t n, const doublereal* const rw) {
        for (size_t j = 0; j < m_ncols; j++) {
            m_data[m_nrows*j + n] = rw[j];
        }
    }

    //! Get the nth row and return it in a vector
    /*!
     *   @param n    Index of the row to be returned.
     *   @param rw   Return Vector  for the operation.
     *               Must have a length of m_ncols.
     */
    void getRow(size_t n, doublereal* const rw) {
        for (size_t j = 0; j < m_ncols; j++) {
            rw[j] = m_data[m_nrows*j + n];
        }
    }

    //! Set the values in column m to those in array col
    /*!
     *  A(i,m) = col(i)
     *
     *  @param m    Column to set
     *  @param col  pointer to a col vector. Vector
     *              must have a length of m_nrows.
     */
    void setColumn(size_t m, doublereal* const col) {
        for (size_t i = 0; i < m_nrows; i++) {
            m_data[m_nrows*m + i] = col[i];
        }
    }

    //! Get the values in column m
    /*!
     * col(i) =  A(i,m)
     *
     *  @param m    Column to set
     *  @param col  pointer to a col vector that will be returned
     */
    void getColumn(size_t m, doublereal* const col) {
        for (size_t i = 0; i < m_nrows; i++) {
            col[i] = m_data[m_nrows*m + i];
        }
    }

    /**
     * Destructor. Does nothing, since no memory allocated on the
     * heap.
     */
    virtual ~Array2D() {}

    //! Evaluate z = a*x + y.
    /*!
     *  This function evaluates the AXPY operation, and stores
     *  the result in the object's Array2D object.
     *  It's assumed that all 3 objects have the same dimensions,
     *  but no error checking is done.
     *
     *  @param a  scalar to multiply x with
     *  @param x  First Array2D object to be used
     *  @param y  Second Array2D object to be used
     *
     */
    void axpy(doublereal a, const Array2D& x, const Array2D& y) {
        iterator b = begin();
        const_iterator xb = x.begin();
        const_iterator yb = y.begin();
        for (; b != end(); ++b, ++xb, ++yb) {
            *b = a*(*xb) + *yb;
        }
    }

    //! Set all of the entries to zero
    inline void zero() {
        size_t nn = m_nrows * m_ncols;
        if (nn > 0) {
            /*
             * Using memset is the fastest way to zero a contiguous
             * section of memory.
             */
            (void) memset((void*) &m_data[0], 0, nn * sizeof(doublereal));
        }
    }

    //! Allows setting elements using the syntax A(i,j) = x.
    /*!
     *  @param  i            row index
     *  @param  j            column index.
     *
     *  @return Returns a reference to A(i,j) which may be assigned.
     */
    doublereal& operator()(size_t i, size_t j) {
        return value(i,j);
    }


    //! Allows retrieving elements using the syntax x = A(i,j).
    /*!
     *   @param i   Index for the row to be retrieved
     *   @param j   Index for the column to be retrieved.
     *
     *   @return    Returns the value of the matrix entry
     */
    doublereal operator()(size_t i, size_t j) const {
        return value(i,j);
    }

    //! Returns a changeable reference to position in the matrix
    /*!
     * This is a key entry. Returns a reference to the matrixes (i,j)
     * element. This may be used as an L value.
     *
     * @param i   The row index
     * @param j   The column index
     *
     * @return  Returns a changeable reference to the matrix entry
     */
    doublereal& value(size_t i, size_t j) {
        return m_data[m_nrows*j + i];
    }

    //! Returns the value of a single matrix entry
    /*!
     * This is a key entry. Returns the value of the  matrix position (i,j)
     * element.
     *
     * @param i   The row index
     * @param j   The column index
     */
    doublereal value(size_t i, size_t j) const {
        return m_data[m_nrows*j + i];
    }

    /// Number of rows
    size_t nRows() const {
        return m_nrows;
    }

    /// Number of columns
    size_t nColumns() const {
        return m_ncols;
    }

    /// Return an iterator pointing to the first element
    iterator begin() {
        return m_data.begin();
    }

    /// Return an iterator pointing past the last element
    iterator end() {
        return m_data.end();
    }

    /// Return a const iterator pointing to the first element
    const_iterator begin() const {
        return m_data.begin();
    }

    /// Return a const iterator pointing to past the last element
    const_iterator end() const {
        return m_data.end();
    }

    /// Return a reference to the data vector
    vector_fp& data() {
        return m_data;
    }

    /// Return a const reference to the data vector
    const vector_fp& data() const {
        return m_data;
    }

    //! Return a pointer to the top of column j, columns are contiguous
    //! in memory
    /*!
     *  @param j   Value of the column
     *
     *  @return  Returns a pointer to the top of the column
     */
    doublereal* ptrColumn(size_t j) {
        return &(m_data[m_nrows*j]);
    }

    //! Return a const pointer to the top of column j, columns are contiguous
    //! in memory
    /*!
     *  @param j   Value of the column
     *
     *  @return  Returns a const pointer to the top of the column
     */
    const doublereal* ptrColumn(size_t j) const {
        return &(m_data[m_nrows*j]);
    }

protected:

    //! Data stored in a single array
    vector_fp m_data;

    //! Number of rows
    size_t m_nrows;

    //! Number of columns
    size_t m_ncols;
};

//! Output the current contents of the Array2D object
/*!
 *  Example of usage:
 *        s << m << endl;
 *
 *  @param s   Reference to the ostream to write to
 *  @param m   Object of type Array2D that you are querying
 *
 *  @return    Returns a reference to the ostream.
 */
inline std::ostream& operator<<(std::ostream& s, const Array2D& m)
{
    size_t nr = m.nRows();
    size_t nc = m.nColumns();
    size_t i,j;
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
inline void operator*=(Array2D& m, doublereal a)
{
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
inline void operator+=(Array2D& x, const Array2D& y)
{
    sum_each(x.begin(), x.end(), y.begin());
}

}

#endif
