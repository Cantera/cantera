/**
 *  @file Array.h Header file for class Cantera::Array2D
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_ARRAY_H
#define CT_ARRAY_H

#include "ct_defs.h"
#include <iostream>

namespace Cantera
{

//!  A class for 2D arrays stored in column-major (Fortran-compatible) form.
/*!
 * In this form, the data entry for an n row, m col matrix is
 *
 *       index = i + (n-1) * j
 *
 * where
 *
 *     J(i,j) = data_start + index
 *         i = row
 *         j = column
 */
class Array2D
{
public:
    //! Type definition for the iterator class that is can be used by Array2D
    //! types.
    /*!
     * This is just equal to vector_fp::iterator.
     */
    typedef vector_fp::iterator iterator;

    //! Type definition for the const_iterator class that is can be used by
    //! Array2D types.
    /*!
     * This is just equal to vector_fp::const_iterator.
     */
    typedef vector_fp::const_iterator const_iterator;

    /**
     * Default constructor. Create an empty array.
     */
    Array2D();

    //! Constructor.
    /*!
     * Create an \c m by \c n array, and initialize all elements to \c v.
     *
     * @param m   Number of rows
     * @param n   Number of columns
     * @param v   Default fill value. The default is 0.0
     */
    Array2D(const size_t m, const size_t n, const double v=0.0);

    //! Constructor.
    /*!
     * Create an \c m by \c n array, initialized with the contents of the array
     * \c values.
     *
     * @param m   Number of rows
     * @param n   Number of columns
     * @param values Initial values of the array. Must be of length m*n, and
     *     stored in column-major order.
     */
    Array2D(const size_t m, const size_t n, const double* values);

    Array2D(const Array2D& y);

    virtual ~Array2D() {}

    Array2D& operator=(const Array2D& y);

    //! Resize the array, and fill the new entries with 'v'
    /*!
     * @param n  This is the number of rows
     * @param m  This is the number of columns in the new matrix
     * @param v  Default fill value -> defaults to zero.
     */
    void resize(size_t n, size_t m, double v=0.0);

    //! Append a column to the existing matrix using a std vector
    /*!
     * This operation will add a column onto the existing matrix.
     *
     * @param c  This vector is the entries in the column to be added. It must
     *           have a length equal to m_nrows or greater.
     */
    void appendColumn(const vector_fp& c);

    //! Append a column to the existing matrix
    /*!
     * This operation will add a column onto the existing matrix.
     *
     * @param c  This vector of doubles is the entries in the column to be
     *           added. It must have a length equal to m_nrows or greater.
     */
    void appendColumn(const double* const c);

    //! Set the nth row to array rw
    /*!
     * @param  n   Index of the row to be changed
     * @param  rw  Vector for the row. Must have a length of m_ncols.
     */
    void setRow(size_t n, const double* const rw);

    //! Get the nth row and return it in a vector
    /*!
     * @param n    Index of the row to be returned.
     * @param rw   Return Vector for the operation. Must have a length of
     *             m_ncols.
     */
    void getRow(size_t n, double* const rw);

    //! Set the values in column m to those in array col
    /*!
     * A(i,m) = col(i)
     *
     * @param m   Column to set
     * @param col pointer to a col vector. Vector must have a length of m_nrows.
     */
    void setColumn(size_t m, double* const col);

    //! Get the values in column m
    /*!
     * col(i) =  A(i,m)
     *
     * @param m    Column to set
     * @param col  pointer to a col vector that will be returned
     */
    void getColumn(size_t m, double* const col);

    //! Set all of the entries to zero
    void zero() {
        m_data.assign(m_data.size(), 0.0);
    }

    //! Allows setting elements using the syntax A(i,j) = x.
    /*!
     *  @param  i            row index
     *  @param  j            column index.
     *  @returns a reference to A(i,j) which may be assigned.
     */
    doublereal& operator()(size_t i, size_t j) {
        return value(i,j);
    }

    //! Allows retrieving elements using the syntax x = A(i,j).
    /*!
     * @param i   Index for the row to be retrieved
     * @param j   Index for the column to be retrieved.
     * @returns the value of the matrix entry
     */
    doublereal operator()(size_t i, size_t j) const {
        return value(i,j);
    }

    //! Returns a changeable reference to position in the matrix
    /*!
     * Returns a reference to the matrix's (i,j) element. This may be used as an
     * L value.
     *
     * @param i   The row index
     * @param j   The column index
     * @returns a changeable reference to the matrix entry
     */
    doublereal& value(size_t i, size_t j) {
        return m_data[m_nrows*j + i];
    }

    //! Returns the value of a single matrix entry
    /*!
     * Returns the value of the matrix position (i,j) element.
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
     * @param j   Value of the column
     * @returns a pointer to the top of the column
     */
    doublereal* ptrColumn(size_t j) {
        return &m_data[m_nrows*j];
    }

    //! Return a const pointer to the top of column j, columns are contiguous
    //! in memory
    /*!
     * @param j   Value of the column
     * @returns a const pointer to the top of the column
     */
    const doublereal* ptrColumn(size_t j) const {
        return &m_data[m_nrows*j];
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
 * Example of usage:
 *        s << m << endl;
 *
 * @param s   Reference to the ostream to write to
 * @param m   Object of type Array2D that you are querying
 * @returns a reference to the ostream.
 */
std::ostream& operator<<(std::ostream& s, const Array2D& m);

//! Overload the times equals operator for multiplication of a matrix and a
//! scalar.
/*!
 * Scaled every element of the matrix by the scalar input
 *
 * @param m   Matrix
 * @param a   scalar
 */
void operator*=(Array2D& m, double a);

}

#endif
