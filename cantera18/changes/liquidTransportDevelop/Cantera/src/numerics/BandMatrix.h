/**
 *  @file BandMatrix.h
 *
 *  Banded matrices.
 */

/*
 *  $Author$
 *  $Revision$
 *  $Date$
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_BANDMATRIX_H
#define CT_BANDMATRIX_H

#include "ct_defs.h"
#include "ctlapack.h"
#include "utilities.h"
#include "ctexceptions.h"

namespace Cantera {

  /**
   *  A class for banded matrices. This class has matrix inversion processes.
   *  The class is based upon the LAPACK banded storage matrix format.
   */
  class BandMatrix {

  public:

    //! Base Constructor
    /*!
     * * Create an \c 0 by \c 0 matrix, and initialize all elements to \c 0.
     */
    BandMatrix();

    //! Creates a banded matrix and sets all elements to zero
    /*!
     *   Create an \c n by \c n  banded matrix, and initialize all elements to \c v.
     *
     * @param n   size of the square matrix
     * @param kl  band size on the lower portion of the matrix
     * @param ku  band size on the upper portion of the matrix
     * @param v   intial value of all matrix components.
     */
    BandMatrix(int n, int kl, int ku, doublereal v = 0.0);

    //! Copy constructor
    /*!
     *  @param y  Matrix to be copied
     */
    BandMatrix(const BandMatrix& y);

    //! Destructor. Does nothing.
    virtual ~BandMatrix();

    //! assignment operator
    /*!
     * @param y  reference to the matrix to be copied
     */
    BandMatrix& operator=(const BandMatrix& y);

    //! Resize the matrix problem
    /*!
     *  All data is lost
     *
     * @param  n   size of the square matrix
     * @param kl  band size on the lower portion of the matrix
     * @param ku  band size on the upper portion of the matrix
     * @param v   intial value of all matrix components.
     */
    void resize(int n, int kl, int ku, doublereal v = 0.0);

    //! Fill or zero the matrix
    /*!
     *  @param v  Fill value, defaults to zero.
     */
    void bfill(doublereal v = 0.0);

    //! Index into the (i,j) element
    /*!
     *  @param i  row
     *  @param j  column
     *
     *  Returns a changeable reference to the matrix entry
     */
    doublereal& operator()(int i, int j);


    //! Constant Index into the (i,j) element
    /*!
     *  @param i  row
     *  @param j  column
     *
     *  Returns an unchangeable reference to the matrix entry
     */
    doublereal operator() (int i, int j) const;

    //! Return a changeable reference to element (i,j).
    /*!
     *  Since this method may alter the element value, it may need to be refactored, so
     *  the flag m_factored is set to false.
     *
     *  @param i  row
     *  @param j  column
     *
     *  @return Returns a reference to the value of the matrix entry
     */
    doublereal& value( int i, int j);


    //! Return the value of element (i,j). 
    /*!
     *   This method does not  alter the array.
     *  @param i  row
     *  @param j  column
     *
     *  @return Returns the value of the matrix entry
     */
    doublereal value( int i, int j) const;

    //! Returns the location in the internal 1D array corresponding to the (i,j) element in the banded array
    /*!
     *  @param i  row
     *  @param j  column
     *
     * @return  Returns the index of the matrix entry
     */
    int index(int i, int j) const;

    //! Return the value of the (i,j) element for (i,j) within the bandwidth.
    /*!
     *  For efficiency, this method does not check that (i,j) are within the bandwidth; it is up to the calling
     *  program to insure that this is true.
     *
     *   @param i  row
     *   @param j  column
     * 
     *   @return   Returns the value of the matrix entry
     */
    doublereal _value(int i, int j) const;

    //! Returns the number of rows
    int nRows() const;

    //! Number of columns
    int nColumns() const;

    //! Number of subdiagonals
    int nSubDiagonals() const;

    //! Number of superdiagonals
    int nSuperDiagonals() const;

    //! Return the number of rows of storage needed for the band storage
    int ldim() const;

    //! Return a reference to the pivot vector
    /*!
     *  @return return a reference to the pivot vector
     */
    vector_int& ipiv();

    //! Multiply A*b and write result to prod.
    /*!
     *  @param b    Vector to do the rh multiplcation
     *  @param prod OUTPUT vector to receive the result 
     */
    void mult(const doublereal * const b, doublereal * const prod) const;

    //! Multiply b*A and write result to prod.
    /*!
     *  @param b    Vector to do the lh multiplcation
     *  @param prod OUTPUT vector to receive the result 
     */
    void leftMult(const doublereal * const b, doublereal * const prod) const;

    //! Perform an LU decomposition, the LAPACK routine DGBTRF is used.
    /*!
     *
     * The factorization is saved in ludata.
     *
     * @return Return a success flag. 
     *         0 indicates a success
     *         ~0  Some error occurred, see the LAPACK documentation
     */
    int factor();

 
    //! Solve the matrix problem Ax = b
    /*!
     *  @param n  size of the matrix
     *  @param b  INPUT rhs of the problem
     *  @param x  OUTPUT solution to the problem
     *
     * @return Return a success flag
     *          0 indicates a success
     *         ~0  Some error occurred, see the LAPACK documentation
     */
    int solve(int n, const doublereal * const b, doublereal * const x);

    //! Solve the matrix problem Ax = b
    /*!
     *  @param n  size of the matrix
     *  @param b  INPUT rhs of the problem
     *            OUTPUT solution to the problem
     *
     * @return Return a success flag
     *          0 indicates a success
     *         ~0  Some error occurred, see the LAPACK documentation
     */
    int solve(int n, doublereal * const b);


    //! Returns an iterator for the start of the band storage data
    /*!
     *  Iterator points to the beginning of the data, and it is changeable.
     */
    vector_fp::iterator begin();

    //! Returns an iterator for the end of the band storage data
    /*!
     *  Iterator points to the end of the data, and it is changeable.
     */
    vector_fp::iterator end();

    //! Returns a const iterator for the start of the band storage data
    /*!
     *  Iterator points to the beginning of the data, and it is not changeable.
     */
    vector_fp::const_iterator begin() const;

    //! Returns a const iterator for the end of the band storage data
    /*!
     *  Iterator points to the end of the data, and it is not changeable.
     */
    vector_fp::const_iterator end() const;

  protected:

    //! Matrix data
    vector_fp data;

    //! Factorized data
    vector_fp ludata;

    //! Boolean indicating whether the matrix is factored
    bool m_factored;

    //! Number of rows and columns of the matrix
    int m_n;

    //! Number of subdiagonals of the matrix
    int m_kl;

    //! Number of super diagonals of the matrix
    int m_ku;

    //! value of zero
    doublereal m_zero;

    //! Pivot vector
    vector_int  m_ipiv;

  };

  //! Utility routine to print out the matrix
  /*!
   *  @param s  ostream to print the matrix out to
   *  @param m  Matrix to be printed
   *
   *  @return Returns a reference to the ostream
   */
  std::ostream& operator<<(std::ostream& s, const BandMatrix& m);

}

#endif
