/**
 *  @file DenseMatrix.h
 *  Headers for the %DenseMatrix object, which deals with dense rectangular matrices and
 *  description of the numerics groupings of objects
 *  (see \ref numerics and \link Cantera::DenseMatrix DenseMatrix \endlink) .
 */

/*
 *  $Author$
 *  $Date$
 *  $Revision$
 *
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_DENSEMATRIX_H
#define CT_DENSEMATRIX_H

#include "ct_defs.h"
#include "Array.h"

namespace Cantera { 
  /**
   * @defgroup numerics  Numerical Utilities within Cantera
   *
   *  Cantera contains some capabilities for solving nonlinear equations and
   *  integrating both ODE and DAE equation systems in time. This section describes these
   *  capabilities.
   *
   */

 
  //! A class for full (non-sparse) matrices with Fortran-compatible
  //!  data storage, which adds matrix operations to class Array2D.
  /*!
   *  The dense matrix class adds matrix operations onto the Array2D class.
   *  These matrix operations are carried out by the appropriate BLAS and LAPACK routines
   *
   *  @ingroup numerics
   */
  class DenseMatrix : public Array2D {

  public:

    //! Default Constructor
    DenseMatrix();
 
    //! Constructor.
    /*!
     *  Create an \c n by \c m matrix, and initialize all elements to \c v.
     * 
     *  @param n  New number of rows
     *  @param m  New number of columns
     *  @param v  Default fill value. defaults to zero.
     */
    DenseMatrix(int n, int m, doublereal v = 0.0);

    //! Copy constructor
    /*!
     *   @param y Object to be copied
     */
    DenseMatrix(const DenseMatrix& y);

    //! Assignment operator
    /*!
     *   @param y Object to be copied
     */
    DenseMatrix& operator=(const DenseMatrix& y);

    //! Destructor. Does nothing.
    virtual ~DenseMatrix();

    //! Resize the matrix 
    /*!
     *  Resize the matrix to n rows by m cols.
     *  
     *  @param n  New number of rows
     *  @param m  New number of columns
     *  @param v  Default fill value. defaults to zero.
     */
    void resize(int n, int m, doublereal v = 0.0);

    //! Multiply A*b and write result to \c prod.
    /*!
     * 
     *  @param    b     input      vector b with length N
     *  @param    prod  output     output vector prod length = M
     */
    virtual void mult(const double* b, double* prod) const;

    //! Left-multiply the matrix by transpose(b), and write the result to prod.
    /*!
     *   @param b    left multiply by this vector. The length must be equal to n
     *               the number of rows in the matrix.
     *
     *   @param prod  Resulting vector. This is of length m, the number of columns
     *                in the matrix
     */
    virtual void leftMult(const double * const b, double* const prod) const;

    //! Return a changeable value of the pivot vector
    /*!
     *  @return  Returns a reference to the pivot vector as a vector_int
     */
    vector_int& ipiv();

    //! Return a changeable value of the pivot vector
    /*!
     *  @return  Returns a reference to the pivot vector as a vector_int
     */
    const vector_int& ipiv() const { return m_ipiv; }

  protected:

    //! Vector of pivots. Length is equal to the max of m and n.
    vector_int     m_ipiv;
  };

  //==================================================================================================================

  
  //! Solve Ax = b. Array b is overwritten on exit with x.
  /*!
   *   The solve class uses the LAPACK routine dgetrf to invert the m xy n matrix.
   *
   *   The factorization has the form
   *     A = P * L * U
   *  where P is a permutation matrix, L is lower triangular with unit
   *  diagonal elements (lower trapezoidal if m > n), and U is upper
   *  triangular (upper trapezoidal if m < n).
   *
   *  The system is then solved using the LAPACK routine dgetrs
   *
   *   @param A   Dense matrix to be factored
   *   @param b   rhs to be solved.
   */
  int solve(DenseMatrix& A, double* b);

  //!  Solve Ax = b for multiple right-hand-side vectors. 
  /*!
   *  @param A    Dense matrix to be factored
   *  @param b   Dense matrix of rhs's. Each column is a rhs
   */
  int solve(DenseMatrix& A, DenseMatrix& b);

#ifdef INCL_LEAST_SQUARES
  //! Solve Ax = b in the least squares sense
  /*!
   *  @param A   Matrix to be inverted in the least squares sense
   *  @param b   Vector b to be solved for
   *   @todo fix lwork
   */
  int leastSquares(DenseMatrix& A, double* b);
#endif
  
  //! Multiply \c A*b and return the result in \c prod. Uses BLAS routine DGEMV.
  /*!
   *  \f[
   *          prod_i = sum^N_{j = 1}{A_{ij} b_j}
   *  \f]
   * 
   *  @param    A     input      Dense Matrix A with M rows and N columns
   *  @param    b     input      vector b with length N
   *  @param    prod  output     output vector prod length = M
   */ 
  void multiply(const DenseMatrix& A, const double * const b, double * const prod);

  //! Multiply \c A*b and add it to the result in \c prod. Uses BLAS routine DGEMV.
  /*!
   *  \f[
   *          prod_i += sum^N_{j = 1}{A_{ij} b_j}
   *  \f]
   * 
   *  @param    A     input      Dense Matrix A with M rows and N columns
   *  @param    b     input      vector b with length N
   *  @param    prod  output     output vector prod length = M
   */ 
  void increment(const DenseMatrix& A, const double * const b, double * const prod);

  //! invert A. A is overwritten with A^-1.
  /*!
   *  @param A  Invert the matrix A and store it back in place
   *
   *  @param nn  Size of A. This defaults to -1, which means that the number
   *                        of rows is used as the default size of n
   */
  int invert(DenseMatrix& A, int nn=-1);

}

#endif



