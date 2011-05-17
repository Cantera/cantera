/**
 *  @file SquareMatrix.h
 * Dense, Square (not sparse) matrices.
 */

/*
 *  $Date$
 *  $Revision$
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#ifndef CT_SQUAREMATRIX_H
#define CT_SQUAREMATRIX_H

#include "DenseMatrix.h"

namespace Cantera { 

  /**
   *  A class for full (non-sparse) matrices with Fortran-compatible
   *  data storage. Adds matrix inversion operations to this class from DenseMatrix.
   */
  class SquareMatrix: public DenseMatrix {

  public:


    //! Base Constructor.
    /*!
     * Create an \c 0 by \c 0 matrix, and initialize all elements to \c 0.
     */
    SquareMatrix();

    //! Constructor.
    /*!
     * Create an \c n by \c n matrix, and initialize all elements to \c v.
     *
     * @param n   size of the square matrix
     * @param v   intial value of all matrix components.
     */
    SquareMatrix(int n, doublereal v = 0.0);

    //! Copy Constructor
    /*!
     *  @param right Object to be copied 
     */
    SquareMatrix(const SquareMatrix& right);
   
    //! Assignment operator
    /*!
     *  @param right  Object to be copied
     */
    SquareMatrix& operator=(const SquareMatrix& right);


    //! Destructor. Does nothing.
    virtual ~SquareMatrix();

    
     //! Solves the Ax = b system returning x in the b spot.
    /*!
     *  @param b  Vector for the rhs of the equation system
     */
    int solve(doublereal *b);

    //! Resize the matrix
    /*!
     *  @param n Number of rows
     *  @param m Number of columns
     *  @param v double to fill the new space (defaults to zero)
     */
    void resize(int n, int m, doublereal v = 0.0);


    /**
     * Zero the matrix
     */
    void zero();

    /**
     * Factors the A matrix, overwriting A. We flip m_factored
     * boolean to indicate that the matrix is now A-1.
     */
    int factor();

    //! Factors the A matrix using the QR algorithm, overwriting A
    /*!
     * we set m_factored to 2 to indicate the matrix is now QR factored
     *
     * @return  Returns the info variable from lapack
     */
    int factorQR();

    //! Returns an estimate of the inverse of the condition number for the matrix
    /*!
     *   The matrix must have been previously factored using the QR algorithm
     *
     * @return  returns the inverse of the condition number
     */
    doublereal rcondQR();

    //! Returns an estimate of the inverse of the condition number for the matrix
    /*!
     *   The matrix must have been previously factored using the LU algorithm
     *
     * @param a1norm Norm of the matrix
     *
     * @return  returns the inverse of the condition number
     */
    doublereal rcond(doublereal a1norm);

    //! Solves the linear problem Ax=b using the QR algorithm returning x in the b spot
    /*!
     *  @param b  RHS to be solved.
     */
    int solveQR(doublereal *b); 

    
    //! clear the factored flag
    void clearFactorFlag();
    /**
     * set the factored flag
     */
    void setFactorFlag();

    /*
     * the factor flag
     */
    int m_factored;

    //! Work vector for QR algorithm
    vector_fp tau;

    //! Work vector for QR algorithm
    vector_fp work;

    //! Integer work vector for QR algorithms
    std::vector<int> iwork_;
    
    //! 1-norm of the matrix. This is determined immediately before every factorization
    doublereal a1norm_;

  public:
    //!  Use the QR algorithm to factor and invert the matrix
    int useQR_;
  };
}

#endif



