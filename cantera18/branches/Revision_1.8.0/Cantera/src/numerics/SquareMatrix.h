/**
 *  @file SquareMatrix.h
 * Dense, Square (not sparse) matrices.
 */

/*
 *  $Date: 2009/03/03 21:20:32 $
 *  $Revision: 1.2 $
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

    SquareMatrix():
      DenseMatrix(),
      m_factored(false)
    {
    }

    /** 
     * Constructor. Create an \c n by \c n matrix, and initialize
     * all elements to \c v.
     */
    SquareMatrix(int n, doublereal v = 0.0)  : 
      DenseMatrix(n, n, v),
      m_factored(false)
    {
    }

    /**
     * Copy Constructor
     */
    SquareMatrix(const SquareMatrix&);

    /**
     * Assignment operator
     */
    SquareMatrix& operator=(const SquareMatrix&);


    /// Destructor. Does nothing.
    virtual ~SquareMatrix(){}

    /**
     * Solves the Ax = b system returning x in the b spot.
     */
    int solve(double *b);

    /**
     * Zero the matrix
     */
    void zero();

    /**
     * Factors the A matrix, overwriting A. We flip m_factored
     * boolean to indicate that the matrix is now A-1.
     */
    int factor();
    /**
     * clear the factored flag
     */
    void clearFactorFlag();
    /**
     * set the factored flag
     */
    void setFactorFlag();

    /*
     * the factor flag
     */
    bool m_factored;
  };
}

#endif



