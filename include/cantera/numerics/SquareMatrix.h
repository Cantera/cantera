/**
 *  @file SquareMatrix.h
 * Dense, Square (not sparse) matrices.
 */

/*
 * Copyright 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#ifndef CT_SQUAREMATRIX_H
#define CT_SQUAREMATRIX_H

#include "DenseMatrix.h"
#include "GeneralMatrix.h"

namespace Cantera
{

/**
 *  A class for full (non-sparse) matrices with Fortran-compatible
 *  data storage. Adds matrix inversion operations to this class from DenseMatrix.
 */
class SquareMatrix: public DenseMatrix, public GeneralMatrix
{
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
     * @param v   initial value of all matrix components.
     */
    SquareMatrix(size_t n, doublereal v = 0.0);

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
    int solve(doublereal* b);

    //! Resize the matrix
    /*!
     *  @param n Number of rows
     *  @param m Number of columns
     *  @param v double to fill the new space (defaults to zero)
     */
    void resize(size_t n, size_t m, doublereal v = 0.0);

    /**
     * Zero the matrix
     */
    void zero();

    //! Multiply A*b and write result to prod.
    /*!
     *  @param b    Vector to do the rh multiplication
     *  @param prod OUTPUT vector to receive the result
     */
    virtual void mult(const doublereal* b, doublereal* prod) const;

    //! Multiply A*B and write result to \c prod.
    /*!
     *
     *  @param    b     input      DenseMatrix B of size NxN
     *  @param    prod  output     output DenseMatrix prod size NxN
     */
    virtual void mult(const DenseMatrix& b, DenseMatrix& prod) const;

    //! Multiply b*A and write result to prod.
    /*!
     *  @param b    Vector to do the lh multiplication
     *  @param prod OUTPUT vector to receive the result
     */
    virtual void leftMult(const doublereal* const b, doublereal* const prod) const;

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
    virtual int factorQR();

    //! Returns an estimate of the inverse of the condition number for the matrix
    /*!
     *   The matrix must have been previously factored using the QR algorithm
     *
     * @return  returns the inverse of the condition number
     */
    virtual doublereal rcondQR();

    //! Returns an estimate of the inverse of the condition number for the matrix
    /*!
     *   The matrix must have been previously factored using the LU algorithm
     *
     * @param a1norm Norm of the matrix
     *
     * @return  returns the inverse of the condition number
     */
    virtual doublereal rcond(doublereal a1norm);

    //! Returns the one norm of the matrix
    virtual doublereal oneNorm() const;

    //! Solves the linear problem Ax=b using the QR algorithm returning x in the b spot
    /*!
     *  @param b  RHS to be solved.
     */
    int solveQR(doublereal* b);


    //! clear the factored flag
    virtual void clearFactorFlag();

    //! set the factored flag
    void setFactorFlag();

    //! Report whether the current matrix has been factored.
    virtual bool factored() const;

    //! Change the way the matrix is factored
    /*!
     *  @param fAlgorithm   integer
     *                   0 LU factorization
     *                   1 QR factorization
     */
    virtual void useFactorAlgorithm(int fAlgorithm);

    //! Returns the factor algorithm used
    /*!
     *     0 LU decomposition
     *     1 QR decomposition
     *
     * This routine will always return 0
     */
    virtual int factorAlgorithm() const;

    //! Return a pointer to the top of column j, columns are assumed to be contiguous in memory
    /*!
     *  @param j   Value of the column
     *
     *  @return  Returns a pointer to the top of the column
     */
    virtual doublereal* ptrColumn(size_t j);

    //! Index into the (i,j) element
    /*!
     *  @param i  row
     *  @param j  column
     *
     * (note, tried a using directive here, and it didn't seem to work)
     *
     *  Returns a changeable reference to the matrix entry
     */
    virtual doublereal& operator()(size_t i, size_t j) {
        return Array2D::operator()(i, j);
    }

    //! Copy the data from one array into another without doing any checking
    /*!
     *  This differs from the assignment operator as no resizing is done and memcpy() is used.
     *  @param y Array to be copied
     */
    virtual void copyData(const GeneralMatrix& y);

    //! Constant Index into the (i,j) element
    /*!
     *  @param i  row
     *  @param j  column
     *
     *  Returns an unchangeable reference to the matrix entry
     */
    virtual doublereal operator()(size_t i, size_t j) const {
        return Array2D::operator()(i, j);
    }

    //! Return the number of rows in the matrix
    virtual size_t nRows() const;

    //! Return the size and structure of the matrix
    /*!
     * This is inherited from GeneralMatrix
     *
     * @param iStruct OUTPUT Pointer to a vector of ints that describe the structure of the matrix.
     *    not used
     *
     * @return  returns the number of rows and columns in the matrix.
     */
    size_t nRowsAndStruct(size_t* const iStruct = 0) const;

    //! Duplicate this object
    virtual GeneralMatrix* duplMyselfAsGeneralMatrix() const;


    //! Return an iterator pointing to the first element
    /*!
     */
    virtual  vector_fp::iterator begin();


    //! Return a const iterator pointing to the first element
    virtual vector_fp::const_iterator begin() const;


    //! Return a vector of const pointers to the columns
    /*!
     *  Note the value of the pointers are protected by their being const.
     *  However, the value of the matrix is open to being changed.
     *
     *   @return returns a vector of pointers to the top of the columns
     *           of the matrices.
     */
    virtual doublereal* const* colPts();

    //! Check to see if we have any zero rows in the jacobian
    /*!
     *  This utility routine checks to see if any rows are zero.
     *  The smallest row is returned along with the largest coefficient in that row
     *
     * @param valueSmall  OUTPUT value of the largest coefficient in the smallest row
     *
     * @return index of the row that is most nearly zero
     */
    virtual size_t checkRows(doublereal& valueSmall) const;

    //! Check to see if we have any zero columns in the jacobian
    /*!
     *  This utility routine checks to see if any columns are zero.
     *  The smallest column is returned along with the largest coefficient in that column
     *
     * @param valueSmall  OUTPUT value of the largest coefficient in the smallest column
     *
     * @return index of the column that is most nearly zero
     */
    virtual size_t checkColumns(doublereal& valueSmall) const;

protected:

    //!  the factor flag
    int m_factored;

public:
    //! Work vector for QR algorithm
    vector_fp tau;

    //! Work vector for QR algorithm
    vector_fp work;

    //! Integer work vector for QR algorithms
    std::vector<int> iwork_;
protected:
    //! 1-norm of the matrix. This is determined immediately before every factorization
    doublereal a1norm_;

    //!  Use the QR algorithm to factor and invert the matrix
    int useQR_;
};
}

#endif

