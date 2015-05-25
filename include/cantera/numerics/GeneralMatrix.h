/**
 *  @file GeneralMatrix.h
 *   Declarations for the class GeneralMatrix which is a virtual base class for matrices handled by solvers
 *    (see class \ref numerics and \link Cantera::GeneralMatrix GeneralMatrix\endlink).
 */

/*
 * Copyright 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#ifndef CT_GENERALMATRIX_H
#define CT_GENERALMATRIX_H

#include "cantera/base/ct_defs.h"
#include "cantera/base/ctexceptions.h"

namespace Cantera
{

//! Generic matrix
class GeneralMatrix
{
public:
    //! Base Constructor
    /*!
     *  @param matType  Matrix type
     *       0 full
     *       1 banded
     */
    GeneralMatrix(int matType);

    //! Copy Constructor
    GeneralMatrix(const GeneralMatrix& right);

    //! Assignment operator
    GeneralMatrix& operator=(const GeneralMatrix& right);

    //! Destructor. Does nothing.
    virtual ~GeneralMatrix() {}

    //! Duplicator member function
    /*!
     *  This function will duplicate the matrix given a generic GeneralMatrix
     *
     *  @return Returns a pointer to the malloced object
     */
    virtual GeneralMatrix* duplMyselfAsGeneralMatrix() const = 0;

    //! Zero the matrix elements
    virtual void zero() = 0;

    //! Multiply A*b and write result to prod.
    /*!
     *  @param b    Vector to do the rh multiplication
     *  @param prod OUTPUT vector to receive the result
     */
    virtual void mult(const doublereal* b, doublereal* prod) const = 0;

    //! Multiply b*A and write result to prod.
    /*!
     *  @param b    Vector to do the lh multiplication
     *  @param prod OUTPUT vector to receive the result
     */
    virtual void leftMult(const doublereal* const b, doublereal* const prod) const = 0;

    //! Factors the A matrix, overwriting A.
    /*
     *   We flip m_factored  boolean to indicate that the matrix is now A-1.
     */
    virtual int factor() = 0;

    //! Factors the A matrix using the QR algorithm, overwriting A
    /*!
     * we set m_factored to 2 to indicate the matrix is now QR factored
     *
     * @return  Returns the info variable from LAPACK
     */
    virtual int factorQR() {
        throw NotImplementedError("GeneralMatrix::factorQR");
    }

    //! Returns an estimate of the inverse of the condition number for the matrix
    /*!
     *   The matrix must have been previously factored using the QR algorithm
     *
     * @return  returns the inverse of the condition number
     */
    virtual doublereal rcondQR() {
        throw NotImplementedError("GeneralMatrix::rcondQR");
    }

    //! Returns an estimate of the inverse of the condition number for the matrix
    /*!
     *   The matrix must have been previously factored using the LU algorithm
     *
     * @param a1norm Norm of the matrix
     *
     * @return  returns the inverse of the condition number
     */
    virtual doublereal rcond(doublereal a1norm) = 0;

    //! Change the way the matrix is factored
    /*!
     *  @param fAlgorithm   integer
     *                   0 LU factorization
     *                   1 QR factorization
     */
    virtual void useFactorAlgorithm(int fAlgorithm) {
        throw NotImplementedError("GeneralMatrix::useFactorAlgorithm");
    };

    //! Return the factor algorithm used
    virtual int factorAlgorithm() const = 0;

    //! Calculate the one norm of the matrix
    virtual doublereal oneNorm() const = 0;

    //! Return the number of rows in the matrix
    virtual size_t nRows() const = 0;

    //! Return the size and structure of the matrix
    /*!
     * @param iStruct OUTPUT Pointer to a vector of ints that describe the structure of the matrix.
     *
     * @return  returns the number of rows and columns in the matrix.
     */
    virtual size_t nRowsAndStruct(size_t* const iStruct = 0) const = 0;

    //! clear the factored flag
    virtual void clearFactorFlag() {
        m_factored = 0;
    };

    //! Solves the Ax = b system returning x in the b spot.
    /*!
     *  @param b    Vector for the RHS of the equation system
     *  @param nrhs Number of right-hand sides to solve, default 1
     *  @param ldb  Leading dimension of the right-hand side array.
     *              Defaults to nRows()
     */
    virtual int solve(doublereal* b, size_t nrhs=1, size_t ldb=0) = 0;

    //! true if the current factorization is up to date with the matrix
    virtual bool factored() const {
        return (m_factored != 0);
    }

    //! Return a pointer to the top of column j, columns are assumed to be contiguous in memory
    /*!
     *  @param j   Value of the column
     *
     *  @return  Returns a pointer to the top of the column
     */
    virtual doublereal* ptrColumn(size_t j) = 0;

    //! Index into the (i,j) element
    /*!
     *  @param i  row
     *  @param j  column
     *
     *  Returns a changeable reference to the matrix entry
     */
    virtual doublereal& operator()(size_t i, size_t j) = 0;

    //! Constant Index into the (i,j) element
    /*!
     *  @param i  row
     *  @param j  column
     *
     *  Returns an unchangeable reference to the matrix entry
     */
    virtual doublereal operator()(size_t i, size_t j) const = 0;

    //! Copy the data from one array into another without doing any checking
    /*!
     *  This differs from the assignment operator as no resizing is done and memcpy() is used.
     *  @param y Array to be copied
     *  @deprecated To be removed after Cantera 2.2.
     */
    virtual void copyData(const GeneralMatrix& y) = 0;

    //! Return an iterator pointing to the first element
    /*!
     *  We might drop this later
     */
    virtual vector_fp::iterator begin() = 0;

    //! Return a const iterator pointing to the first element
    /*!
     *  We might drop this later
     */
    virtual vector_fp::const_iterator begin() const = 0;

    //! Return a vector of const pointers to the columns
    /*!
     *  Note the value of the pointers are protected by their being const.
     *  However, the value of the matrix is open to being changed.
     *
     *   @return returns a vector of pointers to the top of the columns
     *           of the matrices.
     */
    virtual doublereal* const* colPts() = 0;

    //! Check to see if we have any zero rows in the Jacobian
    /*!
     *  This utility routine checks to see if any rows are zero.
     *  The smallest row is returned along with the largest coefficient in that row
     *
     * @param valueSmall  OUTPUT value of the largest coefficient in the smallest row
     *
     * @return index of the row that is most nearly zero
     */
    virtual size_t checkRows(doublereal& valueSmall) const = 0;

    //! Check to see if we have any zero columns in the Jacobian
    /*!
     *  This utility routine checks to see if any columns are zero.
     *  The smallest column is returned along with the largest coefficient in that column
     *
     * @param valueSmall  OUTPUT value of the largest coefficient in the smallest column
     *
     * @return index of the column that is most nearly zero
     */
    virtual size_t checkColumns(doublereal& valueSmall) const = 0;

    //! Matrix type
    /*!
     *      0 Square
     *      1 Banded
     */
    int matrixType_;

protected:
    //! Indicates whether the matrix is factored. 0 for unfactored; Non-zero values
    //! indicate a particular factorization (LU=1, QR=2).
    int m_factored;
};

}
#endif
