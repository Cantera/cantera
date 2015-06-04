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
    SquareMatrix(const SquareMatrix& right);

    //! Assignment operator
    SquareMatrix& operator=(const SquareMatrix& right);

    int solve(doublereal* b, size_t nrhs=1, size_t ldb=0);

    void resize(size_t n, size_t m, doublereal v = 0.0);

    //! Zero the matrix
    void zero();

    virtual void mult(const doublereal* b, doublereal* prod) const;
    virtual void mult(const DenseMatrix& b, DenseMatrix& prod) const;
    virtual void leftMult(const doublereal* const b, doublereal* const prod) const;

    int factor();
    virtual int factorQR();

    virtual doublereal rcondQR();
    virtual doublereal rcond(doublereal a1norm);

    virtual doublereal oneNorm() const;

    //! Solves the linear problem Ax=b using the QR algorithm returning x in the b spot
    /*!
     *  @param b  RHS to be solved.
     */
    int solveQR(doublereal* b);

    //! set the factored flag
    void setFactorFlag();

    virtual void useFactorAlgorithm(int fAlgorithm);

    //! Returns the factor algorithm used
    /*!
     *     0 LU decomposition
     *     1 QR decomposition
     *
     * This routine will always return 0
     */
    virtual int factorAlgorithm() const;

    virtual doublereal* ptrColumn(size_t j);

    virtual doublereal& operator()(size_t i, size_t j) {
        return Array2D::operator()(i, j);
    }

    //! @deprecated To be removed after Cantera 2.2.
    virtual void copyData(const GeneralMatrix& y);

    virtual doublereal operator()(size_t i, size_t j) const {
        return Array2D::operator()(i, j);
    }

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

    virtual GeneralMatrix* duplMyselfAsGeneralMatrix() const;

    virtual  vector_fp::iterator begin();
    virtual vector_fp::const_iterator begin() const;

    virtual doublereal* const* colPts();

    virtual size_t checkRows(doublereal& valueSmall) const;
    virtual size_t checkColumns(doublereal& valueSmall) const;

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
