/**
 *  @file DenseMatrix.h
 *  Headers for the DenseMatrix object, which deals with dense rectangular matrices and
 *  description of the numerics groupings of objects
 *  (see \ref numerics and \link Cantera::DenseMatrix DenseMatrix \endlink) .
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_DENSEMATRIX_H
#define CT_DENSEMATRIX_H

#include "cantera/base/ct_defs.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/Array.h"

namespace Cantera
{
/**
 * @defgroup numerics  Numerical Utilities within Cantera
 *
 * Cantera contains some capabilities for solving nonlinear equations and
 * integrating both ODE and DAE equation systems in time. This section describes
 * these capabilities.
 *
 */

//! A class for full (non-sparse) matrices with Fortran-compatible data storage,
//! which adds matrix operations to class Array2D.
/*!
 * The dense matrix class adds matrix operations onto the Array2D class. These
 * matrix operations are carried out by the appropriate BLAS and LAPACK routines
 *
 * Error handling from BLAS and LAPACK are handled via the following
 * formulation. Depending on a variable, a singular matrix or other terminal
 * error condition from LAPACK is handled by either throwing an exception or
 * by returning the error code condition to the calling routine.
 *
 * The int variable, m_useReturnErrorCode, determines which method is used. The
 * default value of zero means that an exception is thrown. A value of 1 means
 * that a return code is used.
 *
 * Reporting of these LAPACK error conditions is handled by the class variable
 * m_printLevel. The default is for no reporting. If m_printLevel is nonzero,
 * the error condition is reported to Cantera's log file.
 *
 * @ingroup numerics
 */
class DenseMatrix : public Array2D
{
public:
    //! Default Constructor
    DenseMatrix();

    //! Constructor.
    /*!
     * Create an \c n by \c m matrix, and initialize all elements to \c v.
     *
     * @param n  New number of rows
     * @param m  New number of columns
     * @param v  Default fill value. defaults to zero.
     */
    DenseMatrix(size_t n, size_t m, doublereal v = 0.0);

    DenseMatrix(const DenseMatrix& y);
    DenseMatrix& operator=(const DenseMatrix& y);

    //! Resize the matrix
    /*!
     * Resize the matrix to n rows by m cols.
     *
     * @param n  New number of rows
     * @param m  New number of columns
     * @param v  Default fill value. defaults to zero.
     */
    void resize(size_t n, size_t m, doublereal v = 0.0);

    virtual doublereal* const* colPts();

    //! Return a const vector of const pointers to the columns
    /*!
     * Note, the Jacobian can not be altered by this routine, and therefore the
     * member function is const.
     *
     * @returns a vector of pointers to the top of the columns of the matrices.
     */
    const doublereal* const* const_colPts() const;

    virtual void mult(const double* b, double* prod) const;

    //! Multiply A*B and write result to \c prod.
    /*!
     * Take this matrix to be of size NxM.
     * @param[in]  b     DenseMatrix B of size MxP
     * @param[out] prod  DenseMatrix prod size NxP
     */
    virtual void mult(const DenseMatrix& b, DenseMatrix& prod) const;

    //! Left-multiply the matrix by transpose(b), and write the result to prod.
    /*!
     * @param b    left multiply by this vector. The length must be equal to n
     *             the number of rows in the matrix.
     * @param prod  Resulting vector. This is of length m, the number of columns
     *              in the matrix
     */
    virtual void leftMult(const double* const b, double* const prod) const;

    //! Return a changeable value of the pivot vector
    /*!
     * @returns a reference to the pivot vector as a vector_int
     */
    vector_int& ipiv();

    //! Return a changeable value of the pivot vector
    /*!
     *  @returns a reference to the pivot vector as a vector_int
     */
    const vector_int& ipiv() const {
        return m_ipiv;
    }

protected:
    //! Vector of pivots. Length is equal to the max of m and n.
    vector_int m_ipiv;

    //! Vector of column pointers
    std::vector<doublereal*> m_colPts;

public:
    //! Error Handling Flag
    /*!
     * The default is to set this to 0. In this case, if a factorization is
     * requested and can't be achieved, a CESingularMatrix exception is
     * triggered. No return code is used, because an exception is thrown. If
     * this is set to 1, then an exception is not thrown. Routines return with
     * an error code, that is up to the calling routine to handle correctly.
     * Negative return codes always throw an exception.
     */
    int m_useReturnErrorCode;

    //! Print Level
    /*!
     * Printing is done to the log file using the routine writelogf().
     *
     * Level of printing that is carried out. Only error conditions are printed
     * out, if this value is nonzero.
     */
    int m_printLevel;

    // Listing of friend functions which are defined below

    friend int solve(DenseMatrix& A, double* b, size_t nrhs, size_t ldb);
    friend int solve(DenseMatrix& A, DenseMatrix& b);
    friend int invert(DenseMatrix& A, int nn);
};


//! Solve Ax = b. Array b is overwritten on exit with x.
/*!
 * The solve function uses the LAPACK routine dgetrf to invert the m xy n matrix.
 *
 * The factorization has the form
 *
 *     A = P * L * U
 *
 * where P is a permutation matrix, L is lower triangular with unit diagonal
 * elements (lower trapezoidal if m > n), and U is upper triangular (upper
 * trapezoidal if m < n).
 *
 * The system is then solved using the LAPACK routine dgetrs
 *
 * @param A    Dense matrix to be factored
 * @param b    RHS(s) to be solved.
 * @param nrhs Number of right hand sides to solve
 * @param ldb  Leading dimension of b, if nrhs > 1
 */
int solve(DenseMatrix& A, double* b, size_t nrhs=1, size_t ldb=0);

//! Solve Ax = b for multiple right-hand-side vectors.
/*!
 * @param A   Dense matrix to be factored
 * @param b   Dense matrix of RHS's. Each column is a RHS
 */
int solve(DenseMatrix& A, DenseMatrix& b);

//! Multiply \c A*b and return the result in \c prod. Uses BLAS routine DGEMV.
/*!
 * \f[
 *     prod_i = sum^N_{j = 1}{A_{ij} b_j}
 * \f]
 *
 * @param[in]  A     Dense Matrix A with M rows and N columns
 * @param[in]  b     vector b with length N
 * @param[out] prod  vector prod length = M
 */
void multiply(const DenseMatrix& A, const double* const b, double* const prod);

//! Multiply \c A*b and add it to the result in \c prod. Uses BLAS routine DGEMV.
/*!
 * \f[
 *     prod_i += sum^N_{j = 1}{A_{ij} b_j}
 * \f]
 *
 * @param[in]  A     Dense Matrix A with M rows and N columns
 * @param[in]  b     vector b with length N
 * @param[out] prod  vector prod length = M
 */
void increment(const DenseMatrix& A, const double* const b, double* const prod);

//! invert A. A is overwritten with A^-1.
/*!
 *  @param A   Invert the matrix A and store it back in place
 *  @param nn  Size of A. This defaults to -1, which means that the number of
 *             rows is used as the default size of n
 */
int invert(DenseMatrix& A, size_t nn=npos);

}

#endif
