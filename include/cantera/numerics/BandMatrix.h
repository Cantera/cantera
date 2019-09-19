/**
 *  @file BandMatrix.h
 *   Declarations for the class BandMatrix
 *   which is a child class of GeneralMatrix for banded matrices handled by solvers
 *    (see class \ref numerics and \link Cantera::BandMatrix BandMatrix\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_BANDMATRIX_H
#define CT_BANDMATRIX_H

#include "GeneralMatrix.h"

namespace Cantera
{

//! A class for banded matrices, involving matrix inversion processes.
//! The class is based upon the LAPACK banded storage matrix format.
/*!
 * An important issue with this class is that it stores both the original data
 * and the LU factorization of the data.  This means that the banded matrix
 * typically will take up twice the room that it is expected to take.
 *
 * QR factorizations of banded matrices are not included in the original LAPACK
 * work. Add-ons are available. However, they are not included here. Instead we
 * just use the stock LU decompositions.
 *
 * This class is a derived class of the base class GeneralMatrix. However,
 * within the oneD directory, the class is used as is, without reference to the
 * GeneralMatrix base type.
 */
class BandMatrix : public GeneralMatrix
{
public:
    //! Base Constructor
    /*!
     * Create an \c 0 by \c 0 matrix, and initialize all elements to \c 0.
     */
    BandMatrix();
    ~BandMatrix();

    //! Creates a banded matrix and sets all elements to zero
    /*!
     * Create an \c n by \c n  banded matrix, and initialize all elements to \c v.
     *
     * @param n   size of the square matrix
     * @param kl  band size on the lower portion of the matrix
     * @param ku  band size on the upper portion of the matrix
     * @param v   initial value of all matrix components.
     */
    BandMatrix(size_t n, size_t kl, size_t ku, doublereal v = 0.0);

    BandMatrix(const BandMatrix& y);
    BandMatrix& operator=(const BandMatrix& y);

    //! Resize the matrix problem
    /*!
     * All data is lost
     *
     * @param n   size of the square matrix
     * @param kl  band size on the lower portion of the matrix
     * @param ku  band size on the upper portion of the matrix
     * @param v   initial value of all matrix components.
     */
    void resize(size_t n, size_t kl, size_t ku, doublereal v = 0.0);

    //! Fill or zero the matrix
    /*!
     *  @param v  Fill value, defaults to zero.
     */
    void bfill(doublereal v = 0.0);

    doublereal& operator()(size_t i, size_t j);
    doublereal operator()(size_t i, size_t j) const;

    //! Return a changeable reference to element (i,j).
    /*!
     * Since this method may alter the element value, it may need to be
     * refactored, so the flag m_factored is set to false.
     *
     * @param i  row
     * @param j  column
     * @returns a reference to the value of the matrix entry
     */
    doublereal& value(size_t i, size_t j);

    //! Return the value of element (i,j).
    /*!
     * This method does not alter the array.
     * @param i  row
     * @param j  column
     * @returns the value of the matrix entry
     */
    doublereal value(size_t i, size_t j) const;

    //! Returns the location in the internal 1D array corresponding to the (i,j)
    //! element in the banded array
    /*!
     * @param i  row
     * @param j  column
     * @returns the index of the matrix entry
     */
    size_t index(size_t i, size_t j) const;

    //! Return the value of the (i,j) element for (i,j) within the bandwidth.
    /*!
     * For efficiency, this method does not check that (i,j) are within the
     * bandwidth; it is up to the calling program to insure that this is true.
     *
     * @param i  row
     * @param j  column
     * @returns the value of the matrix entry
     */
    doublereal _value(size_t i, size_t j) const;

    virtual size_t nRows() const;

    //! Number of columns
    size_t nColumns() const;

    //! Number of subdiagonals
    size_t nSubDiagonals() const;

    //! Number of superdiagonals
    size_t nSuperDiagonals() const;

    //! Return the number of rows of storage needed for the band storage
    size_t ldim() const;

    //! Multiply A*b and write result to \c prod.
    virtual void mult(const doublereal* b, doublereal* prod) const;
    virtual void leftMult(const doublereal* const b, doublereal* const prod) const;

    //! Perform an LU decomposition, the LAPACK routine DGBTRF is used.
    /*!
     * The factorization is saved in ludata.
     *
     * @returns a success flag. 0 indicates a success; ~0 indicates some
     *         error occurred, see the LAPACK documentation
     */
    int factor();

    //! Solve the matrix problem Ax = b
    /*!
     * @param b  INPUT RHS of the problem
     * @param x  OUTPUT solution to the problem
     * @return a success flag. 0 indicates a success; ~0 indicates some error
     *     occurred, see the LAPACK documentation
     */
    int solve(const doublereal* const b, doublereal* const x);

    //! Solve the matrix problem Ax = b
    /*!
     * @param b     INPUT RHS of the problem
     *              OUTPUT solution to the problem
     * @param nrhs  Number of right hand sides to solve
     * @param ldb   Leading dimension of `b`. Default is nColumns()
     * @returns a success flag. 0 indicates a success; ~0 indicates some error
     *     occurred, see the LAPACK documentation
     */
    int solve(doublereal* b, size_t nrhs=1, size_t ldb=0);

    //! Returns an iterator for the start of the band storage data
    /*!
     * Iterator points to the beginning of the data, and it is changeable.
     */
    virtual vector_fp::iterator begin();

    //! Returns an iterator for the end of the band storage data
    /*!
     * Iterator points to the end of the data, and it is changeable.
     */
    vector_fp::iterator end();

    //! Returns a const iterator for the start of the band storage data
    /*!
     * Iterator points to the beginning of the data, and it is not changeable.
     */
    vector_fp::const_iterator begin() const;

    //! Returns a const iterator for the end of the band storage data
    /*!
     * Iterator points to the end of the data, and it is not changeable.
     */
    vector_fp::const_iterator end() const;

    virtual void zero();

    //! Returns an estimate of the inverse of the condition number for the matrix
    /*!
     * The matrix must have been previously factored using the LU algorithm
     *
     * @param a1norm Norm of the matrix
     * @returns the inverse of the condition number
     */
    virtual doublereal rcond(doublereal a1norm);

    //! Returns the factor algorithm used.  This method will always return 0
    //! (LU) for band matrices.
    virtual int factorAlgorithm() const;

    //! Returns the one norm of the matrix
    virtual doublereal oneNorm() const;

    //! Return a pointer to the top of column j
    /*!
     * Column values are assumed to be contiguous in memory (LAPACK band matrix
     * structure)
     *
     * On entry, the matrix A in band storage, in rows KL+1 to 2*KL+KU+1; rows 1
     * to KL of the array need not be set. The j-th column of A is stored in the
     * j-th column of the array AB as follows:
     *
     * AB(KL + KU + 1 + i - j,j) = A(i,j) for max(1, j - KU) <= i <= min(m, j + KL)
     *
     * This routine returns the position of AB(1,j) (fortran-1 indexing) in the
     * above format
     *
     * So to address the (i,j) position, you use the following indexing:
     *
     *     double *colP_j = matrix.ptrColumn(j);
     *     double a_i_j = colP_j[kl + ku + i - j];
     *
     *  @param j   Value of the column
     *  @returns a pointer to the top of the column
     */
    virtual doublereal* ptrColumn(size_t j);

    //! Return a vector of const pointers to the columns
    /*!
     * Note the value of the pointers are protected by their being const.
     * However, the value of the matrix is open to being changed.
     *
     * @returns a vector of pointers to the top of the columns of the matrices.
     */
    virtual doublereal* const* colPts();

    //! Check to see if we have any zero rows in the Jacobian
    /*!
     * This utility routine checks to see if any rows are zero. The smallest row
     * is returned along with the largest coefficient in that row
     *
     * @param valueSmall  OUTPUT value of the largest coefficient in the smallest row
     * @return index of the row that is most nearly zero
     */
    virtual size_t checkRows(doublereal& valueSmall) const;

    //! Check to see if we have any zero columns in the Jacobian
    /*!
     * This utility routine checks to see if any columns are zero. The smallest
     * column is returned along with the largest coefficient in that column
     *
     * @param valueSmall  OUTPUT value of the largest coefficient in the smallest column
     * @return index of the column that is most nearly zero
     */
    virtual size_t checkColumns(doublereal& valueSmall) const;

    //! LAPACK "info" flag after last factor/solve operation
    int info() const { return m_info; };

protected:
    //! Matrix data
    vector_fp data;

    //! Factorized data
    vector_fp ludata;

    //! Number of rows and columns of the matrix
    size_t m_n;

    //! Number of subdiagonals of the matrix
    size_t m_kl;

    //! Number of super diagonals of the matrix
    size_t m_ku;

    //! value of zero
    doublereal m_zero;

    struct PivData; // pImpl wrapper class

    //! Pivot vector
    std::unique_ptr<PivData> m_ipiv;

    //! Vector of column pointers
    std::vector<doublereal*> m_colPtrs;
    std::vector<double*> m_lu_col_ptrs;

    //! Extra work array needed - size = n
    vector_int iwork_;

    //! Extra dp work array needed - size = 3n
    vector_fp work_;

    int m_info;
};

//! Utility routine to print out the matrix
/*!
 *  @param s  ostream to print the matrix out to
 *  @param m  Matrix to be printed
 *  @returns a reference to the ostream
 */
std::ostream& operator<<(std::ostream& s, const BandMatrix& m);

}

#endif
