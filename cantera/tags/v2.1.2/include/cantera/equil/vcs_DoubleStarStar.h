/**
 *  @file  vcs_DoubleStarStar.h
 *
 *  Header file for class DoubleStarStar
 */
#ifndef VCS_DOUBLESTARSTAR_H
#define VCS_DOUBLESTARSTAR_H

#include <vector>

namespace VCSnonideal
{

using std::size_t;

//!  A class for 2D double arrays stored in column-major
//!  (Fortran-compatible) form.
/*!
 *  In this form, the data entry for an `n` row, `m` colum matrix is index =
 *  `i + (n-1) * j` where `Matrix[j][i]` references the element in row `i`,
 *  column `j`.
 *
 * The way this is instantiated is via the constructor,
 * DoubleStarStar Dmatrix(mcol, mrow)`.
 *
 * The way this is referenced is via the notation: `Dmatrix[icol][irow]`.
 */
class DoubleStarStar
{
public:

    //! Default constructor. Create an empty array.
    DoubleStarStar();

    //! Constructor.
    /*!
     *  Create an `nrow` by `mcol` double array, and initialize all elements
     *  to `v`.
     *
     * @param mcol  Number of columns
     * @param nrow  Number of rows
     * @param v     value used to initialize elements
     */
    DoubleStarStar(size_t mcol, size_t nrow, double v = 0.0);

    DoubleStarStar(const DoubleStarStar& y);
    DoubleStarStar& operator=(const DoubleStarStar& y);

    //! Resize the array, and fill the new entries with `v`
    /*!
     * @param mcol  This is the number of columns in the new matrix
     * @param nrow  This is the number of rows
     * @param v     Default fill value -> defaults to zero.
     */
    void resize(size_t mcol, size_t nrow, double v = 0.0);

    //! Pointer to the top of the column
    /*!
     *  @param jcol  This is the jth column
     *
     * @return returns the pointer to the top of the jth column
     */
    double* operator[](size_t jcol);

    //! Returns a const Pointer to the top of the jth column
    /*!
     *  @param jcol  This is the jth column
     *
     * @return returns the pointer to the top of the jth column
     */
    const double* operator[](size_t jcol) const;

    //! Returns a `double**` pointer to the base address
    /*!
     *  This is the second way to get to the data. This returns a `double**`
     *  which can later be used in `Dmatrix[icol][irow]` notation to get to
     *  the data.
     */
    double* const* baseDataAddr();

    //! Returns a `const double**` pointer to the base address
    /*!
     *  This is the second way to get to the data This returns a double **
     *  which can later be used in `Dmatrix[icol][irow]` notation to get to
     *  the data.
     */
    double const* const* constBaseDataAddr() const;

    //! Number of rows
    size_t nRows() const;

    //! Number of columns
    size_t nColumns() const;

private:
    //! Storage area
    std::vector<double> m_data;

    //! Vector of addresses for the top of the columns
    /*!
     *  Length = mcol
     */
    std::vector<double*> m_colAddr;

    //! number of rows
    size_t m_nrows;

    //! number of columns
    size_t m_ncols;
};

}

#endif
