/**
 *  @file  DoubleStarStar.h
 *
 *  Header file for class DoubleStarStar
 */

/*
 *  $Author: hkmoffa $
 *  $Revision: 1.2 $
 *  $Date: 2008/08/01 16:05:38 $
 */

#ifndef VCS_DOUBLESTARSTAR_H
#define VCS_DOUBLESTARSTAR_H

#include <vector>

namespace VCSnonideal {

//!  A class for 2D double arrays storred in column-major
//!  (Fortran-compatible) form.
/*!
 *  In this form, the data entry for an n row, m col
 *  matrix is
 *       index = i + (n-1) * j
 *  where
 *     Matrix[j][i]
 *         i = row
 *         j = column
 *   The way this is instantiated is via the constructor:
 *         DoubleStarStar Dmatrix(mcol, mrow);
 *
 *   The way this is referenced is via the notation:
 *        Dmatrix[icol][irow]
 */
class DoubleStarStar {

public:

  //! Default constructor. Create an empty array.
  DoubleStarStar();

  //! Constructor.
  /*!
   *  Create an \c nrow by \c mcol double array, and initialize
   *  all elements to \c v.
   *
   * @param mcol  Number of columns
   * @param nrow  Number of rows
   */
  DoubleStarStar(int mcol, int nrow, double v = 0.0);
    
  //! copy constructor
  /*!
   * @param y object to be copied
   */
  DoubleStarStar(const DoubleStarStar& y);

  /// assignment operator
  /*!
   * @param y object to be copied
   */
  DoubleStarStar& operator=(const DoubleStarStar& y);

  //! Resize the array, and fill the new entries with 'v'
  /*!
   * @param mrow  This is the number of columns in the new matrix
   * @param ncol  This is the number of rows
   * @param v     Default fill value -> defaults to zero.
   */
  void resize(int mcol, int nrow, double v = 0.0);

  //! Pointer to the top of the column
  /*!
   *  @param jcol  This is the jth column
   *
   * @return returns the pointer to the top of the jth column
   */
  double * const operator[](int jcol);

  //! Returns a const Pointer to the top of the jth column
  /*!
   *  @param jcol  This is the jth column
   *
   * @return returns the pointer to the top of the jth column
   */
  const double * const operator[](int jcol) const;

  //! Returns a double ** pointer to the base address
  /*!
   *  This is the second way to get to the data
   *  This returns a double ** which can later be used in
   *  Dmatrix[icol][irow] notation to get to the data
   */
  double * const * const baseDataAddr();

  //! Returns a const double ** pointer to the base address
  /*!
   *  This is the second way to get to the data
   *  This returns a double ** which can later be used in
   *  Dmatrix[icol][irow] notation to get to the data
   */
  double const * const * const constBaseDataAddr() const;

  //! Number of rows
  int nRows() const;

  //! Number of columns
  int nColumns() const;

private:
  //! Storage area
  std::vector<double> m_data;

  //! Vector of addresses for the top of the columns
  /*!
   *  Length = mcol
   */
  std::vector<double *> m_colAddr; 

  //! number of rows
  int m_nrows;

  //! number of columns
  int m_ncols;
};

}

#endif
  
  
