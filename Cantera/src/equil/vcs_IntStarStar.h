/**
 *  @file  IntStarStar.h
 *
 *  Header file for class IntStarStar
 */

/*
 *  $Author: hkmoffa $
 *  $Revision: 1.1 $
 *  $Date: 2007/12/19 21:19:28 $
 */

#ifndef VCS_INTSTARSTAR_H
#define VCS_INTSTARSTAR_H

#include <vector>

namespace VCSnonideal {

//!  A class for 2D int arrays storred in column-major
//!  (Fortran-compatible) form.
/*!
 *  In this form, the data entry for an n row, m col
 *  matrix is
 *       index = i + (n-1) * j
 *  where
 *     Matrix[j][i]
 *         i = row
 *         j = column
 */
class IntStarStar {

public:

  //! Default constructor. Create an empty array.
  IntStarStar();

  //! Constructor.
  /*!
   *  Create an \c nrow by \c mcol int array, and initialize
   *  all elements to \c v.
   *
   * @param mcol  Number of columns
   * @param nrow  Number of rows
   */
  IntStarStar(int mcol, int nrow, int v = 0);
    
  //! copy constructor
  IntStarStar(const IntStarStar& y);

  /// assignment operator
  IntStarStar& operator=(const IntStarStar& y);


  //! Resize the array, and fill the new entries with 'v'
  /*!
   * @param mcol  This is the number of columns in the new matrix
   * @param nrow  This is the number of rows
   * @param v     Default fill value -> defaults to zero.
   */
  void resize(int mcol, int nrow, int v = 0);

  //! Pointer to the top of the column
  /*!
   *  @param jcol Pointer to the top of the jth column
   */
  int * const operator[](int jcol);

  //! Pointer to the top of the column
  /*!
   *  @param j Pointer to the top of the jth column
   */
  const int * const operator[](int jcol) const;

  //! Returns a int ** pointer to the base address
  /*!
   *  This is the second way to get to the data
   *  This returns a int ** which can later be used in
   *  Imatrix[icol][irow] notation to get to the data
   */
  int * const * const baseDataAddr();

  //! Number of rows
  int nRows() const;

  //! Number of columns
  int nColumns() const;

private:
  //! Storage area
  std::vector<int> m_data;

  std::vector<int *> m_colAddr;

  //! number of rows
  int m_nrows;

  //! number of columns
  int m_ncols;
};

}

#endif
  
  
