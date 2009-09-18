/**
 *  @file DenseMatrix.cpp
 *
 */
/*
 * $Revision: 1.1 $
 * $Date: 2009/01/27 15:59:13 $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "ct_defs.h"
#include "stringUtils.h"
#include "ctlapack.h"
#include "SquareMatrix.h"

#include <iostream>
#include <vector>

#include <cstring>

using namespace std;

namespace Cantera {
  /**
   *
   * copy constructor
   */
  SquareMatrix::SquareMatrix(const SquareMatrix& y) :
    DenseMatrix(y),
    m_factored(y.m_factored)
  {
  }
    
  /**
   * Assignment operator
   */
  SquareMatrix& SquareMatrix::operator=(const SquareMatrix& y) {
    if (&y == this) return *this;
    DenseMatrix::operator=(y);
    m_factored = y.m_factored;
    return *this;
  }
    
  /**
   * Solve Ax = b. Vector b is overwritten on exit with x.
   */
  int SquareMatrix::solve(double* b) 
  {
    int info=0;
    /*
     * Check to see whether the matrix has been factored.
     */
    if (!m_factored) {
      factor();
    }
    /*
     * Solve the factored system
     */
    ct_dgetrs(ctlapack::NoTranspose, static_cast<int>(nRows()),
	      1, &(*(begin())), static_cast<int>(nRows()), 
	      DATA_PTR(ipiv()), b, static_cast<int>(nColumns()), info);
    if (info != 0)
      throw CanteraError("SquareMatrix::solve",
			 "DGETRS returned INFO = "+int2str(info));
    return 0;
  }

  /**
   * Set all entries to zero
   */
  void SquareMatrix::zero() {
    int n = static_cast<int>(nRows());
    if (n > 0) {
      int nn = n * n;
      double *sm = &m_data[0];
      /*
       * Using memset is the fastest way to zero a contiguous
       * section of memory.
       */
      (void) memset((void *) sm, 0, nn * sizeof(double));
    }
  }
   
  /**
   * Factor A. A is overwritten with the LU decomposition of A.
   */
  int SquareMatrix::factor() {
    integer n = static_cast<int>(nRows());
    int info=0;
    m_factored = true;
    ct_dgetrf(n, n, &(*(begin())), static_cast<int>(nRows()),
	      DATA_PTR(ipiv()), info);
    if (info != 0) {
      cout << "Singular matrix, info = " << info << endl;
      throw CanteraError("invert",
			 "DGETRF returned INFO="+int2str(info));
    }
    return 0;
  }
  /*
   * clear the factored flag
   */
  void SquareMatrix::clearFactorFlag() {
    m_factored = false;
  }
  /**
   * set the factored flag
   */
  void SquareMatrix::setFactorFlag() {
    m_factored = true;
  }
}

