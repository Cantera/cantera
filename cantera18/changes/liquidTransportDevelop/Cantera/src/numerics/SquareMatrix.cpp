/**
 *  @file DenseMatrix.cpp
 *
 */
/*
 * $Revision$
 * $Date$
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
#include "global.h"

#include <iostream>
#include <vector>

#include <cstring>

using namespace std;

namespace Cantera {


  //====================================================================================================================
  SquareMatrix::SquareMatrix() :
    DenseMatrix(),
    m_factored(0),
    useQR_(0)
  {
  }

  //====================================================================================================================
  // Constructor.
  /*
   * Create an \c n by \c n matrix, and initialize
   * all elements to \c v.
   *
   * @param n   size of the square matrix
   * @param v   intial value of all matrix components.
   */
  SquareMatrix::SquareMatrix(int n, doublereal v)  : 
    DenseMatrix(n, n, v),
    m_factored(0),
    useQR_(0)
  {
  }
  //====================================================================================================================
  /*
   *
   * copy constructor
   */
  SquareMatrix::SquareMatrix(const SquareMatrix& y) :
    DenseMatrix(y),
    m_factored(y.m_factored),
    useQR_(y.useQR_)
  {
  }
    
  //====================================================================================================================
  /*
   * Assignment operator
   */
  SquareMatrix& SquareMatrix::operator=(const SquareMatrix& y) {
    if (&y == this) return *this;
    DenseMatrix::operator=(y);
    m_factored = y.m_factored;
    return *this;
  }
  //====================================================================================================================
  SquareMatrix::~SquareMatrix() {
  }
  //====================================================================================================================
  /*
   * Solve Ax = b. Vector b is overwritten on exit with x.
   */
  int SquareMatrix::solve(double* b) 
  {
    if (useQR_) {
      return solveQR(b);
    }
    int info=0;
    /*
     * Check to see whether the matrix has been factored.
     */
    if (!m_factored) {
      int retn = factor();
      if (retn) {
	return retn;
      }
    }
    /*
     * Solve the factored system
     */
    ct_dgetrs(ctlapack::NoTranspose, static_cast<int>(nRows()),
	      1, &(*(begin())), static_cast<int>(nRows()), 
	      DATA_PTR(ipiv()), b, static_cast<int>(nColumns()), info);
    if (info != 0) {
      if (m_printLevel) {
	writelogf("SquareMatrix::solve(): DGETRS returned INFO = %d\n", info);
      }
      if (! m_useReturnErrorCode) {
	throw CELapackError("SquareMatrix::solve()", "DGETRS returned INFO = " + int2str(info));
      }
    }
    return info;
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
  //====================================================================================================================
  void SquareMatrix::resize(int n, int m, doublereal v) {
    DenseMatrix::resize(n, m, v);
  } 
  //====================================================================================================================
  /*
   * Factor A. A is overwritten with the LU decomposition of A.
   */
  int SquareMatrix::factor() {
    if (useQR_) {
      return factorQR();
    }
    integer n = static_cast<int>(nRows());
    int info=0;
    m_factored = 1;
    ct_dgetrf(n, n, &(*(begin())), static_cast<int>(nRows()),
	      DATA_PTR(ipiv()), info);
    if (info != 0) {
      if (m_printLevel) {
	writelogf("SquareMatrix::factor(): DGETRS returned INFO = %d\n", info);
      }
      if (! m_useReturnErrorCode) {
        throw CELapackError("SquareMatrix::factor()", "DGETRS returned INFO = "+int2str(info));
      }
    }
    return info;
  }
  //=====================================================================================================================
  /*
   * clear the factored flag
   */
  void SquareMatrix::clearFactorFlag() {
    m_factored = 0;
  }
  //=====================================================================================================================
  /*
   * set the factored flag
   */
  void SquareMatrix::setFactorFlag() {
    m_factored = 1;
  }
  //=====================================================================================================================
  int SquareMatrix::factorQR() {
     if ((int) tau.size() < m_nrows)  {
       tau.resize(m_nrows, 0.0);
       work.resize(8 * m_nrows, 0.0);
     }
     int info;
     m_factored = 2;
     int lwork = work.size(); 
     ct_dgeqrf(m_nrows, m_nrows, &(*(begin())), m_nrows, DATA_PTR(tau), DATA_PTR(work), lwork, info); 
     if (info != 0) {
      if (m_printLevel) {
	writelogf("SquareMatrix::factorQR(): DGEQRF returned INFO = %d\n", info);
      }
      if (! m_useReturnErrorCode) {
        throw CELapackError("SquareMatrix::factorQR()", "DGEQRF returned INFO = " + int2str(info));
      }
     }
     int lworkOpt = work[0];
     if (lworkOpt > lwork) {
       work.resize(lworkOpt);
     }
     return info;
  }
  //=====================================================================================================================
  /*
   * Solve Ax = b. Vector b is overwritten on exit with x.
   */
  int SquareMatrix::solveQR(double* b)
  {
    int info=0;
    /*
     * Check to see whether the matrix has been factored.
     */
    if (!m_factored) {
      int retn = factorQR();
      if (retn) {
        return retn;
      }
    }
    
    int lwork = work.size(); 
    if (lwork < m_nrows) {
      work.resize(8 * m_nrows, 0.0);
      lwork = 8 * m_nrows;
    }

    /*
     * Solve the factored system
     */
    ct_dormqr(ctlapack::Left, ctlapack::Transpose, m_nrows, 1, m_nrows, &(*(begin())), m_nrows, DATA_PTR(tau), b, m_nrows, 
                        DATA_PTR(work), lwork, info);
    if (info != 0) {
      if (m_printLevel) {
        writelogf("SquareMatrix::solveQR(): DORMQR returned INFO = %d\n", info);
      }
      if (! m_useReturnErrorCode) {
        throw CELapackError("SquareMatrix::solveQR()", "DORMQR returned INFO = " + int2str(info));
      }
    }
    int lworkOpt = work[0];
    if (lworkOpt > lwork) {
      work.resize(lworkOpt);
    }

    char dd = 'N';

    ct_dtrtrs(ctlapack::UpperTriangular, ctlapack::NoTranspose, &dd, m_nrows, 1,  &(*(begin())), m_nrows, b,
              m_nrows, info);
    if (info != 0) {
      if (m_printLevel) {
        writelogf("SquareMatrix::solveQR(): DTRTRS returned INFO = %d\n", info);
      }
      if (! m_useReturnErrorCode) {
        throw CELapackError("SquareMatrix::solveQR()", "DTRTRS returned INFO = " + int2str(info));
      }
    }

    return info;
  }
  //=====================================================================================================================




}

