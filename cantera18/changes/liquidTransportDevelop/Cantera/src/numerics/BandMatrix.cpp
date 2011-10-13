/**
 *  @file BandMatrix.cpp
 *
 *  Banded matrices.
 */

// Copyright 2001  California Institute of Technology

#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "BandMatrix.h"
#include "ctlapack.h"
#include "utilities.h"
#include "ctexceptions.h"
#include "stringUtils.h"
#include "global.h"

using namespace std;

namespace Cantera {

  //====================================================================================================================
  BandMatrix::BandMatrix() : 
    m_factored(false),
    m_n(0), 
    m_kl(0), 
    m_ku(0),
    m_zero(0.0)
  { 
    data.clear();
    ludata.clear(); 
  }
  //====================================================================================================================
  BandMatrix::BandMatrix(int n, int kl, int ku, doublereal v)   :
    m_factored(false), 
    m_n(n),
    m_kl(kl),
    m_ku(ku),
    m_zero(0.0)
  {
    data.resize(n*(2*kl + ku + 1));
    ludata.resize(n*(2*kl + ku + 1));
    fill(data.begin(), data.end(), v);
    fill(ludata.begin(), ludata.end(), 0.0);
    m_ipiv.resize(m_n);
  }
  //====================================================================================================================
  BandMatrix::BandMatrix(const BandMatrix& y) :
    m_factored(false),
    m_n(0), 
    m_kl(0), 
    m_ku(0),
    m_zero(0.0)
 {
    m_n = y.m_n;
    m_kl = y.m_kl;
    m_ku = y.m_ku;
    data = y.data;
    ludata = y.ludata;
    m_factored = y.m_factored;
    m_ipiv = y.m_ipiv;
  }
  //====================================================================================================================
  BandMatrix::~BandMatrix() {

  }
  //====================================================================================================================
  BandMatrix& BandMatrix::operator=(const BandMatrix & y) {
    if (&y == this) return *this;
    m_n = y.m_n;
    m_kl = y.m_kl;
    m_ku = y.m_ku;
    m_ipiv = y.m_ipiv;
    data = y.data;
    ludata = y.ludata;
    m_factored = y.m_factored;
    return *this;
  }
  //====================================================================================================================
  void BandMatrix::resize(int n, int kl, int ku, doublereal v) {
    m_n = n;
    m_kl = kl;
    m_ku = ku;
    data.resize(n*(2*kl + ku + 1));
    ludata.resize(n*(2*kl + ku + 1));
    m_ipiv.resize(m_n);
    fill(data.begin(), data.end(), v);
    m_factored = false;
  }  
  //====================================================================================================================
  void BandMatrix::bfill(doublereal v) {
    std::fill(data.begin(), data.end(), v);
    m_factored = false;
  }
  //====================================================================================================================
  doublereal& BandMatrix::operator()(int i, int j) {
    return value(i,j);
  }
  //====================================================================================================================
  doublereal BandMatrix::operator()(int i, int j) const {
    return value(i,j);
  }
  //====================================================================================================================
  doublereal& BandMatrix::value(int i, int j) {
    m_factored = false;
    if (i < j - m_ku || i > j + m_kl) {
      return m_zero;
    }
    return data[index(i,j)];
  }  
  //====================================================================================================================
  doublereal  BandMatrix::value( int i, int j) const {
    if (i < j - m_ku || i > j + m_kl) return 0.0;
    return data[index(i,j)];
  }
  //====================================================================================================================
  int  BandMatrix::index(int i, int j) const {
    int rw = m_kl + m_ku + i - j;
    return (2*m_kl + m_ku + 1)*j + rw;
  }
  //====================================================================================================================
  doublereal BandMatrix::_value(int i, int j) const {
    return data[index(i,j)];
  }
  //====================================================================================================================
  // Number of rows
  int  BandMatrix::nRows() const { 
    return m_n;
  }
  //====================================================================================================================
  // Number of columns
  int BandMatrix::nColumns() const { 
    return m_n;
  }
  //====================================================================================================================
  // Number of subdiagonals
  int BandMatrix::nSubDiagonals() const { 
    return m_kl;
  }
  //====================================================================================================================
  // Number of superdiagonals
  int BandMatrix::nSuperDiagonals() const {
    return m_ku;
  }
  //====================================================================================================================
  int  BandMatrix::ldim() const { 
    return 2*m_kl + m_ku + 1;
  }
  //====================================================================================================================
  vector_int &  BandMatrix::ipiv() { 
    return m_ipiv; 
  }
  //====================================================================================================================            
  /*
   * Multiply A*b and write result to \c prod.
   */
  void BandMatrix::mult(const doublereal * const b, doublereal * const prod) const {
    int nr = nRows();
    doublereal sum = 0.0;
    for (int m = 0; m < nr; m++) {
      sum = 0.0;
      for (int j = m - m_kl; j <= m + m_ku; j++) {
	if (j >= 0 && j < m_n)
	  sum += _value(m,j) * b[j];
      }
      prod[m] = sum;
    }
  }
  //====================================================================================================================
  /*
   * Multiply b*A and write result to \c prod.
   */
  void BandMatrix::leftMult(const doublereal * const b, doublereal * const prod) const {
    int nc = nColumns();
    doublereal sum = 0.0;
    for (int n = 0; n < nc; n++) {
      sum = 0.0;
      for (int i = n - m_ku; i <= n + m_kl; i++) {
	if (i >= 0 && i < m_n)
	  sum += _value(i,n) * b[i];
      }
      prod[n] = sum;
    }
  }
  //====================================================================================================================
  /*
   * Perform an LU decomposition. LAPACK routine DGBTRF is used.
   * The factorization is saved in ludata.
   */
  int BandMatrix::factor() {
    int info=0;
    copy(data.begin(), data.end(), ludata.begin());
    ct_dgbtrf(nRows(), nColumns(), nSubDiagonals(), nSuperDiagonals(), 
	      DATA_PTR(ludata), ldim(), DATA_PTR(ipiv()), info);

    // if info = 0, LU decomp succeeded. 
    if (info == 0) {
      m_factored = true;
    } else {
      m_factored = false;
      ofstream fout("bandmatrix.csv");
      fout << *this << endl;
      fout.close();
    }
    return info;
  }
  //====================================================================================================================
  int BandMatrix::solve(int n, const doublereal * const b, doublereal * const x) {
    copy(b, b+n, x);
    return solve(n, x);
  }
  //====================================================================================================================
  int BandMatrix::solve(int n, doublereal* b) {
    int info = 0;
    if (!m_factored) info = factor();
    if (info == 0)
      ct_dgbtrs(ctlapack::NoTranspose, nColumns(), nSubDiagonals(), 
		nSuperDiagonals(), 1, DATA_PTR(ludata), ldim(), 
		DATA_PTR(ipiv()), b, nColumns(), info);

    // error handling
    if (info != 0) {
      ofstream fout("bandmatrix.csv");
      fout << *this << endl;
      fout.close();
    }
    return info;
  } 
  //====================================================================================================================
  vector_fp::iterator  BandMatrix::begin() {
    m_factored = false;
    return data.begin();
  }
  //====================================================================================================================
  vector_fp::iterator BandMatrix::end() {
    m_factored = false;
    return data.end();
  }
  //====================================================================================================================
  vector_fp::const_iterator BandMatrix::begin() const {
    return data.begin(); 
  }
  //====================================================================================================================
  vector_fp::const_iterator BandMatrix::end() const {
    return data.end();
  }
  //====================================================================================================================
  ostream& operator<<(ostream& s, const BandMatrix& m) {
    int nr = m.nRows();
    int nc = m.nColumns();
    for (int i = 0; i < nr; i++) {
      for (int j = 0; j < nc; j++) {
	s << m(i,j) << ", ";
      }
      s << endl;
    }
    return s;
  }
 //====================================================================================================================
}

