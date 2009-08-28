/**
 *  @file vec_functions.h
 *  Templates for operations on vector-like objects. 
 */
/*
 *  $Date: 2009/01/15 20:08:50 $
 *  $Revision: 1.3 $
 *
 *  Copyright 2001 California Institute of Technology
 *
 */

#ifndef CT_VEC_FUNCTIONS_H
#define CT_VEC_FUNCTIONS_H

#include "ct_defs.h"
#include "utilities.h"
#include <numeric>
#include <functional>
#include <algorithm>

#include <cstring>

namespace Cantera {


  //! Templated function that copies the first n entries from x to y.
  /*!
   *  
   *
   *  The  templated type is the type of x and y
   *
   *  @param n        Number of elements to copy from x to y
   *  @param x        The object x, of templated type const T&
   *  @param y        The object y, of templated type  T&
   */
  template<class T>
  inline void copyn(size_t n, const T& x, T& y) {
    std::copy(x.begin(), x.begin() + n, y.begin());
  }

  /**
   * Divide each element of x by the corresponding element of y.
   * This function replaces x[n] by x[n]/y[n], for 0 <= n < x.size()
   */
  template<class T>   
  inline void divide_each(T& x, const T& y) {
    std::transform(x.begin(), x.end(), y.begin(), 
		   x.begin(), std::divides<TYPENAME_KEYWORD T::value_type>());
  }
    
  /**
   * multiply each element of x by the corresponding element of y.
   * This function replaces x[n] by x[n]*y[n], for 0 <= n < x.size()
   */
  template<class T>
  inline void multiply_each(T& x, const T& y) {
    std::transform(x.begin(), x.end(), y.begin(), 
		   x.begin(), std::multiplies<TYPENAME_KEYWORD T::value_type>());
  }

  /**
   * Multiply each element of x by scale_factor.
   */
  template<class T, class S>
  inline void scale(T& x, S scale_factor) {
    scale(x.begin(), x.end(), x.begin(), scale_factor);
  }

  /**
   * Returns the sum of x[n]*y[n], for 0 <= n < x.size().
   */
  template<class T>
  inline doublereal dot_product(const T& x, const T& y) {
    return std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
  }

  /**
   * Returns the sum of x[n]/y[n], for 0 <= n < x.size().
   */
  template<class T>
  inline doublereal dot_ratio(const T& x, const T& y) {
    return _dot_ratio(x.begin(), x.end(), y.begin(), 0.0);
  }

  /**
   * Replaces x[n] by x[n] + y[n] for 0 <= n < x.size()
   */
  template<class T>
  inline void add_each(T& x, const T& y) {
    std::transform(x.begin(), x.end(), y.begin(), 
		   x.begin(), std::plus<TYPENAME_KEYWORD T::value_type>());
  }


  //! Templated dot ratio class
  /*!
   *  Calculates the quantity:
   *
   *      S += x[n]/y[n]
   *
   *  The first templated type is the iterator type for x[] and y[].
   *  The second templated type is the type of S.
   *
   *  @param x_begin  InputIter type, indicating the address of the
   *                  first element of x
   *  @param x_end    InputIter type, indicating the address of the
   *                  last element of x
   *  @param y_begin  InputIter type, indicating the address of the
   *                  first element of y
   *  @param start_value S type, indicating the type of the
   *                   accumulation result.
   */
  template<class InputIter, class S>
  inline doublereal _dot_ratio(InputIter x_begin, InputIter x_end, 
			       InputIter y_begin, S start_value) {
    for (; x_begin != x_end; ++x_begin, ++y_begin)
      start_value += *x_begin / *y_begin;
    return start_value;
  }

  /**
   * Finds the entry in a vector with maximum absolute
   * value, and return this value.
   */
  template<class T>
  inline T absmax(const std::vector<T>& v) {
    int n = v.size();
    T val;
    T maxval = 0.0;
    for (int i = 0; i < n; i++) {
      val = v[i];
      if (val < 0) val = -val;
      if (val > maxval) maxval = val;
    }
    return maxval;
  }

  //! Copy a vector of doubles in an efficient, fast manner.
  /*!
   *  No checking is done other that to check that len is greater than 0
   *
   *  @param copyTo  Vector to receive the copy
   *  @param copyFrom Vector from which the copy is coming
   *  @param len  Length of the copy
   */
  inline void fbo_copy_dbl_1(doublereal * const copyTo, const doublereal * const copyFrom, 
			     const int len) {
    if (len > 0) {
      (void) memcpy((void *)copyTo, (const void *)copyFrom, len * sizeof(doublereal));
    }
  }

  //! Copy a vector<doubles> in an efficient, fast manner.
  /*!
   *  No checking is done other that to check that len is greater than 0
   *
   *  @param copyTo  Vector to receive the copy
   *  @param copyFrom Vector from which the copy is coming
   *  @param len  Length of the copy
   */
  inline void fvo_copy_dbl_1(std::vector<doublereal> &copyTo, const std::vector<doublereal> &copyFrom, 
			     const int len) {
    if (len > 0) {
      (void) memcpy((void *)(&copyTo[0]), (const void *)(&copyFrom[0]), len * sizeof(doublereal));
    }
  }

  //! Zero a double vector in an efficient, fast manner.
  /*!
   *  No checking is done other that to check that len is greater than 0
   *
   *  @param v    Vector to be zeroed
   *  @param len  Length of the copy
   */
  inline void fbo_zero_dbl_1(doublereal * const v, const int len) {
    if (len > 0) {
      (void) memset((void *)v, 0, len * sizeof(doublereal));
    }
  }

  //! Zero a vector<doubles> in an efficient, fast manner.
  /*!
   *  No checking is done other that to check that len is greater than 0
   *
   *  @param v    Vector to be zeroed
   *  @param len  Length of the copy
   */
  inline void fvo_zero_dbl_1(std::vector<doublereal> &v, const int len) {
    if (len > 0) {
      (void) memset((void *)(&v[0]), 0, len * sizeof(doublereal));
    }
  }

}

#endif
