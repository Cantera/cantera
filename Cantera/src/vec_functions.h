/**
 *
 *  @file vec_functions.h
 *
 *  Templates for operations on vector-like objects. 
 *
 *  $Author$
 *  $Date$
 *  $Revision$
 *
 *  Copyright 2001 California Institute of Technology
 *
 */

#ifndef CT_VEC_FUNCTIONS_H
#define CT_VEC_FUNCTIONS_H

#include "ct_defs.h"
#include <numeric>
#include <functional>
#include <algorithm>

namespace Cantera {

    /**
     * Copy the first n entries from x to y. Both x and y must have
     * size greater than or equal to n.
     */
    template<class T>
    inline void copyn(size_t n, const T& x, T& y) {
        copy(x.begin(), x.begin() + n, y.begin());
    }

    /**
     * Divide each element of x by the corresponding element of y.
     * This function replaces x[n] by x[n]/y[n], for 0 <= n < x.size()
     */
    template<class T>   
    inline void divide_each(T& x, const T& y) {
        transform(x.begin(), x.end(), y.begin(), 
            x.begin(), divides<TYPENAME_KEYWORD T::value_type>());
    }
    
    /**
     * multiply each element of x by the corresponding element of y.
     * This function replaces x[n] by x[n]*y[n], for 0 <= n < x.size()
     */
    template<class T>
    inline void multiply_each(T& x, const T& y) {
        transform(x.begin(), x.end(), y.begin(), 
            x.begin(), multiplies<TYPENAME_KEYWORD T::value_type>());
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
        return inner_product(x.begin(), x.end(), y.begin(), 0.0);
    }

    /**
     * Returns the sum of x[n]/y[n], for 0 <= n < x.size().
     */
    template<class T>
    inline doublereal dot_ratio(const T& x, const T& y) {
        return __dot_ratio(x.begin(), x.end(), y.begin(), 0.0);
    }

    /**
     * Replaces x[n] by x[n] + y[n] for 0 <= n < x.size()
     */
    template<class T>
    inline void add_each(T& x, const T& y) {
        transform(x.begin(), x.end(), y.begin(), 
            x.begin(), plus<TYPENAME_KEYWORD T::value_type>());
    }

    template<class _InputIter, class S>
    inline doublereal __dot_ratio(_InputIter __x_begin, _InputIter __x_end, 
        _InputIter __y_begin, S start_value) {
        for (; __x_begin != __x_end; ++__x_begin, ++__y_begin)
            start_value += *__x_begin / *__y_begin;
        return start_value;
    }

   /**
     * Finds the entry in a vector with maximum absolute
     * value, and return this value.
     */
    template<class T>
    inline T absmax(const vector<T>& v) {
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

}

#endif
