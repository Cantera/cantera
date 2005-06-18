/**
 *  @file utilities.h
 *
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_UTILITIES_H
#define CT_UTILITIES_H

#include "ct_defs.h"

namespace Cantera {

    /**
     * Maximum of i and j. If \a i and \a j have different types, \a j
     * is converted to the type of \a i before the comparison.
     */
    template<class T, class S>
    inline T max(T i, S j) {
        return (i > T(j) ? i : T(j));
    }

    /**
     * Minimum of i and j. If \a i and \a j have different types, \a j
     * is converted to the type of \a i before the comparison.
     */
    template<class T, class S>
    inline T min(T i, S j) {
        return (i < T(j) ? i : T(j));
    }

    /**
     * Inner product of two vectors of length 4. 
     * If either \a x
     * or \a y has length greater than 4, only the first 4 elements
     * will be used.
     */
    template<class V>
    inline doublereal dot4(const V& x, const V& y) {
        return x[0]*y[0] + x[1]*y[1] + x[2]*y[2] + x[3]*y[3];
    }

    /**
     * Inner product of two vectors of length 5.
     * If either \a x
     * or \a y has length greater than 5, only the first 5 elements
     * will be used.
     */
    template<class V>
    inline doublereal dot5(const V& x, const V& y) {
        return x[0]*y[0] + x[1]*y[1] + x[2]*y[2] + x[3]*y[3] +
            x[4]*y[4];
    }

    /**
     * Inner product of two vectors of length 6.
     * If either \a x
     * or \a y has length greater than 6, only the first 6 elements
     * will be used.
     */
    template<class V>
    inline doublereal dot6(const V& x, const V& y) {
        return x[0]*y[0] + x[1]*y[1] + x[2]*y[2] + x[3]*y[3] +
            x[4]*y[4] + x[5]*y[5];
    }

    /**
     * Inner product.
     */
    template<class _InputIter, class _InputIter2>
    inline doublereal dot(_InputIter x_begin, _InputIter x_end, 
        _InputIter2 y_begin) {
        doublereal sum = 0.0;
        for(; x_begin != x_end; ++x_begin, ++y_begin) 
            sum += *x_begin * *y_begin;
        return sum;
    }

    /**
     * Multiply elements of an array by a scale factor.
     * \code
     * vector_fp in(8, 1.0), out(8);
     * scale(in.begin(), in.end(), out.begin(), factor);
     * \endcode 
     */ 
    template<class _InputIter, class _OutputIter, class S>
    inline void scale(_InputIter __begin, _InputIter __end, 
        _OutputIter __out, S scale_factor) {
        for (; __begin != __end; ++__begin, ++__out) 
            *__out = scale_factor * *__begin;
    }

    template<class _InputIter, class _OutputIter, class S>
    inline void increment_scale(_InputIter __begin, _InputIter __end, 
        _OutputIter __out, S scale_factor) {
        for (; __begin != __end; ++__begin, ++__out) 
            *__out += scale_factor * *__begin;
    }

    /**
     * Multiply each entry in x by the corresponding entry in y.
     */
    template<class _InputIter, class _OutputIter>   
    inline void multiply_each(_OutputIter x_begin, _OutputIter x_end, 
        _InputIter y_begin) {
        for(; x_begin != x_end; ++x_begin, ++y_begin) *x_begin *= *y_begin;
    }


    /**
     * Invoke method 'resize' with argument \a m for a sequence of objects.
     */
    template<class _InputIter>
    inline void _resize_each(int m, _InputIter __begin, _InputIter __end) {
        for(; __begin != __end; ++__begin) __begin->resize(m);
    }

    /**
     * The maximum absolute value.
     */
    template<class _InputIter>
    inline doublereal absmax(_InputIter __begin, _InputIter __end) {
        doublereal amax = 0.0;
        for(; __begin != __end; ++__begin) 
            if (fabs(*__begin) > amax) amax = fabs(*__begin);
        return amax;
    }

    /**
     * Normalize the values in a sequence, such that they sum to 1.0.
     */
    template<class _InputIter, class _OutputIter>
    inline void normalize(_InputIter __begin, _InputIter __end, 
        _OutputIter __out) {
        doublereal sum = accumulate(__begin, __end, 0.0);
        for (; __begin != __end; ++__begin, ++__out) *__out = *__begin/sum;
    }

    /**
     * Divide each element of \a x by the corresponding element of \a y.
     */
    template<class _InputIter, class _OutputIter>   
    inline void divide_each(_OutputIter x_begin, _OutputIter x_end, 
        _InputIter y_begin) {
        for(; x_begin != x_end; ++x_begin, ++y_begin) *x_begin /= *y_begin;
    }

    /**
     * Increment each entry in \a x by the corresponding entry in \a y.
     */
    template<class _InputIter, class _OutputIter>   
    inline void sum_each(_OutputIter x_begin, _OutputIter x_end, 
        _InputIter y_begin) {
        for(; x_begin != x_end; ++x_begin, ++y_begin) *x_begin += *y_begin;
    }

    /** Copies a contiguous range in a sequence to indexed
     *  positions in another sequence. Example:
     *
     *  \code
     *  vector<double> x(3), y(20), ;
     *  vector<int> index(3);
     *  index[0] = 9;
     *  index[1] = 2;
     *  index[3] = 16; 
     *  _scatter_copy(x.begin(), x.end(), y.begin(), index.begin());
     *  \endcode
     */
    template<class _InputIter, class _OutputIter, class _IndexIter>
    inline void _scatter_copy(_InputIter __begin, _InputIter __end, 
        _OutputIter __result, _IndexIter __index) {
        for (; __begin != __end; ++__begin, ++__index) {
            *(__result + *__index) = *__begin;
        }
    }

    /**
     * Multiply selected elements in an array by a contiguous 
     * sequence of multipliers. 
     * Example:
     * \code
     * double multipliers[] = {8.9, -2.0, 5.6};
     * int index[] = {7, 4, 13};
     * vector_fp data(20);
     * ...
     * // multiply elements 7, 4, and 13 in data by multipliers
     * _scatter_mult(multipliers, multipliers + 3, data.begin(),
     *               index);
     * \endcode
     */
    
    template<class _InputIter, class _RandAccessIter, class _IndexIter>
    inline void _scatter_mult(_InputIter __mult_begin, _InputIter __mult_end, 
        _RandAccessIter __data, _IndexIter __index) {
	for (; __mult_begin != __mult_end; ++__mult_begin, ++__index) {
            *(__data + *__index) *= *__mult_begin;
	}
    }

    /**
     * Divide selected elements in an array by a contiguous 
     * sequence of divisors.
     * Example: 
     * \code
     * double divisors[] = {8.9, -2.0, 5.6};
     * int index[] = {7, 4, 13};
     * vector_fp data(20);
     * ...
     * // divide elements 7, 4, and 13 in data by divisors
     * _scatter_divide(divisors, divisors + 3, data.begin(),
     *               index);
     * \endcode
     */
    template<class _InputIter, class _OutputIter, class _IndexIter>
    inline void _scatter_divide(_InputIter __begin, _InputIter __end, 
        _OutputIter __result, _IndexIter __index) {
	for (; __begin != __end; ++__begin, ++__index) {
            *(__result + *__index) /= *__begin;
	}
    }

    /**
     * Compute \f[ \sum_k x_k \log x_k. \f]. A small number (1.0E-20)
     * is added before taking the log.
     */ 
    template<class _InputIter>  
    inline doublereal _sum_xlogx(_InputIter __begin, _InputIter __end) {
	doublereal sum = 0.0;
	for (; __begin != __end; ++__begin) {
            sum += (*__begin) * log(*__begin + Tiny);
	}
	return sum;
    }

    /**
     * Compute \f[ \sum_k x_k \log Q_k. \f]. A small number (1.0E-20)
     * is added before taking the log.
     */ 
    template<class _InputIter1, class _InputIter2>  
    inline doublereal _sum_xlogQ(_InputIter1 __begin, _InputIter1 __end,
               _InputIter2 _Q_begin) {
	doublereal sum = 0.0;
	for (; __begin != __end; ++__begin, ++_Q_begin) {
            sum += (*__begin) * log(*_Q_begin + Tiny);
	}
	return sum;
    }

}


#endif

















