/**
 *  @file utilities.h
 *
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_UTILITIES_H
#define CT_UTILITIES_H

#include "ct_defs.h"

#ifdef DARWINNN
#include <Accelerate.h>
#endif

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
    template<class InputIter, class InputIter2>
    inline doublereal dot(InputIter x_begin, InputIter x_end, 
        InputIter2 y_begin) {
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
    template<class InputIter, class OutputIter, class S>
    inline void scale(InputIter begin, InputIter end, 
        OutputIter out, S scale_factor) {
        for (; begin != end; ++begin, ++out) 
            *out = scale_factor * *begin;
    }

    template<class InputIter, class OutputIter, class S>
    inline void increment_scale(InputIter begin, InputIter end, 
        OutputIter out, S scale_factor) {
        for (; begin != end; ++begin, ++out) 
            *out += scale_factor * *begin;
    }

    /**
     * Multiply each entry in x by the corresponding entry in y.
     */
    template<class InputIter, class OutputIter>   
    inline void multiply_each(OutputIter x_begin, OutputIter x_end, 
        InputIter y_begin) {
        for(; x_begin != x_end; ++x_begin, ++y_begin) *x_begin *= *y_begin;
    }


    /**
     * Invoke method 'resize' with argument \a m for a sequence of objects.
     */
    template<class InputIter>
    inline void resize_each(int m, InputIter begin, InputIter end) {
        for(; begin != end; ++begin) begin->resize(m);
    }

    /**
     * The maximum absolute value.
     */
    template<class InputIter>
    inline doublereal absmax(InputIter begin, InputIter end) {
        doublereal amax = 0.0;
        for(; begin != end; ++begin) 
            if (fabs(*begin) > amax) amax = fabs(*begin);
        return amax;
    }

    /**
     * Normalize the values in a sequence, such that they sum to 1.0.
     */
    template<class InputIter, class OutputIter>
    inline void normalize(InputIter begin, InputIter end, 
        OutputIter out) {
        doublereal sum = accumulate(begin, end, 0.0);
        for (; begin != end; ++begin, ++out) *out = *begin/sum;
    }

    /**
     * Divide each element of \a x by the corresponding element of \a y.
     */
    template<class InputIter, class OutputIter>   
    inline void divide_each(OutputIter x_begin, OutputIter x_end, 
        InputIter y_begin) {
        for(; x_begin != x_end; ++x_begin, ++y_begin) *x_begin /= *y_begin;
    }

    /**
     * Increment each entry in \a x by the corresponding entry in \a y.
     */
    template<class InputIter, class OutputIter>   
    inline void sum_each(OutputIter x_begin, OutputIter x_end, 
        InputIter y_begin) {
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
     *  scatter_copy(x.begin(), x.end(), y.begin(), index.begin());
     *  \endcode
     */
    template<class InputIter, class OutputIter, class IndexIter>
    inline void scatter_copy(InputIter begin, InputIter end, 
        OutputIter result, IndexIter index) {
        for (; begin != end; ++begin, ++index) {
            *(result + *index) = *begin;
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
     * scatter_mult(multipliers, multipliers + 3, data.begin(),
     *               index);
     * \endcode
     */
    
    template<class InputIter, class RandAccessIter, class IndexIter>
    inline void scatter_mult(InputIter mult_begin, InputIter mult_end, 
        RandAccessIter data, IndexIter index) {
	for (; mult_begin != mult_end; ++mult_begin, ++index) {
            *(data + *index) *= *mult_begin;
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
     * scatter_divide(divisors, divisors + 3, data.begin(),
     *               index);
     * \endcode
     */
    template<class InputIter, class OutputIter, class IndexIter>
    inline void scatter_divide(InputIter begin, InputIter end, 
        OutputIter result, IndexIter index) {
	for (; begin != end; ++begin, ++index) {
            *(result + *index) /= *begin;
	}
    }

    /**
     * Compute \f[ \sum_k x_k \log x_k. \f]. A small number (1.0E-20)
     * is added before taking the log.
     */ 
    template<class InputIter>  
    inline doublereal sum_xlogx(InputIter begin, InputIter end) {
	doublereal sum = 0.0;
	for (; begin != end; ++begin) {
            sum += (*begin) * log(*begin + Tiny);
	}
	return sum;
    }

    /**
     * Compute \f[ \sum_k x_k \log Q_k. \f]. A small number (1.0E-20)
     * is added before taking the log.
     */ 
    template<class InputIter1, class InputIter2>  
    inline doublereal sum_xlogQ(InputIter1 begin, InputIter1 end,
               InputIter2 Q_begin) {
	doublereal sum = 0.0;
	for (; begin != end; ++begin, ++Q_begin) {
            sum += (*begin) * log(*Q_begin + Tiny);
	}
	return sum;
    }

    template<class OutputIter>
    inline void scale(int N, double alpha, OutputIter x) {
        //#ifdef DARWINNNN
        //cblas_dscal(N, alpha, x, 1);
        //#else
        scale(x, x+N, x, alpha);
        //#endif
    }

}


#endif

















