/**
 *  @file vec_functions.h
 *  Templates for operations on vector-like objects.
 */
/*
 *  Copyright 2001 California Institute of Technology
 */

#ifndef CT_VEC_FUNCTIONS_H
#define CT_VEC_FUNCTIONS_H

#include "ct_defs.h"
#include "utilities.h"
#include <functional>
#include <iostream>
#include <cstring>

namespace Cantera
{
//! Templated function that copies the first n entries from x to y.
/*!
 *
 *
 *  The  templated type is the type of x and y
 *
 *  @param n        Number of elements to copy from x to y
 *  @param x        The object x, of templated type const T&
 *  @param y        The object y, of templated type  T&
 *  @deprecated Unused. To be removed after Cantera 2.2.
 */
template<class T>
inline void copyn(size_t n, const T& x, T& y)
{
    warn_deprecated("copyn", "To be removed after Cantera 2.2.");
    std::copy(x.begin(), x.begin() + n, y.begin());
}

//!  Divide each element of x by the corresponding element of y.
/*!
 * This function replaces x[n] by x[n]/y[n], for 0 <= n < x.size()
 *
 * @param x  Numerator object of the division operation with template type T
 *           At the end of the calculation, it contains the result.
 * @param y  Denominator object of the division template type T
 * @deprecated Unused. To be removed after Cantera 2.2.
 */
template<class T>
inline void divide_each(T& x, const T& y)
{
    warn_deprecated("divide_each", "To be removed after Cantera 2.2.");
    std::transform(x.begin(), x.end(), y.begin(),
                   x.begin(), std::divides<typename T::value_type>());
}

//! Multiply each element of x by the corresponding element of y.
/*!
 * This function replaces x[n] by x[n]*y[n], for 0 <= n < x.size()
 * This is a templated function with just one template type.
 *
 * @param x  First object of the multiplication with template type T
 *           At the end of the calculation, it contains the result.
 * @param y  Second object of the multiplication with template type T
 * @deprecated Unused. To be removed after Cantera 2.2.
 */
template<class T>
inline void multiply_each(T& x, const T& y)
{
    warn_deprecated("multiply_each", "To be removed after Cantera 2.2.");
    std::transform(x.begin(), x.end(), y.begin(),
                   x.begin(), std::multiplies<typename T::value_type>());
}

//! Multiply each element of x by scale_factor.
/*!
 * This function replaces x[n] by x[n]*scale_factor, for 0 <= n < x.size()
 *
 * @param x  First object of the multiplication with template type T
 *           At the end of the calculation, it contains the result.
 * @param scale_factor scale factor with template type S
 * @deprecated Unused. To be removed after Cantera 2.2.
 */
template<class T, class S>
inline void scale(T& x, S scale_factor)
{
    warn_deprecated("scale", "To be removed after Cantera 2.2.");
    scale(x.begin(), x.end(), x.begin(), scale_factor);
}

//! Return the templated dot product of two objects
/*!
 * Returns the sum of x[n]*y[n], for 0 <= n < x.size().
 *
 * @param x  First object of the dot product with template type T
 *           At the end of the calculation, it contains the result.
 * @param y  Second object of the dot product with template type T
 * @deprecated Unused. To be removed after Cantera 2.2.
 */
template<class T>
inline doublereal dot_product(const T& x, const T& y)
{
    warn_deprecated("dot_product", "To be removed after Cantera 2.2.");
    return std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
}

//! Returns the templated dot ratio of two objects
/**
 * Returns the sum of x[n]/y[n], for 0 <= n < x.size().
 *
 * @param x  First object of the dot product with template type T
 *           At the end of the calculation, it contains the result.
 * @param y  Second object of the dot product with template type T
 * @deprecated Unused. To be removed after Cantera 2.2.
 */
template<class T>
inline doublereal dot_ratio(const T& x, const T& y)
{
    warn_deprecated("dot_ratio", "To be removed after Cantera 2.2.");
    return _dot_ratio(x.begin(), x.end(), y.begin(), 0.0);
}

//! Returns a templated addition operation of two objects
/**
 * Replaces x[n] by x[n] + y[n] for 0 <= n < x.size()
 *
 * @param x  First object of the addition with template type T
 *           At the end of the calculation, it contains the result.
 * @param y  Second object of the addition with template type T
 * @deprecated Unused. To be removed after Cantera 2.2.
 */
template<class T>
inline void add_each(T& x, const T& y)
{
    warn_deprecated("add_each", "To be removed after Cantera 2.2.");
    std::transform(x.begin(), x.end(), y.begin(),
                   x.begin(), std::plus<typename T::value_type>());
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
 * @deprecated Unused. To be removed after Cantera 2.2.
 */
template<class InputIter, class S>
inline doublereal _dot_ratio(InputIter x_begin, InputIter x_end,
                             InputIter y_begin, S start_value)
{
    warn_deprecated("_dot_ratio", "To be removed after Cantera 2.2.");
    for (; x_begin != x_end; ++x_begin, ++y_begin) {
        start_value += *x_begin / *y_begin;
    }
    return start_value;
}

//! Finds the entry in a vector with maximum absolute
//! value, and return this value.
/*!
 *  @param v Vector to be queried for maximum value, with template type T
 *
 * @return Returns an object of type T that is the maximum value,
 * @deprecated Unused. To be removed after Cantera 2.2.
 */
template<class T>
inline T absmax(const std::vector<T>& v)
{
    warn_deprecated("absmax", "To be removed after Cantera 2.2.");
    int n = v.size();
    T maxval = 0.0;
    for (int i = 0; i < n; i++) {
        maxval = std::max(std::abs(v[i]), maxval);
    }
    return maxval;
}

//! Write a vector to a stream
template <class T>
inline std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
{
    size_t n = v.size();
    for (size_t i = 0; i < n; i++) {
        os << v[i];
        if (i != n-1) {
            os << ", ";
        }
    }
    return os;
}

}

#endif
