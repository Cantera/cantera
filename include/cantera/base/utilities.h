/**
 *  @file utilities.h
 *  Various templated functions that carry out common vector
 *  operations (see \ref utils).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

/**
 * @defgroup utils Templated Utility Functions
 *
 * These are templates to perform various simple operations on arrays. Note that
 * the compiler will inline these, so using them carries no performance penalty.
 */

#ifndef CT_UTILITIES_H
#define CT_UTILITIES_H

#include "ct_defs.h"
#include <numeric>

namespace Cantera
{
//! Templated Inner product of two vectors of length 4.
/*!
 * If either \a x or \a y has length greater than 4, only the first 4 elements
 * will be used.
 *
 * @param x   first reference to the templated class V
 * @param y   second reference to the templated class V
 * @return This class returns a hard-coded type, doublereal.
 */
template<class V>
inline doublereal dot4(const V& x, const V& y)
{
    return x[0]*y[0] + x[1]*y[1] + x[2]*y[2] + x[3]*y[3];
}

//! Templated Inner product of two vectors of length 5
/*!
 * If either \a x or \a y has length greater than 4, only the first 4 elements
 * will be used.
 *
 * @param x   first reference to the templated class V
 * @param y   second reference to the templated class V
 * @return This class returns a hard-coded type, doublereal.
 */
template<class V>
inline doublereal dot5(const V& x, const V& y)
{
    return x[0]*y[0] + x[1]*y[1] + x[2]*y[2] + x[3]*y[3] +
           x[4]*y[4];
}

//! Function that calculates a templated inner product.
/*!
 * This inner product is templated twice. The output variable is hard coded
 * to return a doublereal.
 *
 * template<class InputIter, class InputIter2>
 *
 * @code
 *     double x[8], y[8];
 *     doublereal dsum = dot<double *,double *>(x, &x+7, y);
 * @endcode
 *
 * @param x_begin  Iterator pointing to the beginning, belonging to the
 *                 iterator class InputIter.
 * @param x_end    Iterator pointing to the end, belonging to the
 *                 iterator class InputIter.
 * @param y_begin Iterator pointing to the beginning of y, belonging to the
 *               iterator class InputIter2.
 * @return The return is hard-coded to return a double.
 */
template<class InputIter, class InputIter2>
inline doublereal dot(InputIter x_begin, InputIter x_end,
                      InputIter2 y_begin)
{
    return std::inner_product(x_begin, x_end, y_begin, 0.0);
}

//! Multiply elements of an array by a scale factor.
/*!
 * \code
 * vector_fp in(8, 1.0), out(8);
 * scale(in.begin(), in.end(), out.begin(), factor);
 * \endcode
 *
 * @param begin  Iterator pointing to the beginning, belonging to the
 *               iterator class InputIter.
 * @param end    Iterator pointing to the end, belonging to the
 *               iterator class InputIter.
 * @param out    Iterator pointing to the beginning of out, belonging to the
 *               iterator class OutputIter. This is the output variable
 *               for this routine.
 * @param scale_factor  input scale factor belonging to the class S.
 */
template<class InputIter, class OutputIter, class S>
inline void scale(InputIter begin, InputIter end,
                  OutputIter out, S scale_factor)
{
    std::transform(begin, end, out,
        [scale_factor](double x) { return x * scale_factor; });
}

//! Templated evaluation of a polynomial of order 6
/*!
 * @param x   Value of the independent variable - First template parameter
 * @param c   Pointer to the polynomial - Second template parameter
 */
template<class D, class R>
R poly6(D x, R* c)
{
    return ((((((c[6]*x + c[5])*x + c[4])*x + c[3])*x +
              c[2])*x + c[1])*x + c[0]);
}

//! Templated evaluation of a polynomial of order 8
/*!
 * @param x   Value of the independent variable - First template parameter
 * @param c   Pointer to the polynomial - Second template parameter
 */
template<class D, class R>
R poly8(D x, R* c)
{
    return ((((((((c[8]*x + c[7])*x + c[6])*x + c[5])*x + c[4])*x + c[3])*x +
              c[2])*x + c[1])*x + c[0]);
}

//! Templated evaluation of a polynomial of order 5
/*!
 * @param x   Value of the independent variable - First template parameter
 * @param c   Pointer to the polynomial - Second template parameter
 */
template<class D, class R>
R poly5(D x, R* c)
{
    return (((((c[5]*x + c[4])*x + c[3])*x +
              c[2])*x + c[1])*x + c[0]);
}

//! Evaluates a polynomial of order 4.
/*!
 * @param x   Value of the independent variable.
 * @param c   Pointer to the polynomial coefficient array.
 */
template<class D, class R>
R poly4(D x, R* c)
{
    return ((((c[4]*x + c[3])*x +
              c[2])*x + c[1])*x + c[0]);
}

//! Templated evaluation of a polynomial of order 3
/*!
 * @param x   Value of the independent variable - First template parameter
 * @param c   Pointer to the polynomial - Second template parameter
 */
template<class D, class R>
R poly3(D x, R* c)
{
    return (((c[3]*x + c[2])*x + c[1])*x + c[0]);
}

//! Check to see that a number is finite (not NaN, +Inf or -Inf)
void checkFinite(const double tmp);

//! Check to see that all elements in an array are finite
/*!
 * Throws an exception if any element is NaN, +Inf, or -Inf
 * @param name    Name to be used in the exception message if the check fails
 * @param values  Array of *N* values to be checked
 * @param N       Number of elements in *values*
 */
void checkFinite(const std::string& name, double* values, size_t N);

//! Const accessor for a value in a std::map.
/*
 * Similar to std::map.at(key), but returns *default_val* if the key is not
 * found instead of throwing an exception.
 */
template <class T, class U>
const U& getValue(const std::map<T, U>& m, const T& key, const U& default_val) {
    typename std::map<T,U>::const_iterator iter = m.find(key);
    return (iter == m.end()) ? default_val : iter->second;
}

//! A macro for generating member function detectors, which can then be used in
//! combination with `std::enable_if` to allow selection of a particular template
//! specialization based on the presence of that member function. See MultiRate for
//! examples of use.
#define CT_DEFINE_HAS_MEMBER(detector_name, func_name)                   \
    template <typename T>                                                \
    struct detector_name {                                               \
        typedef char (& yes)[1], (& no)[2];                              \
        template <typename C> static yes check(decltype(&C::func_name)); \
        template <typename> static no check(...);                        \
        static bool const value = sizeof(check<T>(0)) == sizeof(yes);    \
    };

}

#endif
