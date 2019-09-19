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
#include "global.h"
#include <stdexcept>

#include <numeric>

namespace Cantera
{
//! Unary operator to multiply the argument by a constant.
/*!
 * The form of this operator is designed for use by std::transform.
 * @see @ref  scale().
 *
 * @deprecated To be removed after Cantera 2.5. Replaceable with C++11 lambda.
 */
template<class T> struct timesConstant : public std::unary_function<T, double> {
    //! Constructor
    /*!
     * @param c  Constant of templated type T that will be stored internally
     *           within the object and used in the multiplication operation
     */
    timesConstant(T c) : m_c(c) {
        warn_deprecated("class timesConstant",
            "To be removed after Cantera 2.5. Replaceable with C++11 lambda.");
    }

    //! Parenthesis operator returning a double
    /*!
     * @param x  Variable of templated type T that will be used in the
     *           multiplication operator
     * @returns a value of type double from the internal multiplication
     */
    double operator()(T x) {
        return m_c * x;
    }

    //! Stored constant value of time T
    T m_c;
};

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

//! Multiply each entry in x by the corresponding entry in y.
/*!
 * The template arguments are:  template<class InputIter, class OutputIter>
 *
 * Simple code Equivalent:
 *  \code
 *   double x[10], y[10]
 *   for (n = 0; n < 10; n++) {
 *     x[n] *= y[n];
 *   }
 * \endcode
 * Example of function call usage to implement the simple code example:
 *  \code
 *     double x[10], y[10]
 *     multiply_each(x, x+10, y);
 *  \endcode
 *
 * @param x_begin   Iterator pointing to the beginning of the vector x,
 *                  belonging to the iterator class InputIter.
 * @param x_end     Iterator pointing to the end of the vector x, belonging to
 *                  the iterator class InputIter. The difference between end and
 *                  begin determines the loop length
 * @param y_begin   Iterator pointing to the beginning of the vector y,
 *                  belonging to the iterator class outputIter.
 * @deprecated      Unused. To be removed after Cantera 2.5.
 */
template<class InputIter, class OutputIter>
inline void multiply_each(OutputIter x_begin, OutputIter x_end,
                          InputIter y_begin)
{
    warn_deprecated("multiply_each", "To be removed after Cantera 2.5.");
    for (; x_begin != x_end; ++x_begin, ++y_begin) {
        *x_begin *= *y_begin;
    }
}

//! The maximum absolute value (templated version)
/*!
 * The template arguments are:  template<class InputIter>
 *
 * Simple code Equivalent:
 *  \code
 *   double x[10] amax = 0.0;
 *   for (int n = 0; n < 10; n++) {
 *    if (fabs(x[n]) > amax) amax = fabs(x[10]);
 *   }
 *   return amax;
 * \endcode
 * Example of function call usage to implement the simple code example:
 *  \code
 *   double x[10]
 *   double amax = absmax(x, x+10);
 *  \endcode
 *
 * @param begin     Iterator pointing to the beginning of the x vector,
 *                  belonging to the iterator class InputIter.
 * @param end       Iterator pointing to the end of the x vector, belonging to
 *                  the iterator class InputIter. The difference between end and
 *                  begin determines the loop length
 * @deprecated      Unused. To be removed after Cantera 2.5.
 */
template<class InputIter>
inline doublereal absmax(InputIter begin, InputIter end)
{
    warn_deprecated("absmax", "To be removed after Cantera 2.5.");
    doublereal amax = 0.0;
    for (; begin != end; ++begin) {
        amax = std::max(fabs(*begin), amax);
    }
    return amax;
}

//! Normalize the values in a sequence, such that they sum to 1.0 (templated
//! version)
/*!
 * The template arguments are:  template<class InputIter, class OutputIter>
 *
 * Simple Equivalent:
 *  \code
 *   double x[10], y[10], sum = 0.0;
 *   for (int n = 0; n < 10; n++) {
 *    sum += x[10];
 *   }
 *   for (int n = 0; n < 10; n++) {
 *    y[n] = x[n]/sum;
 *   }
 * \endcode
 * Example of function call usage:
 *  \code
 *   double x[10], y[10];
 *   normalize(x, x+10, y);
 *  \endcode
 *
 * @param begin     Iterator pointing to the beginning of the x vector,
 *                  belonging to the iterator class InputIter.
 * @param end       Iterator pointing to the end of the x vector, belonging to
 *                  the iterator class InputIter. The difference between end and
 *                  begin determines the loop length
 * @param out       Iterator pointing to the beginning of the output vector,
 *                  belonging to the iterator class OutputIter.
 * @deprecated      Unused. To be removed after Cantera 2.5.
 */
template<class InputIter, class OutputIter>
inline void normalize(InputIter begin, InputIter end,
                      OutputIter out)
{
    warn_deprecated("normalize", "To be removed after Cantera 2.5.");
    doublereal sum = std::accumulate(begin, end, 0.0);
    for (; begin != end; ++begin, ++out) {
        *out = *begin/sum;
    }
}

//! Templated divide of each element of \a x by the corresponding element of \a y.
/*!
 * The template arguments are:  template<class InputIter, class OutputIter>
 *
 * Simple Equivalent:
 *  \code
 *   double x[10], y[10];
 *   for (n = 0; n < 10; n++) {
 *     x[n] /= y[n];
 *   }
 * \endcode
 * Example of code usage:
 *  \code
 *   double x[10], y[10];
 *   divide_each(x, x+10, y);
 *  \endcode
 *
 * @param x_begin   Iterator pointing to the beginning of the x vector,
 *                  belonging to the iterator class OutputIter.
 * @param x_end     Iterator pointing to the end of the x vector, belonging to
 *                  the iterator class OutputIter. The difference between end
 *                  and begin determines the number of inner iterations.
 * @param y_begin   Iterator pointing to the beginning of the yvector, belonging
 *                  to the iterator class InputIter.
 * @deprecated      Unused. To be removed after Cantera 2.5.
 */
template<class InputIter, class OutputIter>
inline void divide_each(OutputIter x_begin, OutputIter x_end,
                        InputIter y_begin)
{
    warn_deprecated("divide_each", "To be removed after Cantera 2.5.");
    for (; x_begin != x_end; ++x_begin, ++y_begin) {
        *x_begin /= *y_begin;
    }
}

//! Increment each entry in \a x by the corresponding entry in \a y.
/*!
 * The template arguments are:  template<class InputIter, class OutputIter>
 *
 * @param x_begin   Iterator pointing to the beginning of the x vector,
 *                  belonging to the iterator class OutputIter.
 * @param x_end     Iterator pointing to the end of the x vector, belonging to
 *                  the iterator class OutputIter. The difference between end
 *                  and begin determines the number of inner iterations.
 * @param y_begin   Iterator pointing to the beginning of the yvector, belonging
 *                  to the iterator class InputIter.
 * @deprecated      Unused. To be removed after Cantera 2.5.
 */
template<class InputIter, class OutputIter>
inline void sum_each(OutputIter x_begin, OutputIter x_end,
                     InputIter y_begin)
{
    warn_deprecated("sum_each", "To be removed after Cantera 2.5.");
    for (; x_begin != x_end; ++x_begin, ++y_begin) {
        *x_begin += *y_begin;
    }
}

//!   Copies a contiguous range in a sequence to indexed
//!   positions in another sequence.
/*!
 * The template arguments are:  template<class InputIter, class OutputIter, class IndexIter>
 *
 *  Example:
 *
 *  \code
 *  vector_fp x(3), y(20);
 *  vector_int index(3);
 *  index[0] = 9;
 *  index[1] = 2;
 *  index[3] = 16;
 *  scatter_copy(x.begin(), x.end(), y.begin(), index.begin());
 *  \endcode
 *
 *  This routine is templated 3 times.
 *      InputIter is an iterator for the source vector
 *      OutputIter is an iterator for the destination vector
 *      IndexIter is an iterator for the index into the destination vector.
 *
 * @param begin   Iterator pointing to the beginning of the source vector,
 *                belonging to the iterator class InputIter.
 * @param end     Iterator pointing to the end of the source vector, belonging
 *                to the iterator class InputIter. The difference between end
 *                and begin determines the number of inner iterations.
 * @param result  Iterator pointing to the beginning of the output vector,
 *                belonging to the iterator class outputIter.
 * @param index   Iterator pointing to the beginning of the index vector, belonging to the
 *                iterator class IndexIter.
 * @deprecated      Unused. To be removed after Cantera 2.5.
 */
template<class InputIter, class OutputIter, class IndexIter>
inline void scatter_copy(InputIter begin, InputIter end,
                         OutputIter result, IndexIter index)
{
    warn_deprecated("scatter_copy", "To be removed after Cantera 2.5.");
    for (; begin != end; ++begin, ++index) {
        *(result + *index) = *begin;
    }
}

//! Multiply selected elements in an array by a contiguous sequence of
//! multipliers.
/*!
 * The template arguments are:  template<class InputIter, class RandAccessIter, class IndexIter>
 *
 * Example:
 * \code
 * double multipliers[] = {8.9, -2.0, 5.6};
 * int index[] = {7, 4, 13};
 * vector_fp data(20);
 * ...
 * // Multiply elements 7, 4, and 13 in data by multipliers[0], multipliers[1],and multipliers[2],
 * // respectively
 * scatter_mult(multipliers, multipliers + 3, data.begin(), index);
 * \endcode
 *
 * @param mult_begin   Iterator pointing to the beginning of the multiplier
 *                     vector, belonging to the iterator class InputIter.
 * @param mult_end     Iterator pointing to the end of the multiplier vector,
 *                     belonging to the iterator class InputIter. The difference
 *                     between end and begin determines the number of inner
 *                     iterations.
 * @param data         Iterator pointing to the beginning of the output vector,
 *                     belonging to the iterator class RandAccessIter, that will
 *                     be selectively multiplied.
 * @param index        Iterator pointing to the beginning of the index vector,
 *                     belonging to the iterator class IndexIter.
 * @deprecated      Unused. To be removed after Cantera 2.5.
 */
template<class InputIter, class RandAccessIter, class IndexIter>
inline void scatter_mult(InputIter mult_begin, InputIter mult_end,
                         RandAccessIter data, IndexIter index)
{
    warn_deprecated("scatter_mult", "To be removed after Cantera 2.5.");
    for (; mult_begin != mult_end; ++mult_begin, ++index) {
        *(data + *index) *= *mult_begin;
    }
}

//!  Compute \f[ \sum_k x_k \log x_k. \f].
/*!
 * The template arguments are:  template<class InputIter>
 *
 * A small number (1.0E-20) is added before taking the log. This templated
 * class does the indicated sum. The template must be an iterator.
 *
 * @param begin  Iterator pointing to the beginning, belonging to the
 *               iterator class InputIter.
 * @param end    Iterator pointing to the end, belonging to the
 *               iterator class InputIter.
 * @return The return from this class is a double.
 * @deprecated      Unused. To be removed after Cantera 2.5.
 */
template<class InputIter>
inline doublereal sum_xlogx(InputIter begin, InputIter end)
{
    warn_deprecated("sum_xlogx", "To be removed after Cantera 2.5.");
    doublereal sum = 0.0;
    for (; begin != end; ++begin) {
        sum += (*begin) * std::log(*begin + Tiny);
    }
    return sum;
}

//! Compute \f[ \sum_k x_k \log Q_k. \f].
/*!
 * The template arguments are:  template<class InputIter1, class InputIter2>
 *
 * This class is templated twice. The first template, InputIter1 is the iterator
 * that points to $x_k$. The second iterator InputIter2, point to $Q_k$. A small
 * number (1.0E-20) is added before taking the log.
 *
 * @param begin  Iterator pointing to the beginning, belonging to the
 *               iterator class InputIter1.
 * @param end    Iterator pointing to the end, belonging to the
 *               iterator class InputIter1.
 * @param Q_begin Iterator pointing to the beginning of Q_k, belonging to the
 *               iterator class InputIter2.
 * @return The return from this class is hard coded to a doublereal.
 * @deprecated      Unused. To be removed after Cantera 2.5.
 */
template<class InputIter1, class InputIter2>
inline doublereal sum_xlogQ(InputIter1 begin, InputIter1 end,
                            InputIter2 Q_begin)
{
    warn_deprecated("sum_xlogQ", "To be removed after Cantera 2.5.");
    doublereal sum = 0.0;
    for (; begin != end; ++begin, ++Q_begin) {
        sum += (*begin) * std::log(*Q_begin + Tiny);
    }
    return sum;
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

//@}

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

}

#endif
