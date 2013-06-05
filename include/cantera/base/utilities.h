/**
 *  @file utilities.h
 *  Various templated functions that carry out common vector
 *  operations (see \ref utils).
 */

// Copyright 2001  California Institute of Technology
/**
 * @defgroup utils Templated Utility Functions
 *
 * These are templates to perform various simple operations on arrays.
 * Note that the compiler will inline these, so using them carries no
 * performance penalty.
 */

#ifndef CT_UTILITIES_H
#define CT_UTILITIES_H

#include "ct_defs.h"

#include <algorithm>

namespace Cantera
{
//! Unary operator to multiply the argument by a constant.
/*!
 *  The form of this operator is designed for use by std::transform.
 *  @see @ref  scale().
 */
template<class T> struct timesConstant : public std::unary_function<T, double> {
    //! Constructor
    /*!
     * @param c  Constant of templated type T
     *           that will be stored internally within the object
     *           and used in the multiplication operation
     */
    timesConstant(T c) : m_c(c) {}

    //! Parenthesis operator returning a double
    /*!
     * @param x  Variable of templated type T that will be
     *           used in the multiplication operator
     *
     * @return Returns a value of type double from the internal
     *         multiplication
     */
    double operator()(T x) {
        return m_c * x;
    }

    //! Stored constant value of time T
    T m_c;
};

//!  Templated Inner product of two vectors of length 4.
/*!
 * If either \a x
 * or \a y has length greater than 4, only the first 4 elements
 * will be used.
 *
 * @param x   first reference to the templated class V
 * @param y   second reference to the templated class V
 * @return
 *      This class returns a hard-coded type, doublereal.
 */
template<class V>
inline doublereal dot4(const V& x, const V& y)
{
    return x[0]*y[0] + x[1]*y[1] + x[2]*y[2] + x[3]*y[3];
}


//!  Templated Inner product of two vectors of length 5
/*!
 * If either \a x
 * or \a y has length greater than 4, only the first 4 elements
 * will be used.
 *
 * @param x   first reference to the templated class V
 * @param y   second reference to the templated class V
 * @return
 *      This class returns a hard-coded type, doublereal.
 */
template<class V>
inline doublereal dot5(const V& x, const V& y)
{
    return x[0]*y[0] + x[1]*y[1] + x[2]*y[2] + x[3]*y[3] +
           x[4]*y[4];
}

//!  Templated Inner product of two vectors of length 6
/*!
 * If either \a x
 * or \a y has length greater than 4, only the first 4 elements
 * will be used.
 *
 * @param x   first reference to the templated class V
 * @param y   second reference to the templated class V
 * @return
 *      This class returns a hard-coded type, doublereal.
 */
template<class V>
inline doublereal dot6(const V& x, const V& y)
{
    return x[0]*y[0] + x[1]*y[1] + x[2]*y[2] + x[3]*y[3] +
           x[4]*y[4] + x[5]*y[5];
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
 * @return
 *     The return is hard-coded to return a double.
 */
template<class InputIter, class InputIter2>
inline doublereal dot(InputIter x_begin, InputIter x_end,
                      InputIter2 y_begin)
{
    return inner_product(x_begin, x_end, y_begin, 0.0);
}

//!   Multiply elements of an array by a scale factor.
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
    std::transform(begin, end, out, timesConstant<S>(scale_factor));
}

/*!
 * Multiply elements of an array, y, by a scale factor, f and add the
 * result to an existing array, x. This is essentially a templated daxpy_
 * operation.
 *
 * The template arguments are:  template<class InputIter,
 * class OutputIter, class S>
 *
 *  Simple Code Example of the functionality;
 * @code
 *       double x[10], y[10], f;
 *       for (i = 0; i < n; i++) {
 *         y[i] += f * x[i]
 *       }
 * @endcode
 *  Example of the function call to implement the simple code example
 *   @code
 *      double x[10], y[10], f;
 *      increment_scale(x, x+10, y, f);
 *   @endcode
 *
 * It is templated with three parameters. The first template
 * is the iterator, InputIter, which controls access to y[].
 * The second template is the iterator OutputIter, which controls
 * access to y[]. The third iterator is S, which is f.
 *
 * @param begin InputIter Iterator for beginning of y[]
 * @param end   inputIter Iterator for end of y[]
 * @param out   OutputIter Iterator for beginning of x[]
 * @param scale_factor Scale Factor to multiply y[i] by
 */
template<class InputIter, class OutputIter, class S>
inline void increment_scale(InputIter begin, InputIter end,
                            OutputIter out, S scale_factor)
{
    for (; begin != end; ++begin, ++out) {
        *out += scale_factor * *begin;
    }
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
 * @param x_begin   Iterator pointing to the beginning of the vector x, belonging to the
 *                  iterator class InputIter.
 * @param x_end     Iterator pointing to the end of the vector x, belonging to the
 *                  iterator class InputIter. The difference between end and begin
 *                  determines the loop length
 * @param y_begin   Iterator pointing to the beginning of the vector y, belonging to the
 *                  iterator class outputIter.
 */
template<class InputIter, class OutputIter>
inline void multiply_each(OutputIter x_begin, OutputIter x_end,
                          InputIter y_begin)
{
    for (; x_begin != x_end; ++x_begin, ++y_begin) {
        *x_begin *= *y_begin;
    }
}

//! Invoke method 'resize' with argument \a m for a sequence of objects (templated version)
/*!
 * The template arguments are:  template<class InputIter>
 *
 * Simple code Equivalent:
 *  \code
 *   vector<vector<double> *> VV;
 *   for (n = 0; n < 20; n++) {
 *     vector<double> *vp = VV[n];
 *     vp->resize(m);
 *   }
 * \endcode
 * Example of function call usage to implement the simple code example:
 *  \code
 *    vector<vector<double> *> VV;
 *    resize_each(m, &VV[0], &VV[20]);
 *  \endcode
 *
 * @param m         Integer specifying the size that each object should be resized to.
 * @param begin     Iterator pointing to the beginning of the sequence of object, belonging to the
 *                  iterator class InputIter.
 * @param end       Iterator pointing to the end of the sequence of objects, belonging to the
 *                  iterator class InputIter. The difference between end and begin
 *                  determines the loop length
 *
 * @note This is currently unused.
 */
template<class InputIter>
inline void resize_each(int m, InputIter begin, InputIter end)
{
    for (; begin != end; ++begin) {
        begin->resize(m);
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
 * @param begin     Iterator pointing to the beginning of the x vector, belonging to the
 *                  iterator class InputIter.
 * @param end       Iterator pointing to the end of the x vector, belonging to the
 *                  iterator class InputIter. The difference between end and begin
 *                  determines the loop length
 */
template<class InputIter>
inline doublereal absmax(InputIter begin, InputIter end)
{
    doublereal amax = 0.0;
    for (; begin != end; ++begin)
        if (fabs(*begin) > amax) {
            amax = fabs(*begin);
        }
    return amax;
}

//! Normalize the values in a sequence, such that they sum to 1.0 (templated version)
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
 * @param begin     Iterator pointing to the beginning of the x vector, belonging to the
 *                  iterator class InputIter.
 * @param end       Iterator pointing to the end of the x vector, belonging to the
 *                  iterator class InputIter. The difference between end and begin
 *                  determines the loop length
 * @param out       Iterator pointing to the beginning of the output vector, belonging to the
 *                  iterator class OutputIter.
 */
template<class InputIter, class OutputIter>
inline void normalize(InputIter begin, InputIter end,
                      OutputIter out)
{
    doublereal sum = accumulate(begin, end, 0.0);
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
 * @param x_begin   Iterator pointing to the beginning of the x vector, belonging to the
 *                  iterator class OutputIter.
 * @param x_end     Iterator pointing to the end of the x vector, belonging to the
 *                  iterator class OutputIter. The difference between end and begin
 *                  determines the number of inner iterations.
 * @param y_begin   Iterator pointing to the beginning of the yvector, belonging to the
 *                  iterator class InputIter.
 */
template<class InputIter, class OutputIter>
inline void divide_each(OutputIter x_begin, OutputIter x_end,
                        InputIter y_begin)
{
    for (; x_begin != x_end; ++x_begin, ++y_begin) {
        *x_begin /= *y_begin;
    }
}

//!  Increment each entry in \a x by the corresponding entry in \a y.
/*!
 * The template arguments are:  template<class InputIter, class OutputIter>
 *
 * @param x_begin   Iterator pointing to the beginning of the x vector, belonging to the
 *                  iterator class OutputIter.
 * @param x_end     Iterator pointing to the end of the x vector, belonging to the
 *                  iterator class OutputIter. The difference between end and begin
 *                  determines the number of inner iterations.
 * @param y_begin   Iterator pointing to the beginning of the yvector, belonging to the
 *                  iterator class InputIter.
 */
template<class InputIter, class OutputIter>
inline void sum_each(OutputIter x_begin, OutputIter x_end,
                     InputIter y_begin)
{
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
 *  vector<double> x(3), y(20), ;
 *  vector<int> index(3);
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
 * @param begin   Iterator pointing to the beginning of the source vector, belonging to the
 *                iterator class InputIter.
 * @param end     Iterator pointing to the end of the source vector, belonging to the
 *                iterator class InputIter. The difference between end and begin
 *                determines the number of inner iterations.
 * @param result  Iterator pointing to the beginning of the output vector, belonging to the
 *                iterator class outputIter.
 * @param index   Iterator pointing to the beginning of the index vector, belonging to the
 *                iterator class IndexIter.
 */
template<class InputIter, class OutputIter, class IndexIter>
inline void scatter_copy(InputIter begin, InputIter end,
                         OutputIter result, IndexIter index)
{
    for (; begin != end; ++begin, ++index) {
        *(result + *index) = *begin;
    }
}


//! Multiply selected elements in an array by a contiguous
//! sequence of multipliers.
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
 * @param mult_begin   Iterator pointing to the beginning of the multiplier vector, belonging to the
 *                     iterator class InputIter.
 * @param mult_end     Iterator pointing to the end of the multiplier vector, belonging to the
 *                     iterator class InputIter. The difference between end and begin
 *                     determines the number of inner iterations.
 * @param data         Iterator pointing to the beginning of the output vector, belonging to the
 *                     iterator class RandAccessIter, that will be selectively multiplied.
 * @param index        Iterator pointing to the beginning of the index vector, belonging to the
 *                     iterator class IndexIter.
 */
template<class InputIter, class RandAccessIter, class IndexIter>
inline void scatter_mult(InputIter mult_begin, InputIter mult_end,
                         RandAccessIter data, IndexIter index)
{
    for (; mult_begin != mult_end; ++mult_begin, ++index) {
        *(data + *index) *= *mult_begin;
    }
}


//! Divide selected elements in an array by a contiguous sequence of divisors.
/*!
 * The template arguments are:  template<class InputIter, class OutputIter, class IndexIter>
 *
 * Example:
 * \code
 * double divisors[] = {8.9, -2.0, 5.6};
 * int index[] = {7, 4, 13};
 * vector_fp data(20);
 * ...
 * // divide elements 7, 4, and 13 in data by divisors[7] divisors[4], and divisors[13]
 * // respectively
 * scatter_divide(divisors, divisors + 3, data.begin(), index);
 * \endcode
 *
 * @param begin   Iterator pointing to the beginning of the source vector, belonging to the
 *                iterator class InputIter.
 * @param end     Iterator pointing to the end of the source vector, belonging to the
 *                iterator class InputIter. The difference between end and begin
 *                determines the number of inner iterations.
 * @param result  Iterator pointing to the beginning of the output vector, belonging to the
 *                iterator class outputIter.
 * @param index   Iterator pointing to the beginning of the index vector, belonging to the
 *                iterator class IndexIter.
 */
template<class InputIter, class OutputIter, class IndexIter>
inline void scatter_divide(InputIter begin, InputIter end,
                           OutputIter result, IndexIter index)
{
    for (; begin != end; ++begin, ++index) {
        *(result + *index) /= *begin;
    }
}

//!  Compute \f[ \sum_k x_k \log x_k. \f].
/*!
 * The template arguments are:  template<class InputIter>
 *
 *  A small number (1.0E-20) is added before taking the log. This templated
 *  class does the indicated sun. The template must be an iterator.
 *
 * @param begin  Iterator pointing to the beginning, belonging to the
 *               iterator class InputIter.
 * @param end    Iterator pointing to the end, belonging to the
 *               iterator class InputIter.
 * @return
 *      The return from this class is a double.
 */
template<class InputIter>
inline doublereal sum_xlogx(InputIter begin, InputIter end)
{
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
 * This class is templated twice. The first template, InputIter1
 * is the iterator that points to $x_k$. The second iterator
 * InputIter2, point to $Q_k$.
 *  A small number (1.0E-20) is added before taking the log.
 *
 * @param begin  Iterator pointing to the beginning, belonging to the
 *               iterator class InputIter1.
 * @param end    Iterator pointing to the end, belonging to the
 *               iterator class InputIter1.
 * @param Q_begin Iterator pointing to the beginning of Q_k, belonging to the
 *               iterator class InputIter2.
 * @return
 *      The return from this class is hard coded to a doublereal.
 */
template<class InputIter1, class InputIter2>
inline doublereal sum_xlogQ(InputIter1 begin, InputIter1 end,
                            InputIter2 Q_begin)
{
    doublereal sum = 0.0;
    for (; begin != end; ++begin, ++Q_begin) {
        sum += (*begin) * std::log(*Q_begin + Tiny);
    }
    return sum;
}

//!   Scale a templated vector by a constant factor.
/*!
 *   The template arguments are:  template<class OutputIter>
 *
 * This function is essentially a wrapper around the stl
 * function %scale(). The function is has one template
 * parameter, OutputIter. OutputIter is a templated iterator
 * that points to the vector to be scaled.
 *
 * @param N       Length of the vector
 * @param alpha   scale factor - double
 * @param x       Templated Iterator to the start of the vector
 *                to be scaled.
 */
template<class OutputIter>
inline void scale(int N, double alpha, OutputIter x)
{
    scale(x, x+N, x, alpha);
}


//! Templated evaluation of a polynomial of order 6
/*!
 *  @param x   Value of the independent variable - First template parameter
 *  @param c   Pointer to the polynomial - Second template parameter
 */
template<class D, class R>
R poly6(D x, R* c)
{
    return ((((((c[6]*x + c[5])*x + c[4])*x + c[3])*x +
              c[2])*x + c[1])*x + c[0]);
}

//! Templated evaluation of a polynomial of order 8
/*!
 *  @param x   Value of the independent variable - First template parameter
 *  @param c   Pointer to the polynomial - Second template parameter
 */
template<class D, class R>
R poly8(D x, R* c)
{
    return ((((((((c[8]*x + c[7])*x + c[6])*x + c[5])*x + c[4])*x + c[3])*x +
              c[2])*x + c[1])*x + c[0]);
}

//! Templated evaluation of a polynomial of order 10
/*!
 *  @param x   Value of the independent variable - First template parameter
 *  @param c   Pointer to the polynomial - Second template parameter
 */
template<class D, class R>
R poly10(D x, R* c)
{
    return ((((((((((c[10]*x + c[9])*x + c[8])*x + c[7])*x
                  + c[6])*x + c[5])*x + c[4])*x + c[3])*x
              + c[2])*x + c[1])*x + c[0]);
}

//! Templated evaluation of a polynomial of order 5
/*!
 *  @param x   Value of the independent variable - First template parameter
 *  @param c   Pointer to the polynomial - Second template parameter
 */
template<class D, class R>
R poly5(D x, R* c)
{
    return (((((c[5]*x + c[4])*x + c[3])*x +
              c[2])*x + c[1])*x + c[0]);
}

//! Evaluates a polynomial of order 4.
/*!
 *  @param x   Value of the independent variable.
 *  @param c   Pointer to the polynomial coefficient array.
 */
template<class D, class R>
R poly4(D x, R* c)
{
    return ((((c[4]*x + c[3])*x +
              c[2])*x + c[1])*x + c[0]);
}

//! Templated evaluation of a polynomial of order 3
/*!
 *  @param x   Value of the independent variable - First template parameter
 *  @param c   Pointer to the polynomial - Second template parameter
 */
template<class D, class R>
R poly3(D x, R* c)
{
    return (((c[3]*x + c[2])*x + c[1])*x + c[0]);
}

//! Templated deep copy of a std vector of pointers
/*!
 *  Performs a deep copy of a std vectors of pointers to an object. This template assumes that
 *  that the templated object has a functioning copy constructor.
 *  It also assumes that pointers are zero when they are not malloced.
 *
 *  @param fromVec   Vector of pointers to a templated class. This will be
 *                   deep-copied to the other vector
 *  @param toVec     Vector of pointers to a templated class. This will be
 *                   overwritten and on return will be a copy of the fromVec
 */
template<class D>
void deepStdVectorPointerCopy(const std::vector<D*> &fromVec, std::vector<D*> &toVec)
{
    size_t is = toVec.size();
    for (size_t i = 0; i < is; is++) {
        delete toVec[i];
    }
    is = fromVec.size();
    toVec.resize(is);
    for (size_t i = 0; i < is; is++) {
        toVec[i] = new D(*(fromVec[i]));
    }
}

//@}

//! Check to see that a number is finite (not NaN, +Inf or -Inf)
void checkFinite(const double tmp);
}

#endif
