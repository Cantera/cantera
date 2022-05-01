/**
 *  @file funcs.h Header for a file containing miscellaneous
 *                numerical functions.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_FUNCS_H
#define CT_FUNCS_H

#include "cantera/base/ct_defs.h"
#include "cantera/numerics/eigen_sparse.h"

namespace Cantera
{

//! Linearly interpolate a function defined on a discrete grid.
/*!
 * Vector xpts contains a monotonic sequence of grid points, and vector fpts
 * contains function values defined at these points. The value returned is the
 * linear interpolate at point x. If x is outside the range of xpts, the value
 * of fpts at the nearest end is returned.
 *
 * @param x value of the x coordinate
 * @param xpts value of the grid points
 * @param fpts value of the interpolant at the grid points
 * @returns the value of of the interpolated function at x.
 */
doublereal linearInterp(doublereal x, const vector_fp& xpts,
                        const vector_fp& fpts);

//! Numerical integration of a function using the trapezoidal rule.
/*!
 * Vector x contanins a monotonic sequence of grid points, and
 * Vector f contains function values defined at these points.
 * The size of x and f must be the same.
 *
 * @param  f vector of function value
 * @param  x vector of function coordinate
 */
double trapezoidal(const Eigen::ArrayXd& f, const Eigen::ArrayXd& x);

//! Numerical integration of a function using Simpson's rule
//! with flexibility of taking odd and even number of points.
//! For even number, Simpson's rule is used for the first
//! N-2 intervals with a trapezoidal rule on the last interval.
/*!
 * Vector x contanins a monotonic sequence of grid points, and
 * Vector f contains function values defined at these points.
 * The size of x and f must be the same.
 *
 * @param  f vector of function value
 * @param  x vector of function coordinate
 */
double simpson(const Eigen::ArrayXd& f, const Eigen::ArrayXd& x);

//! Numerical integration of a function.
/*!
 * Vector x contanins a monotonic sequence of grid points, and
 * Vector f contains function values defined at these points.
 * The size of x and f must be the same.
 *
 * @param  method method name
 * @param  f vector of function value
 * @param  x vector of function coordinate
 */
double numericalQuadrature(const std::string& method,
                           const Eigen::ArrayXd& f,
                           const Eigen::ArrayXd& x);
}
#endif
