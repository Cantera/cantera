//! @file polyfit.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_POLYFIT_H
#define CT_POLYFIT_H

#include "cantera/base/ct_defs.h"

namespace Cantera
{

//! Fits a polynomial function to a set of data points
/*!
 * Given a collection of *n* points *x* and a set of values *y* of some function
 * evaluated at those points, this function computes the weighted least-squares
 * polynomial fit of degree *deg*:
 *
 * @f[ f(x) = p[0] + p[1] x + p[2] x^2 + \cdots + p[deg] x^deg @f]
 *
 * @param deg  The degree of the polynomial fit to be computed. `deg <= x.size() - 1`.
 * @param x    Array of points at which the function is evaluated.
 * @param y    Array of function values at the points in *x*.
 * @param w    Array of weights. If empty or `w[0] < 0`, then all the
 *     weights will be set to 1.0.
 * @param[out] p  Array of polynomial coefficients, starting with the constant
 *     term. Length *deg+1*.
 * @returns the root mean squared error of the fit at the input points.
 * @ingroup mathUtils
 */
double polyfit(size_t deg, span<const double> x, span<const double> y,
               span<const double> w, span<double> p);

}
#endif
