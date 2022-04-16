//! @file funcs.cpp file containing miscellaneous numerical functions.

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/numerics/funcs.h"
#include "cantera/numerics/polyfit.h"
#include "cantera/base/ctexceptions.h"

using namespace std;

namespace Cantera
{

doublereal linearInterp(doublereal x, const vector_fp& xpts,
                        const vector_fp& fpts)
{
    if (x <= xpts[0]) {
        return fpts[0];
    }
    if (x >= xpts.back()) {
        return fpts.back();
    }
    auto loc = lower_bound(xpts.begin(), xpts.end(), x);
    int iloc = int(loc - xpts.begin()) - 1;
    doublereal ff = fpts[iloc] +
                    (x - xpts[iloc])*(fpts[iloc + 1]
                                      - fpts[iloc])/(xpts[iloc + 1] - xpts[iloc]);
    return ff;
}

double trapezoidal(const Eigen::ArrayXd& f, const Eigen::ArrayXd& x)
{
    // check length
    if (f.size() != x.size()) {
        throw CanteraError("trapezoidal",
                           "Vector lengths need to be the same.");
    }
    // Vector of f(i+1) + f(i)
    Eigen::VectorXd f_av = f.tail(f.size() - 1) + f.head(f.size() - 1);
    // Vector of x(i+1) - x(i)
    Eigen::VectorXd x_diff = x.tail(x.size() - 1) - x.head(x.size() - 1);
    // check if the coordinate is a monotonically increase vector.
    if ((x_diff.array() <= 0.0).any()) {
        throw CanteraError("trapezoidal",
            "x (coordinate) needs to be the monotonically increasing.");
    }
    return f_av.dot(x_diff) / 2.0;
}

//! Numerical integration of a function using Simpson's rule.
//! Only for odd number of points. This function is used only
//! by calling simpson.
/*!
 * Vector x contanins a monotonic sequence of grid points, and
 * Vector f contains function values defined at these points.
 * The size of x and f must be the same.
 *
 * @param  f vector of function value
 * @param  x vector of function coordinate
 */
double basicSimpson(const Eigen::ArrayXd& f, const Eigen::ArrayXd& x)
{
    if (f.size() < 2) {
        throw CanteraError("basicSimpson",
                           "Vector lengths need to be larger than two.");
    }
    if (f.size()%2 == 0) {
        throw CanteraError("basicSimpson",
                           "Vector lengths need to be an odd number.");
    }

    size_t N = f.size() - 1;
    Eigen::VectorXd h = x.tail(N) - x.head(N);

    double sum = 0.0;
    for (size_t i = 1; i < N; i+=2) {
        double h0 = h[i-1];
        double h1 = h[i];
        double hph = h1 + h0;
        double hdh = h1 / h0;
        double hmh = h1 * h0;
        sum += (hph / 6.0) * (
                    (2.0 - hdh) * f[i - 1] + (pow(hph, 2) / hmh) * f[i] +
                    (2.0 - 1.0 / hdh) * f[i + 1]);
    }
    return sum;
}

double simpson(const Eigen::ArrayXd& f, const Eigen::ArrayXd& x)
{
    Eigen::ArrayXd h = x.tail(x.size() - 1) - x.head(x.size() - 1);
    if ((h <= 0.0).any()) {
        throw CanteraError("simpson",
            "Values of x need to be positive and monotonically increasing.");
    }
    if (f.size() != x.size()) {
        throw CanteraError("simpson", "Vector lengths need to be the same.");
    }

    if (f.size()%2 == 1) {
        return basicSimpson(f, x);
    } else if (f.size() == 2) {
        return 0.5 * h[0] * (f[1] + f[0]);
    } else {
        size_t N = f.size() - 1;
        // pick first N-1 point for simpson
        double headSimps = basicSimpson(f.head(N), x.head(N));
        // Use trapezoidal rules for the last interval
        double tailTrap = 0.5 * h[N-1] * (f[N] + f[N-1]);
        return headSimps + tailTrap;
    }
}

double numericalQuadrature(const std::string& method,
                           const Eigen::ArrayXd& f,
                           const Eigen::ArrayXd& x)
{
    if (method == "simpson") {
        return simpson(f, x);
    } else if (method == "trapezoidal") {
        return trapezoidal(f, x);
    } else {
        throw CanteraError("numericalQuadrature",
                           "Unknown method of numerical quadrature. "
                           "Please use 'simpson' or 'trapezoidal'");
    }
}

}
