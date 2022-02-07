//! @file funcs.cpp file containing miscellaneous numerical functions.

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/numerics/funcs.h"
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

}
