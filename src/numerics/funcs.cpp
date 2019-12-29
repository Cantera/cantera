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

double simpsonQuadrature(const vector_fp& x, const vector_fp& y)
{
    if (x.size() != y.size()) {
        throw CanteraError("Cantera::simpsonQuadrature",
                           "size of x is not equal to size of y");
    }
    size_t N = x.size();
    size_t ns = (N - 1) / 2;
    size_t ms = (N - 1) % 2;
    double sum = 0.0;
    for (size_t i = 0; i < ns; i++) {
        double c[3];
        double xs[3] = {x[2*i],x[2*i+1],x[2*i+2]};
        double ys[3] = {y[2*i],y[2*i+1],y[2*i+2]};
        polyfit(3, 2, xs, ys, nullptr, c);
        sum += c[0] * x[2*i+2] +
               c[1] * 0.5 * x[2*i+2] * x[2*i+2] +
               c[2] * 1./3. * x[2*i+2] * x[2*i+2] * x[2*i+2] -
               c[0] * x[2*i] -
               c[1] * 0.5 * x[2*i] * x[2*i] -
               c[2] * 1./3. * x[2*i] * x[2*i] * x[2*i];
    }
    if (ms == 1) {
        sum += 0.5 * (x[N-1] - x[N-2]) * (y[N-1] + y[N-2]);
    }
    return sum;
}

}
