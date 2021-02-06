//! @file funcs.cpp file containing miscellaneous numerical functions.

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/numerics/funcs.h"

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

}
