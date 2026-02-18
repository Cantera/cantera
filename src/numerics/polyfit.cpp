//! @file polyfit.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/numerics/polyfit.h"
#include "cantera/numerics/eigen_dense.h"
#include "cantera/base/global.h"
#include "cantera/base/ctexceptions.h"

namespace Cantera
{

double polyfit(size_t deg, span<const double> x, span<const double> y,
               span<const double> w, span<double> p)
{
    size_t n = x.size();
    checkArraySize("polyfit", y.size(), n);
    checkArraySize("polyfit", p.size(), deg + 1);
    if (!w.empty()) {
        checkArraySize("polyfit", w.size(), n);
    }

    if (deg >= n) {
        throw CanteraError("polyfit", "Polynomial degree ({}) must be less "
            "than number of input data points ({})", deg, n);
    }

    ConstMappedVector xx(x.data(), n);
    Eigen::VectorXd yy = ConstMappedVector(y.data(), n);
    MappedVector pp(p.data(), deg + 1);

    // Construct A such that each row i of A has the elements
    // 1, x[i], x[i]^2, x[i]^3 ... + x[i]^deg
    Eigen::MatrixXd A(n, deg+1);
    A.col(0).setConstant(1.0);

    if (deg > 0) {
        A.col(1) = xx;
    }
    for (size_t i = 1; i < deg; i++) {
        A.col(i+1) = A.col(i).array() * xx.array();
    }

    if (!w.empty() && w[0] > 0) {
        // For compatibility with old Fortran dpolft, input weights are the
        // squares of the weight vector used in this algorithm
        Eigen::VectorXd ww = ConstMappedVector(w.data(), n).cwiseSqrt().eval();

        // Multiply by the weights on both sides
        A = ww.asDiagonal() * A;
        yy.array() *= ww.array();
    }

    // Solve W*A*p = W*y to find the polynomial coefficients
    pp = A.colPivHouseholderQr().solve(yy);

    // Evaluate the computed polynomial at the input x coordinates to compute
    // the RMS error as the return value
    return (A * pp - yy).eval().norm() / sqrt(n);
}

}
