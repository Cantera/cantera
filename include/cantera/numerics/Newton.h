//! @file Newton.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_NEWTON_H
#define CT_NEWTON_H

#include "FuncEval.h"
#include "DenseMatrix.h"

namespace Cantera
{

struct Configuration {
    vector_fp rtol;
    vector_fp atol;
    double convtol;
    double dt;
    size_t jac_maxage;
    double jac_rtol;
    double jac_atol;
};

/**
 * A Newton solver.
 */
class Newton
{
public:
    Newton(FuncEval& func);
    virtual ~Newton() {};
    Newton(const Newton&) = delete;
    Newton& operator=(const Newton&) = delete;

    size_t size() {
        return m_nv;
    }

    //! Compute the undamped Newton step. The residual function is evaluated
    //! at `x`, but the Jacobian is not recomputed.
    void step(doublereal* x, doublereal* step);

    //! Compute the weighted 2-norm of `step`.
    doublereal weightedNorm(const doublereal* x, const doublereal* step) const;

    int hybridSolve();

    int solve(double* x);

    //TODO: implement get methods
    //nice implementation for steady vs transient below
    //! Relative tolerance of the nth component.
    // doublereal rtol(size_t n) {
    //     return (m_rdt == 0.0 ? m_rtol_ss[n] : m_rtol_ts[n]);
    // }

    //! Set upper and lower bounds on the nth component
    void setBounds(size_t n, doublereal lower, doublereal upper) {
        m_lower_bounds[n] = lower;
        m_upper_bounds[n] = upper;
    }

    void setConstants(vector_int constantComponents) {
        m_constantComponents = constantComponents;
    }

    void evalJacobian(doublereal* x, doublereal* xdot);

    void getSolution(double* x) {
        for (size_t i = 0; i < m_nv; i++) {
            x[i] = m_x[i];
        }
    }

protected:
    FuncEval* m_residfunc;

    //! number of variables
    size_t m_nv;

    DenseMatrix m_jacobian, m_jacFactored;
    size_t m_jacAge;

    //! work arrays of size #m_nv used in solve().
    vector_fp m_x, m_x1, m_stp, m_stp1;

    vector_fp m_upper_bounds, m_lower_bounds;

    vector_fp m_xlast, m_xsave;

    //! the indexes of any constant variables
    vector_int m_constantComponents;

    //! current timestep reciprocal
    doublereal m_rdt = 0;

    Configuration m_directsolve_config;
    Configuration m_timestep_config;
    Configuration* m_config;
};

// //! Returns the weighted Root Mean Square Deviation given a vector of residuals and
// //  vectors of the corresponding weights and absolute tolerances for each component.
// double weightedRMS(vector_fp residuals, vector_fp weights, vector_fp atols) {
//     size_t n = residuals.size();
//     double square = 0.0;
//     for (size_t i = 0; i < n; i++) {
//         square += pow(residuals[i]/(weights[i] + atols[i]), 2);
//     }
//     return sqrt(square/n);
// }

}

#endif
