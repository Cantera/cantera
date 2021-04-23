//! @file Newton.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_NEWTON_H
#define CT_NEWTON_H

#include "FuncEval.h"
#include "DenseMatrix.h"

namespace Cantera
{

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
    void step(doublereal* x, doublereal* step, int loglevel);

    /**
     * Return the factor by which the undamped Newton step 'step0'
     * must be multiplied in order to keep all solution components in
     * all domains between their specified lower and upper bounds.
     */
    doublereal boundStep(const doublereal* x0, const doublereal* step0, int loglevel);

    /**
     * On entry, step0 must contain an undamped Newton step for the solution x0.
     * This method attempts to find a damping coefficient such that the next
     * undamped step would have a norm smaller than that of step0. If
     * successful, the new solution after taking the damped step is returned in
     * x1, and the undamped step at x1 is returned in step1.
     */
    int dampStep(const doublereal* x0, const doublereal* step0,
                 doublereal* x1, doublereal* step1, doublereal& s1,
                 int loglevel, bool writetitle);

    //! Compute the weighted 2-norm of `step`.
    doublereal weightedNorm(const doublereal* x, const doublereal* step) const;

    int hybridSolve();

    int timestep();

    /**
     * Find the solution to F(X) = 0 by damped Newton iteration.
     */
    int solve(int loglevel=0);

    /// Set options.
    void setOptions(int maxJacAge = 5) {
        m_jacMaxAge = maxJacAge;
    }

    //TODO: implement get methods
    //nice implementation for steady vs transient below
    //! Relative tolerance of the nth component.
    // doublereal rtol(size_t n) {
    //     return (m_rdt == 0.0 ? m_rtol_ss[n] : m_rtol_ts[n]);
    // }

    //! Set upper and lower bounds on the nth component
    void setBounds(size_t n, doublereal lower, doublereal upper) {
        m_min[n] = lower;
        m_max[n] = upper;
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

    //! solution converged if [weightedNorm(sol, step) < m_convergenceThreshold]
    doublereal m_convergenceThreshold;

    DenseMatrix m_jacobian, m_jacFactored;
    int m_jacAge, m_jacMaxAge;
    doublereal m_jacRtol, m_jacAtol;


    //! work arrays of size #m_nv used in solve().
    vector_fp m_x, m_x1, m_stp, m_stp1;

    vector_fp m_max, m_min;
    vector_fp m_rtol_ss, m_rtol_ts;
    vector_fp m_atol_ss, m_atol_ts;

    vector_fp m_xlast, m_xsave;

    //! the indexes of any constant variables
    vector_int m_constantComponents;

    //! current timestep reciprocal
    doublereal m_rdt = 0;
};
}

#endif
