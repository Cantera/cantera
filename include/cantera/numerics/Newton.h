//! @file Newton.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_NEWTON_H
#define CT_NEWTON_H

#include "Jacobian.h"
#include "FuncEval.h"

namespace Cantera
{

/**
 * A Newton solver.
 */
class Newton
{
public:
    Newton(FuncEval& func, Jacobian& jac);
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

    /**
     * Find the solution to F(X) = 0 by damped Newton iteration.
     */
    int solve(int loglevel=0);

    /// Set options.
    void setOptions(int maxJacAge = 5) {
        m_maxAge = maxJacAge;
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

protected:
    doublereal m_rdt = 0.0;

    FuncEval* m_residfunc;
    Jacobian* m_jac;

    //! Work arrays of size #m_nv used in solve().
    vector_fp m_x, m_x1, m_stp, m_stp1;

    vector_fp m_max, m_min;
    vector_fp m_rtol_ss, m_rtol_ts;
    vector_fp m_atol_ss, m_atol_ts;

    int m_maxAge;

    //! number of variables
    size_t m_nv;

    doublereal m_elapsed;
};
}

#endif
