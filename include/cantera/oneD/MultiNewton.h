//! @file MultiNewton.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_MULTINEWTON_H
#define CT_MULTINEWTON_H

#include "MultiJac.h"

namespace Cantera
{

//! @defgroup onedUtilsGroup Utilities
//! Utility classes and functions for one-dimensional problems.
//! @ingroup onedGroup

/**
 * Newton iterator for multi-domain, one-dimensional problems.
 * Used by class OneDim.
 * @ingroup onedUtilsGroup
 */
class MultiNewton
{
public:
    MultiNewton(int sz);
    virtual ~MultiNewton() {};
    MultiNewton(const MultiNewton&) = delete;
    MultiNewton& operator=(const MultiNewton&) = delete;

    size_t size() {
        return m_n;
    }

    //! Compute the undamped Newton step.  The residual function is evaluated
    //! at `x`, but the Jacobian is not recomputed.
    void step(double* x, double* step, OneDim& r, MultiJac& jac, int loglevel);

    /**
     * Return the factor by which the undamped Newton step 'step0'
     * must be multiplied in order to keep all solution components in
     * all domains between their specified lower and upper bounds.
     */
    double boundStep(const double* x0, const double* step0,
                     const OneDim& r, int loglevel);

    /**
     * Performs a damped Newton step to solve the system of nonlinear equations.
     *
     * On entry, `step0` must contain an undamped Newton step for the solution `x0`.
     * This method attempts to find a damping coefficient `alpha_k` such that the next
     * undamped step would have a norm smaller than that of `step0`. If successful,
     * the new solution after taking the damped step is returned in `x1`, and the
     * undamped step at `x1` is returned in `step1`.
     *
     * This uses the method outlined in @cite kee2003.
     *
     * The system of equations can be written in the form:
     * @f[
     *   F(x) = 0
     * @f]
     *
     * Where `F` is the system of nonlinear equations, `x` is the solution vector.
     *
     * For the damped Newton method we are solving:
     *
     * @f[
     *   x_{k+1} - x_k = \inc x_k = -\alpha_k J^(-1)(x_k) F(x_k)
     * @f]
     *
     * Where `J` is the Jacobian matrix of `F` with respect to `x`, and `alpha_k` is
     * the damping factor, and @f$ \inc x_k @f$ is the Newton step at `x_k`, sometimes
     * called the correction vector. In the equations here, k is just an iteration
     * variable.
     *
     * In this method, the Jacobian does not update, even when the solution vector is
     * evaluated at different points.
     *
     * The general algorithm is described below.
     *
     * We want to solve the equation:
     * @f[
     *   x_{k+1} = x_k + \alpha_k \inc x_k
     * @f]
     *
     * Pick @f$ \alpha_k @f$ such that @f$ \norm{\inc x_{k+1}} < \norm{\inc x_k} @f$.
     * Where @f$ \inc x_k = J^{-1}(x_k) F(x_k) @f$, and
     * @f$ \inc x_{k+1} = J^{-1}(x_{k}) F(x_{k+1}) @f$.
     *
     * @param x0 initial solution about which a Newton step will be taken
     * @param step0 initial Newton step without any damping
     * @param x1 solution after taking the damped Newton step
     * @param step1 Newton step after taking the damped Newton step
     * @param loglevel controls amount of printed diagnostics
     * @param writetitle controls if logging title is printed
     *
     * Returns:
     * - int : Status code
     *   - `1` if a damping coefficient is found and the solution converges.
     *   - `0` if a damping coefficient is found but the solution does not converge.
     *   - `-2` if no suitable damping coefficient is found within the maximum iterations.
     *   - `-3` if the current solution `x0` is too close to the boundary and the step points out of the allowed domain.
     */
    int dampStep(const double* x0, const double* step0, double* x1, double* step1,
                 double& s1, OneDim& r, MultiJac& jac, int loglevel, bool writetitle);

    //! Compute the weighted 2-norm of `step`.
    double norm2(const double* x, const double* step, OneDim& r) const;

    /**
     * Find the solution to F(X) = 0 by damped Newton iteration. On entry, x0
     * contains an initial estimate of the solution. On successful return, x1
     * contains the converged solution.
     */
    int solve(double* x0, double* x1, OneDim& r, MultiJac& jac, int loglevel);

    //! Set options.
    void setOptions(int maxJacAge = 5) {
        m_maxAge = maxJacAge;
    }

    //! Change the problem size.
    void resize(size_t points);

protected:
    //! Work arrays of size #m_n used in solve().
    vector<double> m_x, m_stp, m_stp1, temp_x0, temp_stp0;

    int m_maxAge = 5;

    //! number of variables
    size_t m_n;

    double m_elapsed = 0.0;
};
}

#endif
