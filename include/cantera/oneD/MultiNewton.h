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
    //! Constructor
    //! @param sz  Number of variables in the system
    MultiNewton(int sz);
    virtual ~MultiNewton() {};
    MultiNewton(const MultiNewton&) = delete;
    MultiNewton& operator=(const MultiNewton&) = delete;

    //! Get the number of variables in the system.
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
     * This uses the method outlined in Kee et al. @cite kee2003.
     *
     * The system of equations can be written in the form:
     * @f[
     *   F(x) = 0
     * @f]
     *
     * Where @f$ F @f$ is the system of nonlinear equations and @f$ x @f$ is the
     * solution vector.
     *
     * For the damped Newton method we are solving:
     *
     * @f[
     *   x_{k+1} - x_k = \Delta x_k = -\alpha_k J^{-1}(x_k) F(x_k)
     * @f]
     *
     * Where @f$ J @f$ is the Jacobian matrix of @f$ F @f$ with respect to @f$ x @f$,
     * and @f$ \alpha_k @f$ is the damping factor, and @f$ \Delta x_k @f$ is the Newton
     * step at @f$ x_k @f$, sometimes called the correction vector. In the equations
     * here, @f$ k @f$ is just an iteration variable.
     *
     * In this method, the Jacobian does not update, even when the solution vector is
     * evaluated at different points.
     *
     * The general algorithm is described below.
     *
     * We want to solve the equation:
     * @f[
     *   x_{k+1} = x_k + \alpha_k \Delta x_k
     * @f]
     *
     * Pick @f$ \alpha_k @f$ such that
     * @f$ \Vert \Delta x_{k+1} \Vert < \Vert \Delta x_k \Vert @f$
     * where @f$ \Delta x_k = J^{-1}(x_k) F(x_k) @f$, and
     * @f$ \Delta x_{k+1} = J^{-1}(x_{k}) F(x_{k+1}) @f$.
     *
     * @param[in] x0  initial solution about which a Newton step will be taken
     * @param[in] step0  initial undamped Newton step
     * @param[out] x1  solution after taking the damped Newton step
     * @param[out] step1  Newton step after taking the damped Newton step
     * @param[out] s1  norm of the subsequent Newton step after taking the damped Newton
     *     step
     * @param[in] r  domain object, used for evaluating residuals over all domains
     * @param[in] jac  Jacobian evaluator
     * @param[in] loglevel  controls amount of printed diagnostics
     * @param[in] writetitle  controls if logging title is printed
     *
     * @returns
     *   - `1` a damping coefficient was found and the solution converges.
     *   - `0` a damping coefficient was found, but the solution has not converged yet.
     *   - `-2` no suitable damping coefficient was found within the maximum iterations.
     *   - `-3` the current solution `x0` is too close to the solution bounds and the
     *          step would exceed the bounds on one or more components.
     */
    int dampStep(const double* x0, const double* step0, double* x1, double* step1,
                 double& s1, OneDim& r, MultiJac& jac, int loglevel, bool writetitle);

    //! Compute the weighted 2-norm of `step`.
    double norm2(const double* x, const double* step, OneDim& r) const;

    /**
     * Find the solution to F(x) = 0 by damped Newton iteration. On entry, x0
     * contains an initial estimate of the solution. On successful return, x1
     * contains the converged solution. If failure occurs, x1 will contain the
     * value of x0 i.e. no change in solution.
     *
     * The convergence criteria is when the 2-norm of the Newton step is less than one.
     *
     * @returns
     *   - `1`  a converged solution was found.
     *   - `-2` no suitable damping coefficient was found within the maximum iterations.
     *   - `-3` the current solution `x0` is too close to the solution bounds and the
     *          step would exceed the bounds on one or more components.
     */
    int solve(double* x0, double* x1, OneDim& r, MultiJac& jac, int loglevel);

    //! Set options.
    //! @param maxJacAge  Maximum number of steps that can be taken before requiring
    //!                   a Jacobian update
    void setOptions(int maxJacAge = 5) {
        m_maxAge = maxJacAge;
    }

    //! Change the problem size.
    void resize(size_t points);

protected:
    //! Work array holding the system state after the last successful step. Size #m_n.
    vector<double> m_x;

    //! Work array holding the undamped Newton step or the system residual. Size #m_n.
    vector<double> m_stp;

    //! Work array holding the damped Newton step. Size #m_n.
    vector<double> m_stp1;

    //! Maximum allowable Jacobian age before it is recomputed.
    int m_maxAge = 5;

    //! Factor by which the damping coefficient is reduced in each iteration
    double m_dampFactor = sqrt(2.0);

    //! Maximum number of damping iterations
    size_t m_maxDampIter = 7;

    //! number of variables
    size_t m_n;

    //! Elapsed CPU time spent computing the Jacobian.
    double m_elapsed = 0.0;
};
}

#endif
