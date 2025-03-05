//! @file SteadyStateSystem.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_STEADYSTATESYSTEM_H
#define CT_STEADYSTATESYSTEM_H

#include "cantera/base/ct_defs.h"
#include "SystemJacobian.h"

namespace Cantera
{

class Func1;
class MultiNewton;

//! Base class for representing a system of differential-algebraic equations and solving
//! for its steady-state response.
class SteadyStateSystem
{
public:
    SteadyStateSystem();
    virtual ~SteadyStateSystem();
    SteadyStateSystem(const SteadyStateSystem&) = delete;
    SteadyStateSystem& operator=(const SteadyStateSystem&) = delete;

    //! Evaluate the residual function
    //!
    //! @param[in] x  State vector
    //! @param[out] r  On return, contains the residual vector
    //! @param rdt  Reciprocal of the time step. if omitted, then the internally stored
    //!     value (accessible using the rdt() method) is used.
    //! @param count   Set to zero to omit this call from the statistics
    virtual void eval(double* x, double* r, double rdt=-1.0, int count=1) = 0;

    //! Evaluates the Jacobian at x0 using finite differences.
    //!
    //! The Jacobian is computed by perturbing each component of `x0`, evaluating the
    //! residual function, and then estimating the partial derivatives numerically using
    //! finite differences to determine the corresponding column of the Jacobian.
    //!
    //! @param x0  State vector at which to evaluate the Jacobian
    virtual void evalJacobian(double* x0) = 0;

    //! Compute the weighted norm of `step`.
    //!
    //! The purpose is to measure the "size" of the step vector @f$ \Delta x @f$ in a
    //! way that takes into account the relative importance or scale of different
    //! solution components. Each component of the step vector is normalized by a weight
    //! that depends on both the current magnitude of the solution vector and specified
    //! tolerances. This makes the norm dimensionless and scaled appropriately, avoiding
    //! issues where some components dominate due to differences in their scales. See
    //! OneDim::weightedNorm() for a representative implementation.
    virtual double weightedNorm(const double* step) const = 0;

    //! Solve @f$ F(x) = 0 @f$, where @f$ F(x) @f$ is the residual function.
    //!
    //! @param[in] x0  Starting estimate of solution.
    //! @param[out] x1  Final solution satisfying @f$ F(x1) = 0 @f$.
    //! @param loglevel   Controls amount of diagnostic output.
    //!
    //! @returns
    //! - 1 for success
    //! - -2 failure (maximum number of damping steps was reached)
    //! - -3 failure (solution was up against the bounds
    //!
    //! @deprecated To be removed after Cantera 3.2. Use setInitialGuess(), solve(int),
    //!     and getState() instead.
    int solve(double* x0, double* x1, int loglevel);

    //! Set the initial guess. Should be called before solve().
    void setInitialGuess(const double* x);

    //! Solve the steady-state problem, taking internal timesteps as necessary until
    //! the Newton solver can converge for the steady problem.
    //! @param loglevel   Controls amount of diagnostic output.
    void solve(int loglevel=0);

    //! Get the converged steady-state solution after calling solve().
    void getState(double* x) const;

    //! Steady-state max norm (infinity norm) of the residual evaluated using solution
    //! x. On return, array r contains the steady-state residual values.
    double ssnorm(double* x, double* r);

    //! Total solution vector length;
    size_t size() const {
        return m_size;
    }

    //! Call to set the size of internal data structures after first defining the system
    //! or if the problem size changes, for example after grid refinement.
    //!
    //! The #m_size variable should be updated before calling this method
    virtual void resize();

    //! Jacobian bandwidth.
    size_t bandwidth() const {
        return m_bw;
    }

    //! Get the name of the i-th component of the state vector
    virtual string componentName(size_t i) const { return fmt::format("{}", i); }

    //! Get header lines describing the column names included in a component label.
    //! Headings should be aligned with items formatted in the output of
    //! componentTableLabel(). Two lines are allowed. If only one is needed, the first
    //! line should be blank.
    virtual pair<string, string> componentTableHeader() const {
        return {"", "component name"};
    }

    //! Get elements of the component name, aligned with the column headings given
    //! by componentTableHeader().
    virtual string componentTableLabel(size_t i) const { return componentName(i); }

    //! Get the upper bound for global component *i* in the state vector
    virtual double upperBound(size_t i) const = 0;

    //! Get the lower bound for global component *i* in the state vector
    virtual double lowerBound(size_t i) const = 0;

    //! Return a reference to the Newton iterator.
    MultiNewton& newton();

    //! Set the linear solver used to hold the Jacobian matrix and solve linear systems
    //! as part of each Newton iteration.
    void setLinearSolver(shared_ptr<SystemJacobian> solver);

    //! Get the the linear solver being used to hold the Jacobian matrix and solve
    //! linear systems as part of each Newton iteration.
    shared_ptr<SystemJacobian> linearSolver() const { return m_jac; }

    //! Reciprocal of the time step.
    double rdt() const {
        return m_rdt;
    }

    //! True if transient mode.
    bool transient() const {
        return (m_rdt != 0.0);
    }

    //! True if steady mode.
    bool steady() const {
        return (m_rdt == 0.0);
    }

    //! Prepare for time stepping beginning with solution *x* and timestep *dt*.
    virtual void initTimeInteg(double dt, double* x);

    //! Prepare to solve the steady-state problem. After invoking this method,
    //! subsequent calls to solve() will solve the steady-state problem. Sets the
    //! reciprocal of the time step to zero, and, if it was previously non-zero, signals
    //! that a new Jacobian will be needed.
    virtual void setSteadyMode();

    //! Access the vector indicating which equations contain a transient term.
    //! Elements are 1 for equations with a transient terms and 0 otherwise.
    vector<int>& transientMask() {
        return m_mask;
    }

    //! Take time steps using Backward Euler.
    //!
    //! @param nsteps  number of steps
    //! @param dt  initial step size
    //! @param x  current solution vector
    //! @param r  solution vector after time stepping
    //! @param loglevel  controls amount of printed diagnostics
    //! @returns size of last timestep taken
    double timeStep(int nsteps, double dt, double* x, double* r, int loglevel);

    //! Reset values such as negative species concentrations
    virtual void resetBadValues(double* x) {}

    //! @name Options
    //! @{

    //! Set the number of time steps to try when the steady Newton solver is
    //! unsuccessful.
    //! @param stepsize  Initial time step size [s]
    //! @param n  Length of `tsteps` array
    //! @param tsteps  A sequence of time step counts to take after subsequent failures
    //!     of the steady-state solver. The last value in `tsteps` will be used again
    //!     after further unsuccessful solution attempts.
    void setTimeStep(double stepsize, size_t n, const int* tsteps);

    //! Set the minimum time step allowed during time stepping
    void setMinTimeStep(double tmin) {
        m_tmin = tmin;
    }

    //! Set the maximum time step allowed during time stepping
    void setMaxTimeStep(double tmax) {
        m_tmax = tmax;
    }

    //! Sets a factor by which the time step is reduced if the time stepping fails.
    //! The default value is 0.5.
    //!
    //! @param tfactor  factor time step is multiplied by if time stepping fails
    void setTimeStepFactor(double tfactor) {
        m_tfactor = tfactor;
    }

    //! Set the maximum number of timeteps allowed before successful steady-state solve
    void setMaxTimeStepCount(int nmax) {
        m_nsteps_max = nmax;
    }

    //! Get the maximum number of timeteps allowed before successful steady-state solve
    int maxTimeStepCount() const {
        return m_nsteps_max;
    }
    //! @}

    //! Set the maximum number of steps that can be taken using the same Jacobian
    //! before it must be re-evaluated.
    //! @param ss_age  Age limit during steady-state mode
    //! @param ts_age  Age limit during time stepping mode. If not specified, the
    //!     steady-state age is also used during time stepping.
    void setJacAge(int ss_age, int ts_age=-1);

    //! Set a function that will be called every time #eval is called.
    //! Can be used to provide keyboard interrupt support in the high-level
    //! language interfaces.
    void setInterrupt(Func1* interrupt) {
        m_interrupt = interrupt;
    }

    //! Set a function that will be called after each successful timestep. The
    //! function will be called with the size of the timestep as the argument.
    //! Intended to be used for observing solver progress for debugging
    //! purposes.
    void setTimeStepCallback(Func1* callback) {
        m_time_step_callback = callback;
    }

    //! Configure perturbations used to evaluate finite difference Jacobian
    //! @param relative  Relative perturbation (multiplied by the absolute value of
    //!     each component). Default `1.0e-5`.
    //! @param absolute  Absolute perturbation (independent of component value).
    //!     Default `1.0e-10`.
    //! @param threshold  Threshold below which to exclude elements from the Jacobian
    //!     Default `0.0`.
    void setJacobianPerturbation(double relative, double absolute, double threshold) {
        m_jacobianRelPerturb = relative;
        m_jacobianAbsPerturb = absolute;
        m_jacobianThreshold = threshold;
    }

    //! Write solver debugging based on the specified log level.
    //!
    //! @see Sim1D::writeDebugInfo for a specific implementation of this capability.
    virtual void writeDebugInfo(const string& header_suffix, const string& message,
                                int loglevel, int attempt_counter) {}

    //! Delete debug output file that may be created by writeDebugInfo() when solving
    //! with high `loglevel`.
    virtual void clearDebugFile() {}

protected:
    //! Evaluate the steady-state Jacobian, accessible via linearSolver()
    //! @param[in] x  Current state vector, length size()
    //! @param[out] rsd  Storage for the residual, length size()
    void evalSSJacobian(double* x, double* rsd);

    //! Array of number of steps to take after each unsuccessful steady-state solve
    //! before re-attempting the steady-state solution. For subsequent time stepping
    //! calls, the final value is reused. See setTimeStep().
    vector<int> m_steps = { 10 };

    double m_tstep = 1.0e-5; //!< Initial timestep
    double m_tmin = 1e-16; //!< minimum timestep size
    double m_tmax = 1e+08; //!< maximum timestep size

    //! factor time step is multiplied by  if time stepping fails ( < 1 )
    double m_tfactor = 0.5;

    shared_ptr<vector<double>> m_state; //!< Solution vector

    //! Work array used to hold the residual or the new solution
    vector<double> m_xnew;

    //! State vector after the last successful set of time steps
    vector<double> m_xlast_ts;

    unique_ptr<MultiNewton> m_newt; //!< Newton iterator
    double m_rdt = 0.0; //!< reciprocal of time step

    shared_ptr<SystemJacobian> m_jac; //!< Jacobian evaluator
    bool m_jac_ok = false; //!< if true, Jacobian is current

    size_t m_bw = 0; //!< Jacobian bandwidth
    size_t m_size = 0; //!< solution vector size

    //! Work arrays used during Jacobian evaluation
    vector<double> m_work1, m_work2;

    vector<int> m_mask; //!< Transient mask. See transientMask()
    int m_ss_jac_age = 20; //!< Maximum age of the Jacobian in steady-state mode.
    int m_ts_jac_age = 20; //!< Maximum age of the Jacobian in time-stepping mode.

    //! Counter used to manage the number of states stored in the debug log file
    //! generated by writeDebugInfo()
    int m_attempt_counter = 0;

    //! Constant that determines the maximum number of states stored in the debug log
    //! file generated by writeDebugInfo()
    int m_max_history = 10;

    //! Function called at the start of every call to #eval.
    Func1* m_interrupt = nullptr;

    //! User-supplied function called after each successful timestep.
    Func1* m_time_step_callback = nullptr;

    int m_nsteps = 0; //!< Number of time steps taken in the current call to solve()
    int m_nsteps_max = 500; //!< Maximum number of timesteps allowed per call to solve()

    //! Threshold for ignoring small elements in Jacobian
    double m_jacobianThreshold = 0.0;
    //! Relative perturbation of each component in finite difference Jacobian
    double m_jacobianRelPerturb = 1e-5;
    //! Absolute perturbation of each component in finite difference Jacobian
    double m_jacobianAbsPerturb = 1e-10;
};

}

#endif
