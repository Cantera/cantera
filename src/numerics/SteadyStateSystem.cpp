//! @file SteadyStateSystem.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/numerics/SteadyStateSystem.h"
#include "cantera/oneD/MultiNewton.h"
#include "cantera/numerics/Func1.h"

using namespace std;

namespace Cantera
{

SteadyStateSystem::SteadyStateSystem()
{
    m_state = make_shared<vector<double>>();
    m_newt = make_unique<MultiNewton>(1);
}

SteadyStateSystem::~SteadyStateSystem()
{
}

void SteadyStateSystem::setInitialGuess(span<const double> x)
{
    clearDebugFile();
    m_attempt_counter = 0;
    m_state->assign(x.begin(), x.end());
}

void SteadyStateSystem::getState(span<double> x) const
{
    copy(m_xnew.begin(), m_xnew.end(), x.begin());
}

void SteadyStateSystem::solve(int loglevel)
{
    size_t istep = 0;
    int nsteps = m_steps[istep];
    m_nsteps = 0;
    double dt = m_tstep;

    while (true) {
        // Keep the attempt_counter in the range of [1, max_history]
        m_attempt_counter = (m_attempt_counter % m_max_history) + 1;

        // Attempt to solve the steady problem
        setSteadyMode();
        newton().setOptions(m_ss_jac_age);
        debuglog("\nAttempt Newton solution of steady-state problem.", loglevel);
        if (!m_jac_ok) {
            evalJacobian(*m_state);
            try {
                m_jac->updateTransient(m_rdt, m_mask);
                m_jac_ok = true;
            } catch (CanteraError& err) {
                // Allow solver to continue after failure to factorize the steady-state
                // Jacobian by returning to time stepping mode
                if (m_rdt == 0.0) {
                    if (loglevel > 1) {
                        writelog("\nSteady Jacobian factorization failed:"
                                "\n  {}:\n    {}",
                                err.getMethod(), err.getMessage());
                    }
                } else {
                    throw;
                }
            }
        }
        int m = -1;
        if (m_jac_ok) {
            m = newton().solve(*m_state, m_xnew, *this, loglevel);
        }
        if (m == 1) {
            *m_state = m_xnew;
            writeDebugInfo("NewtonSuccess", "After successful Newton solve",
                           loglevel, m_attempt_counter);

            return;
        } else {
            debuglog("\nNewton steady-state solve failed.\n", loglevel);
            writeDebugInfo("NewtonFail", "After unsuccessful Newton solve",
                           loglevel, m_attempt_counter);

            if (loglevel > 0) {
                writelog("\nAttempt {} timesteps.", nsteps);
            }

            dt = timeStep(nsteps, dt, *m_state, m_xnew, loglevel-1);
            m_xlast_ts = *m_state;
            writeDebugInfo("Timestepping", "After timestepping", loglevel,
                           m_attempt_counter);

            // Repeat the last timestep's data for logging purposes
            if (loglevel == 1) {
                writelog("\nFinal timestep info: dt= {:<10.4g} log(ss)= {:<10.4g}\n", dt,
                         log10(ssnorm(*m_state, m_xnew)));
            }
            istep++;
            if (istep >= m_steps.size()) {
                nsteps = m_steps.back();
            } else {
                nsteps = m_steps[istep];
            }
            dt = std::min(dt, m_tmax);
        }
    }
}

double SteadyStateSystem::timeStepIncreaseFactor(span<const double> x_before,
                                                 span<const double> x_after)
{
    int heuristic = 2; // Can be 1-4
    m_work1.resize(m_size);
    const double grow_factor = 1.25;

    if (heuristic == 1) {
        // Steady-state residual gate. If the steady-state residual decreased, allow
        // increase in timestep.
        double ss_before = ssnorm(x_before, m_work1);
        double ss_after = ssnorm(x_after, m_work1);
        return (ss_after < ss_before) ? grow_factor : 1.0;
    } else if (heuristic == 2) {
        // Transient residual gate. If the transient residual decreased, allow
        // increase in timestep.
        double ts_before = tsnorm(x_before, m_work1);
        double ts_after = tsnorm(x_after, m_work1);
        return (ts_after < ts_before) ? grow_factor : 1.0;
    } else if (heuristic == 3) {
        // Residual-ratio scaling gate (transient residual). If the transient residual
        // decreased significantly, allow larger increase in timestep.
        double ts_before = tsnorm(x_before, m_work1);
        double ts_after = tsnorm(x_after, m_work1);
        if (!(ts_after > 0.0) || !(ts_before > ts_after)) {
            return 1.0;
        }
        const double max_growth = 1.5;
        const double exponent = 0.2;
        double ratio = ts_before / ts_after;
        double factor = std::pow(ratio, exponent);
        return std::min(max_growth, std::max(1.0, factor));
    } else if (heuristic == 4) {
        // Newton-iteration gate. If the last Newton solve required only a few
        // iterations, allow increase in timestep.
        const int max_iters_for_growth = 3;
        return (newton().lastIterations() <= max_iters_for_growth) ? grow_factor : 1.0;
    }
    return 1.0;
}

double SteadyStateSystem::timeStep(int nsteps, double dt, span<double> x,
                                   span<double> r, int loglevel)
{
    // set the Jacobian age parameter to the transient value
    newton().setOptions(m_ts_jac_age);

    int n = 0;
    int successiveFailures = 0;

    // Only output this if nothing else under this function call will be output
    if (loglevel == 1) {
        writelog("\n============================");
        writelog("\n{:<5s}  {:<11s}   {:<7s}\n", "step", "dt (s)", "log(ss)");
        writelog("============================");
    }
    while (n < nsteps) {
        if (loglevel >= 1) {
            double ss_before = ssnorm(x, r);
            if (loglevel == 1) { // At level 1, output concise information
                writelog("\n{:<5d}  {:<6.4e}   {:>7.4f}", n, dt, log10(ss_before));
            } else {
                writelog("\nTimestep ({}) dt= {:<11.4e}  log(ss)= {:<7.4f}", n, dt,
                         log10(ss_before));
            }
        }

        // set up for time stepping with stepsize dt
        initTimeInteg(dt, x);

        int j0 = m_jac->nEvals(); // Store the current number of Jacobian evaluations

        // solve the transient problem
        int status = newton().solve(x, r, *this, loglevel);

        // successful time step. Copy the new solution in r to
        // the current solution in x.
        if (status >= 0) {
            if (loglevel > 1) {
                writelog("\nTimestep ({}) succeeded", n);
            }
            successiveFailures = 0;
            m_nsteps++;
            n += 1;
            double grow_factor = 1.0;
            if (m_jac->nEvals() == j0) {
                grow_factor = timeStepIncreaseFactor(x, r);
            }
            copy(r.begin(), r.end(), x.begin());
            dt *= grow_factor;
            if (m_time_step_callback) {
                m_time_step_callback->eval(dt);
            }
            dt = std::min(dt, m_tmax);
            if (m_nsteps >= m_nsteps_max) {
                throw CanteraError("SteadyStateSystem::timeStep",
                    "Took maximum number of timesteps allowed ({}) without "
                    "reaching steady-state solution.", m_nsteps_max);
            }
        } else {
            // No solution could be found with this time step.
            // Decrease the stepsize and try again.
            successiveFailures++;
            if (loglevel == 1) {
                writelog("\nTimestep failed");
            } else if (loglevel > 1) {
                writelog("\nTimestep ({}) failed", n);
            }
            if (successiveFailures > 2) {
                debuglog("--> Resetting negative species concentrations", loglevel);
                resetBadValues(x);
                successiveFailures = 0;
            } else {
                debuglog("--> Reducing timestep", loglevel);
                dt *= m_tfactor;
                if (dt < m_tmin) {
                    string err_msg = fmt::format(
                        "Time integration failed. Minimum timestep ({}) reached.", m_tmin);
                    throw CanteraError("SteadyStateSystem::timeStep", err_msg);
                }
            }
        }
    }

    // Write the final step to the log
    if (loglevel == 1) {
        double ss = ssnorm(x, r);
        writelog("\n{:<5d}  {:<6.4e}   {:>7.4f}", n, dt, log10(ss));
        writelog("\n============================");
    } else if (loglevel > 1) {
        double ss = ssnorm(x, r);
        writelog("\nTimestep ({}) dt= {:<11.4e} log10(ss)= {:<7.4f}\n", n, dt, log10(ss));
    }

    // return the value of the last stepsize, which may be smaller
    // than the initial stepsize
    return dt;
}

double SteadyStateSystem::ssnorm(span<const double> x, span<double> r)
{
    eval(x, r, 0.0, 0);
    double ss = 0.0;
    for (size_t i = 0; i < m_size; i++) {
        ss = std::max(fabs(r[i]),ss);
    }
    return ss;
}

double SteadyStateSystem::tsnorm(span<const double> x, span<double> r)
{
    eval(x, r, m_rdt, 0);
    double ts = 0.0;
    for (size_t i = 0; i < m_size; i++) {
        ts = std::max(fabs(r[i]), ts);
    }
    return ts;
}

void SteadyStateSystem::setTimeStep(double stepsize, span<const int> tsteps)
{
    m_tstep = stepsize;
    m_steps.assign(tsteps.begin(), tsteps.end());
}

void SteadyStateSystem::resize()
{
    m_state->resize(size());
    m_xnew.resize(size());
    m_newt->resize(size());
    m_mask.resize(size());
    if (!m_jac) {
        throw CanteraError("SteadyStateSystem::resize",
            "Jacobian evaluator must be instantiated before calling resize()");
    }
    m_jac->initialize(size());
    m_jac->setBandwidth(bandwidth());
    m_jac->clearStats();
    m_jac_ok = false;
}

void SteadyStateSystem::setLinearSolver(shared_ptr<SystemJacobian> solver)
{
    m_jac = solver;
    m_jac->initialize(size());
    m_jac->setBandwidth(bandwidth());
    m_jac->clearStats();
    m_jac_ok = false;
}

void SteadyStateSystem::evalSSJacobian(span<const double> x)
{
    double rdt_save = m_rdt;
    m_jac_ok = false;
    setSteadyMode();
    evalJacobian(x);
    m_rdt = rdt_save;
}

void SteadyStateSystem::initTimeInteg(double dt, span<const double> x)
{
    double rdt_old = m_rdt;
    m_rdt = 1.0/dt;

    // if the stepsize has changed, then update the transient part of the Jacobian
    if (fabs(rdt_old - m_rdt) > Tiny && m_jac_ok) {
        m_jac->updateTransient(m_rdt, m_mask);
    }
}

void SteadyStateSystem::setSteadyMode()
{
    if (m_rdt == 0) {
        return;
    }

    m_rdt = 0.0;
    if (m_jac_ok) {
        m_jac->updateTransient(m_rdt, m_mask);
    }
}

void SteadyStateSystem::setJacAge(int ss_age, int ts_age)
{
    m_ss_jac_age = ss_age;
    if (ts_age > 0) {
        m_ts_jac_age = ts_age;
    } else {
        m_ts_jac_age = m_ss_jac_age;
    }
}

MultiNewton& SteadyStateSystem::newton()
{
    return *m_newt;
}

}
