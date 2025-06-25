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

void SteadyStateSystem::setInitialGuess(const double* x)
{
    clearDebugFile();
    m_attempt_counter = 0;
    m_state->assign(x, x + size());
}

void SteadyStateSystem::getState(double* x) const
{
    copy(m_xnew.begin(), m_xnew.end(), x);
}

int SteadyStateSystem::solve(double* x, double* xnew, int loglevel)
{
    warn_deprecated("SteadyStateSystem::solve(double*, double*, int)",
        "To be removed after Cantera 3.2. Use setInitialGuess, solve(int), and "
        "getState instead.");
    setInitialGuess(x);
    solve(loglevel);
    getState(xnew);
    return 1;
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
            evalJacobian(m_state->data());
            try {
                m_jac->updateTransient(m_rdt, m_mask.data());
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
            m = newton().solve(m_state->data(), m_xnew.data(), *this, loglevel);
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

            dt = timeStep(nsteps, dt, m_state->data(), m_xnew.data(), loglevel-1);
            m_xlast_ts = *m_state;
            writeDebugInfo("Timestepping", "After timestepping", loglevel,
                           m_attempt_counter);

            // Repeat the last timestep's data for logging purposes
            if (loglevel == 1) {
                writelog("\nFinal timestep info: dt= {:<10.4g} log(ss)= {:<10.4g}\n", dt,
                         log10(ssnorm(m_state->data(), m_xnew.data())));
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

double SteadyStateSystem::timeStep(int nsteps, double dt, double* x, double* r, int loglevel)
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
        if (loglevel == 1) { // At level 1, output concise information
            double ss = ssnorm(x, r);
            writelog("\n{:<5d}  {:<6.4e}   {:>7.4f}", n, dt, log10(ss));
        } else if (loglevel > 1) {
            double ss = ssnorm(x, r);
            writelog("\nTimestep ({}) dt= {:<11.4e}  log(ss)= {:<7.4f}", n, dt, log10(ss));
        }

        // set up for time stepping with stepsize dt
        initTimeInteg(dt,x);

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
            copy(r, r + m_size, x);
            // No Jacobian evaluations were performed, so a larger timestep can be used
            if (m_jac->nEvals() == j0) {
                dt *= 1.5;
            }
            if (m_time_step_callback) {
                m_time_step_callback->eval(dt);
            }
            dt = std::min(dt, m_tmax);
            if (m_nsteps >= m_nsteps_max) {
                throw CanteraError("OneDim::timeStep",
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
                    throw CanteraError("OneDim::timeStep", err_msg);
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

double SteadyStateSystem::ssnorm(double* x, double* r)
{
    eval(x, r, 0.0, 0);
    double ss = 0.0;
    for (size_t i = 0; i < m_size; i++) {
        ss = std::max(fabs(r[i]),ss);
    }
    return ss;
}

void SteadyStateSystem::setTimeStep(double stepsize, size_t n, const int* tsteps)
{
    m_tstep = stepsize;
    m_steps.resize(n);
    for (size_t i = 0; i < n; i++) {
        m_steps[i] = tsteps[i];
    }
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

void SteadyStateSystem::evalSSJacobian(double* x, double* rsd)
{
    double rdt_save = m_rdt;
    m_jac_ok = false;
    setSteadyMode();
    evalJacobian(x);
    m_rdt = rdt_save;
}

void SteadyStateSystem::initTimeInteg(double dt, double* x)
{
    double rdt_old = m_rdt;
    m_rdt = 1.0/dt;

    // if the stepsize has changed, then update the transient part of the Jacobian
    if (fabs(rdt_old - m_rdt) > Tiny) {
        m_jac->updateTransient(m_rdt, m_mask.data());
    }
}

void SteadyStateSystem::setSteadyMode()
{
    if (m_rdt == 0) {
        return;
    }

    m_rdt = 0.0;
    m_jac->updateTransient(m_rdt, m_mask.data());
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
