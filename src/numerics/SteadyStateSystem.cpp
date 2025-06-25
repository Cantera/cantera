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
    m_newt = make_unique<MultiNewton>(1);
}

SteadyStateSystem::~SteadyStateSystem()
{
}

int SteadyStateSystem::solve(double* x, double* xnew, int loglevel)
{
    if (!m_jac_ok) {
        evalJacobian(x);
        m_jac->updateTransient(m_rdt, m_mask.data());
        m_jac_ok = true;
    }
    return m_newt->solve(x, xnew, *this, loglevel);
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
        int status = solve(x, r, loglevel);

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
