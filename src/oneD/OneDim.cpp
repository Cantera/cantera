//! @file OneDim.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/oneD/OneDim.h"
#include "cantera/numerics/Func1.h"
#include "cantera/oneD/MultiNewton.h"
#include "cantera/base/AnyMap.h"
#include "cantera/numerics/SystemJacobianFactory.h"

#include <fstream>
#include <ctime>

using namespace std;

namespace Cantera
{

OneDim::OneDim()
{
    m_newt = make_unique<MultiNewton>(1);
}

OneDim::OneDim(vector<shared_ptr<Domain1D>>& domains)
{
    // create a Newton iterator, and add each domain.
    m_newt = make_unique<MultiNewton>(1);
    m_state = make_shared<vector<double>>();
    for (auto& dom : domains) {
        addDomain(dom);
    }
    init();
    resize();
}

OneDim::~OneDim()
{
}

size_t OneDim::domainIndex(const string& name)
{
    for (size_t n = 0; n < m_dom.size(); n++) {
        if (domain(n).id() == name) {
            return n;
        }
    }
    throw CanteraError("OneDim::domainIndex","no domain named >>"+name+"<<");
}

std::tuple<string, size_t, string> OneDim::component(size_t i) {
    size_t n;
    for (n = nDomains()-1; n != npos; n--) {
        if (i >= start(n)) {
            break;
        }
    }
    Domain1D& dom = domain(n);
    size_t offset = i - start(n);
    size_t pt = offset / dom.nComponents();
    size_t comp = offset - pt*dom.nComponents();
    return make_tuple(dom.id(), pt, dom.componentName(comp));
}

void OneDim::addDomain(shared_ptr<Domain1D> d)
{
    // if 'd' is not the first domain, link it to the last domain
    // added (the rightmost one)
    size_t n = m_dom.size();
    if (n > 0) {
        m_dom.back()->append(d.get());
    }

    // every other domain is a connector
    if (n % 2 == 0) {
        m_connect.push_back(d);
    } else {
        m_bulk.push_back(d);
    }

    // add it also to the global domain list, and set its container and position
    m_dom.push_back(d);
    d->setData(m_state);
    d->setContainer(this, m_dom.size()-1);
    resize();
}

MultiJac& OneDim::jacobian()
{
    warn_deprecated("OneDim::jacobian",
                    "Replaced by getJacobian. To be removed after Cantera 3.2.");
    auto multijac = dynamic_pointer_cast<MultiJac>(m_jac);
    if (multijac) {
        return *multijac;
    } else {
        throw CanteraError("OneDim::jacobian", "Active Jacobian is not a MultiJac");
    }
}

MultiNewton& OneDim::newton()
{
    return *m_newt;
}

void OneDim::setLinearSolver(shared_ptr<SystemJacobian> solver)
{
    m_jac = solver;
    m_jac->initialize(size());
    m_jac->setBandwidth(bandwidth());
    m_jac->clearStats();
    m_jac_ok = false;
}

void OneDim::setJacAge(int ss_age, int ts_age)
{
    m_ss_jac_age = ss_age;
    if (ts_age > 0) {
        m_ts_jac_age = ts_age;
    } else {
        m_ts_jac_age = m_ss_jac_age;
    }
}

void OneDim::writeStats(int printTime)
{
    saveStats();
    writelog("\nStatistics:\n\n Grid   Timesteps  Functions      Time  Jacobians      Time\n");
    size_t n = m_gridpts.size();
    for (size_t i = 0; i < n; i++) {
        if (printTime) {
            writelog("{:5d}       {:5d}     {:6d} {:9.4f}      {:5d} {:9.4f}\n",
                     m_gridpts[i], m_timeSteps[i], m_funcEvals[i], m_funcElapsed[i],
                     m_jacEvals[i], m_jacElapsed[i]);
        } else {
            writelog("{:5d}       {:5d}     {:6d}        NA      {:5d}        NA\n",
                     m_gridpts[i], m_timeSteps[i], m_funcEvals[i], m_jacEvals[i]);
        }
    }
}

void OneDim::saveStats()
{
    if (m_jac) {
        int nev = m_jac->nEvals();
        if (nev > 0 && m_nevals > 0) {
            m_gridpts.push_back(m_pts);
            m_jacEvals.push_back(m_jac->nEvals());
            m_jacElapsed.push_back(m_jac->elapsedTime());
            m_funcEvals.push_back(m_nevals);
            m_nevals = 0;
            m_funcElapsed.push_back(m_evaltime);
            m_evaltime = 0.0;
            m_timeSteps.push_back(m_nsteps);
            m_nsteps = 0;
        }
    }
}

void OneDim::clearStats()
{
    m_gridpts.clear();
    m_jacEvals.clear();
    m_jacElapsed.clear();
    m_funcEvals.clear();
    m_funcElapsed.clear();
    m_timeSteps.clear();
    m_nevals = 0;
    m_evaltime = 0.0;
    m_nsteps = 0;
}

void OneDim::resize()
{
    m_bw = 0;
    m_nvars.clear();
    m_loc.clear();
    size_t lc = 0;

    // save the statistics for the last grid
    saveStats();
    m_pts = 0;
    for (size_t i = 0; i < nDomains(); i++) {
        const auto& d = m_dom[i];

        size_t np = d->nPoints();
        size_t nv = d->nComponents();
        for (size_t n = 0; n < np; n++) {
            m_nvars.push_back(nv);
            m_loc.push_back(lc);
            lc += nv;
            m_pts++;
        }

        // update the Jacobian bandwidth

        // bandwidth of the local block
        size_t bw1 = d->bandwidth();
        if (bw1 == npos) {
            bw1 = std::max<size_t>(2*d->nComponents(), 1) - 1;
        }
        m_bw = std::max(m_bw, bw1);

        // bandwidth of the block coupling the first point of this
        // domain to the last point of the previous domain
        if (i > 0) {
            size_t bw2 = m_dom[i-1]->bandwidth();
            if (bw2 == npos) {
                bw2 = m_dom[i-1]->nComponents();
            }
            bw2 += d->nComponents() - 1;
            m_bw = std::max(m_bw, bw2);
        }
        m_size = d->loc() + d->size();
    }

    m_state->resize(size());

    m_newt->resize(size());
    m_mask.resize(size());

    if (!m_jac) {
        m_jac = newSystemJacobian("banded-direct");
    }
    m_jac->initialize(size());
    m_jac->setBandwidth(bandwidth());
    m_jac->clearStats();
    m_jac_ok = false;
}

int OneDim::solve(double* x, double* xnew, int loglevel)
{
    if (!m_jac_ok) {
        evalJacobian(x);
        m_jac->updateTransient(m_rdt, m_mask.data());
        m_jac_ok = true;
    }
    return m_newt->solve(x, xnew, *this, loglevel);
}

void OneDim::evalSSJacobian(double* x, double* rsd)
{
    double rdt_save = m_rdt;
    m_jac_ok = false;
    setSteadyMode();
    evalJacobian(x);
    m_rdt = rdt_save;
}

Domain1D* OneDim::pointDomain(size_t i)
{
    Domain1D* d = right();
    while (d) {
        if (d->loc() <= i) {
            return d;
        }
        d = d->left();
    }
    return 0;
}

void OneDim::eval(size_t j, double* x, double* r, double rdt, int count)
{
    clock_t t0 = clock();
    if (m_interrupt) {
        m_interrupt->eval(m_nevals);
    }
    fill(r, r + m_size, 0.0);
    if (j == npos) {
        fill(m_mask.begin(), m_mask.end(), 0);
    }
    if (rdt < 0.0) {
        rdt = m_rdt;
    }

    // iterate over the bulk domains first
    for (const auto& d : m_bulk) {
        d->eval(j, x, r, m_mask.data(), rdt);
    }

    // then over the connector domains
    for (const auto& d : m_connect) {
        d->eval(j, x, r, m_mask.data(), rdt);
    }

    // increment counter and time
    if (count) {
        clock_t t1 = clock();
        m_evaltime += double(t1 - t0)/CLOCKS_PER_SEC;
        m_nevals++;
    }
}

void OneDim::evalJacobian(double* x0)
{
    m_jac->reset();
    clock_t t0 = clock();
    m_work1.resize(size());
    m_work2.resize(size());
    eval(npos, x0, m_work1.data(), 0.0, 0);
    size_t ipt = 0;
    for (size_t j = 0; j < points(); j++) {
        size_t nv = nVars(j);
        for (size_t n = 0; n < nv; n++) {
            // perturb x(n); preserve sign(x(n))
            double xsave = x0[ipt];
            double dx = fabs(xsave) * m_jacobianRelPerturb + m_jacobianAbsPerturb;
            if (xsave < 0) {
                dx = -dx;
            }
            x0[ipt] = xsave + dx;
            double rdx = 1.0 / (x0[ipt] - xsave);

            // calculate perturbed residual
            eval(j, x0, m_work2.data(), 0.0, 0);

            // compute nth column of Jacobian
            for (size_t i = j - 1; i != j+2; i++) {
                if (i != npos && i < points()) {
                    size_t mv = nVars(i);
                    size_t iloc = loc(i);
                    for (size_t m = 0; m < mv; m++) {
                        double delta = m_work2[m+iloc] - m_work1[m+iloc];
                        if (std::abs(delta) > m_jacobianThreshold || m+iloc == ipt) {
                            m_jac->setValue(m + iloc, ipt, delta * rdx);
                        }
                    }
                }
            }
            x0[ipt] = xsave;
            ipt++;
        }
    }

    m_jac->updateElapsed(double(clock() - t0) / CLOCKS_PER_SEC);
    m_jac->incrementEvals();
    m_jac->setAge(0);
}

double OneDim::ssnorm(double* x, double* r)
{
    eval(npos, x, r, 0.0, 0);
    double ss = 0.0;
    for (size_t i = 0; i < m_size; i++) {
        ss = std::max(fabs(r[i]),ss);
    }
    return ss;
}

void OneDim::initTimeInteg(double dt, double* x)
{
    double rdt_old = m_rdt;
    m_rdt = 1.0/dt;

    // if the stepsize has changed, then update the transient part of the
    // Jacobian
    if (fabs(rdt_old - m_rdt) > Tiny) {
        m_jac->updateTransient(m_rdt, m_mask.data());
    }

    // iterate over all domains, preparing each one to begin time stepping
    Domain1D* d = left();
    while (d) {
        d->initTimeInteg(dt, x);
        d = d->right();
    }
}

void OneDim::setSteadyMode()
{
    if (m_rdt == 0) {
        return;
    }

    m_rdt = 0.0;
    m_jac->updateTransient(m_rdt, m_mask.data());

    // iterate over all domains, preparing them for steady-state solution
    Domain1D* d = left();
    while (d) {
        d->setSteadyMode();
        d = d->right();
    }
}

void OneDim::init()
{
    if (!m_init) {
        Domain1D* d = left();
        while (d) {
            d->init();
            d = d->right();
        }
    }
    m_init = true;
}

double OneDim::timeStep(int nsteps, double dt, double* x, double* r, int loglevel)
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

void OneDim::resetBadValues(double* x)
{
    for (auto dom : m_dom) {
        dom->resetBadValues(x);
    }
}

}
