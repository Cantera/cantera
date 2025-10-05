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

OneDim::OneDim(vector<shared_ptr<Domain1D>>& domains)
{
    for (auto& dom : domains) {
        addDomain(dom);
    }
    init();
    resize();
}

size_t OneDim::domainIndex(const string& name) const
{
    for (size_t n = 0; n < m_dom.size(); n++) {
        if (domain(n).id() == name) {
            return n;
        }
    }
    throw CanteraError("OneDim::domainIndex","no domain named >>"+name+"<<");
}

std::tuple<string, size_t, string> OneDim::component(size_t i) const {
    const auto& [n, j, k] = m_componentInfo[i];
    Domain1D& dom = domain(n);
    return make_tuple(dom.id(), j, dom.componentName(k));
}

string OneDim::componentName(size_t i) const {
    const auto& [dom, pt, comp] = component(i);
    return fmt::format("domain {}, component {} at point {}", dom, comp, pt);
}

pair<string, string> OneDim::componentTableHeader() const
{
    return {"", "Domain   Pt. Component"};
}

string OneDim::componentTableLabel(size_t i) const
{
    const auto& [dom, pt, comp] = component(i);
    return fmt::format("{:8s} {:3d} {:<12s}", dom, pt, comp);
}

double OneDim::upperBound(size_t i) const
{
    const auto& [n, j, k] = m_componentInfo[i];
    Domain1D& dom = domain(n);
    return dom.upperBound(k);
}

double OneDim::lowerBound(size_t i) const
{
    const auto& [n, j, k] = m_componentInfo[i];
    Domain1D& dom = domain(n);
    return dom.lowerBound(k);
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

double OneDim::weightedNorm(const double* step) const
{
    double sum = 0.0;
    const double* x = m_state->data();
    size_t nd = nDomains();
    for (size_t n = 0; n < nd; n++) {
        Domain1D& dom = domain(n);
        double d_sum = 0.0;
        size_t nv = dom.nComponents();
        size_t np = dom.nPoints();
        size_t dstart = start(n);

        for (size_t i = 0; i < nv; i++) {
            double esum = 0.0;
            for (size_t j = 0; j < np; j++) {
                esum += fabs(x[dstart + nv*j + i]);
            }
            double ewt = dom.rtol(i)*esum/np + dom.atol(i);
            for (size_t j = 0; j < np; j++) {
                double f = step[dstart + nv*j + i]/ewt;
                d_sum += f*f;
            }
        }
        sum += d_sum;
    }
    return sqrt(sum / size());
}

MultiJac& OneDim::jacobian()
{
    warn_deprecated("OneDim::jacobian",
                    "Replaced by linearSolver(). To be removed after Cantera 3.2.");
    auto multijac = dynamic_pointer_cast<MultiJac>(m_jac);
    if (multijac) {
        return *multijac;
    } else {
        throw CanteraError("OneDim::jacobian", "Active Jacobian is not a MultiJac");
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
    m_componentInfo.clear();
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
            for (size_t k = 0; k < nv; k++) {
                m_componentInfo.emplace_back(i, n, k);
            }
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
    if (!m_jac) {
        m_jac = newSystemJacobian("banded-direct");
    }
    SteadyStateSystem::resize();
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

void OneDim::initTimeInteg(double dt, double* x)
{
    SteadyStateSystem::initTimeInteg(dt, x);
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
    SteadyStateSystem::setSteadyMode();
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

void OneDim::resetBadValues(double* x)
{
    for (auto dom : m_dom) {
        dom->resetBadValues(x);
    }
}

}
