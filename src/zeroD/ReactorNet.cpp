//! @file ReactorNet.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/ReactorNet.h"
#include "cantera/zeroD/FlowDevice.h"
#include "cantera/zeroD/Wall.h"
#include "cantera/base/utilities.h"
#include "cantera/base/Array.h"
#include "cantera/numerics/Integrator.h"

namespace Cantera
{

ReactorNet::ReactorNet() :
    m_integ(newIntegrator("CVODE"))
{
    suppressErrors(true);

    // use backward differencing, with a full Jacobian computed
    // numerically, and use a Newton linear iterator
    m_integ->setMethod(BDF_Method);
    m_integ->setLinearSolverType("DENSE");
}

ReactorNet::~ReactorNet()
{
    // Defined in .cpp to limit dependence on Integrator.h
}

void ReactorNet::setInitialTime(double time)
{
    m_time = time;
    m_integrator_init = false;
}

void ReactorNet::setMaxTimeStep(double maxstep)
{
    m_maxstep = maxstep;
    m_init = false;
}

void ReactorNet::setMaxErrTestFails(int nmax)
{
    m_maxErrTestFails = nmax;
    m_init = false;
}

void ReactorNet::setTolerances(double rtol, double atol)
{
    if (rtol >= 0.0) {
        m_rtol = rtol;
    }
    if (atol >= 0.0) {
        m_atols = atol;
    }
    m_init = false;
}

void ReactorNet::setSensitivityTolerances(double rtol, double atol)
{
    if (rtol >= 0.0) {
        m_rtolsens = rtol;
    }
    if (atol >= 0.0) {
        m_atolsens = atol;
    }
    m_init = false;
}

void ReactorNet::initialize()
{
    m_nv = 0;
    debuglog("Initializing reactor network.\n", m_verbose);
    if (m_reactors.empty()) {
        throw CanteraError("ReactorNet::initialize",
                           "no reactors in network!");
    }
    m_start.assign(1, 0);
    for (size_t n = 0; n < m_reactors.size(); n++) {
        Reactor& r = *m_reactors[n];
        r.initialize(m_time);
        size_t nv = r.neq();
        m_nv += nv;
        m_start.push_back(m_nv);

        if (m_verbose) {
            writelog("Reactor {:d}: {:d} variables.\n", n, nv);
            writelog("              {:d} sensitivity params.\n", r.nSensParams());
        }
        if (r.type() == "FlowReactor" && m_reactors.size() > 1) {
            throw CanteraError("ReactorNet::initialize",
                               "FlowReactors must be used alone.");
        }
    }

    m_ydot.resize(m_nv,0.0);
    m_yest.resize(m_nv,0.0);
    m_advancelimits.resize(m_nv,-1.0);
    m_atol.resize(neq());
    fill(m_atol.begin(), m_atol.end(), m_atols);
    m_integ->setTolerances(m_rtol, neq(), m_atol.data());
    m_integ->setSensitivityTolerances(m_rtolsens, m_atolsens);
    m_integ->setMaxStepSize(m_maxstep);
    m_integ->setMaxErrTestFails(m_maxErrTestFails);
    if (m_verbose) {
        writelog("Number of equations: {:d}\n", neq());
        writelog("Maximum time step:   {:14.6g}\n", m_maxstep);
    }
    m_integ->initialize(m_time, *this);
    if (m_integ->preconditionerSide() != PreconditionerSide::NO_PRECONDITION) {
        checkPreconditionerSupported();
    }
    m_integrator_init = true;
    m_init = true;
}

void ReactorNet::reinitialize()
{
    if (m_init) {
        debuglog("Re-initializing reactor network.\n", m_verbose);
        m_integ->reinitialize(m_time, *this);
        if (m_integ->preconditionerSide() != PreconditionerSide::NO_PRECONDITION) {
            checkPreconditionerSupported();
        }
        m_integrator_init = true;
    } else {
        initialize();
    }
}

void ReactorNet::setLinearSolverType(const std::string& linSolverType)
{
    m_integ->setLinearSolverType(linSolverType);
    m_integrator_init = false;
}

void ReactorNet::setPreconditioner(shared_ptr<PreconditionerBase> preconditioner)
{
    m_integ->setPreconditioner(preconditioner);
    m_integrator_init = false;
}

void ReactorNet::setMaxSteps(int nmax)
{
    m_integ->setMaxSteps(nmax);
}

int ReactorNet::maxSteps()
{
    return m_integ->maxSteps();
}

void ReactorNet::advance(doublereal time)
{
    if (!m_init) {
        initialize();
    } else if (!m_integrator_init) {
        reinitialize();
    }
    m_integ->integrate(time);
    m_time = time;
    updateState(m_integ->solution());
}

double ReactorNet::advance(double time, bool applylimit)
{
    if (!m_init) {
        initialize();
    } else if (!m_integrator_init) {
        reinitialize();
    }

    if (!applylimit) {
        // take full step
        advance(time);
        return time;
    }

    if (!hasAdvanceLimits()) {
        // take full step
        advance(time);
        return time;
    }

    getAdvanceLimits(m_advancelimits.data());

    // ensure that gradient is available
    while (lastOrder() < 1) {
        step();
    }

    int k = lastOrder();
    double t = time, delta;
    double* y = m_integ->solution();

    // reduce time step if limits are exceeded
    while (true) {
        bool exceeded = false;
        getEstimate(t, k, &m_yest[0]);
        for (size_t j = 0; j < m_nv; j++) {
            delta = abs(m_yest[j] - y[j]);
            if ( (m_advancelimits[j] > 0.) && ( delta > m_advancelimits[j]) ) {
                exceeded = true;
                if (m_verbose) {
                    writelog("    Limiting global state vector component {:d} (dt = {:9.4g}):"
                             "{:11.6g} > {:9.4g}\n",
                             j, t - m_time, delta, m_advancelimits[j]);
                }
            }
        }
        if (!exceeded) {
            break;
        }
        t = .5 * (m_time + t);
    }
    advance(t);
    return t;
}

double ReactorNet::step()
{
    if (!m_init) {
        initialize();
    } else if (!m_integrator_init) {
        reinitialize();
    }
    m_time = m_integ->step(m_time + 1.0);
    updateState(m_integ->solution());
    return m_time;
}

void ReactorNet::getEstimate(double time, int k, double* yest)
{
    // initialize
    double* cvode_dky = m_integ->solution();
    for (size_t j = 0; j < m_nv; j++) {
        yest[j] = cvode_dky[j];
    }

    // Taylor expansion
    double factor = 1.;
    double deltat = time - m_time;
    for (int n = 1; n <= k; n++) {
        factor *= deltat / n;
        cvode_dky = m_integ->derivative(m_time, n);
        for (size_t j = 0; j < m_nv; j++) {
            yest[j] += factor * cvode_dky[j];
        }
    }
}

int ReactorNet::lastOrder()
{
    return m_integ->lastOrder();
}

void ReactorNet::addReactor(Reactor& r)
{
    r.setNetwork(this);
    m_reactors.push_back(&r);
}

void ReactorNet::eval(doublereal t, doublereal* y,
                      doublereal* ydot, doublereal* p)
{
    m_time = t;
    updateState(y);
    m_LHS.assign(m_nv, 1);
    m_RHS.assign(m_nv, 0);
    for (size_t n = 0; n < m_reactors.size(); n++) {
        m_reactors[n]->applySensitivity(p);
        m_reactors[n]->eval(t, m_LHS.data() + m_start[n], m_RHS.data() + m_start[n]);
        size_t yEnd = 0;
        if (n == m_reactors.size() - 1) {
            yEnd = m_RHS.size();
        } else {
            yEnd = m_start[n + 1];
        }
        for (size_t i = m_start[n]; i < yEnd; i++) {
            ydot[i] = m_RHS[i] / m_LHS[i];
        }
        m_reactors[n]->resetSensitivity(p);
    }
    checkFinite("ydot", ydot, m_nv);
}

double ReactorNet::sensitivity(size_t k, size_t p)
{
    if (!m_init) {
        initialize();
    }
    if (p >= m_sens_params.size()) {
        throw IndexError("ReactorNet::sensitivity",
                         "m_sens_params", p, m_sens_params.size()-1);
    }
    double denom = m_integ->solution(k);
    if (denom == 0.0) {
        denom = SmallNumber;
    }
    return m_integ->sensitivity(k, p) / denom;
}

void ReactorNet::evalJacobian(doublereal t, doublereal* y,
                              doublereal* ydot, doublereal* p, Array2D* j)
{
    //evaluate the unperturbed ydot
    eval(t, y, ydot, p);
    for (size_t n = 0; n < m_nv; n++) {
        // perturb x(n)
        double ysave = y[n];
        double dy = m_atol[n] + fabs(ysave)*m_rtol;
        y[n] = ysave + dy;
        dy = y[n] - ysave;

        // calculate perturbed residual
        eval(t, y, m_ydot.data(), p);

        // compute nth column of Jacobian
        for (size_t m = 0; m < m_nv; m++) {
            j->value(m,n) = (m_ydot[m] - ydot[m])/dy;
        }
        y[n] = ysave;
    }
}

void ReactorNet::updateState(doublereal* y)
{
    checkFinite("y", y, m_nv);
    for (size_t n = 0; n < m_reactors.size(); n++) {
        m_reactors[n]->updateState(y + m_start[n]);
    }
}

void ReactorNet::getState(double* y)
{
    for (size_t n = 0; n < m_reactors.size(); n++) {
        m_reactors[n]->getState(y + m_start[n]);
    }
}

void ReactorNet::getDerivative(int k, double* dky)
{
    double* cvode_dky = m_integ->derivative(m_time, k);
    for (size_t j = 0; j < m_nv; j++) {
        dky[j] = cvode_dky[j];
    }
}

void ReactorNet::setAdvanceLimits(const double *limits)
{
    if (!m_init) {
        initialize();
    }
    for (size_t n = 0; n < m_reactors.size(); n++) {
        m_reactors[n]->setAdvanceLimits(limits + m_start[n]);
    }
}

bool ReactorNet::hasAdvanceLimits()
{
    bool has_limit = false;
    for (size_t n = 0; n < m_reactors.size(); n++) {
        has_limit |= m_reactors[n]->hasAdvanceLimits();
    }
    return has_limit;
}

bool ReactorNet::getAdvanceLimits(double *limits)
{
    bool has_limit = false;
    for (size_t n = 0; n < m_reactors.size(); n++) {
        has_limit |= m_reactors[n]->getAdvanceLimits(limits + m_start[n]);
    }
    return has_limit;
}

size_t ReactorNet::globalComponentIndex(const string& component, size_t reactor)
{
    if (!m_init) {
        initialize();
    }
    return m_start[reactor] + m_reactors[reactor]->componentIndex(component);
}

std::string ReactorNet::componentName(size_t i) const
{
    for (auto r : m_reactors) {
        if (i < r->neq()) {
            return r->name() + ": " + r->componentName(i);
        } else {
            i -= r->neq();
        }
    }
    throw CanteraError("ReactorNet::componentName", "Index out of bounds");
}

size_t ReactorNet::registerSensitivityParameter(
    const std::string& name, double value, double scale)
{
    if (m_integrator_init) {
        throw CanteraError("ReactorNet::registerSensitivityParameter",
                           "Sensitivity parameters cannot be added after the"
                           "integrator has been initialized.");
    }
    m_paramNames.push_back(name);
    m_sens_params.push_back(value);
    m_paramScales.push_back(scale);
    return m_sens_params.size() - 1;
}

void ReactorNet::setDerivativeSettings(AnyMap& settings)
{
    // Apply given settings to all reactors
    for (size_t i = 0; i < m_reactors.size(); i++) {
        m_reactors[i]->setDerivativeSettings(settings);
    }
}

AnyMap ReactorNet::solverStats() const
{
    return m_integ->solverStats();
}

std::string ReactorNet::linearSolverType() const
{
    return m_integ->linearSolverType();
}

void ReactorNet::preconditionerSolve(double* rhs, double* output)
{
    m_integ->preconditionerSolve(m_nv, rhs, output);
}

void ReactorNet::preconditionerSetup(double t, double* y, double gamma)
{
    // ensure state is up to date.
    updateState(y);
    // get the preconditioner
    auto precon = m_integ->preconditioner();
    // Reset preconditioner
    precon->reset();
    // Set gamma value for M =I - gamma*J
    precon->setGamma(gamma);
    // Make a copy of state to adjust it for preconditioner
    vector_fp yCopy(m_nv);
    // Get state of reactor
    getState(yCopy.data());
    // transform state based on preconditioner rules
    precon->stateAdjustment(yCopy);
    // update network with adjusted state
    updateState(yCopy.data());
    // Get jacobians and give elements to preconditioners
    for (size_t i = 0; i < m_reactors.size(); i++) {
        Eigen::SparseMatrix<double> rJac = m_reactors[i]->jacobian();
        for (int k=0; k<rJac.outerSize(); ++k) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(rJac, k); it; ++it) {
                precon->setValue(it.row() + m_start[i], it.col() + m_start[i],
                    it.value());
            }
        }
    }
    // post reactor setup operations
    precon->setup();
}

void ReactorNet::checkPreconditionerSupported()
{
    // preconditioner currently not supported for surfaces
    for (size_t i = 0; i < m_reactors.size(); i++) {
        if (m_reactors[i]->nSurfs() > 0) {
            throw CanteraError("ReactorNet::checkPreconditionerSupported",
                "Preconditioning is not supported for networks with surfaces.");
        }
    }
}

void ReactorNet::updatePreconditioner(double gamma)
{
    auto precon = m_integ->preconditioner();
    precon->setGamma(gamma);
    precon->updatePreconditioner();
}

}
