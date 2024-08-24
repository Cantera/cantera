//! @file ReactorNet.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/ReactorNet.h"
#include "cantera/zeroD/FlowDevice.h"
#include "cantera/zeroD/ReactorSurface.h"
#include "cantera/zeroD/Wall.h"
#include "cantera/base/utilities.h"
#include "cantera/base/Array.h"
#include "cantera/numerics/Integrator.h"
#include "cantera/zeroD/FlowReactor.h"

#include <cstdio>

namespace Cantera
{

ReactorNet::ReactorNet()
{
    suppressErrors(true);
}

ReactorNet::~ReactorNet()
{
}

void ReactorNet::setInitialTime(double time)
{
    m_time = time;
    m_initial_time = time;
    m_integrator_init = false;
}

void ReactorNet::setMaxTimeStep(double maxstep)
{
    m_maxstep = maxstep;
    integrator().setMaxStepSize(m_maxstep);
}

void ReactorNet::setMaxErrTestFails(int nmax)
{
    integrator().setMaxErrTestFails(nmax);
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

double ReactorNet::time() {
    if (m_timeIsIndependent) {
        return m_time;
    } else {
        throw CanteraError("ReactorNet::time", "Time is not the independent variable"
            " for this reactor network.");
    }
}

double ReactorNet::distance() {
    if (!m_timeIsIndependent) {
        return m_time;
    } else {
        throw CanteraError("ReactorNet::distance", "Distance is not the independent"
            " variable for this reactor network.");
    }
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
    if (!m_linearSolverType.empty()) {
        m_integ->setLinearSolverType(m_linearSolverType);
    }
    if (m_precon) {
        m_integ->setPreconditioner(m_precon);
    }
    m_integ->initialize(m_time, *this);
    if (m_verbose) {
        writelog("Number of equations: {:d}\n", neq());
        writelog("Maximum time step:   {:14.6g}\n", m_maxstep);
    }
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

void ReactorNet::setLinearSolverType(const string& linSolverType)
{
    m_linearSolverType = linSolverType;
    m_integrator_init = false;
}

void ReactorNet::setPreconditioner(shared_ptr<PreconditionerBase> preconditioner)
{
    m_precon = preconditioner;
    m_integrator_init = false;
}

void ReactorNet::setMaxSteps(int nmax)
{
    integrator().setMaxSteps(nmax);
}

int ReactorNet::maxSteps()
{
    return integrator().maxSteps();
}

void ReactorNet::advance(double time)
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
    if (!m_init) {
        initialize();
    }
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

int ReactorNet::lastOrder() const
{
    if (m_integ) {
        return m_integ->lastOrder();
    } else {
        return 0;
    }
}

void ReactorNet::addReactor(Reactor& r)
{
    warn_deprecated("ReactorNet::addReactor",
        "To be removed after Cantera 3.1. Replaced by version using shared pointer.");
    for (auto current : m_reactors) {
        if (current->isOde() != r.isOde()) {
            throw CanteraError("ReactorNet::addReactor",
                "Cannot mix Reactor types using both ODEs and DAEs ({} and {})",
                current->type(), r.type());
        }
        if (current->timeIsIndependent() != r.timeIsIndependent()) {
            throw CanteraError("ReactorNet::addReactor",
                "Cannot mix Reactor types using time and space as independent variables"
                "\n({} and {})", current->type(), r.type());
        }
    }
    m_timeIsIndependent = r.timeIsIndependent();
    r.setNetwork(this);
    m_reactors.push_back(&r);
    if (!m_integ) {
        m_integ.reset(newIntegrator(r.isOde() ? "CVODE" : "IDA"));
        // use backward differencing, with a full Jacobian computed
        // numerically, and use a Newton linear iterator
        m_integ->setMethod(BDF_Method);
        m_integ->setLinearSolverType("DENSE");
    }
    updateNames(r);
}

void ReactorNet::addReactor(shared_ptr<ReactorNode> node)
{
    auto reactor = std::dynamic_pointer_cast<Reactor>(node);
    if (!reactor) {
        throw CanteraError("ReactorNet::addReactor",
            "Invalid reactor type {}", node->type());
    }
    for (auto current : m_reactors) {
        if (current->isOde() != reactor->isOde()) {
            throw CanteraError("ReactorNet::addReactor",
                "Cannot mix Reactor types using both ODEs and DAEs ({} and {})",
                current->type(), reactor->type());
        }
        if (current->timeIsIndependent() != reactor->timeIsIndependent()) {
            throw CanteraError("ReactorNet::addReactor",
                "Cannot mix Reactor types using time and space as independent variables"
                "\n({} and {})", current->type(), reactor->type());
        }
    }
    m_timeIsIndependent = reactor->timeIsIndependent();
    reactor->setNetwork(this);
    m_reactors.push_back(reactor.get());
    if (!m_integ) {
        m_integ.reset(newIntegrator(reactor->isOde() ? "CVODE" : "IDA"));
        // use backward differencing, with a full Jacobian computed
        // numerically, and use a Newton linear iterator
        m_integ->setMethod(BDF_Method);
        m_integ->setLinearSolverType("DENSE");
    }
    updateNames(*reactor);
}

void ReactorNet::updateNames(Reactor& r)
{
    // ensure that reactors and components have reproducible names
    r.setDefaultName(m_counts);

    for (size_t i=0; i<r.nWalls(); i++) {
        auto& w = r.wall(i);
        w.setDefaultName(m_counts);
        if (w.left().type() == "Reservoir") {
            w.left().setDefaultName(m_counts);
        }
        if (w.right().type() == "Reservoir") {
            w.right().setDefaultName(m_counts);
        }
    }

    for (size_t i=0; i<r.nInlets(); i++) {
        auto& in = r.inlet(i);
        in.setDefaultName(m_counts);
        if (in.in().type() == "Reservoir") {
            in.in().setDefaultName(m_counts);
        }
    }

    for (size_t i=0; i<r.nOutlets(); i++) {
        auto& out = r.outlet(i);
        out.setDefaultName(m_counts);
        if (out.out().type() == "Reservoir") {
            out.out().setDefaultName(m_counts);
        }
    }

    for (size_t i=0; i<r.nSurfs(); i++) {
        r.surface(i)->setDefaultName(m_counts);
    }
}

Integrator& ReactorNet::integrator() {
    if (m_integ == nullptr) {
        throw CanteraError("ReactorNet::integrator",
            "Integrator has not been instantiated. Add one or more reactors first.");
    }
    return *m_integ;
}

void ReactorNet::eval(double t, double* y, double* ydot, double* p)
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

void ReactorNet::evalDae(double t, double* y, double* ydot, double* p, double* residual)
{
    m_time = t;
    updateState(y);
    for (size_t n = 0; n < m_reactors.size(); n++) {
        m_reactors[n]->applySensitivity(p);
        m_reactors[n]->evalDae(t, y, ydot, residual);
        m_reactors[n]->resetSensitivity(p);
    }
    checkFinite("ydot", ydot, m_nv);
}

void ReactorNet::getConstraints(double* constraints)
{
    for (size_t n = 0; n < m_reactors.size(); n++) {
        m_reactors[n]->getConstraints(constraints + m_start[n]);
    }
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

void ReactorNet::evalJacobian(double t, double* y, double* ydot, double* p, Array2D* j)
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

void ReactorNet::updateState(double* y)
{
    checkFinite("y", y, m_nv);
    for (size_t n = 0; n < m_reactors.size(); n++) {
        m_reactors[n]->updateState(y + m_start[n]);
    }
}

void ReactorNet::getDerivative(int k, double* dky)
{
    if (!m_init) {
        initialize();
    }
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

bool ReactorNet::hasAdvanceLimits() const
{
    bool has_limit = false;
    for (size_t n = 0; n < m_reactors.size(); n++) {
        has_limit |= m_reactors[n]->hasAdvanceLimits();
    }
    return has_limit;
}

bool ReactorNet::getAdvanceLimits(double *limits) const
{
    bool has_limit = false;
    for (size_t n = 0; n < m_reactors.size(); n++) {
        has_limit |= m_reactors[n]->getAdvanceLimits(limits + m_start[n]);
    }
    return has_limit;
}

void ReactorNet::getState(double* y)
{
    for (size_t n = 0; n < m_reactors.size(); n++) {
        m_reactors[n]->getState(y + m_start[n]);
    }
}

void ReactorNet::getStateDae(double* y, double* ydot)
{
    for (size_t n = 0; n < m_reactors.size(); n++) {
        m_reactors[n]->getStateDae(y + m_start[n], ydot + m_start[n]);
    }
}

size_t ReactorNet::globalComponentIndex(const string& component, size_t reactor)
{
    if (!m_init) {
        initialize();
    }
    return m_start[reactor] + m_reactors[reactor]->componentIndex(component);
}

string ReactorNet::componentName(size_t i) const
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
    const string& name, double value, double scale)
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
    if (m_integ) {
        return m_integ->solverStats();
    } else {
        return AnyMap();
    }
}

string ReactorNet::linearSolverType() const
{
    if (m_integ) {
        return m_integ->linearSolverType();
    } else {
        return "";
    }
}

void ReactorNet::preconditionerSolve(double* rhs, double* output)
{
    if (!m_integ) {
        throw CanteraError("ReactorNet::preconditionerSolve",
                           "Must only be called after ReactorNet is initialized.");
    }
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
    vector<double> yCopy(m_nv);
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

void ReactorNet::updatePreconditioner(double gamma)
{
    if (!m_integ) {
        throw CanteraError("ReactorNet::updatePreconditioner",
                           "Must only be called after ReactorNet is initialized.");
    }
    auto precon = m_integ->preconditioner();
    precon->setGamma(gamma);
    precon->updatePreconditioner();
}

void ReactorNet::checkPreconditionerSupported() const {
    // check for non-mole-based reactors and throw an error otherwise
    for (auto reactor : m_reactors) {
        if (!reactor->preconditionerSupported()) {
            throw CanteraError("ReactorNet::checkPreconditionerSupported",
                "Preconditioning is only supported for type *MoleReactor,\n"
                "Reactor type given: '{}'.",
                reactor->type());
        }
    }
}

}
