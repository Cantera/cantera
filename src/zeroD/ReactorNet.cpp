//! @file ReactorNet.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/ReactorNet.h"
#include "cantera/zeroD/FlowDevice.h"
#include "cantera/zeroD/ReactorSurface.h"
#include "cantera/zeroD/Wall.h"
#include "cantera/zeroD/Reservoir.h"
#include "cantera/base/utilities.h"
#include "cantera/base/Array.h"
#include "cantera/base/Solution.h"
#include "cantera/numerics/Integrator.h"
#include "cantera/zeroD/FlowReactor.h"
#include "cantera/numerics/SystemJacobianFactory.h"
#include "cantera/numerics/EigenSparseJacobian.h"
#include "cantera/oneD/MultiJac.h"
#include "cantera/oneD/MultiNewton.h"

#include <cstdio>
#include <deque>

namespace Cantera
{

ReactorNet::~ReactorNet()
{
}

ReactorNet::ReactorNet(shared_ptr<ReactorBase> reactor)
    : ReactorNet(span<shared_ptr<ReactorBase>>{&reactor, 1})
{
}

ReactorNet::ReactorNet(span<shared_ptr<ReactorBase>> reactors)
{
    if (reactors.empty()) {
        throw CanteraError("ReactorNet::ReactorNet", "No reactors in network!");
    }

    suppressErrors(true);
    bool isOde = true;

    // Names of Reactors and ReactorSurfaces using each Solution; should be only one
    // reactor per Solution.
    map<Solution*, set<string>> solutions;

    // All ReactorBase objects connected to this network, starting with the ones
    // explicitly included.
    std::deque<shared_ptr<ReactorBase>> reactorQueue(reactors.begin(), reactors.end());
    std::set<shared_ptr<ReactorBase>> visited;

    while (!reactorQueue.empty()) {
        auto R = reactorQueue.front();
        reactorQueue.pop_front();

        if (visited.find(R) != visited.end()) {
            continue;
        }
        visited.insert(R);

        if (auto bulk = std::dynamic_pointer_cast<Reactor>(R)) {
            m_bulkReactors.push_back(bulk);
            m_reactors.push_back(R);
        } else if (auto surf = std::dynamic_pointer_cast<ReactorSurface>(R)) {
            m_surfaces.push_back(surf);
            m_reactors.push_back(R);
            for (size_t i = 0; i < surf->nAdjacent(); i++) {
                reactorQueue.push_back(surf->adjacent(i));
            }
        } else if (auto res = std::dynamic_pointer_cast<Reservoir>(R)) {
            m_reservoirs.push_back(res);
        }

        // Discover reservoirs, flow devices, and walls connected to this reactor
        for (size_t i = 0; i < R->nInlets(); i++) {
            auto& inlet = R->inlet(i);
            m_flowDevices.insert(&inlet);
            reactorQueue.push_back(inlet.in().shared_from_this());
        }
        for (size_t i = 0; i < R->nOutlets(); i++) {
            auto& outlet = R->outlet(i);
            m_flowDevices.insert(&outlet);
            reactorQueue.push_back(outlet.out().shared_from_this());
        }
        for (size_t i = 0; i < R->nWalls(); i++) {
            auto& wall = R->wall(i);
            m_walls.insert(&wall);
            reactorQueue.push_back(wall.left().shared_from_this());
            reactorQueue.push_back(wall.right().shared_from_this());
        }

        auto phase = R->phase();
        for (size_t i = 0; i < R->nSurfs(); i++) {
            if (R->surface(i)->phase()->adjacent(phase->name()) != phase) {
                throw CanteraError("ReactorNet::initialize",
                    "Bulk phase '{}' used by interface '{}' must be the same object\n"
                    "as the contents of the adjacent reactor '{}'.",
                    phase->name(), R->surface(i)->name(), R->name());
            }
            reactorQueue.push_back(R->surface(i)->shared_from_this());
        }
        R->setNetwork(this);
        updateNames(*R);
        solutions[phase.get()].insert(R->name());
    }

    for (auto& R : m_bulkReactors) {
        if (R->type() == "FlowReactor" && m_bulkReactors.size() > 1) {
            throw CanteraError("ReactorNet::initialize",
                               "FlowReactors must be used alone.");
        }
        m_timeIsIndependent = R->timeIsIndependent();
        R->initialize(m_time);
        isOde = R->isOde();
        R->setOffset(m_nv);
        size_t nv = R->neq();
        m_nv += nv;

        for (auto current : m_bulkReactors) {
            if (current->isOde() != R->isOde()) {
                throw CanteraError("ReactorNet::ReactorNet",
                    "Cannot mix Reactor types using both ODEs and DAEs ({} and {})",
                    current->type(), R->type());
            }
            if (current->timeIsIndependent() != R->timeIsIndependent()) {
                throw CanteraError("ReactorNet::ReactorNet",
                    "Cannot mix Reactor types using time and space as independent variables"
                    "\n({} and {})", current->type(), R->type());
            }
        }
    }

    for (auto surf : m_surfaces) {
        surf->setOffset(m_nv);
        surf->initialize(m_time);
        m_nv += surf->neq();
        solutions[surf->phase().get()].insert(surf->name());
    }

    for (auto& [soln, reactors] : solutions) {
        if (reactors.size() > 1) {
            string shared;
            for (auto r : reactors) {
                if (r != *reactors.begin()) {
                    shared += ", ";
                }
                shared += fmt::format("'{}'", r);
            }
            throw CanteraError("ReactorNet::initialize", "The following reactors /"
                " reactor surfaces are using the same Solution object: {}. Use"
                " independent Solution objects or set the 'clone' argument to 'true'"
                " when creating the Reactor or ReactorSurface objects.", shared);
        }
    }

    m_ydot.resize(m_nv, 0.0);
    m_yest.resize(m_nv, 0.0);
    m_advancelimits.resize(m_nv, -1.0);
    m_atol.resize(m_nv);

    m_integ.reset(newIntegrator(isOde ? "CVODE" : "IDA"));
    // use backward differencing, with a full Jacobian computed
    // numerically, and use a Newton linear iterator
    m_integ->setMethod(BDF_Method);
    m_integ->setLinearSolverType("DENSE");
    setTolerances(-1.0, -1.0);
}

void ReactorNet::setInitialTime(double time)
{
    m_time = time;
    m_initial_time = time;
    m_needIntegratorInit = true;
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
    fill(m_atol.begin(), m_atol.end(), m_atols);
    m_integ->setTolerances(m_rtol, neq(), m_atol.data());
}

void ReactorNet::setSensitivityTolerances(double rtol, double atol)
{
    if (rtol >= 0.0) {
        m_rtolsens = rtol;
    }
    if (atol >= 0.0) {
        m_atolsens = atol;
    }
    m_integ->setSensitivityTolerances(m_rtolsens, m_atolsens);
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
    if (m_integratorInitialized) {
        if (m_needIntegratorInit) {
            reinitialize();
        }
        return;
    }
    if (!m_linearSolverType.empty()) {
        m_integ->setLinearSolverType(m_linearSolverType);
    }
    if (m_precon) {
        m_integ->setPreconditioner(m_precon);
    }
    for (auto& reactor : m_reactors) {
        reactor->updateConnected(true);
    }
    m_integ->initialize(m_time, *this);
    if (m_verbose) {
        writelog("Number of equations: {:d}\n", neq());
        writelog("Maximum time step:   {:14.6g}\n", m_maxstep);
    }
    if (m_integ->preconditionerSide() != PreconditionerSide::NO_PRECONDITION) {
        checkPreconditionerSupported();
    }
    m_needIntegratorInit = false;
    m_integratorInitialized = true;
}

void ReactorNet::reinitialize()
{
    if (m_integratorInitialized) {
        debuglog("Re-initializing reactor network.\n", m_verbose);
        for (auto& reactor : m_reactors) {
            reactor->updateConnected(true);
        }
        m_integ->reinitialize(m_time, *this);
        if (m_integ->preconditionerSide() != PreconditionerSide::NO_PRECONDITION) {
            checkPreconditionerSupported();
        }
        m_needIntegratorInit = false;
    } else {
        initialize();
    }
}

void ReactorNet::setLinearSolverType(const string& linSolverType)
{
    m_linearSolverType = linSolverType;
    m_needIntegratorInit = true;
}

void ReactorNet::setPreconditioner(shared_ptr<SystemJacobian> preconditioner)
{
    m_precon = preconditioner;
    m_needIntegratorInit = true;
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
    initialize();
    m_integ->integrate(time);
    m_time = m_integ->currentTime();
    updateState(span<const double>(m_integ->solution(), m_nv));
}

double ReactorNet::advance(double time, bool applylimit)
{
    initialize();
    if (!applylimit || !hasAdvanceLimits()) {
        // No limit enforcement requested; integrate to requested time
        advance(time);
        return time;
    }

    // Enable root-based limit detection and set the base state to the current state
    m_ybase.assign(m_nv, 0.0);
    getState(m_ybase);
    m_ybase_time = m_time;
    m_limit_check_active = true;
    m_integ->setRootFunctionCount(nRootFunctions());

    // Integrate toward the requested time; integrator will return early if a limit is
    // reached (CV_ROOT_RETURN). The try/catch ensures the temporary root-finding state
    // is cleared even when CVODE throws so subsequent calls start clean.
    try {
        m_integ->integrate(time);
    } catch (...) {
        m_limit_check_active = false;
        m_integ->setRootFunctionCount(nRootFunctions());
        throw;
    }
    m_time = m_integ->currentTime();

    // Update reactor states to match the integrator solution at the time reached
    // (which may be earlier than 'time' if a limit was triggered)
    updateState(span<const double>(m_integ->solution(), m_nv));

    // Disable limit checking after this call
    m_limit_check_active = false;
    m_integ->setRootFunctionCount(nRootFunctions());

    // When a root event stopped integration before reaching the requested time, report
    // the most limiting component and details about the step.
    if (m_verbose && m_time < time) {
        // Ensure limits are available
        if (m_advancelimits.size() != m_nv) {
            m_advancelimits.assign(m_nv, -1.0);
        }
        getAdvanceLimits(m_advancelimits);
        double* ycurr = m_integ->solution();
        size_t jmax = npos;
        double max_ratio = -1.0;
        double best_limit = 0.0;
        for (size_t j = 0; j < m_nv; j++) {
            double lim = m_advancelimits[j];
            if (lim > 0.0) {
                double delta = std::abs(ycurr[j] - m_ybase[j]);
                double ratio = delta / lim;
                if (ratio > max_ratio) {
                    max_ratio = ratio;
                    jmax = j;
                    best_limit = lim;
                }
            }
        }
        if (jmax != npos) {
            double dt = m_time - m_ybase_time;
            double y_start = m_ybase[jmax];
            double y_end = ycurr[jmax];
            double delta = y_end - y_start;
            writelog("    Advance limit triggered for component {:d} (dt = {:9.4g}):"
                     " y_start = {:11.6g}, y_end = {:11.6g},"
                     " delta = {:11.6g}, limit = {:9.4g}\n",
                     jmax, dt, y_start, y_end, delta, best_limit);
        }
    }

    // m_time is tracked via callbacks during integration
    return m_time;
}

double ReactorNet::step()
{
    initialize();
    m_time = m_integ->step(m_time + 1.0);
    updateState(span<const double>(m_integ->solution(), m_nv));
    return m_time;
}

void ReactorNet::solveSteady(int loglevel)
{
    initialize();
    vector<double> y(neq());
    getState(y);
    SteadyReactorSolver solver(this, y.data());
    solver.setMaxTimeStepCount(maxSteps());
    solver.solve(loglevel);
    solver.getState(y.data());
    updateState(y);
}

Eigen::SparseMatrix<double> ReactorNet::steadyJacobian(double rdt)
{
    initialize();
    vector<double> y0(neq());
    vector<double> y1(neq());
    getState(y0);
    SteadyReactorSolver solver(this, y0.data());
    solver.evalJacobian(y0.data());
    if (rdt) {
        solver.linearSolver()->updateTransient(rdt, solver.transientMask().data());
    }
    return std::dynamic_pointer_cast<EigenSparseJacobian>(solver.linearSolver())->jacobian();
}

void ReactorNet::getEstimate(double time, int k, span<double> yest)
{
    initialize();
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

size_t ReactorNet::nRootFunctions() const
{
    return (m_limit_check_active && hasAdvanceLimits()) ? 1 : 0;
}

void ReactorNet::evalRootFunctions(double t, span<const double> y, span<double> gout)
{
    // Default: no root detected
    double g = 1.0;

    if (m_limit_check_active) {
        // Ensure limits vector is current
        if (m_advancelimits.size() != m_nv) {
            m_advancelimits.assign(m_nv, -1.0);
        }
        getAdvanceLimits(m_advancelimits);

        double max_ratio = 0.0;
        for (size_t i = 0; i < m_nv; i++) {
            double lim = m_advancelimits[i];
            if (lim > 0.0) {
                double delta = std::abs(y[i] - m_ybase[i]);
                double ratio = delta / lim;
                if (ratio > max_ratio) {
                    max_ratio = ratio;
                }
            }
        }
        g = 1.0 - max_ratio; // root at g = 0 when any component reaches its limit
    }

    gout[0] = g;
}

void ReactorNet::updateNames(ReactorBase& r)
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

void ReactorNet::eval(double t, span<const double> y, span<double> ydot,
                      span<const double> p)
{
    m_time = t;
    updateState(y);
    m_LHS.assign(m_nv, 1);
    m_RHS.assign(m_nv, 0);
    for (auto& R : m_reactors) {
        size_t offset = R->offset();
        R->applySensitivity(p);
        R->eval(t, span<double>(m_LHS.data() + offset, R->neq()),
                span<double>(m_RHS.data() + offset, R->neq()));
        for (size_t i = offset; i < offset + R->neq(); i++) {
            ydot[i] = m_RHS[i] / m_LHS[i];
        }
        R->resetSensitivity(p);
    }
    checkFinite("ydot", ydot);
}

void ReactorNet::evalSteady(span<const double> y, span<double> residual)
{
    updateState(y);
    m_LHS.assign(m_nv, 1);
    m_RHS.assign(m_nv, 0);
    for (auto& R : m_reactors) {
        size_t offset = R->offset();
        R->evalSteady(m_time, span<double>(m_LHS.data() + offset, R->neq()),
                      span<double>(m_RHS.data() + offset, R->neq()));
        for (size_t i = offset; i < offset + R->neq(); i++) {
            residual[i] = m_RHS[i] / m_LHS[i];
        }
    }
    checkFinite("residual", residual);
}

void ReactorNet::evalDae(double t, span<const double> y, span<const double> ydot,
                         span<const double> p, span<double> residual)
{
    m_time = t;
    updateState(y);
    for (auto& R : m_reactors) {
        size_t offset = R->offset();
        R->applySensitivity(p);
        R->evalDae(t, y.subspan(offset, R->neq()),
                   ydot.subspan(offset, R->neq()),
                   residual.subspan(offset, R->neq()));
        R->resetSensitivity(p);
    }
    checkFinite("ydot", ydot);
}

void ReactorNet::getConstraints(span<double> constraints)
{
    for (auto& R : m_reactors) {
        R->getConstraints(constraints.subspan(R->offset(), R->neq()));
    }
}

double ReactorNet::sensitivity(size_t k, size_t p)
{
    initialize();
    if (p >= m_sens_params.size()) {
        throw IndexError("ReactorNet::sensitivity",
                         "m_sens_params", p, m_sens_params.size());
    }
    double denom = m_integ->solution(k);
    if (denom == 0.0) {
        denom = SmallNumber;
    }
    return m_integ->sensitivity(k, p) / denom;
}

void ReactorNet::evalJacobian(double t, span<double> y, span<double> ydot,
                              span<const double> p, Array2D* j)
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
        eval(t, y, m_ydot, p);

        // compute nth column of Jacobian
        for (size_t m = 0; m < m_nv; m++) {
            j->value(m,n) = (m_ydot[m] - ydot[m])/dy;
        }
        y[n] = ysave;
    }
}

void ReactorNet::updateState(span<const double> y)
{
    checkFinite("y", y);
    for (auto& R : m_reactors) {
        R->updateState(y.subspan(R->offset(), R->neq()));
    }
}

void ReactorNet::getDerivative(int k, span<double> dky)
{
    initialize();
    double* cvode_dky = m_integ->derivative(m_time, k);
    for (size_t j = 0; j < m_nv; j++) {
        dky[j] = cvode_dky[j];
    }
}

void ReactorNet::setAdvanceLimits(span<const double> limits)
{
    initialize();
    for (auto& R : m_bulkReactors) {
        R->setAdvanceLimits(limits.subspan(R->offset(), R->neq()));
    }
}

bool ReactorNet::hasAdvanceLimits() const
{
    bool has_limit = false;
    for (size_t n = 0; n < m_bulkReactors.size(); n++) {
        has_limit |= m_bulkReactors[n]->hasAdvanceLimits();
    }
    return has_limit;
}

bool ReactorNet::getAdvanceLimits(span<double> limits) const
{
    bool has_limit = false;
    for (auto& R : m_bulkReactors) {
        has_limit |= R->getAdvanceLimits(limits.subspan(R->offset(), R->neq()));
    }
    return has_limit;
}

void ReactorNet::getState(span<double> y)
{
    for (auto& R : m_reactors) {
        R->getState(y.subspan(R->offset(), R->neq()));
    }
}

void ReactorNet::getStateDae(span<double> y, span<double> ydot)
{
    // Iterate in reverse order so that surfaces will be handled first and up-to-date
    // values of the surface production rates of bulk species will be available when
    // bulk reactors are processed.
    // TODO: Replace with view::reverse once Cantera requires C++20
    for (size_t n = m_reactors.size(); n != 0 ; n--) {
        auto& R = m_reactors[n-1];
        R->getStateDae(y.subspan(R->offset(), R->neq()),
                       ydot.subspan(R->offset(), R->neq()));
    }
}

size_t ReactorNet::globalComponentIndex(const string& component, size_t reactor)
{
    initialize();
    return m_reactors[reactor]->offset()
           + m_reactors[reactor]->componentIndex(component);
}

string ReactorNet::componentName(size_t i) const
{
    size_t iTot = 0;
    size_t i0 = i;
    for (auto r : m_reactors) {
        if (i < r->neq()) {
            return r->name() + ": " + r->componentName(i);
        } else {
            i -= r->neq();
        }
        iTot += r->neq();
    }
    throw IndexError("ReactorNet::componentName", "component", i0, iTot);
}

double ReactorNet::upperBound(size_t i) const
{
    size_t iTot = 0;
    size_t i0 = i;
    for (auto r : m_reactors) {
        if (i < r->neq()) {
            return r->upperBound(i);
        } else {
            i -= r->neq();
        }
        iTot += r->neq();
    }
    throw IndexError("ReactorNet::upperBound", "upperBound", i0, iTot);
}

double ReactorNet::lowerBound(size_t i) const
{
    size_t iTot = 0;
    size_t i0 = i;
    for (auto r : m_reactors) {
        if (i < r->neq()) {
            return r->lowerBound(i);
        } else {
            i -= r->neq();
        }
        iTot += r->neq();
    }
    throw IndexError("ReactorNet::lowerBound", "lowerBound", i0, iTot);
}

void ReactorNet::resetBadValues(span<double> y) {
    for (auto& R : m_reactors) {
        R->resetBadValues(y.subspan(R->offset(), R->neq()));
    }
}

size_t ReactorNet::registerSensitivityParameter(
    const string& name, double value, double scale)
{
    if (m_integratorInitialized) {
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
    for (auto& R : m_bulkReactors) {
        R->setDerivativeSettings(settings);
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

void ReactorNet::preconditionerSolve(span<const double> rhs, span<double> output)
{
    if (!m_integ) {
        throw CanteraError("ReactorNet::preconditionerSolve",
                           "Must only be called after ReactorNet is initialized.");
    }
    m_integ->preconditionerSolve(m_nv, const_cast<double*>(rhs.data()), output.data());
}

void ReactorNet::preconditionerSetup(double t, span<const double> y, double gamma)
{
    // ensure state is up to date.
    updateState(y);
    // get the preconditioner
    auto precon = std::dynamic_pointer_cast<EigenSparseJacobian>(
        m_integ->preconditioner());
    // Reset preconditioner
    precon->reset();
    // Set gamma value for M =I - gamma*J
    precon->setGamma(gamma);
    // Make a copy of state to adjust it for preconditioner
    vector<double> yCopy(m_nv);
    // Get state of reactor
    getState(yCopy);
    // transform state based on preconditioner rules
    precon->stateAdjustment(yCopy);
    // update network with adjusted state
    updateState(yCopy);
    // Get jacobians and give elements to preconditioners
    vector<Eigen::Triplet<double>> trips;
    for (auto& R : m_reactors) {
        R->getJacobianElements(trips);
    }
    precon->setFromTriplets(trips);
    // post reactor setup operations
    precon->updatePreconditioner();
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
    for (auto reactor : m_bulkReactors) {
        if (!reactor->preconditionerSupported()) {
            throw CanteraError("ReactorNet::checkPreconditionerSupported",
                "Preconditioning is only supported for type *MoleReactor,\n"
                "Reactor type given: '{}'.",
                reactor->type());
        }
    }
}

SteadyReactorSolver::SteadyReactorSolver(ReactorNet* net, double* x0)
    : m_net(net)
{
    m_size = m_net->neq();
    m_jac = newSystemJacobian("eigen-sparse-direct");
    SteadyStateSystem::resize();
    m_initialState.assign(x0, x0 + m_size);
    setInitialGuess(x0);
    setJacobianPerturbation(m_jacobianRelPerturb, 1000 * m_net->atol(),
                            m_jacobianThreshold);
    m_mask.assign(m_size, 1);
    size_t start = 0;
    for (size_t i = 0; i < net->nReactors(); i++) {
        auto& R = net->reactor(i);
        R.updateState(span<const double>(x0 + start, R.neq()));
        auto algebraic = R.initializeSteady();
        for (auto& m : algebraic) {
            m_mask[start + m] = 0;
        }
        start += R.neq();
    }
}

void SteadyReactorSolver::eval(double* x, double* r, double rdt, int count)
{
    if (rdt < 0.0) {
        rdt = m_rdt;
    }
    vector<double> xv(x, x + size());
    m_net->evalSteady(span<const double>(x, size()),
                      span<double>(r, size()));
    for (size_t i = 0; i < size(); i++) {
        r[i] -= (x[i] - m_initialState[i]) * rdt * m_mask[i];
    }
}

void SteadyReactorSolver::initTimeInteg(double dt, double* x)
{
    SteadyStateSystem::initTimeInteg(dt, x);
    m_initialState.assign(x, x + size());
}

void SteadyReactorSolver::evalJacobian(double* x0)
{
    m_jac->reset();
    clock_t t0 = clock();
    m_work1.resize(size());
    m_work2.resize(size());
    eval(x0, m_work1.data(), 0.0, 0);
    for (size_t j = 0; j < size(); j++) {
        // perturb x(n); preserve sign(x(n))
        double xsave = x0[j];
        double dx = fabs(xsave) * m_jacobianRelPerturb + m_jacobianAbsPerturb;
        if (xsave < 0) {
            dx = -dx;
        }
        x0[j] = xsave + dx;
        double rdx = 1.0 / (x0[j] - xsave);

        // calculate perturbed residual
        eval(x0, m_work2.data(), 0.0, 0);

        // compute nth column of Jacobian
        for (size_t i = 0; i < size(); i++) {
            double delta = m_work2[i] - m_work1[i];
            if (std::abs(delta) > m_jacobianThreshold || i == j) {
                m_jac->setValue(i, j, delta * rdx);
            }
        }
        x0[j] = xsave;
    }
    // Restore system to unperturbed state
    m_net->updateState(span<const double>(x0, m_net->neq()));

    m_jac->updateElapsed(double(clock() - t0) / CLOCKS_PER_SEC);
    m_jac->incrementEvals();
    m_jac->setAge(0);
}

double SteadyReactorSolver::weightedNorm(const double* step) const
{
    double sum = 0.0;
    const double* x = m_state->data();
    for (size_t i = 0; i < size(); i++) {
        double ewt = m_net->rtol()*fabs(x[i]) + m_net->atol();
        double f = step[i] / ewt;
        sum += f*f;
    }
    return sqrt(sum / size());
}

string SteadyReactorSolver::componentName(size_t i) const
{
    return m_net->componentName(i);
}

double SteadyReactorSolver::upperBound(size_t i) const
{
    return m_net->upperBound(i);
}

double SteadyReactorSolver::lowerBound(size_t i) const
{
    return m_net->lowerBound(i);
}

void SteadyReactorSolver::resetBadValues(double* x)
{
    m_net->resetBadValues(span<double>(x, m_net->neq()));
}

void SteadyReactorSolver::writeDebugInfo(const string& header_suffix,
    const string& message, int loglevel, int attempt_counter)
{
    if (loglevel >= 6 && !m_state->empty()) {
        const auto& state = *m_state;
        writelog("Current state ({}):\n[", header_suffix);
        for (size_t i = 0; i < state.size() - 1; i++) {
            writelog("{}, ", state[i]);
        }
        writelog("{}]\n", state.back());
    }
    if (loglevel >= 7 && !m_xnew.empty()) {
        writelog("Current residual ({}):\n[", header_suffix);
        for (size_t i = 0; i < m_xnew.size() - 1; i++) {
            writelog("{}, ", m_xnew[i]);
        }
        writelog("{}]\n", m_xnew.back());
    }
}

shared_ptr<ReactorNet> newReactorNet(vector<shared_ptr<ReactorBase>>& reactors)
{
    return make_shared<ReactorNet>(reactors);
}

}
