//! @file ReactorNet.cpp
#include "cantera/zeroD/ReactorNet.h"
#include "cantera/zeroD/FlowDevice.h"
#include "cantera/zeroD/Wall.h"

#include <cstdio>

using namespace std;

namespace Cantera
{

ReactorNet::ReactorNet() :
    m_integ(0), m_time(0.0), m_init(false), m_integrator_init(false),
    m_nv(0), m_rtol(1.0e-9), m_rtolsens(1.0e-4),
    m_atols(1.0e-15), m_atolsens(1.0e-4),
    m_maxstep(0.0), m_maxErrTestFails(0),
    m_verbose(false), m_ntotpar(0)
{
    m_integ = newIntegrator("CVODE");

    // use backward differencing, with a full Jacobian computed
    // numerically, and use a Newton linear iterator
    m_integ->setMethod(BDF_Method);
    m_integ->setProblemType(DENSE + NOJAC);
    m_integ->setIterator(Newton_Iter);
}

ReactorNet::~ReactorNet()
{
    delete m_integ;
}

void ReactorNet::initialize()
{
    size_t n, nv;
    m_nv = 0;
    debuglog("Initializing reactor network.\n", m_verbose);
    if (m_reactors.empty()) {
        throw CanteraError("ReactorNet::initialize",
                           "no reactors in network!");
    }
    size_t sensParamNumber = 0;
    m_start.assign(1, 0);
    for (n = 0; n < m_reactors.size(); n++) {
        Reactor& r = *m_reactors[n];
        r.initialize(m_time);
        nv = r.neq();
        m_nparams.push_back(r.nSensParams());
        for (const auto& sens_obj : r.getSensitivityOrder()) {
            for (const auto& order : m_sensOrder[sens_obj]) {
                m_sensIndex.resize(std::max(order.second + 1, m_sensIndex.size()));
                m_sensIndex[order.second] = sensParamNumber++;
            }
        }
        m_nv += nv;
        m_start.push_back(m_nv);

        if (m_verbose) {
            writelog("Reactor {:d}: {:d} variables.\n", n, nv);
            writelog("              {:d} sensitivity params.\n", r.nSensParams());
        }
        if (r.type() == FlowReactorType && m_reactors.size() > 1) {
            throw CanteraError("ReactorNet::initialize",
                               "FlowReactors must be used alone.");
        }
    }

    m_ydot.resize(m_nv,0.0);
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
    m_integrator_init = true;
    m_init = true;
}

void ReactorNet::reinitialize()
{
    if (m_init) {
        debuglog("Re-initializing reactor network.\n", m_verbose);
        m_integ->reinitialize(m_time, *this);
        m_integrator_init = true;
    } else {
        initialize();
    }
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

double ReactorNet::step(doublereal time)
{
    if (time != -999) {
        warn_deprecated("ReactorNet::step(t)", "The argument to this function"
            " is deprecated and will be removed after Cantera 2.3.");
    }
    if (!m_init) {
        initialize();
    } else if (!m_integrator_init) {
        reinitialize();
    }
    m_time = m_integ->step(m_time + 1.0);
    updateState(m_integ->solution());
    return m_time;
}

void ReactorNet::addReactor(Reactor& r)
{
    r.setNetwork(this);
    m_reactors.push_back(&r);
}

void ReactorNet::eval(doublereal t, doublereal* y,
                      doublereal* ydot, doublereal* p)
{
    size_t n;
    size_t pstart = 0;
    updateState(y);
    for (n = 0; n < m_reactors.size(); n++) {
        m_reactors[n]->evalEqs(t, y + m_start[n],
                               ydot + m_start[n], p + pstart);
        pstart += m_nparams[n];
    }
    checkFinite("ydot", ydot, m_nv);
}

void ReactorNet::evalJacobian(doublereal t, doublereal* y,
                              doublereal* ydot, doublereal* p, Array2D* j)
{
    doublereal ysave, dy;
    Array2D& jac = *j;

    //evaluate the unperturbed ydot
    eval(t, y, ydot, p);
    for (size_t n = 0; n < m_nv; n++) {
        // perturb x(n)
        ysave = y[n];
        dy = m_atol[n] + fabs(ysave)*m_rtol;
        y[n] = ysave + dy;
        dy = y[n] - ysave;

        // calculate perturbed residual
        eval(t, y, m_ydot.data(), p);

        // compute nth column of Jacobian
        for (size_t m = 0; m < m_nv; m++) {
            jac(m,n) = (m_ydot[m] - ydot[m])/dy;
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

void ReactorNet::getInitialConditions(doublereal t0,
                                      size_t leny, doublereal* y)
{
    for (size_t n = 0; n < m_reactors.size(); n++) {
        m_reactors[n]->getInitialConditions(t0, m_start[n+1]-m_start[n],
                                            y + m_start[n]);
    }
}

size_t ReactorNet::globalComponentIndex(const string& component, size_t reactor)
{
    if (!m_init) {
        initialize();
    }
    return m_start[reactor] + m_reactors[reactor]->componentIndex(component);
}

void ReactorNet::registerSensitivityReaction(void* reactor,
        size_t reactionIndex, const std::string& name, int leftright)
{
    if (m_integrator_init) {
        throw CanteraError("ReactorNet::registerSensitivityReaction",
                           "Sensitivity reactions cannot be added after the"
                           "integrator has been initialized.");
    }
    std::pair<void*, int> R = {reactor, leftright};
    if (m_sensOrder.count(R) &&
            m_sensOrder[R].count(reactionIndex)) {
        throw CanteraError("ReactorNet::registerSensitivityReaction",
                           "Attempted to register duplicate sensitivity reaction");
    }
    m_paramNames.push_back(name);
    m_sensOrder[R][reactionIndex] = m_ntotpar;
    m_ntotpar++;
}

}
