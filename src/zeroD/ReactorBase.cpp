//! @file ReactorBase.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/ReactorBase.h"
#include "cantera/zeroD/FlowDevice.h"
#include "cantera/zeroD/ReactorNet.h"
#include "cantera/zeroD/ReactorSurface.h"
#include "cantera/zeroD/Wall.h"
#include "cantera/base/Solution.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/kinetics/Kinetics.h"

namespace Cantera
{

ReactorBase::ReactorBase(const string& name) : m_name(name)
{
}

ReactorBase::ReactorBase(shared_ptr<Solution> sol, const string& name)
    : ReactorBase(sol, true, name)
{
}

ReactorBase::ReactorBase(shared_ptr<Solution> sol, bool clone, const string& name)
    : ReactorBase(name)
{
    if (!sol || !(sol->thermo())) {
        throw CanteraError("ReactorBase::ReactorBase",
                           "Missing or incomplete Solution object.");
    }
    if (clone) {
        m_solution = sol->clone({}, true, false);
    } else {
        m_solution = sol;
    }
    m_solution->thermo()->addSpeciesLock();
    m_thermo = m_solution->thermo().get();
    m_nsp = m_thermo->nSpecies();
    m_enthalpy = m_thermo->enthalpy_mass(); // Needed for flow and wall interactions
    m_pressure = m_thermo->pressure(); // Needed for flow and wall interactions
}

ReactorBase::~ReactorBase()
{
    if (m_solution) {
        m_solution->thermo()->removeSpeciesLock();
    }
}

bool ReactorBase::setDefaultName(map<string, int>& counts)
{
    if (m_defaultNameSet) {
        return false;
    }
    m_defaultNameSet = true;
    if (m_name == "(none)" || m_name == "") {
        m_name = fmt::format("{}_{}", type(), counts[type()]);
    }
    counts[type()]++;
    return true;
}

void ReactorBase::addInlet(FlowDevice& inlet)
{
    if (m_net) {
        throw CanteraError("ReactorBase::addInlet",
                           "Cannot add an inlet to reactor '{}' after it has been"
                           " added to a ReactorNet.", m_name);
    }
    m_inlet.push_back(&inlet);
}

void ReactorBase::addOutlet(FlowDevice& outlet)
{
    if (m_net) {
        throw CanteraError("ReactorBase::addOutlet",
                           "Cannot add an outlet to reactor '{}' after it has been"
                           " added to a ReactorNet.", m_name);
    }
    m_outlet.push_back(&outlet);
}

void ReactorBase::addWall(WallBase& w, int lr)
{
    if (m_net) {
        throw CanteraError("ReactorBase::addWall",
                           "Cannot add a wall to reactor '{}' after it has been"
                           " added to a ReactorNet.", m_name);
    }
    m_wall.push_back(&w);
    if (lr == 0) {
        m_lr.push_back(0);
    } else {
        m_lr.push_back(1);
    }
}

WallBase& ReactorBase::wall(size_t n)
{
    return *m_wall[n];
}

void ReactorBase::addSurface(ReactorSurface* surf)
{
    if (m_net) {
        throw CanteraError("ReactorBase::addSurface",
                           "Cannot add a surface to reactor '{}' after it has been"
                           " added to a ReactorNet.", m_name);
    }
    if (find(m_surfaces.begin(), m_surfaces.end(), surf) == m_surfaces.end()) {
        m_surfaces.push_back(surf);
    }
}

ReactorSurface* ReactorBase::surface(size_t n)
{
    return m_surfaces[n];
}

Eigen::SparseMatrix<double> ReactorBase::jacobian()
{
    vector<Eigen::Triplet<double>> trips;
    getJacobianElements(trips);
    for (auto& trip : trips) {
        trip = Eigen::Triplet<double>(trip.row() - m_offset, trip.col() - m_offset, trip.value());
    }
    Eigen::SparseMatrix<double> J(m_nv, m_nv);
    J.setFromTriplets(trips.begin(), trips.end());
    return J;
}

void ReactorBase::syncState()
{
    warn_deprecated("ReactorBase::syncState",
        "To be removed after Cantera 4.0. Use ReactorNet::reinitialize to indicate "
        "a change in state that requires integrator reinitialization.");
    if (m_net) {
        m_net->setNeedsReinit();
    }
}

void ReactorBase::updateConnected(bool updatePressure) {
    // save parameters needed by other connected reactors
    m_enthalpy = m_thermo->enthalpy_mass();
    if (updatePressure) {
        m_pressure = m_thermo->pressure();
    }

    // Update the mass flow rate of connected flow devices
    double time = 0.0;
    if (m_net != nullptr) {
        time = (timeIsIndependent()) ? m_net->time() : m_net->distance();
    }
    for (size_t i = 0; i < m_outlet.size(); i++) {
        m_outlet[i]->setSimTime(time);
        m_outlet[i]->updateMassFlowRate(time);
    }
    for (size_t i = 0; i < m_inlet.size(); i++) {
        m_inlet[i]->setSimTime(time);
        m_inlet[i]->updateMassFlowRate(time);
    }
    for (size_t i = 0; i < m_wall.size(); i++) {
        m_wall[i]->setSimTime(time);
    }
}

ReactorNet& ReactorBase::network()
{
    if (m_net) {
        return *m_net;
    } else {
        throw CanteraError("ReactorBase::network",
                           "Reactor is not part of a ReactorNet");
    }
}

void ReactorBase::setNetwork(ReactorNet* net)
{
    if (m_net) {
        throw CanteraError("ReactorBase::setNetwork",
                           "Reactor {} is already part of a ReactorNet", m_name);
    }
    m_net = net;
}

double ReactorBase::residenceTime()
{
    double mout = 0.0;
    for (size_t i = 0; i < m_outlet.size(); i++) {
        mout += m_outlet[i]->massFlowRate();
    }
    return mass()/mout;
}

double ReactorBase::density() const
{
    return m_thermo->density();
}

double ReactorBase::temperature() const
{
    return m_thermo->temperature();
}

const double* ReactorBase::massFractions() const
{
    return m_thermo->massFractions().data();
}

double ReactorBase::massFraction(size_t k) const
{
    return m_thermo->massFraction(k);
}

FlowDevice& ReactorBase::inlet(size_t n)
{
    return *m_inlet[n];
}
FlowDevice& ReactorBase::outlet(size_t n)
{
    return *m_outlet[n];
}

}
