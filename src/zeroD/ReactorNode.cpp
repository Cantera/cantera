//! @file ReactorNode.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/ReactorNode.h"
#include "cantera/zeroD/FlowDevice.h"
#include "cantera/zeroD/ReactorNet.h"
#include "cantera/zeroD/ReactorSurface.h"
#include "cantera/base/Solution.h"
#include "cantera/thermo/ThermoPhase.h"

namespace Cantera
{

ReactorNode::ReactorNode(const string& name) : m_name(name)
{
}

ReactorNode::ReactorNode(shared_ptr<Solution> sol, const string& name)
    : ReactorNode(name)
{
    if (!sol || !(sol->thermo())) {
        throw CanteraError("ReactorNode::ReactorNode",
                           "Missing or incomplete Solution object.");
    }
    m_solution = sol;
    setThermo(*sol->thermo());
    try {
        setKinetics(*sol->kinetics());
    } catch (NotImplementedError&) {
        // kinetics not used (example: Reservoir)
    }
    m_solution->thermo()->addSpeciesLock();
}

ReactorNode::~ReactorNode()
{
    if (m_solution) {
        m_solution->thermo()->removeSpeciesLock();
    }
}

bool ReactorNode::setDefaultName(map<string, int>& counts)
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

void ReactorNode::setSolution(shared_ptr<Solution> sol)
{
    warn_deprecated("ReactorNode::setSolution",
        "After Cantera 3.2, a change of reactor contents after instantiation "
        "will be disabled.");
    if (!sol || !(sol->thermo())) {
        throw CanteraError("ReactorNode::setSolution",
            "Missing or incomplete Solution object.");
    }
    if (m_solution) {
        m_solution->thermo()->removeSpeciesLock();
    }
    m_solution = sol;
    setThermo(*sol->thermo());
    try {
        setKinetics(*sol->kinetics());
    } catch (NotImplementedError&) {
        // kinetics not used (example: Reservoir)
    }
    m_solution->thermo()->addSpeciesLock();
}

void ReactorNode::setThermo(ThermoPhase& thermo)
{
    m_thermo = &thermo;
    m_nsp = m_thermo->nSpecies();
    m_thermo->saveState(m_state);
    m_enthalpy = m_thermo->enthalpy_mass();
    m_intEnergy = m_thermo->intEnergy_mass();
    m_pressure = m_thermo->pressure();
}

void ReactorNode::syncState()
{
    m_thermo->saveState(m_state);
    m_enthalpy = m_thermo->enthalpy_mass();
    m_intEnergy = m_thermo->intEnergy_mass();
    m_pressure = m_thermo->pressure();
    if (m_net) {
        m_net->setNeedsReinit();
    }
}

void ReactorNode::addInlet(FlowDevice& inlet)
{
    m_inlet.push_back(&inlet);
}

void ReactorNode::addOutlet(FlowDevice& outlet)
{
    m_outlet.push_back(&outlet);
}

void ReactorNode::addWall(WallBase& w, int lr)
{
    m_wall.push_back(&w);
    if (lr == 0) {
        m_lr.push_back(0);
    } else {
        m_lr.push_back(1);
    }
}

WallBase& ReactorNode::wall(size_t n)
{
    return *m_wall[n];
}

void ReactorNode::addSurface(ReactorSurface* surf)
{
    if (find(m_surfaces.begin(), m_surfaces.end(), surf) == m_surfaces.end()) {
        m_surfaces.push_back(surf);
        surf->setReactor(this);
    }
}

ReactorSurface* ReactorNode::surface(size_t n)
{
    return m_surfaces[n];
}

void ReactorNode::restoreState() {
    if (!m_thermo) {
        throw CanteraError("ReactorNode::restoreState", "No phase defined.");
    }
    m_thermo->restoreState(m_state);
}

ReactorNet& ReactorNode::network()
{
    if (m_net) {
        return *m_net;
    } else {
        throw CanteraError("ReactorNode::network",
                           "Reactor is not part of a ReactorNet");
    }
}

void ReactorNode::setNetwork(ReactorNet* net)
{
    m_net = net;
}

double ReactorNode::residenceTime()
{
    double mout = 0.0;
    for (size_t i = 0; i < m_outlet.size(); i++) {
        mout += m_outlet[i]->massFlowRate();
    }
    return mass()/mout;
}

FlowDevice& ReactorNode::inlet(size_t n)
{
    return *m_inlet[n];
}
FlowDevice& ReactorNode::outlet(size_t n)
{
    return *m_outlet[n];
}

}
