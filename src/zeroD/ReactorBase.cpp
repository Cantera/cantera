//! @file ReactorBase.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/ReactorBase.h"
#include "cantera/zeroD/FlowDevice.h"
#include "cantera/zeroD/ReactorNet.h"
#include "cantera/zeroD/ReactorSurface.h"
#include "cantera/base/Solution.h"
#include "cantera/thermo/ThermoPhase.h"

namespace Cantera
{

ReactorBase::ReactorBase(const string& name) : m_name(name)
{
}

// TODO: After Cantera 3.2, this method can be replaced by delegating to
// ReactorBase::ReactorBase(sol, true, name)
ReactorBase::ReactorBase(shared_ptr<Solution> sol, const string& name)
    : ReactorBase(name)
{
    warn_deprecated("ReactorBase::ReactorBase", "`clone` argument not specified; "
        "Default behavior will change from `clone=False` to `clone=True` after "
        "Cantera 3.2.");
    if (!sol || !(sol->thermo())) {
        throw CanteraError("ReactorBase::ReactorBase",
                           "Missing or incomplete Solution object.");
    }
    m_solution = sol->clone({}, true, false);
    m_solution->thermo()->addSpeciesLock();
    m_thermo = m_solution->thermo().get();
    m_nsp = m_thermo->nSpecies();
    m_thermo->saveState(m_state);
    m_enthalpy = m_thermo->enthalpy_mass(); // Needed for flow and wall interactions
    m_pressure = m_thermo->pressure(); // Needed for flow and wall interactions
    try {
        m_intEnergy = m_thermo->intEnergy_mass();
    } catch (NotImplementedError&) {
        // some ThermoPhase objects do not implement intEnergy_mass()
    }
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
    m_thermo->saveState(m_state);
    m_enthalpy = m_thermo->enthalpy_mass(); // Needed for flow and wall interactions
    m_pressure = m_thermo->pressure(); // Needed for flow and wall interactions
    try {
        m_intEnergy = m_thermo->intEnergy_mass();
    } catch (NotImplementedError&) {
        // some ThermoPhase objects do not implement intEnergy_mass()
    }
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

void ReactorBase::setSolution(shared_ptr<Solution> sol)
{
    warn_deprecated("ReactorBase::setSolution",
        "After Cantera 3.2, a change of reactor contents after instantiation "
        "will be disabled.");
    if (!sol || !(sol->thermo())) {
        throw CanteraError("ReactorBase::setSolution",
            "Missing or incomplete Solution object.");
    }
    if (m_solution) {
        m_solution->thermo()->removeSpeciesLock();
    }
    m_solution = sol;
    m_solution->thermo()->addSpeciesLock();
    setThermo(*m_solution->thermo());
    try {
        setKinetics(*m_solution->kinetics());
    } catch (NotImplementedError&) {
        // kinetics not used (example: Reservoir)
    }
}

void ReactorBase::setThermo(ThermoPhase& thermo)
{
    warn_deprecated("ReactorBase::setThermo",
        "After Cantera 3.2, a change of reactor contents after instantiation "
        "will be disabled.");
    m_thermo = &thermo;
    m_nsp = m_thermo->nSpecies();
    m_thermo->saveState(m_state);
    m_enthalpy = m_thermo->enthalpy_mass(); // Needed for flow and wall interactions
    m_pressure = m_thermo->pressure(); // Needed for flow and wall interactions
    try {
        m_intEnergy = m_thermo->intEnergy_mass();
    } catch (NotImplementedError&) {
        // some ThermoPhase objects do not implement intEnergy_mass()
    }
}

void ReactorBase::addInlet(FlowDevice& inlet)
{
    m_inlet.push_back(&inlet);
}

void ReactorBase::addOutlet(FlowDevice& outlet)
{
    m_outlet.push_back(&outlet);
}

void ReactorBase::addWall(WallBase& w, int lr)
{
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
    if (find(m_surfaces.begin(), m_surfaces.end(), surf) == m_surfaces.end()) {
        m_surfaces.push_back(surf);
        surf->setReactor(this);
    }
}

void ReactorBase::addSurface(shared_ptr<ReactorBase> surf)
{
    auto r = std::dynamic_pointer_cast<ReactorSurface>(surf);
    if (!r) {
        throw CanteraError("ReactorBase::addSurface",
                           "Invalid reactor type '{}'.", surf->type());
    }
    addSurface(r.get());
}

ReactorSurface* ReactorBase::surface(size_t n)
{
    return m_surfaces[n];
}

void ReactorBase::restoreState() {
    if (!m_thermo) {
        throw CanteraError("ReactorBase::restoreState", "No phase defined.");
    }
    m_thermo->restoreState(m_state);
}

void ReactorBase::syncState()
{
    m_thermo->saveState(m_state);
    m_enthalpy = m_thermo->enthalpy_mass();
    try {
        m_intEnergy = m_thermo->intEnergy_mass();
    } catch (NotImplementedError&) {
        m_intEnergy = NAN;
    }
    m_pressure = m_thermo->pressure();
    m_mass = m_thermo->density() * m_vol;
    if (m_net) {
        m_net->setNeedsReinit();
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

FlowDevice& ReactorBase::inlet(size_t n)
{
    return *m_inlet[n];
}
FlowDevice& ReactorBase::outlet(size_t n)
{
    return *m_outlet[n];
}

}
