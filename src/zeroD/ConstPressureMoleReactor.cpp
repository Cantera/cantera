//! @file ConstPressureMoleReactor.cpp A constant pressure
//! zero-dimensional reactor with moles as the state

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/Wall.h"
#include "cantera/zeroD/ReactorSurface.h"
#include "cantera/zeroD/FlowDevice.h"
#include "cantera/zeroD/ConstPressureMoleReactor.h"
#include "cantera/base/utilities.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/thermo/PlasmaPhase.h"

namespace Cantera
{

ConstPressureMoleReactor::ConstPressureMoleReactor(shared_ptr<Solution> sol,
                                                   const string& name)
    : ConstPressureMoleReactor(sol, true, name)
{
}

ConstPressureMoleReactor::ConstPressureMoleReactor(shared_ptr<Solution> sol, bool clone,
                                                   const string& name)
    : MoleReactor(sol, clone, name)
{
    m_nv = 1 + m_nsp; // enthalpy and moles of each species
}

void ConstPressureMoleReactor::getState(double* y)
{
    // set mass to be used in getMoles function
    m_mass = m_thermo->density() * m_vol;
    // set the first array element to enthalpy
    y[0] = m_thermo->enthalpy_mass() * m_thermo->density() * m_vol;
    // get moles of species in remaining state
    getMoles(y + m_sidx);
}

void ConstPressureMoleReactor::updateState(double* y)
{
    // the components of y are: [0] the enthalpy, [1...K+1) are the
    // moles of each species, and [K+1...] are the moles of surface
    // species on each wall.
    setMassFromMoles(y + m_sidx);
    m_thermo->setMolesNoTruncate(y + m_sidx);
    if (m_energy) {
        m_thermo->setState_HP(y[0] / m_mass, m_pressure);
    } else {
        m_thermo->setPressure(m_pressure);
    }
    m_vol = m_mass / m_thermo->density();
    updateConnected(false);
}

void ConstPressureMoleReactor::eval(double time, double* LHS, double* RHS)
{
    double* dndt = RHS + m_sidx; // kmol per s

    evalWalls(time);
    updateSurfaceProductionRates();

    const vector<double>& imw = m_thermo->inverseMolecularWeights();

    if (m_chem) {
        m_kin->getNetProductionRates(&m_wdot[0]); // "omega dot"
    }

    // external heat transfer
    double dHdt = m_Qdot;

    if (auto* plasma = dynamic_cast<PlasmaPhase*>(m_thermo)) {
        const double qJ = plasma->jouleHeatingPower(); // σE^2 [W/m^3]
        const double qElastic = plasma->elasticPowerLoss(); // elastic transfer [W/m^3]
        const double q_total = (qJ + qElastic) * m_vol; // total power [W]
        if (std::isfinite(q_total)) {
            dHdt += q_total;
        }
    }

    for (size_t n = 0; n < m_nsp; n++) {
        // production in gas phase and from surfaces
        dndt[n] = m_wdot[n] * m_vol + m_sdot[n];
    }

    // add terms for outlets
    for (auto outlet : m_outlet) {
        // determine enthalpy contribution
        dHdt -= outlet->massFlowRate() * m_enthalpy;
        // flow of species into system and dilution by other species
        for (size_t n = 0; n < m_nsp; n++) {
            dndt[n] -= outlet->outletSpeciesMassFlowRate(n) * imw[n];
        }
    }

    // add terms for inlets
    for (auto inlet : m_inlet) {
        // enthalpy contribution from inlets
        dHdt += inlet->enthalpy_mass() * inlet->massFlowRate();
        // flow of species into system and dilution by other species
        for (size_t n = 0; n < m_nsp; n++) {
            dndt[n] += inlet->outletSpeciesMassFlowRate(n) * imw[n];
        }
    }

    if (m_energy) {
        RHS[0] = dHdt;
    } else {
        RHS[0] = 0.0;
    }
}

size_t ConstPressureMoleReactor::componentIndex(const string& nm) const
{
    if (nm == "enthalpy") {
        return 0;
    }
    try {
        return m_thermo->speciesIndex(nm) + m_sidx;
    } catch (const CanteraError&) {
        throw CanteraError("ConstPressureMoleReactor::componentIndex",
            "Component '{}' not found", nm);
    }
}

string ConstPressureMoleReactor::componentName(size_t k) {
    if (k == 0) {
        return "enthalpy";
    } else if (k >= m_sidx && k < neq()) {
        return m_thermo->speciesName(k - m_sidx);
    } else {
        throw IndexError("ConstPressureMoleReactor::componentName",
                         "component", k, m_nv);
    }
}

double ConstPressureMoleReactor::upperBound(size_t k) const {
    // Component is either enthalpy or moles of a bulk or surface species
    return BigNumber;
}

double ConstPressureMoleReactor::lowerBound(size_t k) const {
    if (k == 0) {
        return -BigNumber; // enthalpy
    } else if (k >= 1 && k < m_nv) {
        return -Tiny; // moles of bulk or surface species
    } else {
        throw CanteraError("ConstPressureMoleReactor::lowerBound", "Index {} is out of bounds.", k);
    }
}

void ConstPressureMoleReactor::resetBadValues(double* y) {
    for (size_t k = m_sidx; k < m_nv; k++) {
        y[k] = std::max(y[k], 0.0);
    }
}

}
