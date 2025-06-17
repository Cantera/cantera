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

namespace Cantera
{

void ConstPressureMoleReactor::getState(double* y)
{
    if (m_thermo == 0) {
        throw CanteraError("ConstPressureMoleReactor::getState",
                           "Error: reactor is empty.");
    }
    m_thermo->restoreState(m_state);
    // set mass to be used in getMoles function
    m_mass = m_thermo->density() * m_vol;
    // set the first array element to enthalpy
    y[0] = m_thermo->enthalpy_mass() * m_thermo->density() * m_vol;
    // get moles of species in remaining state
    getMoles(y + m_sidx);
    // set the remaining components to the surface species moles on the walls
    getSurfaceInitialConditions(y+m_nsp+m_sidx);
}

void ConstPressureMoleReactor::initialize(double t0)
{
    MoleReactor::initialize(t0);
    m_nv -= 1; // const pressure system loses 1 more variable from MoleReactor
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
    updateSurfaceState(y + m_nsp + m_sidx);
}

void ConstPressureMoleReactor::eval(double time, double* LHS, double* RHS)
{
    double* dndt = RHS + m_sidx; // kmol per s

    evalWalls(time);

    m_thermo->restoreState(m_state);

    const vector<double>& imw = m_thermo->inverseMolecularWeights();

    if (m_chem) {
        m_kin->getNetProductionRates(&m_wdot[0]); // "omega dot"
    }

    // evaluate reactor surfaces
    evalSurfaces(LHS + m_nsp + m_sidx, RHS + m_nsp + m_sidx, m_sdot.data());

    // external heat transfer
    double dHdt = m_Qdot;

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
    size_t k = speciesIndex(nm);
    if (k != npos) {
        return k + m_sidx;
    } else if (nm == "enthalpy") {
        return 0;
    } else {
        return npos;
    }
}

string ConstPressureMoleReactor::componentName(size_t k) {
    if (k == 0) {
        return "enthalpy";
    } else if (k >= m_sidx && k < neq()) {
        k -= m_sidx;
        if (k < m_thermo->nSpecies()) {
            return m_thermo->speciesName(k);
        } else {
            k -= m_thermo->nSpecies();
        }
        for (auto& S : m_surfaces) {
            ThermoPhase* th = S->thermo();
            if (k < th->nSpecies()) {
                return th->speciesName(k);
            } else {
                k -= th->nSpecies();
            }
        }
    }
    throw CanteraError("ConstPressureMoleReactor::componentName",
                       "Index is out of bounds.");
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
