//! @file ConstPressureReactor.cpp A constant pressure zero-dimensional reactor

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/ConstPressureReactor.h"
#include "cantera/zeroD/FlowDevice.h"
#include "cantera/zeroD/Wall.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/base/utilities.h"

using namespace std;

namespace Cantera
{

void ConstPressureReactor::getState(double* y)
{
    if (m_thermo == 0) {
        throw CanteraError("ConstPressureReactor::getState",
                           "Error: reactor is empty.");
    }
    m_thermo->restoreState(m_state);

    // set the first component to the total mass
    y[0] = m_thermo->density() * m_vol;

    // set the second component to the total enthalpy
    y[1] = m_thermo->enthalpy_mass() * m_thermo->density() * m_vol;

    // set components y+2 ... y+K+1 to the mass fractions Y_k of each species
    m_thermo->getMassFractions(y+2);

    // set the remaining components to the surface species
    // coverages on the walls
    getSurfaceInitialConditions(y + m_nsp + 2);
}

void ConstPressureReactor::initialize(doublereal t0)
{
    Reactor::initialize(t0);
    m_nv -= 1; // Constant pressure reactor has one fewer state variable
}

void ConstPressureReactor::updateState(doublereal* y)
{
    // The components of y are [0] the total mass, [1] the total enthalpy,
    // [2...K+2) are the mass fractions of each species, and [K+2...] are the
    // coverages of surface species on each wall.
    m_mass = y[0];
    m_thermo->setMassFractions_NoNorm(y+2);
    if (m_energy) {
        m_thermo->setState_HP(y[1]/m_mass, m_pressure);
    } else {
        m_thermo->setPressure(m_pressure);
    }
    m_vol = m_mass / m_thermo->density();
    updateConnected(false);
    updateSurfaceState(y + m_nsp + 2);
}

void ConstPressureReactor::eval(double time, double* LHS, double* RHS)
{
    double& dmdt = RHS[0];
    double* mdYdt = RHS + 2; // mass * dY/dt

    dmdt = 0.0;

    evalWalls(time);

    m_thermo->restoreState(m_state);
    const vector_fp& mw = m_thermo->molecularWeights();
    const doublereal* Y = m_thermo->massFractions();

    evalSurfaces(LHS + m_nsp + 2, RHS + m_nsp + 2, m_sdot.data());
    double mdot_surf = dot(m_sdot.begin(), m_sdot.end(), mw.begin());
    dmdt += mdot_surf;

    if (m_chem) {
        m_kin->getNetProductionRates(&m_wdot[0]); // "omega dot"
    }

    for (size_t k = 0; k < m_nsp; k++) {
        // production in gas phase and from surfaces
        mdYdt[k] = (m_wdot[k] * m_vol + m_sdot[k]) * mw[k];
        // dilution by net surface mass flux
        mdYdt[k] -= Y[k] * mdot_surf;
        //Assign left-hand side of dYdt ODE as total mass
        LHS[k+2] = m_mass;
    }

    // external heat transfer
    double dHdt = m_Qdot;

    // add terms for outlets
    for (auto outlet : m_outlet) {
        double mdot = outlet->massFlowRate();
        dmdt -= mdot;
        dHdt -= mdot * m_enthalpy;
    }

    // add terms for inlets
    for (auto inlet : m_inlet) {
        double mdot = inlet->massFlowRate();
        dmdt += mdot; // mass flow into system
        for (size_t n = 0; n < m_nsp; n++) {
            double mdot_spec = inlet->outletSpeciesMassFlowRate(n);
            // flow of species into system and dilution by other species
            mdYdt[n] += mdot_spec - mdot * Y[n];
        }
        dHdt += mdot * inlet->enthalpy_mass();
    }

    if (m_energy) {
        RHS[1] = dHdt;
    } else {
        RHS[1] = 0.0;
    }
}

size_t ConstPressureReactor::componentIndex(const string& nm) const
{
    size_t k = speciesIndex(nm);
    if (k != npos) {
        return k + 2;
    } else if (nm == "mass") {
        return 0;
    } else if (nm == "enthalpy") {
        return 1;
    } else {
        return npos;
    }
}

std::string ConstPressureReactor::componentName(size_t k) {
    if (k == 0) {
        return "mass";
    } else if (k == 1) {
        return "enthalpy";
    } else if (k >= 2 && k < neq()) {
        k -= 2;
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
    throw CanteraError("ConstPressureReactor::componentName",
                       "Index is out of bounds.");
}

}
