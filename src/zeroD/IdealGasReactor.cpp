//! @file IdealGasReactor.cpp A zero-dimensional reactor

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/IdealGasReactor.h"
#include "cantera/zeroD/FlowDevice.h"
#include "cantera/zeroD/Wall.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/base/utilities.h"

namespace Cantera
{

void IdealGasReactor::getState(double* y)
{
    // set the first component to the total mass
    m_mass = m_thermo->density() * m_vol;
    y[0] = m_mass;

    // set the second component to the total volume
    y[1] = m_vol;

    // Set the third component to the temperature
    y[2] = m_thermo->temperature();

    // set components y+3 ... y+K+2 to the mass fractions of each species
    m_thermo->getMassFractions(y+3);
}

void IdealGasReactor::initialize(double t0)
{
    //! @todo: Add a method to ThermoPhase that indicates whether a given
    //! subclass is compatible with this reactor model
    if (m_thermo->type() != "ideal-gas") {
        throw CanteraError("IdealGasReactor::initialize",
                           "Incompatible phase type '{}' provided", m_thermo->type());
    }
    Reactor::initialize(t0);
    m_uk.resize(m_nsp, 0.0);
}

void IdealGasReactor::updateState(double* y)
{
    // The components of y are [0] the total mass, [1] the total volume,
    // [2] the temperature, [3...K+3] are the mass fractions of each species,
    // and [K+3...] are the coverages of surface species on each wall.
    m_mass = y[0];
    m_vol = y[1];
    m_thermo->setMassFractions_NoNorm(y+3);
    m_thermo->setState_TD(y[2], m_mass / m_vol);
    updateConnected(true);
}

void IdealGasReactor::eval(double time, double* LHS, double* RHS)
{
    double& dmdt = RHS[0]; // dm/dt (gas phase)
    double& mcvdTdt = RHS[2]; // m * c_v * dT/dt
    double* mdYdt = RHS + 3; // mass * dY/dt

    evalWalls(time);
    updateSurfaceProductionRates();
    m_thermo->getPartialMolarIntEnergies(&m_uk[0]);
    const vector<double>& mw = m_thermo->molecularWeights();
    const double* Y = m_thermo->massFractions();

    if (m_chem) {
        m_kin->getNetProductionRates(&m_wdot[0]); // "omega dot"
    }

    double mdot_surf = dot(m_sdot.begin(), m_sdot.end(), mw.begin());
    dmdt += mdot_surf;

    // compression work and external heat transfer
    mcvdTdt += - m_pressure * m_vdot + m_Qdot;

    for (size_t n = 0; n < m_nsp; n++) {
        // heat release from gas phase and surface reactions
        mcvdTdt -= m_wdot[n] * m_uk[n] * m_vol;
        mcvdTdt -= m_sdot[n] * m_uk[n];
        // production in gas phase and from surfaces
        mdYdt[n] = (m_wdot[n] * m_vol + m_sdot[n]) * mw[n];
        // dilution by net surface mass flux
        mdYdt[n] -= Y[n] * mdot_surf;
        //Assign left-hand side of dYdt ODE as total mass
        LHS[n+3] = m_mass;
    }

    // add terms for outlets
    for (auto outlet : m_outlet) {
        double mdot = outlet->massFlowRate();
        dmdt -= mdot; // mass flow out of system
        mcvdTdt -= mdot * m_pressure * m_vol / m_mass; // flow work
    }

    // add terms for inlets
    for (auto inlet : m_inlet) {
        double mdot = inlet->massFlowRate();
        dmdt += mdot; // mass flow into system
        mcvdTdt += inlet->enthalpy_mass() * mdot;
        for (size_t n = 0; n < m_nsp; n++) {
            double mdot_spec = inlet->outletSpeciesMassFlowRate(n);
            // flow of species into system and dilution by other species
            mdYdt[n] += mdot_spec - mdot * Y[n];

            // In combination with h_in*mdot_in, flow work plus thermal
            // energy carried with the species
            mcvdTdt -= m_uk[n] / mw[n] * mdot_spec;
        }
    }

    RHS[1] = m_vdot;
    if (m_energy) {
        LHS[2] = m_mass * m_thermo->cv_mass();
    } else {
        RHS[2] = 0;
    }
}

vector<size_t> IdealGasReactor::steadyConstraints() const
{
    if (nSurfs() != 0) {
        throw CanteraError("IdealGasReactor::steadyConstraints",
            "Steady state solver cannot currently be used with IdealGasReactor"
            " when reactor surfaces are present.\n"
            "See https://github.com/Cantera/enhancements/issues/234");
    }
    if (energyEnabled()) {
        return {1}; // volume
    } else {
        return {1, 2}; // volume and temperature
    }
}

size_t IdealGasReactor::componentIndex(const string& nm) const
{
    if (nm == "mass") {
        return 0;
    }
    if (nm == "volume") {
        return 1;
    }
    if (nm == "temperature") {
        return 2;
    }
    try {
        return m_thermo->speciesIndex(nm) + 3;
    } catch (const CanteraError&) {
        throw CanteraError("IdealGasReactor::componentIndex",
            "Component '{}' not found", nm);
    }
}

string IdealGasReactor::componentName(size_t k) {
    if (k == 2) {
        return "temperature";
    } else {
        return Reactor::componentName(k);
    }
}

double IdealGasReactor::upperBound(size_t k) const {
    if (k == 2) {
        //@todo: Revise pending resolution of https://github.com/Cantera/enhancements/issues/229
        return 1.5 * m_thermo->maxTemp();
    } else {
        return Reactor::upperBound(k);
    }
}

double IdealGasReactor::lowerBound(size_t k) const {
    if (k == 2) {
        //@todo: Revise pending resolution of https://github.com/Cantera/enhancements/issues/229
        return 0.5 * m_thermo->minTemp();
    } else {
        return Reactor::lowerBound(k);
    }
}

}
