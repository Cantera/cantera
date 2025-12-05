//! @file ConstPressureReactor.cpp A constant pressure zero-dimensional reactor

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/IdealGasConstPressureReactor.h"
#include "cantera/zeroD/FlowDevice.h"
#include "cantera/zeroD/ReactorNet.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/base/utilities.h"

namespace Cantera
{

void IdealGasConstPressureReactor::getState(double* y)
{
    // set the first component to the total mass
    y[0] = m_thermo->density() * m_vol;

    // set the second component to the temperature
    y[1] = m_thermo->temperature();

    // set components y+2 ... y+K+1 to the mass fractions Y_k of each species
    m_thermo->getMassFractions(y+2);

    // set the remaining components to the surface species
    // coverages on the walls
    getSurfaceInitialConditions(y + m_nsp + 2);
}

void IdealGasConstPressureReactor::initialize(double t0)
{
    //! @todo: Add a method to ThermoPhase that indicates whether a given
    //! subclass is compatible with this reactor model
    if (m_thermo->type() != "ideal-gas") {
        throw CanteraError("IdealGasConstPressureReactor::initialize",
                           "Incompatible phase type '{}' provided", m_thermo->type());
    }    ConstPressureReactor::initialize(t0);
    m_hk.resize(m_nsp, 0.0);
}

void IdealGasConstPressureReactor::updateState(double* y)
{
    // The components of y are [0] the total mass, [1] the temperature,
    // [2...K+2) are the mass fractions of each species, and [K+2...] are the
    // coverages of surface species on each wall.
    m_mass = y[0];
    m_thermo->setMassFractions_NoNorm(y+2);
    m_thermo->setState_TP(y[1], m_pressure);
    m_vol = m_mass / m_thermo->density();
    updateConnected(false);
    updateSurfaceState(y + m_nsp + 2);
}

void IdealGasConstPressureReactor::eval(double time, double* LHS, double* RHS)
{
    double& dmdt = RHS[0]; // dm/dt (gas phase)
    double& mcpdTdt = RHS[1]; // m * c_p * dT/dt
    double* mdYdt = RHS + 2; // mass * dY/dt

    dmdt = 0.0;
    mcpdTdt = 0.0;

    evalWalls(time);
    const vector<double>& mw = m_thermo->molecularWeights();
    const double* Y = m_thermo->massFractions();

    evalSurfaces(LHS + m_nsp + 2, RHS + m_nsp + 2, m_sdot.data());
    double mdot_surf = dot(m_sdot.begin(), m_sdot.end(), mw.begin());
    dmdt += mdot_surf;

    m_thermo->getPartialMolarEnthalpies(&m_hk[0]);

    if (m_chem) {
        m_kin->getNetProductionRates(&m_wdot[0]); // "omega dot"
    }

    // external heat transfer
    mcpdTdt += m_Qdot;

    for (size_t n = 0; n < m_nsp; n++) {
        // heat release from gas phase and surface reactions
        mcpdTdt -= m_wdot[n] * m_hk[n] * m_vol;
        mcpdTdt -= m_sdot[n] * m_hk[n];
        // production in gas phase and from surfaces
        mdYdt[n] = (m_wdot[n] * m_vol + m_sdot[n]) * mw[n];
        // dilution by net surface mass flux
        mdYdt[n] -= Y[n] * mdot_surf;
        //Assign left-hand side of dYdt ODE as total mass
        LHS[n+2] = m_mass;
    }

    // add terms for outlets
    for (auto outlet : m_outlet) {
        dmdt -= outlet->massFlowRate(); // mass flow out of system
    }

    // add terms for inlets
    for (auto inlet : m_inlet) {
        double mdot = inlet->massFlowRate();
        dmdt += mdot; // mass flow into system
        mcpdTdt += inlet->enthalpy_mass() * mdot;
        for (size_t n = 0; n < m_nsp; n++) {
            double mdot_spec = inlet->outletSpeciesMassFlowRate(n);
            // flow of species into system and dilution by other species
            mdYdt[n] += mdot_spec - mdot * Y[n];
            mcpdTdt -= m_hk[n] / mw[n] * mdot_spec;
        }
    }

    if (m_energy) {
        LHS[1] = m_mass * m_thermo->cp_mass();
    } else {
        RHS[1] = 0.0;
    }
}

vector<size_t> IdealGasConstPressureReactor::steadyConstraints() const
{
    if (nSurfs() != 0) {
        throw CanteraError("IdealGasConstPressureReactor::steadyConstraints",
            "Steady state solver cannot currently be used with "
            " IdealGasConstPressureReactor when reactor surfaces are present.\n"
            "See https://github.com/Cantera/enhancements/issues/234");
    }
    if (energyEnabled()) {
        return {0}; // mass
    } else {
        return {0, 1}; // mass and temperature
    }
}

size_t IdealGasConstPressureReactor::componentIndex(const string& nm) const
{
    if (nm == "mass") {
        return 0;
    }
    if (nm == "temperature") {
        return 1;
    }
    try {
        return speciesIndex(nm) + 2;
    } catch (const CanteraError&) {
        throw CanteraError("IdealGasConstPressureReactor::componentIndex",
            "Component '{}' not found", nm);
    }
}

string IdealGasConstPressureReactor::componentName(size_t k) {
    if (k == 1) {
        return "temperature";
    } else {
        return ConstPressureReactor::componentName(k);
    }
}

double IdealGasConstPressureReactor::upperBound(size_t k) const
{
    if (k == 1) {
        //@todo: Revise pending resolution of https://github.com/Cantera/enhancements/issues/229
        return 1.5 * m_thermo->maxTemp();
    } else {
        return ConstPressureReactor::upperBound(k);
    }
}

double IdealGasConstPressureReactor::lowerBound(size_t k) const
{
    if (k == 1) {
        //@todo: Revise pending resolution of https://github.com/Cantera/enhancements/issues/229
        return 0.5 * m_thermo->minTemp();
    } else {
        return ConstPressureReactor::lowerBound(k);
    }
}

}
