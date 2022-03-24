//! @file ConstPressureReactor.cpp A constant pressure zero-dimensional reactor

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/IdealGasConstPressureReactor.h"
#include "cantera/zeroD/FlowDevice.h"
#include "cantera/zeroD/ReactorNet.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/base/utilities.h"

using namespace std;

namespace Cantera
{

void IdealGasConstPressureReactor::setThermoMgr(ThermoPhase& thermo)
{
    //! @todo: Add a method to ThermoPhase that indicates whether a given
    //! subclass is compatible with this reactor model
    if (thermo.type() != "IdealGas") {
        throw CanteraError("IdealGasConstPressureReactor::setThermoMgr",
                           "Incompatible phase type provided");
    }
    Reactor::setThermoMgr(thermo);
}

void IdealGasConstPressureReactor::getState(double* y)
{
    if (m_thermo == 0) {
        throw CanteraError("IdealGasConstPressureReactor::getState",
                           "Error: reactor is empty.");
    }
    m_thermo->restoreState(m_state);

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

void IdealGasConstPressureReactor::initialize(doublereal t0)
{
    ConstPressureReactor::initialize(t0);
    m_hk.resize(m_nsp, 0.0);
}

void IdealGasConstPressureReactor::updateState(doublereal* y)
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

    m_thermo->restoreState(m_state);
    const vector_fp& mw = m_thermo->molecularWeights();
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

size_t IdealGasConstPressureReactor::componentIndex(const string& nm) const
{
    size_t k = speciesIndex(nm);
    if (k != npos) {
        return k + 2;
    } else if (nm == "mass") {
        return 0;
    } else if (nm == "temperature") {
        return 1;
    } else {
        return npos;
    }
}

std::string IdealGasConstPressureReactor::componentName(size_t k) {
    if (k == 1) {
        return "temperature";
    } else {
        return ConstPressureReactor::componentName(k);
    }
}

}
