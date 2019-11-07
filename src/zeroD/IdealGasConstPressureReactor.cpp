//! @file ConstPressureReactor.cpp A constant pressure zero-dimensional reactor

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/IdealGasConstPressureReactor.h"
#include "cantera/zeroD/FlowDevice.h"

using namespace std;

namespace Cantera
{

void IdealGasConstPressureReactor::setThermoMgr(ThermoPhase& thermo)
{
    //! @TODO: Add a method to ThermoPhase that indicates whether a given
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
    updateSurfaceState(y + m_nsp + 2);

    // save parameters needed by other connected reactors
    m_enthalpy = m_thermo->enthalpy_mass();
    m_intEnergy = m_thermo->intEnergy_mass();
    m_thermo->saveState(m_state);
}

void IdealGasConstPressureReactor::evalEqs(doublereal time, doublereal* y,
                                   doublereal* ydot, doublereal* params)
{
    double dmdt = 0.0; // dm/dt (gas phase)
    double mcpdTdt = 0.0; // m * c_p * dT/dt
    double* dYdt = ydot + 2;

    evalFlowDevices(time);
    applySensitivity(params);
    evalWalls(time);

    m_thermo->restoreState(m_state);
    double mdot_surf = evalSurfaces(time, ydot + m_nsp + 2);
    dmdt += mdot_surf;

    m_thermo->getPartialMolarEnthalpies(&m_hk[0]);
    const vector_fp& mw = m_thermo->molecularWeights();
    const doublereal* Y = m_thermo->massFractions();

    if (m_chem) {
        m_kin->getNetProductionRates(&m_wdot[0]); // "omega dot"
    }

    // external heat transfer
    mcpdTdt -= m_Q;

    for (size_t n = 0; n < m_nsp; n++) {
        // heat release from gas phase and surface reactions
        mcpdTdt -= m_wdot[n] * m_hk[n] * m_vol;
        mcpdTdt -= m_sdot[n] * m_hk[n];
        // production in gas phase and from surfaces
        dYdt[n] = (m_wdot[n] * m_vol + m_sdot[n]) * mw[n] / m_mass;
        // dilution by net surface mass flux
        dYdt[n] -= Y[n] * mdot_surf / m_mass;
    }

    // add terms for outlets
    for (size_t i = 0; i < m_outlet.size(); i++) {
        dmdt -= m_mdot_out[i]; // mass flow out of system
    }

    // add terms for inlets
    for (size_t i = 0; i < m_inlet.size(); i++) {
        dmdt += m_mdot_in[i]; // mass flow into system
        mcpdTdt += m_inlet[i]->enthalpy_mass() * m_mdot_in[i];
        for (size_t n = 0; n < m_nsp; n++) {
            double mdot_spec = m_inlet[i]->outletSpeciesMassFlowRate(n);
            // flow of species into system and dilution by other species
            dYdt[n] += (mdot_spec - m_mdot_in[i] * Y[n]) / m_mass;
            mcpdTdt -= m_hk[n] / mw[n] * mdot_spec;
        }
    }

    ydot[0] = dmdt;
    if (m_energy) {
        ydot[1] = mcpdTdt / (m_mass * m_thermo->cp_mass());
    } else {
        ydot[1] = 0.0;
    }

    resetSensitivity(params);
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
