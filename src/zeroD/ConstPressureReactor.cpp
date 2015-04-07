/**
 *  @file ConstPressureReactor.cpp A constant pressure zero-dimensional
 *      reactor
 */

// Copyright 2001  California Institute of Technology

#include "cantera/zeroD/ConstPressureReactor.h"
#include "cantera/zeroD/FlowDevice.h"

using namespace std;

namespace Cantera
{

ConstPressureReactor::ConstPressureReactor() : Reactor() {}

void ConstPressureReactor::getInitialConditions(double t0, size_t leny, double* y)
{
    m_init = true;
    if (m_thermo == 0) {
        throw CanteraError("getInitialConditions",
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
        m_thermo->setState_HP(y[1]/m_mass, m_pressure, 1.0e-4);
    } else {
        m_thermo->setPressure(m_pressure);
    }
    m_vol = m_mass / m_thermo->density();
    updateSurfaceState(y + m_nsp + 2);

    // save parameters needed by other connected reactors
    m_enthalpy = m_thermo->enthalpy_mass();
    m_intEnergy = m_thermo->intEnergy_mass();
    m_thermo->saveState(m_state);
}

void ConstPressureReactor::evalEqs(doublereal time, doublereal* y,
                                   doublereal* ydot, doublereal* params)
{
    double dmdt = 0.0; // dm/dt (gas phase)
    double* dYdt = ydot + 2;

    m_thermo->restoreState(m_state);
    applySensitivity(params);
    evalWalls(time);
    double mdot_surf = evalSurfaces(time, ydot + m_nsp + 2);
    dmdt += mdot_surf;

    const vector_fp& mw = m_thermo->molecularWeights();
    const doublereal* Y = m_thermo->massFractions();

    if (m_chem) {
        m_kin->getNetProductionRates(&m_wdot[0]); // "omega dot"
    }

    for (size_t k = 0; k < m_nsp; k++) {
        // production in gas phase and from surfaces
        dYdt[k] = (m_wdot[k] * m_vol + m_sdot[k]) * mw[k] / m_mass;
        // dilution by net surface mass flux
        dYdt[k] -= Y[k] * mdot_surf / m_mass;
    }

    // external heat transfer
    double dHdt = - m_Q;

    // add terms for outlets
    for (size_t i = 0; i < m_outlet.size(); i++) {
        double mdot_out = m_outlet[i]->massFlowRate(time); // mass flow out of system
        dmdt -= mdot_out;
        dHdt -= mdot_out * m_enthalpy;
    }

    // add terms for inlets
    for (size_t i = 0; i < m_inlet.size(); i++) {
        double mdot_in = m_inlet[i]->massFlowRate(time);
        dmdt += mdot_in; // mass flow into system
        for (size_t n = 0; n < m_nsp; n++) {
            double mdot_spec = m_inlet[i]->outletSpeciesMassFlowRate(n);
            // flow of species into system and dilution by other species
            dYdt[n] += (mdot_spec - mdot_in * Y[n]) / m_mass;
        }
        dHdt += mdot_in * m_inlet[i]->enthalpy_mass();
    }

    ydot[0] = dmdt;
    if (m_energy) {
        ydot[1] = dHdt;
    } else {
        ydot[1] = 0.0;
    }

    // reset sensitivity parameters
    resetSensitivity(params);
}

size_t ConstPressureReactor::componentIndex(const string& nm) const
{
    size_t k = speciesIndex(nm);
    if (k != npos) {
        return k + 2;
    } else if (nm == "m" || nm == "mass") {
        return 0;
    } else if (nm == "H" || nm == "enthalpy") {
        return 1;
    } else {
        return npos;
    }
}

}
