//! @file FlowReactor.cpp A steady-state plug flow reactor

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/FlowReactor.h"
#include "cantera/base/global.h"

using namespace std;

namespace Cantera
{

FlowReactor::FlowReactor() :
    m_speed(0.0),
    m_dist(0.0),
    m_T(0.0),
    m_fctr(1.0e10),
    m_rho0(0.0),
    m_speed0(0.0),
    m_P0(0.0),
    m_h0(0.0)
{
}

void FlowReactor::getState(double* y)
{
    if (m_thermo == 0) {
        throw CanteraError("FlowReactor::getState",
                           "Error: reactor is empty.");
    }
    m_thermo->restoreState(m_state);
    m_thermo->getMassFractions(y+2);
    y[0] = 0.0; // distance

    // set the second component to the initial speed
    y[1] = m_speed0;
}

void FlowReactor::initialize(doublereal t0)
{
    m_thermo->restoreState(m_state);
    m_nv = m_nsp + 2;
}

void FlowReactor::updateState(doublereal* y)
{
    // Set the mass fractions and density of the mixture.
    m_dist = y[0];
    m_speed = y[1];
    doublereal* mss = y + 2;
    m_thermo->setMassFractions(mss);
    doublereal rho = m_rho0 * m_speed0/m_speed;

    // assumes frictionless
    doublereal pmom = m_P0 - rho*m_speed*m_speed;

    doublereal hmom;
    // assumes adiabatic
    if (m_energy) {
        hmom = m_h0 - 0.5*m_speed*m_speed;
        m_thermo->setState_HP(hmom, pmom);
    } else {
        m_thermo->setState_TP(m_T, pmom);
    }
    m_thermo->saveState(m_state);
}

void FlowReactor::evalEqs(doublereal time, doublereal* y,
                          doublereal* ydot, doublereal* params)
{
    m_thermo->restoreState(m_state);
    applySensitivity(params);

    // distance equation
    ydot[0] = m_speed;

    // speed equation. Set m_fctr to a large value, so that rho*u is held fixed
    ydot[1] = m_fctr*(m_speed0 - m_thermo->density()*m_speed/m_rho0);

    // species equations //
    const vector_fp& mw = m_thermo->molecularWeights();

    if (m_chem) {
        m_kin->getNetProductionRates(ydot+2); // "omega dot"
    } else {
        fill(ydot + 2, ydot + 2 + m_nsp, 0.0);
    }
    doublereal rrho = 1.0/m_thermo->density();
    for (size_t n = 0; n < m_nsp; n++) {
        ydot[n+2] *= mw[n]*rrho;
    }
    resetSensitivity(params);
}

size_t FlowReactor::componentIndex(const string& nm) const
{
    // check for a gas species name
    size_t k = m_thermo->speciesIndex(nm);
    if (k != npos) {
        return k + 2;
    } else if (nm == "X" || nm == "distance") {
        return 0;
    } else if (nm == "U" || nm == "velocity") {
        return 1;
    } else {
        return npos;
    }
}

}
