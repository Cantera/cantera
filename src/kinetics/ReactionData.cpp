//! @file ReactionData.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/ReactionData.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/global.h"

namespace Cantera
{

void ReactionData::update(double T, double extra)
{
    throw NotImplementedError("ReactionData::update",
        "ReactionData type does not use extra scalar argument.");
}

void ReactionData::update(double T, const vector_fp& extra)
{
    throw NotImplementedError("ReactionData::update",
        "ReactionData type does not use extra vector argument.");
}

void ReactionData::perturbTemperature(double deltaT)
{
    if (m_temperature_buf > 0.) {
        throw CanteraError("ReactionData::perturbTemperature",
            "Cannot apply another perturbation as state is already perturbed.");
    }
    m_temperature_buf = temperature;
    ReactionData::update(temperature * (1. + deltaT));
}

void ReactionData::restore()
{
    // only restore if there is a valid buffered value
    if (m_temperature_buf < 0.) {
        return;
    }
    ReactionData::update(m_temperature_buf);
    m_temperature_buf = -1.;
}

BlowersMaselData::BlowersMaselData()
    : ready(false)
    , density(NAN)
    , m_state_mf_number(-1)
{
}

void BlowersMaselData::update(double T) {
    warn_user("BlowersMaselData::update",
        "This method does not update the change of reaction enthalpy.");
    ReactionData::update(T);
}

bool BlowersMaselData::update(const ThermoPhase& phase, const Kinetics& kin)
{
    double rho = phase.density();
    int mf = phase.stateMFNumber();
    double T = phase.temperature();
    bool changed = false;
    if (T != temperature) {
        ReactionData::update(T);
        changed = true;
    }
    if (changed || rho != density || mf != m_state_mf_number) {
        density = rho;
        m_state_mf_number = mf;
        phase.getPartialMolarEnthalpies(partialMolarEnthalpies.data());
        changed = true;
    }
    return changed;
}

void PlogData::update(double T)
{
    throw CanteraError("PlogData::update",
        "Missing state information: 'PlogData' requires pressure.");
}

bool PlogData::update(const ThermoPhase& phase, const Kinetics& kin)
{
    double T = phase.temperature();
    double P = phase.pressure();
    if (P != pressure || T != temperature) {
        update(T, P);
        return true;
    }
    return false;
}

void PlogData::perturbPressure(double deltaP)
{
    if (m_pressure_buf > 0.) {
        throw CanteraError("PlogData::perturbPressure",
            "Cannot apply another perturbation as state is already perturbed.");
    }
    m_pressure_buf = pressure;
    update(temperature, pressure * (1. + deltaP));
}

void PlogData::restore()
{
    ReactionData::restore();
    // only restore if there is a valid buffered value
    if (m_pressure_buf < 0.) {
        return;
    }
    update(temperature, m_pressure_buf);
    m_pressure_buf = -1.;
}

void ChebyshevData::update(double T)
{
    throw CanteraError("ChebyshevData::update",
        "Missing state information: 'ChebyshevData' requires pressure.");
}

bool ChebyshevData::update(const ThermoPhase& phase, const Kinetics& kin)
{
    double T = phase.temperature();
    double P = phase.pressure();
    if (P != pressure || T != temperature) {
        update(T, P);
        return true;
    }
    return false;
}

void ChebyshevData::perturbPressure(double deltaP)
{
    if (m_pressure_buf > 0.) {
        throw CanteraError("ChebyshevData::perturbPressure",
            "Cannot apply another perturbation as state is already perturbed.");
    }
    m_pressure_buf = pressure;
    update(temperature, pressure * (1. + deltaP));
}

void ChebyshevData::restore()
{
    ReactionData::restore();
    // only restore if there is a valid buffered value
    if (m_pressure_buf < 0.) {
        return;
    }
    update(temperature, m_pressure_buf);
    m_pressure_buf = -1.;
}

}
