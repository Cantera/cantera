//! @file ReactionData.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/ReactionData.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/base/ctexceptions.h"

namespace Cantera
{

void ReactionData::update(double T, double P, double extra)
{
    throw NotImplementedError("ReactionData::update",
        "Update with extra argument not implemented.");
}

bool ArrheniusData::update(const ThermoPhase& bulk, const Kinetics& kin)
{
    double T = bulk.temperature();
    pressure = bulk.pressure();
    if (T == temperature) {
        return false;
    }
    update(T);
    return true;
}

BlowersMaselData::BlowersMaselData()
    : ready(false)
    , density(NAN)
    , m_state_mf_number(-1)
{
    dH.resize(1, NAN);
}

void BlowersMaselData::update(double T)
{
    throw CanteraError("BlowersMaselData::update",
        "Missing state information: reaction type requires temperature, "
        "pressure\nand enthalpy change.");
}

void BlowersMaselData::update(double T, double P)
{
    throw CanteraError("BlowersMaselData::update",
        "Missing state information: reaction type requires temperature, "
        "pressure\nand enthalpy change.");
}

void BlowersMaselData::update(double T, double P, double deltaH)
{
    if (ready) {
        throw CanteraError("BlowersMaselData::update", "This should not be reachable.");
    }
    ReactionData::update(T, P);
    dH[0] = deltaH;
}

bool BlowersMaselData::update(const ThermoPhase& bulk, const Kinetics& kin)
{
    double rho = bulk.density();
    int mf = bulk.stateMFNumber();
    double T = bulk.temperature();
    bool changed = false;
    if (T != temperature) {
        ReactionData::update(T);
        changed = true;
    }
    if (changed || rho != density || mf != m_state_mf_number) {
        density = rho;
        m_state_mf_number = mf;
        bulk.getPartialMolarEnthalpies(m_grt.data());
        kin.getReactionDelta(m_grt.data(), dH.data());
        changed = true;
    }
    pressure = bulk.pressure();
    return changed;
}

FalloffData::FalloffData()
    : ready(false)
    , molar_density(NAN)
    , m_state_mf_number(-1)
{
    conc_3b.resize(1, NAN);
}

void FalloffData::update(double T)
{
    throw CanteraError("FalloffData::update",
        "Missing state information: reaction type requires temperature, "
        "pressure\nand third-body concentration.");
}

void FalloffData::update(double T, double P)
{
    throw CanteraError("FalloffData::update",
        "Missing state information: reaction type requires temperature, "
        "pressure\nand third-body concentration.");
}

void FalloffData::update(double T, double P, double M)
{
    if (ready) {
        throw CanteraError("FalloffData::update", "This should not be reachable.");
    }
    ReactionData::update(T, P);
    conc_3b[0] = M;
}

bool FalloffData::update(const ThermoPhase& bulk, const Kinetics& kin)
{
    double rho_m = bulk.molarDensity();
    int mf = bulk.stateMFNumber();
    double T = bulk.temperature();
    bool changed = false;
    if (T != temperature) {
        ReactionData::update(T);
        changed = true;
    }
    if (rho_m != molar_density || mf != m_state_mf_number) {
        molar_density = rho_m;
        m_state_mf_number = mf;
        kin.getThirdBodyConcentrations(conc_3b.data());
        changed = true;
    }
    pressure = bulk.pressure();
    return changed;
}

void PlogData::update(double T)
{
    throw CanteraError("PlogData::update",
        "Missing state information: reaction type requires pressure.");
}

bool PlogData::update(const ThermoPhase& bulk, const Kinetics& kin)
{
    double T = bulk.temperature();
    double P = bulk.pressure();
    if (P != pressure || T != temperature) {
        update(T, P);
        return true;
    }
    return false;
}

void ChebyshevData::update(double T)
{
    throw CanteraError("ChebyshevData::update",
        "Missing state information: reaction type requires pressure.");
}

bool ChebyshevData::update(const ThermoPhase& bulk, const Kinetics& kin)
{
    double T = bulk.temperature();
    double P = bulk.pressure();
    if (P != pressure || T != temperature) {
        update(T, P);
        return true;
    }
    return false;
}

}
