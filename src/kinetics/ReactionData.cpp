//! @file ReactionData.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/ReactionData.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/base/ctexceptions.h"

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

bool ArrheniusData::update(const ThermoPhase& phase, const Kinetics& kin)
{
    double T = phase.temperature();
    if (T == temperature) {
        return false;
    }
    update(T);
    return true;
}

bool TwoTempPlasmaData::update(const ThermoPhase& phase, const Kinetics& kin)
{
    double T = phase.temperature();
    double Te = phase.electronTemperature();
    bool changed = false;
    if (T != temperature) {
        ReactionData::update(T);
        changed = true;
    }
    if (Te != electronTemp) {
        updateTe(Te);
        changed = true;
    }
    return changed;
}

void TwoTempPlasmaData::update(double T)
{
    throw CanteraError("TwoTempPlasmaData::update",
        "Missing state information: 'TwoTempPlasmaData' requires electron temperature.");
}

void TwoTempPlasmaData::update(double T, double Te)
{
    ReactionData::update(T);
    updateTe(Te);
}

void TwoTempPlasmaData::updateTe(double Te)
{
    electronTemp = Te;
    logTe = std::log(Te);
    recipTe = 1./Te;
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
        phase.getPartialMolarEnthalpies(grt.data());
        changed = true;
    }
    return changed;
}

FalloffData::FalloffData()
    : ready(false)
    , molar_density(NAN)
    , m_state_mf_number(-1)
    , m_perturbed(false)
{
    conc_3b.resize(1, NAN);
    m_conc_3b_buf.resize(1, NAN);
}

void FalloffData::update(double T)
{
    throw CanteraError("FalloffData::update",
        "Missing state information: 'FalloffData' requires third-body concentration.");
}

void FalloffData::update(double T, double M)
{
    ReactionData::update(T);
    conc_3b[0] = M;
}

bool FalloffData::update(const ThermoPhase& phase, const Kinetics& kin)
{
    double rho_m = phase.molarDensity();
    int mf = phase.stateMFNumber();
    double T = phase.temperature();
    bool changed = false;
    if (T != temperature) {
        ReactionData::update(T);
        changed = true;
    }
    if (rho_m != molar_density || mf != m_state_mf_number) {
        molar_density = rho_m;
        m_state_mf_number = mf;
        conc_3b = kin.thirdBodyConcentrations();
        changed = true;
    }
    return changed;
}

void FalloffData::perturbThirdBodies(double deltaM)
{
    if (m_perturbed) {
        throw CanteraError("FalloffData::perturbThirdBodies",
            "Cannot apply another perturbation as state is already perturbed.");
    }
    m_conc_3b_buf = conc_3b;
    for (auto& c3b : conc_3b) {
        c3b *= 1. + deltaM;
    }
    m_perturbed = true;
}

void FalloffData::restore()
{
    ReactionData::restore();
    // only restore if there is a valid buffered value
    if (!m_perturbed) {
        return;
    }
    conc_3b = m_conc_3b_buf;
    m_perturbed = false;
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

CoverageData::CoverageData()
    : sqrtT(NAN)
{
}

void CoverageData::update(double T)
{
    throw CanteraError("CoverageData::update",
        "Missing state information: 'CoverageData' requires species coverages.");
}

void CoverageData::update(double T, const vector_fp& values)
{
    warn_user("CoverageData::update",
        "This method does not update the site density.");
    ReactionData::update(T);
    sqrtT = sqrt(T);
    if (coverages.size() == 0) {
        coverages = values;
        logCoverages.resize(values.size());
    } else if (values.size() == coverages.size()) {
        std::copy(values.begin(), values.end(), coverages.begin());
    } else {
        throw CanteraError("CoverageData::update",
            "Incompatible lengths of coverage arrays: received {} elements while "
            "{} are required.", values.size(), coverages.size());
    }
    for (size_t n = 0; n < coverages.size(); n++) {
        logCoverages[n] = std::log(std::max(coverages[n], Tiny));
    }
}

bool CoverageData::update(const ThermoPhase& phase, const Kinetics& kin)
{
    int mf = phase.stateMFNumber();
    double T = phase.temperature();
    bool changed = false;
    const auto& surf = dynamic_cast<const SurfPhase&>(
        kin.thermo(kin.surfacePhaseIndex()));
    double site_density = surf.siteDensity();
    if (density != site_density) {
        density = surf.siteDensity();
        changed = true;
    }
    if (T != temperature) {
        ReactionData::update(T);
        sqrtT = sqrt(T);
        changed = true;
    }
    if (changed || mf != m_state_mf_number) {
        phase.getActivityConcentrations(coverages.data());
        surf.getCoverages(coverages.data());
        for (size_t n = 0; n < coverages.size(); n++) {
            logCoverages[n] = std::log(std::max(coverages[n], Tiny));
        }
        for (size_t n = 0; n < kin.nPhases(); n++) {
            kin.thermo(n).getPartialMolarEnthalpies(grt.data() + kin.startIndex(n));
        }
        m_state_mf_number = mf;
        changed = true;
    }
    return changed;
}

void CoverageData::perturbTemperature(double deltaT)
{
    throw NotImplementedError("CoverageData::perturbTemperature");
}

}
