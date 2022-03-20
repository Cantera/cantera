//! @file BlowersMaselRate.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/BlowersMaselRate.h"
#include "cantera/kinetics/Reaction.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/thermo/ThermoPhase.h"

namespace Cantera
{

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

BlowersMaselRate::BlowersMaselRate()
    : m_deltaH_R(0.)
{
    m_Ea_str = "Ea0";
    m_E4_str = "w";
}

BlowersMaselRate::BlowersMaselRate(double A, double b, double Ea0, double w)
    : ArrheniusBase(A, b, Ea0)
    , m_deltaH_R(0.)
{
    m_Ea_str = "Ea0";
    m_E4_str = "w";
    m_E4_R = w / GasConstant;
}

double BlowersMaselRate::ddTScaledFromStruct(const BlowersMaselData& shared_data) const
{
    warn_user("BlowersMaselRate::ddTScaledFromStruct",
        "Temperature derivative does not consider changes of reaction enthalpy.");
    double Ea_R = effectiveActivationEnergy_R(m_deltaH_R);
    return (Ea_R * shared_data.recipT + m_b) * shared_data.recipT;
}

void BlowersMaselRate::setContext(const Reaction& rxn, const Kinetics& kin)
{
    m_stoich_coeffs.clear();
    for (const auto& sp : rxn.reactants) {
        m_stoich_coeffs.emplace_back(kin.kineticsSpeciesIndex(sp.first), -sp.second);
    }
    for (const auto& sp : rxn.products) {
        m_stoich_coeffs.emplace_back(kin.kineticsSpeciesIndex(sp.first), sp.second);
    }
}

}
