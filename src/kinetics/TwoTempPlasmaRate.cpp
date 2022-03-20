//! @file TwoTempPlasmaRate.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/TwoTempPlasmaRate.h"
#include "cantera/kinetics/Reaction.h"
#include "cantera/thermo/ThermoPhase.h"

namespace Cantera
{

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

TwoTempPlasmaRate::TwoTempPlasmaRate()
    : ArrheniusBase()
{
    m_Ea_str = "Ea-gas";
    m_E4_str = "Ea-electron";
}

TwoTempPlasmaRate::TwoTempPlasmaRate(double A, double b, double Ea, double EE)
    : ArrheniusBase(A, b, Ea)
{
    m_Ea_str = "Ea-gas";
    m_E4_str = "Ea-electron";
    m_E4_R = EE / GasConstant;
}

double TwoTempPlasmaRate::ddTScaledFromStruct(const TwoTempPlasmaData& shared_data) const
{
    warn_user("TwoTempPlasmaRate::ddTScaledFromStruct",
        "Temperature derivative does not consider changes of electron temperature.");
    return (m_Ea_R - m_E4_R) * shared_data.recipT * shared_data.recipT;
}

void TwoTempPlasmaRate::setContext(const Reaction& rxn, const Kinetics& kin)
{
    // TwoTempPlasmaReaction is for a non-equilibrium plasma, and the reverse rate
    // cannot be calculated from the conventional thermochemistry.
    // @todo implement the reversible rate for non-equilibrium plasma
    if (rxn.reversible) {
        throw InputFileError("TwoTempPlasmaRate::setContext", rxn.input,
            "TwoTempPlasmaRate does not support reversible reactions");
    }
}

}
