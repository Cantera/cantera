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
{
    m_Ea_str = "Ea-gas";
    m_E4_str = "Ea-electron";
}

TwoTempPlasmaRate::TwoTempPlasmaRate(double A, double b, double Ea, double EE,
                                     double bg, double Tinv)
    : ArrheniusBase(A, b, Ea)
{
    m_Ea_str = "Ea-gas";
    m_E4_str = "Ea-electron";
    m_E4_R = EE / GasConstant;
    m_bg = bg;
    if (Tinv != 0) {
        m_recip_Tinv = 1.0 / Tinv;
    } // no else required here because m_recip_Tinv is otherwise initialized at 0.
}

TwoTempPlasmaRate::TwoTempPlasmaRate(const AnyMap& node, const UnitStack& rate_units)
    : TwoTempPlasmaRate()
{
    setParameters(node, rate_units);
}

void TwoTempPlasmaRate::getParameters(AnyMap& node) const
{
    ArrheniusBase::getParameters(node);

    if (!node.hasKey("rate-constant")) {
        return;
    }

    auto& rateNode = node["rate-constant"].as<AnyMap>();

    if (m_bg != 0.0) {
        rateNode["b-gas"] = m_bg;
    }

    if (m_recip_Tinv != 0.0) {
        rateNode["T-inv"] = 1.0/m_recip_Tinv;
    }

    rateNode.setFlowStyle();
}

void TwoTempPlasmaRate::setParameters(const AnyMap& node, const UnitStack& rate_units)
{
    // First, set back to zero these two first parameters to avoid any unwanted
    // remanence should setParameters be called several times and not contain b-gas or
    // T-inv in one of the calls.
    m_bg = 0.0;
    m_recip_Tinv = 0.0;

    // Option 1: there is no rate constant provided
    if (!node.hasKey("rate-constant")) {
        ArrheniusBase::setParameters(node, rate_units);
        return;
    }

    // Option 2: the rate constant is an AnyMap
    if (node["rate-constant"].is<AnyMap>()) {
        ArrheniusBase::setParameters(node, rate_units);

        const auto& rate = node["rate-constant"].as<AnyMap>();

        m_bg = rate.getDouble("b-gas", 0.0);

        double Tinv = 0.0;
        if (rate.hasKey("T-inv")) {
            Tinv = rate.convert("T-inv", "K");
        }

        if (Tinv != 0) {
            m_recip_Tinv = 1.0 / Tinv;
        }
        return;
    }

    // If both other cases fail, the rate constant is
    // a vector and is treated below.
    const auto& rate = node["rate-constant"].asVector<AnyValue>(2, 6);

    // Option 3a: the vector is classical Arrhenius
    if (rate.size() <= 4) {
        ArrheniusBase::setParameters(node, rate_units);
        return;
    }

    ReactionRate::setParameters(node, rate_units);
    m_negativeA_ok = node.getBool("negative-A", false);

    AnyValue baseRate = node["rate-constant"];
    baseRate.asVector<AnyValue>().resize(4);

    setRateParameters(baseRate, node.units(), rate_units);

    // Option 3b: the vector has the b-gas parameter
    m_bg = rate[4].asDouble();

    // Option 3c: the vector has the T-inv parameter
    if (rate.size() == 6) {
        double Tinv = node.units().convert(rate[5], "K");
        m_recip_Tinv = Tinv != 0.0 ? 1.0 / Tinv : 0.0;
    }
}

double TwoTempPlasmaRate::ddTScaledFromStruct(const TwoTempPlasmaData& shared_data) const
{
    warn_user("TwoTempPlasmaRate::ddTScaledFromStruct",
        "Temperature derivative does not consider changes of electron temperature.");
    return m_bg * shared_data.recipT + (m_Ea_R - m_E4_R)
                    * shared_data.recipT * shared_data.recipT - m_recip_Tinv;
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
