//! @file Arrhenius.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/Arrhenius.h"
#include "cantera/thermo/ThermoPhase.h"

namespace Cantera
{

ArrheniusBase::ArrheniusBase(double A, double b, double Ea)
    : m_A(A)
    , m_b(b)
    , m_Ea_R(Ea / GasConstant)
{
    if (m_A > 0.0) {
        m_logA = std::log(m_A);
    }
    m_valid = true;
}

ArrheniusBase::ArrheniusBase(const AnyValue& rate, const UnitSystem& units,
                             const UnitStack& rate_units)
{
    setRateUnits(rate_units);
    setRateParameters(rate, units, rate_units);
}

ArrheniusBase::ArrheniusBase(const AnyMap& node, const UnitStack& rate_units)
{
    setParameters(node, rate_units);
}

void ArrheniusBase::setRateParameters(
    const AnyValue& rate, const UnitSystem& units, const UnitStack& rate_units)
{
    m_Ea_R = 0.; // assume zero if not provided
    m_E4_R = 0.; // assume zero if not provided
    if (rate.empty()) {
        m_A = NAN;
        m_b = NAN;
        m_logA = NAN;
        setRateUnits(Units(0.));
        return;
    }

    if (rate.is<AnyMap>()) {

        auto& rate_map = rate.as<AnyMap>();
        m_A = units.convertRateCoeff(rate_map[m_A_str], conversionUnits());
        m_b = rate_map[m_b_str].asDouble();
        if (rate_map.hasKey(m_Ea_str)) {
            m_Ea_R = units.convertActivationEnergy(rate_map[m_Ea_str], "K");
        }
        if (rate_map.hasKey(m_E4_str)) {
            m_E4_R = units.convertActivationEnergy(rate_map[m_E4_str], "K");
        }
    } else {
        auto& rate_vec = rate.asVector<AnyValue>(2, 4);
        m_A = units.convertRateCoeff(rate_vec[0], conversionUnits());
        m_b = rate_vec[1].asDouble();
        if (rate_vec.size() > 2) {
            m_Ea_R = units.convertActivationEnergy(rate_vec[2], "K");
        }
        if (rate_vec.size() > 3) {
            m_E4_R = units.convertActivationEnergy(rate_vec[3], "K");
        }
    }
    if (m_A > 0.0) {
        m_logA = std::log(m_A);
    }
    m_valid = true;
}

void ArrheniusBase::getRateParameters(AnyMap& node) const
{
    if (!valid()) {
        // Return empty/unmodified AnyMap
        return;
    }

    if (conversionUnits().factor() != 0.0) {
        node[m_A_str].setQuantity(m_A, conversionUnits());
    } else {
        node[m_A_str] = m_A;
        // This can't be converted to a different unit system because the dimensions of
        // the rate constant were not set. Can occur if the reaction was created outside
        // the context of a Kinetics object and never added to a Kinetics object.
        node["__unconvertible__"] = true;
    }
    node[m_b_str] = m_b;
    node[m_Ea_str].setQuantity(m_Ea_R, "K", true);
    if (m_E4_str != "") {
        node[m_E4_str].setQuantity(m_E4_R, "K", true);
    }
    node.setFlowStyle();
}

void ArrheniusBase::setParameters(const AnyMap& node, const UnitStack& rate_units)
{
    ReactionRate::setParameters(node, rate_units);
    m_negativeA_ok = node.getBool("negative-A", false);
    if (!node.hasKey("rate-constant")) {
        setRateParameters(AnyValue(), node.units(), rate_units);
        return;
    }
    setRateParameters(node["rate-constant"], node.units(), rate_units);
}

void ArrheniusBase::getParameters(AnyMap& node) const {
    if (m_negativeA_ok) {
        node["negative-A"] = true;
    }
    AnyMap rateNode;
    getRateParameters(rateNode);
    if (!rateNode.empty()) {
        // RateType object is configured
        node["rate-constant"] = std::move(rateNode);
    }
}

void ArrheniusBase::check(const string& equation)
{
    if (!m_negativeA_ok && m_A < 0) {
        if (equation == "") {
            throw CanteraError("ArrheniusBase::check",
                "Detected negative pre-exponential factor (A={}).\n"
                "Enable 'allowNegativePreExponentialFactor' to suppress "
                "this message.", m_A);
        }
        throw InputFileError("ArrheniusBase::check", m_input,
            "Undeclared negative pre-exponential factor found in reaction '{}'",
            equation);
    }
}

void ArrheniusBase::validate(const string& equation, const Kinetics& kin)
{
    if (!valid()) {
        throw InputFileError("ArrheniusBase::validate", m_input,
            "Rate object for reaction '{}' is not configured.", equation);
    }
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

}
