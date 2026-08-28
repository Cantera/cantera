//! @file VibrationalRelaxationRate.cpp
//! @since New in %Cantera 4.0

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/VibrationalRelaxationRate.h"

namespace Cantera
{

namespace
{

const string WhereSetContext =
    "VibrationalRelaxationRate::setContext";

} // namespace

void VibrationalRelaxationData::update(double T)
{
    ReactionData::update(T);
    recipT13 = std::cbrt(recipT);
}

bool VibrationalRelaxationData::update(const ThermoPhase& phase, const Kinetics& kin)
{
    const double T = phase.temperature();

    if (T == temperature) {
        return false;
    }

    update(T);
    return true;
}

// Default constructor
VibrationalRelaxationRate::VibrationalRelaxationRate() = default;

VibrationalRelaxationRate::VibrationalRelaxationRate(
    double A, double b, double C0, double C13,
    double Cm, double m, double Cn, double n)
    : m_A(A), m_b(b) , m_C0(C0), m_C13(C13), m_Cm(Cm),
    m_m(m), m_Cn(Cn), m_n(n)
{
    m_valid = true;
}

void VibrationalRelaxationRate::setParameters(const AnyMap& node,
                                              const UnitStack& rate_units)
{
    m_valid = false; // m_valid will be set to true by the overriding class's setParameters() function.

    ReactionRate::setParameters(node, rate_units);

    m_negativeA_ok = node.getBool("negative-A", false);
}

void VibrationalRelaxationRate::getParameters(AnyMap& node) const
{
    if (!valid()) {
        return;
    }

    if (m_negativeA_ok) {
        node["negative-A"] = true;
    }
}

void VibrationalRelaxationRate::getPreExponentialFactor(
    AnyMap& rateNode) const
{
    if (conversionUnits().factor() != 0.0) {
        rateNode["A"].setQuantity(m_A, conversionUnits());
    } else {
        rateNode["A"] = m_A;
        rateNode["__unconvertible__"] = true;
    }
}

double VibrationalRelaxationRate::ddTScaledFromStruct(
    const VibrationalRelaxationData& shared_data) const
{

    return m_b * shared_data.recipT
           - (m_C13 / 3.0) * shared_data.recipT13 * shared_data.recipT
           - m_m * m_Cm * std::pow(shared_data.recipT, m_m) * shared_data.recipT
           - m_n * m_Cn * std::pow(shared_data.recipT, m_n) * shared_data.recipT;
}

void VibrationalRelaxationRate::setContext(const Reaction& rxn, const Kinetics& kin)
{
    if (rxn.reversible) {
        throw InputFileError(WhereSetContext, rxn.input,
            "Vibrational relaxation rates do not support reversible "
            "reactions.");
    }
}

void VibrationalRelaxationRate::check(const string& equation)
{
    if (!m_negativeA_ok && m_A < 0.0) {
        if (equation.empty()) {
            throw CanteraError(
                "VibrationalRelaxationRate::check",
                "Detected negative pre-exponential constant A={}. "
                "Set 'negative-A: true' to allow it.",
                m_A);
        }

        throw InputFileError(
            "VibrationalRelaxationRate::check", m_input,
            "Undeclared negative leading pre-exponential constant found "
            "in reaction '{}'.",
            equation);
    }
}

void VibrationalRelaxationRate::validate(
    const string& equation, const Kinetics& kin)
{
    if (!valid()) {
        throw InputFileError(
            "VibrationalRelaxationRate::validate", m_input,
            "Rate object for reaction '{}' is not configured.",
            equation);
    }
}

void VibrationalRelaxationRate::requireKeys(const AnyMap& node,
    const string& rateType, const string& where,
    std::initializer_list<string> keys)
{
    for (const auto& key : keys) {
        if (!node.hasKey(key)) {
            throw InputFileError(
                where,
                node,
                "Missing required key '{}' for reaction rate type '{}'.",
                key,
                rateType);
        }
    }
}


void VibrationalRelaxationRate::forbidKeys(const AnyMap& node,
    const string& rateType, const string& where,
    std::initializer_list<string> keys)
{
    for (const auto& key : keys) {
        if (node.hasKey(key)) {
            throw InputFileError(
                where,
                node,
                "Key '{}' is not allowed for reaction rate type '{}'.",
                key,
                rateType);
        }
    }
}

} // namespace Cantera