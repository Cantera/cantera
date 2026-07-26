//! @file VibrationalRelaxationRate.cpp
//! @since New in %Cantera 4.0

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/VibrationalRelaxationRate.h"

namespace Cantera
{

namespace
{

const string WhereSetParameters = "VibrationalRelaxationRate::setParameters";
const string WhereSetContext = "VibrationalRelaxationRate::setContext";

void requireKeys(const AnyMap& node, const string& model,
                 const string& where, std::initializer_list<string> keys)
{
    for (const auto& key : keys) {
        if (!node.hasKey(key)) {
            throw InputFileError(where, node,
                "Missing required key '{}' for vibration-model '{}'.", key, model);
        }
    }
}

void forbidKeys(const AnyMap& node, const string& model,
                const string& where, std::initializer_list<string> keys)
{
    for (const auto& key : keys) {
        if (node.hasKey(key)) {
            throw InputFileError(where, node,
                "Key '{}' is not allowed for vibration-model '{}'.", key, model);
        }
    }
}
} // end of namespace where all the helpers and safety functions are defined.

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

// Constructor
VibrationalRelaxationRate::VibrationalRelaxationRate(
    double A, double B, double C, double D, double b)
    : m_A(A) , m_b(b) , m_B(B) , m_C(C), m_D(D)
{
    m_valid = true;
}

// AnyMap constructor
VibrationalRelaxationRate::VibrationalRelaxationRate(
    const AnyMap& node, const UnitStack& rate_units)
    : VibrationalRelaxationRate()
{
    setParameters(node, rate_units);
}

void VibrationalRelaxationRate::setParameters(const AnyMap& node,
                                              const UnitStack& rate_units)
{
    // The reaction rate type is always:
    //
    //   type: vibrational-relaxation
    //
    // The physical model is selected separately by:
    //
    //   vibration-model: constant
    //   vibration-model: multi-state-resolved
    //   vibration-model: Starikovskiy
    //   vibration-model: Castela
    //
    // The default model is multi-state-resolved.

    m_valid = false;
    ReactionRate::setParameters(node, rate_units);

    m_negativeA_ok = node.getBool("negative-A", false);
    m_vibration_model = node.getString("vibration-model", "multi-state-resolved");
    const auto& rateMap = node["rate-constant"].as<AnyMap>();

    if (m_vibration_model == "constant") {
        setConstantParameters(node, rateMap);
    } else if (m_vibration_model == "multi-state-resolved") {
        setMultiStateParameters(node, rateMap);
    } else if (m_vibration_model == "Starikovskiy") {
        setStarikovskiyParameters(node, rateMap);
    } else if (m_vibration_model == "Castela") {
        setCastelaParameters(node, rateMap);
    } else {
        throw InputFileError(WhereSetParameters, node,
            "Unrecognized vibration-model '{}'. Expected "
            "'multi-state-resolved', 'Starikovskiy', 'Castela', "
            "or 'constant'.",
            m_vibration_model);
    }

    m_valid = true;
}

void VibrationalRelaxationRate::setConstantParameters(
    const AnyMap& node, const AnyMap& rateMap)
{
    // Constant model:
    //
    //   k(T) = A
    requireKeys(rateMap, m_vibration_model, WhereSetParameters, {"A"});

    forbidKeys(rateMap, m_vibration_model, WhereSetParameters,
        {"b", "n", "B", "C", "D", "m",
         "E", "z", "Ea", "K", "a", "reference-pressure"});

    m_A = node.units().convertRateCoeff(rateMap["A"], conversionUnits());
    m_b = 0.0;
    m_B = 0.0;
    m_C = 0.0;
    m_D = 0.0;
    m_m = 2.0 / 3.0;
    m_E = 0.0;
    m_z = 1.0;
}

void VibrationalRelaxationRate::setMultiStateParameters(
    const AnyMap& node, const AnyMap& rateMap)
{
    // Detailed VV/VT relaxation model:
    //
    //   k(T) = A * exp(
    //       b * log(T)
    //       + B
    //       + C * T^(-1/3)
    //       + D * T^(-2/3)
    //   )
    requireKeys(rateMap, m_vibration_model, WhereSetParameters, {"A", "b"});

    forbidKeys(rateMap, m_vibration_model, WhereSetParameters,
        {"n", "m", "E", "z", "Ea", "K", "a", "reference-pressure"});

    m_A = node.units().convertRateCoeff(rateMap["A"], conversionUnits());
    m_b = rateMap["b"].asDouble();
    m_B = rateMap.getDouble("B", 0.0);
    m_C = rateMap.getDouble("C", 0.0);
    m_D = rateMap.getDouble("D", 0.0);
    m_m = 2.0 / 3.0;
    m_E = 0.0;
    m_z = 1.0;
}

void VibrationalRelaxationRate::setStarikovskiyParameters(
    const AnyMap& node, const AnyMap& rateMap)
{
    // User-facing formula:
    //
    //   k(T) = A * T^n * exp(
    //       K
    //       + B * T^(-1/3)
    //       + C * T^(-m)
    //       + D * T^(-z)
    //   )
    //
    // B, C, and D are signed coefficients read directly from YAML.
    requireKeys(rateMap, m_vibration_model, WhereSetParameters, {"A"});

    forbidKeys(rateMap, m_vibration_model, WhereSetParameters, {"b", "Ea",
        "a", "reference-pressure"});

    const double m = rateMap.getDouble("m", 1.0);
    const double z = rateMap.getDouble("z", 1.0);

    if (m <= 0.0 || z <= 0.0) {
        throw InputFileError(WhereSetParameters, node,
            "The Starikovskiy exponents 'm' and 'z' must be positive.");
    }

    m_A = node.units().convertRateCoeff(rateMap["A"], conversionUnits());
    m_b = rateMap.getDouble("n", 0.0);
    m_B = rateMap.getDouble("K", 0.0);
    m_C = rateMap.getDouble("B", 0.0);
    m_D = rateMap.getDouble("C", 0.0);
    m_m = m;
    m_E = rateMap.getDouble("D", 0.0);
    m_z = z;
}

void VibrationalRelaxationRate::setCastelaParameters(
    const AnyMap& node, const AnyMap& rateMap)
{
    // Castela model:
    //
    // Original relaxation time:
    //
    //   tau_k = p0 / p_k
    //           * exp[a_k * (T^(-1/3) - b_k) - 18.42]
    //
    // Equivalent bimolecular rate coefficient:
    //
    //   k_k(T) = R T / p0
    //            * exp[18.42 + a_k b_k - a_k T^(-1/3)]
    //
    // Internal mapping:
    //
    //   A = R / p0
    //   b = 1
    //   B = 18.42 + a_k b_k
    //   C = -a_k
    //   D = 0
    //   E = 0
    requireKeys(rateMap, m_vibration_model, WhereSetParameters, {"a", "b"});

    forbidKeys(rateMap, m_vibration_model, WhereSetParameters,
        {"A", "n", "K", "B", "C", "D", "m",
         "E", "z", "Ea"});

    m_castela_a = rateMap["a"].asDouble();
    m_castela_b = rateMap["b"].asDouble();

    if (rateMap.hasKey("reference-pressure")) {
        m_referencePressure = rateMap.convert("reference-pressure", "Pa");
    } else {
        m_referencePressure = OneAtm;
    }

    if (m_referencePressure <= 0.0) {
        throw InputFileError(WhereSetParameters, node,
            "Castela reference-pressure must be positive.");
    }

    m_A = GasConstant / m_referencePressure;
    m_b = 1.0;
    m_B = 18.42 + m_castela_a * m_castela_b;
    m_C = -m_castela_a;
    m_D = 0.0;
    m_m = 2.0 / 3.0;
    m_E = 0.0;
    m_z = 1.0;
}

void VibrationalRelaxationRate::getParameters(AnyMap& node) const
{
    if (!valid()) {
        return;
    }

    if (m_negativeA_ok) {
        node["negative-A"] = true;
    }

    node["vibration-model"] = m_vibration_model;

    AnyMap rateNode;

    if (m_vibration_model == "constant") {
        getConstantParameters(rateNode);
    } else if (m_vibration_model == "multi-state-resolved") {
        getMultiStateParameters(rateNode);
    } else if (m_vibration_model == "Starikovskiy") {
        getStarikovskiyParameters(rateNode);
    } else if (m_vibration_model == "Castela") {
        getCastelaParameters(rateNode);
    }

    rateNode.setFlowStyle();
    node["rate-constant"] = std::move(rateNode);
}

void VibrationalRelaxationRate::getConstantParameters(AnyMap& rateNode) const
{
    if (conversionUnits().factor() != 0.0) {
        rateNode["A"].setQuantity(m_A, conversionUnits());
    } else {
        rateNode["A"] = m_A;
        rateNode["__unconvertible__"] = true;
    }
}

void VibrationalRelaxationRate::getMultiStateParameters(
    AnyMap& rateNode) const
{
    if (conversionUnits().factor() != 0.0) {
        rateNode["A"].setQuantity(m_A, conversionUnits());
    } else {
        rateNode["A"] = m_A;
        rateNode["__unconvertible__"] = true;
    }

    rateNode["b"] = m_b;
    rateNode["B"] = m_B;
    rateNode["C"] = m_C;
    rateNode["D"] = m_D;
}

void VibrationalRelaxationRate::getStarikovskiyParameters(
    AnyMap& rateNode) const
{
    if (conversionUnits().factor() != 0.0) {
        rateNode["A"].setQuantity(m_A, conversionUnits());
    } else {
        rateNode["A"] = m_A;
        rateNode["__unconvertible__"] = true;
    }

    rateNode["n"] = m_b;
    rateNode["K"] = m_B;
    rateNode["B"] = m_C;
    rateNode["C"] = m_D;
    rateNode["m"] = m_m;
    rateNode["D"] = m_E;
    rateNode["z"] = m_z;
}

void VibrationalRelaxationRate::getCastelaParameters(AnyMap& rateNode) const
{
    rateNode["a"] = m_castela_a;
    rateNode["b"] = m_castela_b;
    rateNode["reference-pressure"].setQuantity(m_referencePressure, "Pa");
}

double VibrationalRelaxationRate::ddTScaledFromStruct(
    const VibrationalRelaxationData& shared_data) const
{

    return m_b * shared_data.recipT
           - (m_C / 3.0) * shared_data.recipT13 * shared_data.recipT
           - m_m * m_D * std::pow(shared_data.recipT, m_m) * shared_data.recipT
           - m_z * m_E * std::pow(shared_data.recipT, m_z) * shared_data.recipT;
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

} // namespace Cantera