//! @file CastelaVibrationalRelaxationRate.cpp
//! @since New in %Cantera 4.0

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/CastelaVibrationalRelaxationRate.h"

namespace Cantera
{

namespace
{

const string WhereSetParameters =
    "CastelaVibrationalRelaxationRate::setParameters";

} // namespace

CastelaVibrationalRelaxationRate::
    CastelaVibrationalRelaxationRate(
        const AnyMap& node, const UnitStack& rate_units){
    setParameters(node, rate_units);
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

void CastelaVibrationalRelaxationRate::setParameters(
    const AnyMap& node, const UnitStack& rate_units){
    VibrationalRelaxationRate::setParameters(node, rate_units);

    const auto& rateMap = node["rate-constant"].as<AnyMap>();

    requireKeys(rateMap, type(), WhereSetParameters,{"a", "b"});

    forbidKeys(rateMap, type(), WhereSetParameters,
        {"A", "n", "K", "B", "C", "D",
        "m", "E", "z", "Ea"});

    m_castela_a = rateMap["a"].asDouble();

    m_castela_b = rateMap["b"].asDouble();

    if (rateMap.hasKey("reference-pressure")) {
        m_referencePressure =
            rateMap.convert(
                "reference-pressure", "Pa");
    } else {
        m_referencePressure = OneAtm;
    }

    if (m_referencePressure <= 0.0) {
        throw InputFileError(
            WhereSetParameters,
            node,
            "Castela reference-pressure must be positive.");
    }

    // Castela relaxation time:
    //
    // tau_k = p0 / p_k
    //         * exp[a_k * (T^(-1/3) - b_k) - 18.42]
    //
    // Equivalent rate:
    //
    // k_k(T) = R T / p0
    //          * exp[18.42 + a_k b_k
    //                - a_k T^(-1/3)]
    //
    // Generic mapping:
    //
    // A = R / p0
    // b = 1
    // B = 18.42 + a_k b_k
    // C = -a_k
    // D = 0
    // E = 0

    m_A = GasConstant / m_referencePressure;

    m_b = 1.0;

    m_B = 18.42 + m_castela_a * m_castela_b;

    m_C = -m_castela_a;

    m_D = 0.0;
    m_m = 2.0 / 3.0;

    m_E = 0.0;
    m_z = 1.0;

    m_valid = true;
}


void CastelaVibrationalRelaxationRate::getParameters(
    AnyMap& node) const
{
    if (!valid()) {
        return;
    }

    VibrationalRelaxationRate::getParameters(node);

    AnyMap rateNode;

    rateNode["a"] = m_castela_a;

    rateNode["b"] = m_castela_b;

    rateNode["reference-pressure"].setQuantity(
        m_referencePressure, "Pa");

    rateNode.setFlowStyle();
    node["rate-constant"] = std::move(rateNode);
}

} // namespace Cantera