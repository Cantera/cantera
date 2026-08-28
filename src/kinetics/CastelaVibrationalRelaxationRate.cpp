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

CastelaVibrationalRelaxationRate::
    CastelaVibrationalRelaxationRate(double a, double b,
        double referencePressure){

            if (referencePressure <= 0.0) {
                throw CanteraError(
                    "CastelaVibrationalRelaxationRate",
                    "Reference pressure must be positive.");
            }

            m_A = GasConstant / referencePressure;
            m_b = 1.0;
            m_C0 = 18.42 + a * b;
            m_C13 = -a;
            m_Cm = 0.0;
            m_m = 2.0 / 3.0;
            m_Cn = 0.0;
            m_n = 1.0;

            m_valid = true;
            m_castela_a = a;
            m_castela_b = b;
            m_referencePressure = referencePressure;
}

void CastelaVibrationalRelaxationRate::setParameters(
    const AnyMap& node, const UnitStack& rate_units){
    VibrationalRelaxationRate::setParameters(node, rate_units);

    const auto& rateMap = node["rate-constant"].as<AnyMap>();

    requireKeys(rateMap, type(), WhereSetParameters,{"a", "b"});

    forbidKeys(rateMap, type(), WhereSetParameters,
        {"A", "n", "C0", "C13", "C23", "Cn",
        "m", "Cm", "Ea"});

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
    // C0 = 18.42 + a_k b_k
    // C13 = -a_k
    // Cm = 0
    // Cn = 0

    m_A = GasConstant / m_referencePressure;

    m_b = 1.0;

    m_C0 = 18.42 + m_castela_a * m_castela_b;

    m_C13 = -m_castela_a;

    m_Cm = 0.0;
    m_m = 2.0 / 3.0;

    m_Cn = 0.0;
    m_n = 1.0;

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