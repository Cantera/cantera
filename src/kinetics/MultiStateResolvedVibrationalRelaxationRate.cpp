//! @file MultiStateResolvedVibrationalRelaxationRate.cpp
//! @since New in %Cantera 4.0

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/MultiStateResolvedVibrationalRelaxationRate.h"

namespace Cantera
{

namespace
{

const string WhereSetParameters =
    "MultiStateResolvedVibrationalRelaxationRate::setParameters";

} // namespace

MultiStateResolvedVibrationalRelaxationRate::
    MultiStateResolvedVibrationalRelaxationRate(
        double A, double B, double C, double D, double b)
    : VibrationalRelaxationRate(A, B, C, D, b)
{
}

MultiStateResolvedVibrationalRelaxationRate::
    MultiStateResolvedVibrationalRelaxationRate(
        const AnyMap& node,
        const UnitStack& rate_units)
{
    setParameters(node, rate_units);
}

void MultiStateResolvedVibrationalRelaxationRate::setParameters(
    const AnyMap& node,
    const UnitStack& rate_units)
{
    VibrationalRelaxationRate::setParameters(node, rate_units);

    const auto& rateMap =
        node["rate-constant"].as<AnyMap>();

    requireKeys(
        rateMap, type(), WhereSetParameters, {"A", "b"});

    forbidKeys(
        rateMap, type(), WhereSetParameters,
        {"n", "m", "E", "z", "Ea", "K",
        "a", "reference-pressure"});

    m_A = node.units().convertRateCoeff(
        rateMap["A"], conversionUnits());

    m_b = rateMap["b"].asDouble();

    m_B = rateMap.getDouble("B", 0.0);
    m_C = rateMap.getDouble("C", 0.0);
    m_D = rateMap.getDouble("D", 0.0);

    m_m = 2.0 / 3.0;

    m_E = 0.0;
    m_z = 1.0;

    m_valid = true;
}

void MultiStateResolvedVibrationalRelaxationRate::getParameters(
    AnyMap& node) const
{
    if (!valid()) {
        return;
    }

    VibrationalRelaxationRate::getParameters(node);

    AnyMap rateNode;

    getPreExponentialFactor(rateNode);

    rateNode["b"] = m_b;
    rateNode["B"] = m_B;
    rateNode["C"] = m_C;
    rateNode["D"] = m_D;

    rateNode.setFlowStyle();
    node["rate-constant"] = std::move(rateNode);
}

} // namespace Cantera