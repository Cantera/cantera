//! @file ConstantVibrationalRelaxationRate.cpp
//! @since New in %Cantera 4.0

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/ConstantVibrationalRelaxationRate.h"

namespace Cantera
{

namespace
{

const string WhereSetParameters =
    "ConstantVibrationalRelaxationRate::setParameters";

} // namespace

ConstantVibrationalRelaxationRate::
    ConstantVibrationalRelaxationRate(
        const AnyMap& node, const UnitStack& rate_units){
    setParameters(node, rate_units);
}

ConstantVibrationalRelaxationRate::
    ConstantVibrationalRelaxationRate(double A)
    : VibrationalRelaxationRate(A, 0.0, 0.0, 0.0, 0.0,
        2.0 / 3.0, 0.0, 1.0)
{
}

void ConstantVibrationalRelaxationRate::setParameters(
    const AnyMap& node, const UnitStack& rate_units){
    VibrationalRelaxationRate::setParameters(node, rate_units);

    const auto& rateMap =
        node["rate-constant"].as<AnyMap>();

    requireKeys(
        rateMap, type(), WhereSetParameters,
        {"A"});

    forbidKeys(
        rateMap, type(), WhereSetParameters,
        {"b", "C23", "C0", "C13", "Cm", "m",
        "Cn", "n", "Ea", "K", "a",
        "reference-pressure"});

    m_A = node.units().convertRateCoeff(
        rateMap["A"], conversionUnits());

    m_b = 0.0;
    m_C0 = 0.0;
    m_C13 = 0.0;
    m_Cm = 0.0;
    m_m = 2.0 / 3.0;
    m_Cn = 0.0;
    m_n = 1.0;

    m_valid = true;
}

void ConstantVibrationalRelaxationRate::getParameters(
    AnyMap& node) const
{
    if (!valid()) {
        return;
    }

    VibrationalRelaxationRate::getParameters(node);

    AnyMap rateNode;

    getPreExponentialFactor(rateNode);

    rateNode.setFlowStyle();
    node["rate-constant"] = std::move(rateNode);
}

} // namespace Cantera