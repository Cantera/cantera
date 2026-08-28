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
        double A, double b, double C0, double C13, double C23)
    : VibrationalRelaxationRate(
        A, b, C0, C13, C23, 2.0 / 3.0, 0.0, 1.0)
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
        {"n", "m", "Cn", "Cm", "Ea",
        "a", "reference-pressure"});

    m_A = node.units().convertRateCoeff(
        rateMap["A"], conversionUnits());

    m_b = rateMap["b"].asDouble();

    m_C0 = rateMap.getDouble("C0", 0.0);
    m_C13 = rateMap.getDouble("C13", 0.0);
    m_Cm = rateMap.getDouble("C23", 0.0);

    m_m = 2.0 / 3.0;

    m_Cn = 0.0;
    m_n = 1.0;

    m_valid = true;
}

void MultiStateResolvedVibrationalRelaxationRate::getRateParameters(
    AnyMap& rateNode) const
{
    getPreExponentialFactor(rateNode);

    rateNode["b"] = m_b;
    rateNode["C0"] = m_C0;
    rateNode["C13"] = m_C13;
    rateNode["C23"] = m_Cm;
}

} // namespace Cantera
