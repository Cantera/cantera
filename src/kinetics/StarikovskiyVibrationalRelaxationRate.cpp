//! @file StarikovskiyVibrationalRelaxationRate.cpp
//! @since New in %Cantera 4.0

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/StarikovskiyVibrationalRelaxationRate.h"

namespace Cantera
{

namespace
{

const string WhereSetParameters =
    "StarikovskiyVibrationalRelaxationRate::setParameters";

} // namespace

StarikovskiyVibrationalRelaxationRate::
    StarikovskiyVibrationalRelaxationRate(
        const AnyMap& node, const UnitStack& rate_units){
    setParameters(node, rate_units);
}

StarikovskiyVibrationalRelaxationRate::
    StarikovskiyVibrationalRelaxationRate(
        double A,
        double n,
        double K,
        double B,
        double C,
        double m,
        double D,
        double z)
    : VibrationalRelaxationRate(
        A,
        n,  // internal b
        K,  // internal B
        B,  // internal C
        C,  // internal D
        m,
        D,  // internal E
        z)
{
    if (m <= 0.0 || z <= 0.0) {
        throw CanteraError(
            "StarikovskiyVibrationalRelaxationRate",
            "The exponents 'm' and 'z' must be positive.");
    }
}

void StarikovskiyVibrationalRelaxationRate::setParameters(
    const AnyMap& node, const UnitStack& rate_units)
{
    VibrationalRelaxationRate::setParameters(node, rate_units);

    const auto& rateMap = node["rate-constant"].as<AnyMap>();

    requireKeys(rateMap, type(), WhereSetParameters, {"A"});

    forbidKeys(
        rateMap, type(), WhereSetParameters,
        {"b", "Ea", "a", "reference-pressure"});

    const double m = rateMap.getDouble("m", 1.0);

    const double z = rateMap.getDouble("z", 1.0);

    if (m <= 0.0 || z <= 0.0) {
        throw InputFileError(
            WhereSetParameters,
            node,
            "The Starikovskiy exponents 'm' and 'z' "
            "must be positive.");
    }

    m_A = node.units().convertRateCoeff(
        rateMap["A"], conversionUnits());

    m_b = rateMap.getDouble("n", 0.0);

    m_B = rateMap.getDouble("K", 0.0);
    m_C = rateMap.getDouble("B", 0.0);
    m_D = rateMap.getDouble("C", 0.0);
    m_m = m;

    m_E = rateMap.getDouble("D", 0.0);
    m_z = z;

    m_valid = true;
}


void StarikovskiyVibrationalRelaxationRate::getParameters(
    AnyMap& node) const
{
    if (!valid()) {
        return;
    }

    VibrationalRelaxationRate::getParameters(node);

    AnyMap rateNode;

    getPreExponentialFactor(rateNode);

    rateNode["n"] = m_b;

    rateNode["K"] = m_B;
    rateNode["B"] = m_C;

    rateNode["C"] = m_D;
    rateNode["m"] = m_m;

    rateNode["D"] = m_E;
    rateNode["z"] = m_z;

    rateNode.setFlowStyle();
    node["rate-constant"] = std::move(rateNode);
}

} // namespace Cantera