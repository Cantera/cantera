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
        const AnyMap& node, const UnitStack& rate_units)
{
    setParameters(node, rate_units);
}

StarikovskiyVibrationalRelaxationRate::
    StarikovskiyVibrationalRelaxationRate(
        double A, double b, double C0, double C13, double Cm, double m,
        double Cn, double n)
    : VibrationalRelaxationRate(A, b, C0, C13, Cm, m, Cn, n)
{
    if (m <= 0.0 || n <= 0.0) {
        throw CanteraError(
            "StarikovskiyVibrationalRelaxationRate",
            "The exponents 'm' and 'n' must be positive.");
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
        {"C23", "Ea", "a", "reference-pressure"});

    const double m = rateMap.getDouble("m", 1.0);

    const double n = rateMap.getDouble("n", 1.0);

    if (m <= 0.0 || n <= 0.0) {
        throw InputFileError(
            WhereSetParameters,
            node,
            "The Starikovskiy exponents 'm' and 'n' "
            "must be positive.");
    }

    m_A = node.units().convertRateCoeff(
        rateMap["A"], conversionUnits());

    m_b = rateMap.getDouble("b", 0.0);

    m_C0 = rateMap.getDouble("C0", 0.0);
    m_C13 = rateMap.getDouble("C13", 0.0);
    m_Cm = rateMap.getDouble("Cm", 0.0);
    m_m = m;

    m_Cn = rateMap.getDouble("Cn", 0.0);
    m_n = n;

    m_valid = true;
}

void StarikovskiyVibrationalRelaxationRate::getRateParameters(AnyMap& rateNode) const
{
    getPreExponentialFactor(rateNode);

    rateNode["b"] = m_b;

    rateNode["C0"] = m_C0;
    rateNode["C13"] = m_C13;

    rateNode["Cm"] = m_Cm;
    rateNode["m"] = m_m;

    rateNode["Cn"] = m_Cn;
    rateNode["n"] = m_n;
}

} // namespace Cantera
