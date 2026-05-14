//! @file DetailedVVVTRate.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/DetailedVVVTRate.h"
#include "cantera/kinetics/Reaction.h"
#include "cantera/thermo/ThermoPhase.h"

#include <cmath>
#include <utility>

namespace Cantera
{

bool DetailedVibData::update(const ThermoPhase& phase, const Kinetics& kin)
{
    double T = phase.temperature();

    if (T == temperature) {
        return false;
    }

    ReactionData::update(T);
    return true;
}


DetailedVVVTRate::DetailedVVVTRate()
{
}


DetailedVVVTRate::DetailedVVVTRate(double A, double B, double C, double D,
                                   double b, double scaling)
    : ArrheniusBase(A, b, 0.0)
    , m_B(B)
    , m_C(C)
    , m_D(D)
    , m_scaling(scaling)
{
}


DetailedVVVTRate::DetailedVVVTRate(const AnyMap& node,
                                   const UnitStack& rate_units)
    : DetailedVVVTRate()
{
    setParameters(node, rate_units);
}


void DetailedVVVTRate::setParameters(const AnyMap& node,
                                     const UnitStack& rate_units)
{
    // Let ArrheniusBase handle:
    //
    // - storage of the input AnyMap
    // - negative-A handling
    // - conversion of A using the normal unit system
    // - reading of b
    //
    // The Ea parameter, if present, is parsed by ArrheniusBase but is not used
    // by DetailedVVVTRate. The rate expression implemented here uses B, C, D
    // instead of an Arrhenius activation energy.
    ArrheniusBase::setParameters(node, rate_units);

    if (!node.hasKey("rate-constant")) {
        return;
    }

    const auto& rate = node["rate-constant"];

    if (!rate.is<AnyMap>()) {
        throw InputFileError("DetailedVVVTRate::setParameters", node,
            "The 'rate-constant' field for a 'detailed-vv-vt' reaction "
            "must be a mapping containing A, b, B, C, D, and optionally "
            "scaling.");
    }

    const auto& rate_map = rate.as<AnyMap>();

    // B is dimensionless.
    if (rate_map.hasKey(m_B_str)) {
        m_B = rate_map[m_B_str].asDouble();
    }

    // C is interpreted as K^(1/3), assuming T is in K.
    // It is intentionally read as a raw number and is not converted by the
    // Cantera unit system.
    if (rate_map.hasKey(m_C_str)) {
        m_C = rate_map[m_C_str].asDouble();
    }

    // D is interpreted as K^(2/3), assuming T is in K.
    // It is intentionally read as a raw number and is not converted by the
    // Cantera unit system.
    if (rate_map.hasKey(m_D_str)) {
        m_D = rate_map[m_D_str].asDouble();
    }

    // scaling is dimensionless.
    if (rate_map.hasKey(m_scaling_str)) {
        m_scaling = rate_map[m_scaling_str].asDouble();
    }
}


void DetailedVVVTRate::getParameters(AnyMap& node) const
{
    if (!valid()) {
        return;
    }

    if (allowNegativePreExponentialFactor()) {
        node["negative-A"] = true;
    }

    AnyMap rateNode;

    // Store A using the same convention as ArrheniusBase::getRateParameters.
    // When the reaction has been associated with a Kinetics object, Cantera
    // knows the conversion units for the leading rate coefficient.
    if (conversionUnits().factor() != 0.0) {
        rateNode[m_A_str].setQuantity(m_A, conversionUnits());
    } else {
        // This case can occur when the rate was created outside the context of
        // a Kinetics object and therefore the reaction-order-dependent units
        // are not known.
        rateNode[m_A_str] = m_A;
        rateNode["__unconvertible__"] = true;
    }

    // b is dimensionless.
    rateNode[m_b_str] = m_b;

    // Custom DetailedVVVTRate parameters.
    //
    // B is dimensionless.
    // C is interpreted as K^(1/3), assuming T is in K.
    // D is interpreted as K^(2/3), assuming T is in K.
    // scaling is dimensionless.
    //
    // These values are deliberately serialized as raw floating-point numbers.
    rateNode[m_B_str] = m_B;
    rateNode[m_C_str] = m_C;
    rateNode[m_D_str] = m_D;
    rateNode[m_scaling_str] = m_scaling;

    rateNode.setFlowStyle();

    node["rate-constant"] = std::move(rateNode);
}


double DetailedVVVTRate::ddTScaledFromStruct(
    const DetailedVibData& shared_data) const
{
    const double invT = shared_data.recipT;
    const double invT13 = std::cbrt(invT);

    // For:
    //
    //   k = scaling * A * exp(b*log(T) + B + C*T^(-1/3) + D*T^(-2/3))
    //
    // the logarithmic derivative is:
    //
    //   (1/k) * dk/dT
    //     = b/T
    //       - (C/3)   * T^(-4/3)
    //       - (2D/3)  * T^(-5/3)
    //
    // Using invT = 1/T and invT13 = T^(-1/3):
    //
    //   T^(-4/3) = invT13 * invT
    //   T^(-5/3) = invT13 * invT13 * invT
    return m_b * invT
           - (m_C / 3.0) * invT13 * invT
           - (2.0 * m_D / 3.0) * invT13 * invT13 * invT;
}


void DetailedVVVTRate::setContext(const Reaction& rxn, const Kinetics& kin)
{
    // DetailedVVVTRate is intended for non-equilibrium plasma kinetics.
    // The reverse rate cannot be calculated from conventional thermochemistry
    // without an additional non-equilibrium model.
    if (rxn.reversible) {
        throw InputFileError("DetailedVVVTRate::setContext", rxn.input,
            "DetailedVVVTRate does not support reversible reactions.");
    }
}

}