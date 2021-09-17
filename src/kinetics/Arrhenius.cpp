//! @file Arrhenius.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/Arrhenius.h"
#include "cantera/base/AnyMap.h"

namespace Cantera
{

ArrheniusBase::ArrheniusBase()
    : allow_negative_pre_exponential_factor(false)
    , rate_index(npos)
    , m_A(NAN)
    , m_b(NAN)
    , m_Ea_R(NAN)
{
}

ArrheniusBase::ArrheniusBase(double A, double b, double Ea)
    : allow_negative_pre_exponential_factor(false)
    , rate_index(npos)
    , m_A(A)
    , m_b(b)
    , m_Ea_R(Ea / GasConstant)
{
}

void ArrheniusBase::setParameters(
    const AnyValue& rate, const UnitSystem& units, const Units& rate_units)
{
    if (rate.empty()) {
        m_A = NAN;
        m_b = NAN;
        m_Ea_R = NAN;
    } else if (rate.is<AnyMap>()) {

        auto& rate_map = rate.as<AnyMap>();
        if (rate_units.factor() == 0) {
            // A zero rate units factor is used as a sentinel to detect
            // stand-alone reaction rate objects
            if (rate_map["A"].is<std::string>()) {
                throw InputFileError("Arrhenius::setParameters", rate_map,
                    "Specification of units is not supported for pre-exponential "
                    "factor when\ncreating a standalone 'ReactionRate' object.");
            }
            m_A = rate_map["A"].asDouble();
        } else {
            m_A = units.convert(rate_map["A"], rate_units);
        }
        m_b = rate_map["b"].asDouble();
        if (rate_map.hasKey("Ea")) {
            m_Ea_R = units.convertActivationEnergy(rate_map["Ea"], "K");
        } else {
            m_Ea_R = NAN;
        }
    } else {
        auto& rate_vec = rate.asVector<AnyValue>(3);
        m_A = units.convert(rate_vec[0], rate_units);
        m_b = rate_vec[1].asDouble();
        m_Ea_R = units.convertActivationEnergy(rate_vec[2], "K");
    }
}

void ArrheniusBase::getParameters(AnyMap& node, const Units& rate_units) const
{
    if (std::isnan(m_A)) {
        // Return empty/unmodified AnyMap
        return;
    } else if (rate_units.factor() != 0.0) {
        node["A"].setQuantity(m_A, rate_units);
    } else {
        node["A"] = m_A;
        // This can't be converted to a different unit system because the dimensions of
        // the rate constant were not set. Can occur if the reaction was created outside
        // the context of a Kinetics object and never added to a Kinetics object.
        node["__unconvertible__"] = true;
    }

    node["b"] = m_b;
    node["Ea"].setQuantity(m_Ea_R, "K", true);
    node.setFlowStyle();
}

void ArrheniusBase::validate(const std::string& equation)
{
    if (!allow_negative_pre_exponential_factor && m_A < 0) {
        throw CanteraError("ArrheniusBase::validate",
            "Undeclared negative pre-exponential factor found in reaction '{}'",
            equation);
    }
}

}
