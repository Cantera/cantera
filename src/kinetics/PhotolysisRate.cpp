//! @file PhotolysisRate.cpp
//! @since  New in Cantera 3.0.

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/PhotolysisRate.h"
#include "cantera/thermo/ThermoPhase.h"

namespace Cantera
{

PhotolysisRate::PhotolysisRate(double J)
    : m_J(J)
{
    m_valid = true;
    m_rate_evaluated = true;
}

PhotolysisRate::PhotolysisRate(double l, double m, double n)
    : m_l(l)
    , m_m(m)
    , m_n(n)
{
    m_valid = true;
}

PhotolysisRate::PhotolysisRate(const AnyValue& rate, const UnitSystem& units,
                             const UnitStack& rate_units)
{
    setRateUnits(rate_units);
    setRateParameters(rate, units, rate_units);
}

PhotolysisRate::PhotolysisRate(const AnyMap& node, const UnitStack& rate_units)
{
    setParameters(node, rate_units);
}

void PhotolysisRate::setRateParameters(
    const AnyValue& rate, const UnitSystem& units, const UnitStack& rate_units)
{
    if (rate.empty()) {
        m_l = NAN;
        m_m = NAN;
        m_n = NAN;
        setRateUnits(Units(0.));
        return;
    }

    if (rate.is<AnyMap>()) {
        auto& rate_map = rate.as<AnyMap>();
        // try to use constant J first then decide to use parameterizations
        if (rate_map.hasKey(m_J_str)) {
            m_J = units.convertRateCoeff(rate_map[m_J_str], conversionUnits());
            m_rate_evaluated = true;
        } else {
            // The rate has units which are assigned to l as it is a scalar
            m_l = units.convertRateCoeff(rate_map[m_l_str], conversionUnits());
            m_m = rate_map[m_m_str].asDouble();
            m_n = rate_map[m_n_str].asDouble();
            m_rate_evaluated = false;
        }
    } else {
        auto& rate_vec = rate.asVector<AnyValue>(1, 3);
        if (rate_vec.size() < 2) {
            m_J = units.convertRateCoeff(rate_vec[0], conversionUnits());
            m_rate_evaluated = true;
        } else {
            m_l = units.convertRateCoeff(rate_vec[0], conversionUnits());
            m_m = rate_vec[1].asDouble();
            m_n = rate_vec[2].asDouble();
            m_rate_evaluated = false;
        }
    }
    m_valid = true;
}

void PhotolysisRate::getRateParameters(AnyMap& node) const
{
    if (!valid()) {
        // Return empty/unmodified AnyMap
        return;
    }
    // Use constant J value if it is not 0 else use parameters
    if (m_J != 0) {
        if (conversionUnits().factor() != 0.0) {
            node[m_J_str].setQuantity(m_J, conversionUnits());
        } else {
            node[m_J_str] = m_J;
            // This can't be converted to a different unit system because the dimensions of
            // the rate constant were not set. Can occur if the reaction was created outside
            // the context of a Kinetics object and never added to a Kinetics object.
            node["__unconvertible__"] = true;
        }
    } else {
        if (conversionUnits().factor() != 0.0) {
            node[m_l_str].setQuantity(m_l, conversionUnits());
        } else {
            node[m_l_str] = m_l;
        }
        node[m_m_str] = m_m;
        node[m_n_str] = m_n;
    }
    node.setFlowStyle();
}

void PhotolysisRate::setParameters(const AnyMap& node, const UnitStack& rate_units)
{
    ReactionRate::setParameters(node, rate_units);
    if (!node.hasKey("rate-constant")) {
        setRateParameters(AnyValue(), node.units(), rate_units);
        return;
    }
    setRateParameters(node["rate-constant"], node.units(), rate_units);
}

void PhotolysisRate::getParameters(AnyMap& node) const {
    AnyMap rateNode;
    getRateParameters(rateNode);
    if (!rateNode.empty()) {
        // RateType object is configured
        node["rate-constant"] = std::move(rateNode);
    }
}

void PhotolysisRate::check(const std::string& equation)
{
    //TODO: add check if necessary
}

void PhotolysisRate::validate(const std::string& equation, const Kinetics& kin)
{
    if (!valid()) {
        throw InputFileError("PhotolysisRate::validate", m_input,
            "Rate object for reaction '{}' is not configured.", equation);
    }
}

bool PhotolysisData::update(const ThermoPhase& phase, const Kinetics& kin)
{
    // set a zenith angle if there is one
    X = phase.zenithAngle();
    return true;
}

}
