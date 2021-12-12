//! @file Arrhenius.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/Arrhenius.h"
#include "cantera/base/AnyMap.h"

namespace Cantera
{

ArrheniusBase::ArrheniusBase()
    : m_negativeA_ok(false)
    , m_A(NAN)
    , m_b(NAN)
    , m_Ea_R(NAN)
    , m_order(NAN)
    , m_rate_units(Units(0.))
{
}

ArrheniusBase::ArrheniusBase(double A, double b, double Ea)
    : m_negativeA_ok(false)
    , m_A(A)
    , m_b(b)
    , m_Ea_R(Ea / GasConstant)
    , m_order(NAN)
    , m_rate_units(Units(0.))
{
}

void ArrheniusBase::setRateParameters(
    const AnyValue& rate, const UnitSystem& units, const UnitStack& rate_units)
{
    if (rate.empty()) {
        m_A = NAN;
        m_b = NAN;
        m_Ea_R = NAN;
        m_order = NAN;
        m_rate_units = Units(0.);
        return;
    }

    setRateUnits(rate_units);
    if (rate.is<AnyMap>()) {

        auto& rate_map = rate.as<AnyMap>();
        if (m_rate_units.factor() == 0) {
            // A zero rate units factor is used as a sentinel to detect
            // stand-alone reaction rate objects
            if (rate_map["A"].is<std::string>()) {
                throw InputFileError("Arrhenius::setRateParameters", rate_map,
                    "Specification of units is not supported for pre-exponential "
                    "factor when\ncreating a standalone 'ReactionRate' object.");
            }
            m_A = rate_map["A"].asDouble();
        } else {
            m_A = units.convert(rate_map["A"], m_rate_units);
        }
        m_b = rate_map["b"].asDouble();
        if (rate_map.hasKey("Ea")) {
            m_Ea_R = units.convertActivationEnergy(rate_map["Ea"], "K");
        } else {
            m_Ea_R = NAN;
        }
    } else {
        auto& rate_vec = rate.asVector<AnyValue>(3);
        m_A = units.convert(rate_vec[0], m_rate_units);
        m_b = rate_vec[1].asDouble();
        m_Ea_R = units.convertActivationEnergy(rate_vec[2], "K");
    }
}

void ArrheniusBase::getParameters(AnyMap& node) const
{
    getParameters(node, m_rate_units);
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

void ArrheniusBase::check(const std::string& equation, const AnyMap& node)
{
    if (!m_negativeA_ok && m_A < 0) {
        if (equation == "") {
            //
            throw CanteraError("ArrheniusBase::check",
                "Detected negative pre-exponential factor (A={}).\n"
                "Enable 'allowNegativePreExponentialFactor' to suppress "
                "this message.", m_A);
        }
        throw InputFileError("ArrheniusBase::check", node,
            "Undeclared negative pre-exponential factor found in reaction '{}'",
            equation);
    }
}

void ArrheniusRate::setParameters(const AnyMap& node, const UnitStack& rate_units)
{
    m_negativeA_ok = node.getBool("negative-A", false);
    if (!node.hasKey("rate-constant")) {
        ArrheniusBase::setRateParameters(AnyValue(), node.units(), rate_units);
        return;
    }

    ArrheniusBase::setRateParameters(node["rate-constant"], node.units(), rate_units);
}

void ArrheniusRate::getParameters(AnyMap& rateNode) const
{
    if (m_negativeA_ok) {
        rateNode["negative-A"] = true;
    }
    AnyMap node;
    ArrheniusBase::getParameters(node);
    if (!node.empty()) {
        // Arrhenius object is configured
        rateNode["rate-constant"] = std::move(node);
    }
}

BlowersMaselRate::BlowersMaselRate()
    : m_w_R(NAN)
{
}

BlowersMaselRate::BlowersMaselRate(double A, double b, double Ea0, double w)
    : ArrheniusBase(A, b, Ea0)
    , m_w_R(w / GasConstant)
{
}

void BlowersMaselRate::setRateParameters(
    const AnyValue& rate, const UnitSystem& units, const UnitStack& rate_units)
{
    if (rate.empty()) {
        m_A = NAN;
        m_b = NAN;
        m_Ea_R = NAN;
        m_w_R = NAN;
        m_order = NAN;
        m_rate_units = Units(0.);
        return;
    }

    if (rate.is<AnyMap>()) {
        ArrheniusBase::setRateParameters(rate, units, rate_units);
        auto& rate_map = rate.as<AnyMap>();
        m_Ea_R = units.convertActivationEnergy(rate_map["Ea0"], "K");
        m_w_R = units.convertActivationEnergy(rate_map["w"], "K");
    } else {
        setRateUnits(rate_units);
        auto& rate_vec = rate.asVector<AnyValue>(4);
        m_A = units.convert(rate_vec[0], m_rate_units);
        m_b = rate_vec[1].asDouble();
        m_Ea_R = units.convertActivationEnergy(rate_vec[2], "K");
        m_w_R = units.convertActivationEnergy(rate_vec[3], "K");
    }
}

void BlowersMaselRate::setParameters(const AnyMap& node, const UnitStack& rate_units)
{
    m_negativeA_ok = node.getBool("negative-A", false);
    if (!node.hasKey("rate-constant")) {
        BlowersMaselRate::setRateParameters(AnyValue(), node.units(), rate_units);
        return;
    }
    BlowersMaselRate::setRateParameters(node["rate-constant"], node.units(), rate_units);
}

void BlowersMaselRate::getParameters(AnyMap& rateNode) const
{
    if (m_negativeA_ok) {
        rateNode["negative-A"] = true;
    }
    AnyMap node;
    ArrheniusBase::getParameters(node);
    if (!node.empty()) {
        // object is configured
        node.erase("Ea");
        node["Ea0"].setQuantity(m_Ea_R, "K", true);
        node["w"].setQuantity(m_w_R, "K", true);
        rateNode["rate-constant"] = std::move(node);
    }
    rateNode["type"] = type();
}

}
