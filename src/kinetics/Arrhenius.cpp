//! @file Arrhenius.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/Arrhenius.h"
#include "cantera/base/AnyMap.h"
#include "cantera/base/global.h"

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
            if (rate_map[m_A_str].is<std::string>()) {
                throw InputFileError("Arrhenius::setRateParameters", rate_map,
                    "Specification of units is not supported for pre-exponential "
                    "factor when\ncreating a standalone 'ReactionRate' object.");
            }
            m_A = rate_map[m_A_str].asDouble();
        } else {
            m_A = units.convert(rate_map[m_A_str], m_rate_units);
        }
        m_b = rate_map[m_b_str].asDouble();
        if (rate_map.hasKey(m_Ea_str)) {
            m_Ea_R = units.convertActivationEnergy(rate_map[m_Ea_str], "K");
        } else {
            m_Ea_R = NAN;
        }
    } else {
        auto& rate_vec = rate.asVector<AnyValue>(2, 3);
        m_A = units.convert(rate_vec[0], m_rate_units);
        m_b = rate_vec[1].asDouble();
        if (rate_vec.size() == 3) {
            m_Ea_R = units.convertActivationEnergy(rate_vec[2], "K");
        } else {
            m_Ea_R = NAN;
        }
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
        node[m_A_str].setQuantity(m_A, rate_units);
    } else {
        node[m_A_str] = m_A;
        // This can't be converted to a different unit system because the dimensions of
        // the rate constant were not set. Can occur if the reaction was created outside
        // the context of a Kinetics object and never added to a Kinetics object.
        node["__unconvertible__"] = true;
    }
    node[m_b_str] = m_b;
    node[m_Ea_str].setQuantity(m_Ea_R, "K", true);

    node.setFlowStyle();
}

void ArrheniusBase::check(const std::string& equation, const AnyMap& node)
{
    if (!m_negativeA_ok && m_A < 0) {
        if (equation == "") {
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

TwoTempPlasmaRate::TwoTempPlasmaRate()
    : ArrheniusBase()
    , m_EE_R(0.0)
{
    m_Ea_R = 0.0;
    m_Ea_str = "Ea-gas";
}

TwoTempPlasmaRate::TwoTempPlasmaRate(double A, double b, double Ea, double EE)
    : ArrheniusBase(A, b, Ea)
    , m_EE_R(EE  / GasConstant)
{
    m_Ea_str = "Ea-gas";
}

void TwoTempPlasmaRate::setRateParameters(
    const AnyValue& rate, const UnitSystem& units, const UnitStack& rate_units)
{
    if (rate.empty()) {
        m_A = NAN;
        m_b = NAN;
        m_Ea_R = NAN;
        m_EE_R = NAN;
        m_order = NAN;
        m_rate_units = Units(0.);
        return;
    }

    if (rate.is<AnyMap>()) {
        ArrheniusBase::setRateParameters(rate, units, rate_units);
        auto& rate_map = rate.as<AnyMap>();
        // check Ea-gas value. Set to zero in the case of absent.
        if (isnan(m_Ea_R)) {
            m_Ea_R = 0.0;
        }
        // Get Ea-electron value. If not set it to zero.
        if (rate_map.hasKey("Ea-electron")) {
            m_EE_R = units.convertActivationEnergy(rate_map["Ea-electron"], "K");
        } else {
            m_EE_R = 0.0;
        }
    } else {
        setRateUnits(rate_units);
        auto& rate_vec = rate.asVector<AnyValue>(2,4);
        m_A = units.convert(rate_vec[0], m_rate_units);
        m_b = rate_vec[1].asDouble();
        if (rate_vec.size() == 4) {
            m_EE_R = units.convertActivationEnergy(rate_vec[3], "K");
            m_Ea_R = units.convertActivationEnergy(rate_vec[2], "K");
        } else if (rate_vec.size() == 3) {
            m_EE_R = 0.0;
            m_Ea_R = units.convertActivationEnergy(rate_vec[2], "K");
        } else {
            m_Ea_R = 0.0;
            m_EE_R = 0.0;
        }
    }
}

void TwoTempPlasmaRate::setParameters(const AnyMap& node, const UnitStack& rate_units)
{
    m_negativeA_ok = node.getBool("negative-A", false);
    if (!node.hasKey("rate-constant")) {
        TwoTempPlasmaRate::setRateParameters(AnyValue(), node.units(), rate_units);
        return;
    }
    TwoTempPlasmaRate::setRateParameters(node["rate-constant"], node.units(), rate_units);
}

void TwoTempPlasmaRate::getParameters(AnyMap& rateNode) const
{
    if (m_negativeA_ok) {
        rateNode["negative-A"] = true;
    }
    AnyMap node;
    ArrheniusBase::getParameters(node);
    if (!node.empty()) {
        // object is configured
        node["Ea-electron"].setQuantity(m_EE_R, "K", true);
        rateNode["rate-constant"] = std::move(node);
    }
    rateNode["type"] = type();
}

BlowersMaselRate::BlowersMaselRate()
    : m_w_R(NAN)
{
    m_Ea_str = "Ea0";
}

BlowersMaselRate::BlowersMaselRate(double A, double b, double Ea0, double w)
    : ArrheniusBase(A, b, Ea0)
    , m_w_R(w / GasConstant)
{
    m_Ea_str = "Ea0";
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

double BlowersMaselRate::ddTScaledFromStruct(const BlowersMaselData& shared_data) const
{
    warn_user("BlowersMaselRate::ddTScaledFromStruct",
        "Temperature derivative does not consider changes of reaction enthalpy.");
    double deltaH_R;
    if (shared_data.ready) {
        deltaH_R = shared_data.dH[m_rate_index] / GasConstant;
    } else {
        deltaH_R = shared_data.dH[0] / GasConstant;
    }
    double Ea_R = activationEnergy_R(deltaH_R);
    return m_A * std::exp(m_b * shared_data.logT - Ea_R * shared_data.recipT);
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
        node["w"].setQuantity(m_w_R, "K", true);
        rateNode["rate-constant"] = std::move(node);
    }
    rateNode["type"] = type();
}

}
