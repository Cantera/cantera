//! @file Arrhenius.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/Arrhenius.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/kinetics/Reaction.h"

namespace Cantera
{

ArrheniusBase::ArrheniusBase()
    : m_negativeA_ok(false)
    , m_A(NAN)
    , m_b(NAN)
    , m_Ea_R(0.)
    , m_E4_R(0.)
    , m_logA(NAN)
    , m_order(NAN)
    , m_rate_units(Units(0.))
{
}

ArrheniusBase::ArrheniusBase(double A, double b, double Ea)
    : m_negativeA_ok(false)
    , m_A(A)
    , m_b(b)
    , m_Ea_R(Ea / GasConstant)
    , m_E4_R(0.)
    , m_order(NAN)
    , m_rate_units(Units(0.))
{
    if (m_A > 0.0) {
        m_logA = std::log(m_A);
    }
}

void ArrheniusBase::setRateParameters(
    const AnyValue& rate, const UnitSystem& units, const UnitStack& rate_units)
{
    m_Ea_R = 0.; // assume zero if not provided
    m_E4_R = 0.; // assume zero if not provided
    if (rate.empty()) {
        m_A = NAN;
        m_b = NAN;
        m_logA = NAN;
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
        }
        if (rate_map.hasKey(m_E4_str)) {
            m_E4_R = units.convertActivationEnergy(rate_map[m_E4_str], "K");
        }
    } else {
        auto& rate_vec = rate.asVector<AnyValue>(2, 4);
        m_A = units.convert(rate_vec[0], m_rate_units);
        m_b = rate_vec[1].asDouble();
        if (rate_vec.size() > 2) {
            m_Ea_R = units.convertActivationEnergy(rate_vec[2], "K");
        }
        if (rate_vec.size() > 3) {
            m_E4_R = units.convertActivationEnergy(rate_vec[3], "K");
        }
    }
    if (m_A > 0.0) {
        m_logA = std::log(m_A);
    }
}

void ArrheniusBase::getRateParameters(AnyMap& node) const
{
    if (std::isnan(m_A)) {
        // Return empty/unmodified AnyMap
        return;
    } else if (m_rate_units.factor() != 0.0) {
        node[m_A_str].setQuantity(m_A, m_rate_units);
    } else {
        node[m_A_str] = m_A;
        // This can't be converted to a different unit system because the dimensions of
        // the rate constant were not set. Can occur if the reaction was created outside
        // the context of a Kinetics object and never added to a Kinetics object.
        node["__unconvertible__"] = true;
    }
    node[m_b_str] = m_b;
    node[m_Ea_str].setQuantity(m_Ea_R, "K", true);
    if (m_E4_str != "") {
        node[m_E4_str].setQuantity(m_E4_R, "K", true);
    }
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

TwoTempPlasma::TwoTempPlasma()
    : ArrheniusBase()
{
    m_Ea_str = "Ea-gas";
    m_E4_str = "Ea-electron";
}

TwoTempPlasma::TwoTempPlasma(double A, double b, double Ea, double EE)
    : ArrheniusBase(A, b, Ea)
{
    m_Ea_str = "Ea-gas";
    m_E4_str = "Ea-electron";
    m_E4_R = EE / GasConstant;
}

double TwoTempPlasma::ddTScaledFromStruct(const TwoTempPlasmaData& shared_data) const
{
    warn_user("TwoTempPlasma::ddTScaledFromStruct",
        "Temperature derivative does not consider changes of electron temperature.");
    return (m_Ea_R - m_E4_R) * shared_data.recipT * shared_data.recipT;
}

void TwoTempPlasma::setContext(const Reaction& rxn, const Kinetics& kin)
{
    // TwoTempPlasmaReaction is for a non-equilibrium plasma, and the reverse rate
    // cannot be calculated from the conventional thermochemistry.
    // @todo implement the reversible rate for non-equilibrium plasma
    if (rxn.reversible) {
        throw InputFileError("TwoTempPlasma::setContext", rxn.input,
            "TwoTempPlasmaRate does not support reversible reactions");
    }
}

BlowersMasel::BlowersMasel()
    : m_deltaH_R(0.)
{
    m_Ea_str = "Ea0";
    m_E4_str = "w";
}

BlowersMasel::BlowersMasel(double A, double b, double Ea0, double w)
    : ArrheniusBase(A, b, Ea0)
    , m_deltaH_R(0.)
{
    m_Ea_str = "Ea0";
    m_E4_str = "w";
    m_E4_R = w / GasConstant;
}

double BlowersMasel::ddTScaledFromStruct(const BlowersMaselData& shared_data) const
{
    warn_user("BlowersMasel::ddTScaledFromStruct",
        "Temperature derivative does not consider changes of reaction enthalpy.");
    double Ea_R = effectiveActivationEnergy_R(m_deltaH_R);
    return (Ea_R * shared_data.recipT + m_b) * shared_data.recipT;
}

void BlowersMasel::setContext(const Reaction& rxn, const Kinetics& kin)
{
    m_stoich_coeffs.clear();
    for (const auto& sp : rxn.reactants) {
        m_stoich_coeffs.emplace_back(kin.kineticsSpeciesIndex(sp.first), -sp.second);
    }
    for (const auto& sp : rxn.products) {
        m_stoich_coeffs.emplace_back(kin.kineticsSpeciesIndex(sp.first), sp.second);
    }
}

}
