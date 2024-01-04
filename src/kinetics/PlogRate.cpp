//! @file PlogRate.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/PlogRate.h"
#include "cantera/thermo/ThermoPhase.h"

namespace Cantera
{

void PlogData::update(double T)
{
    throw CanteraError("PlogData::update",
        "Missing state information: 'PlogData' requires pressure.");
}

bool PlogData::update(const ThermoPhase& phase, const Kinetics& kin)
{
    double T = phase.temperature();
    double P = phase.pressure();
    if (P != pressure || T != temperature) {
        update(T, P);
        return true;
    }
    return false;
}

void PlogData::perturbPressure(double deltaP)
{
    if (m_pressure_buf > 0.) {
        throw CanteraError("PlogData::perturbPressure",
            "Cannot apply another perturbation as state is already perturbed.");
    }
    m_pressure_buf = pressure;
    update(temperature, pressure * (1. + deltaP));
}

void PlogData::restore()
{
    ReactionData::restore();
    // only restore if there is a valid buffered value
    if (m_pressure_buf < 0.) {
        return;
    }
    update(temperature, m_pressure_buf);
    m_pressure_buf = -1.;
}

// Methods of class PlogRate

PlogRate::PlogRate(const std::multimap<double, ArrheniusRate>& rates)
{
    setRates(rates);
}

PlogRate::PlogRate(const AnyMap& node, const UnitStack& rate_units)
{
    setParameters(node, rate_units);
}

void PlogRate::setParameters(const AnyMap& node, const UnitStack& rate_units)
{
    ReactionRate::setParameters(node, rate_units);
    std::multimap<double, ArrheniusRate> multi_rates;
    if (node.hasKey("rate-constants")) {
        auto& rates = node["rate-constants"].asVector<AnyMap>();
        for (const auto& rate : rates) {
            multi_rates.insert({rate.convert("P", "Pa"),
                ArrheniusRate(AnyValue(rate), node.units(), rate_units)});
        }
    }
    setRates(multi_rates);
}

void PlogRate::getParameters(AnyMap& rateNode, const Units& rate_units) const
{
    vector<AnyMap> rateList;
    if (!valid()) {
        // object not fully set up
        return;
    }
    for (const auto& [pressure, rate] : getRates()) {
        AnyMap rateNode_;
        rateNode_["P"].setQuantity(pressure, "Pa");
        rate.getRateParameters(rateNode_);
        rateList.push_back(std::move(rateNode_));
    }
    rateNode["rate-constants"] = std::move(rateList);
}

void PlogRate::setRates(const std::multimap<double, ArrheniusRate>& rates)
{
    size_t j = 0;
    rates_.clear();
    pressures_.clear();
    m_valid = !rates.empty();
    rates_.reserve(rates.size());
    // Insert intermediate pressures
    for (const auto& [pressure, rate] : rates) {
        double logp = std::log(pressure);
        if (pressures_.empty() || pressures_.rbegin()->first != logp) {
            // starting a new group
            pressures_[logp] = {j, j+1};
        } else {
            // another rate expression at the same pressure
            pressures_[logp].second = j+1;
        }

        j++;
        rates_.push_back(rate);
    }
    if (!m_valid) {
        // ensure that reaction rate can be evaluated (but returns NaN)
        rates_.reserve(1);
        pressures_[std::log(OneBar)] = {0, 0};
        rates_.push_back(ArrheniusRate());
    }

    // Duplicate the first and last groups to handle P < P_0 and P > P_N
    pressures_.insert({-1000.0, pressures_.begin()->second});
    pressures_.insert({1000.0, pressures_.rbegin()->second});
}

void PlogRate::validate(const string& equation, const Kinetics& kin)
{
    if (!valid()) {
        throw InputFileError("PlogRate::validate", m_input,
            "Rate object for reaction '{}' is not configured.", equation);
    }

    fmt::memory_buffer err_reactions;
    double T[] = {300.0, 500.0, 1000.0, 2000.0, 5000.0, 10000.0};
    PlogData data;

    for (auto iter = ++pressures_.begin(); iter->first < 1000; iter++) {
        data.update(T[0], exp(iter->first)); // iter->first contains log(p)
        updateFromStruct(data);
        for (size_t i=0; i < 6; i++) {
            double k = 0;
            for (size_t p = ilow1_; p < ilow2_; p++) {
                k += rates_.at(p).evalRate(log(T[i]), 1.0 / T[i]);
            }
            if (!(k > 0)) {
                fmt_append(err_reactions,
                    "at P = {:.5g}, T = {:.1f}\n", std::exp(iter->first), T[i]);
            }
        }
    }
    if (err_reactions.size()) {
        throw InputFileError("PlogRate::validate", m_input,
            "\nInvalid rate coefficient for reaction '{}'\n{}\n"
            "To fix this error, remove this reaction or contact the author of the\n"
            "reaction/mechanism in question, because the rate expression is\n"
            "mathematically unsound at the temperatures and pressures noted above.\n",
            equation, to_string(err_reactions));
    }
}

std::multimap<double, ArrheniusRate> PlogRate::getRates() const
{
    std::multimap<double, ArrheniusRate> rateMap;
    // initial preincrement to skip rate for P --> 0
    for (auto iter = ++pressures_.begin();
            iter->first < 1000; // skip rates for (P --> infinity)
            ++iter) {
        for (size_t i = iter->second.first; i < iter->second.second; i++) {
            rateMap.insert({std::exp(iter->first), rates_[i]});
        }
    }
    return rateMap;
}

}
