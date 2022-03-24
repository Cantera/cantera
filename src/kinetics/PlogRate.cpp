//! @file PlogRate.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/PlogRate.h"
#include "cantera/kinetics/RxnRates.h"
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

PlogRate::PlogRate()
    : logP_(-1000)
    , logP1_(1000)
    , logP2_(-1000)
    , rDeltaP_(-1.0)
{
}

PlogRate::PlogRate(const std::multimap<double, ArrheniusRate>& rates)
    : PlogRate()
{
    setRates(rates);
}

PlogRate::PlogRate(const std::multimap<double, Arrhenius2>& rates)
    : PlogRate()
{
    setup(rates);
}

void PlogRate::setParameters(const AnyMap& node, const UnitStack& units)
{
    std::multimap<double, ArrheniusRate> multi_rates;
    if (!node.hasKey("rate-constants")) {
         // ensure that reaction rate can be evaluated (but returns NaN)
        multi_rates.insert({1.e-7, ArrheniusRate(NAN, NAN, NAN)});
        multi_rates.insert({1.e14, ArrheniusRate(NAN, NAN, NAN)});
    } else {
        auto& rates = node["rate-constants"].asVector<AnyMap>();
        for (const auto& rate : rates) {
            multi_rates.insert({rate.convert("P", "Pa"),
                ArrheniusRate(AnyValue(rate), node.units(), units)});
        }
    }
    setRates(multi_rates);
}

void PlogRate::getParameters(AnyMap& rateNode, const Units& rate_units) const
{
    std::vector<AnyMap> rateList;
    rateNode["type"] = type();
    if (!rates_.size()
        || (rates_.size() > 1 && std::isnan(rates_[1].preExponentialFactor())))
    {
        // object not fully set up
        return;
    }
    for (const auto& r : getRates()) {
        AnyMap rateNode_;
        rateNode_["P"].setQuantity(r.first, "Pa");
        if (rate_units.factor() == 0) {
            r.second.getRateParameters(rateNode_);
        } else {
            // legacy rate
            Arrhenius2(r.second).getParameters(rateNode_, rate_units);
        }
        rateList.push_back(std::move(rateNode_));
    }
    rateNode["rate-constants"] = std::move(rateList);
}

void PlogRate::setup(const std::multimap<double, Arrhenius2>& rates)
{
    std::multimap<double, ArrheniusRate> rates2;
    for (const auto& item : rates) {
        rates2.emplace(item.first, item.second);
    }
    setRates(rates2);
}

void PlogRate::setRates(const std::multimap<double, ArrheniusRate>& rates)
{
    size_t j = 0;
    rates_.clear();
    pressures_.clear();
    rates_.reserve(rates.size());
    // Insert intermediate pressures
    for (const auto& rate : rates) {
        double logp = std::log(rate.first);
        if (pressures_.empty() || pressures_.rbegin()->first != logp) {
            // starting a new group
            pressures_[logp] = {j, j+1};
        } else {
            // another rate expression at the same pressure
            pressures_[logp].second = j+1;
        }

        j++;
        rates_.push_back(rate.second);
    }

    // Duplicate the first and last groups to handle P < P_0 and P > P_N
    pressures_.insert({-1000.0, pressures_.begin()->second});
    pressures_.insert({1000.0, pressures_.rbegin()->second});
}

void PlogRate::validate(const std::string& equation)
{
    fmt::memory_buffer err_reactions;
    double T[] = {200.0, 500.0, 1000.0, 2000.0, 5000.0, 10000.0};

    for (auto iter = ++pressures_.begin(); iter->first < 1000; iter++) {
        update_C(&iter->first);
        for (size_t i=0; i < 6; i++) {
            double k = 0;
            for (size_t p = ilow1_; p < ilow2_; p++) {
                k += rates_.at(p).evalRate(log(T[i]), 1.0 / T[i]);
            }
            if (!(k > 0)) {
                fmt_append(err_reactions,
                           "\nInvalid rate coefficient for reaction '{}'\n"
                           "at P = {:.5g}, T = {:.1f}\n",
                           equation, std::exp(iter->first), T[i]);
            }
        }
    }
    if (err_reactions.size()) {
        throw CanteraError("PlogRate::validate", to_string(err_reactions));
    }
}

std::vector<std::pair<double, Arrhenius2> > PlogRate::rates() const
{
    auto rateMap = getRates();
    std::vector<std::pair<double, Arrhenius2>> out;
    for (const auto& item : rateMap) {
        out.emplace_back(item.first, Arrhenius2(item.second));
    }
    return out;
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
