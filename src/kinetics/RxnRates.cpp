//! @file RxnRates.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/RxnRates.h"
#include "cantera/base/AnyMap.h"

namespace Cantera
{
Arrhenius2::Arrhenius2()
    : Arrhenius3()
{
    m_b = 0.0;
    m_A = 0.0;
    m_logA = -1.0E300;
}

Arrhenius2::Arrhenius2(doublereal A, doublereal b, doublereal E)
    : Arrhenius3(A, b, E * GasConstant)
{
    if (m_A  <= 0.0) {
        m_logA = -1.0E300;
    }
}

Arrhenius2::Arrhenius2(const AnyValue& rate,
                       const UnitSystem& units, const Units& rate_units)
{
    setRateParameters(rate, units, rate_units);
}

Arrhenius2::Arrhenius2(const Arrhenius3& other)
    : Arrhenius3(other.preExponentialFactor(),
                 other.temperatureExponent(),
                 other.activationEnergy())
{
}

void Arrhenius2::setRateParameters(const AnyValue& rate,
                                   const UnitSystem& units, const Units& rate_units)
{
    UnitStack units_stack(rate_units);
    Arrhenius3::setRateParameters(rate, units, units_stack);
    if (m_A <= 0.0) {
        m_logA = -1.0E300;
    }
}

void Arrhenius2::getParameters(AnyMap& node, const Units& rate_units) const
{
    if (rate_units.factor() != 0.0) {
        node["A"].setQuantity(m_A, rate_units);
    } else {
        node["A"] = preExponentialFactor();
        // This can't be converted to a different unit system because the dimensions of
        // the rate constant were not set. Can occur if the reaction was created outside
        // the context of a Kinetics object and never added to a Kinetics object.
        node["__unconvertible__"] = true;
    }
    node["b"] = m_b;
    node["Ea"].setQuantity(m_Ea_R, "K", true);
    node.setFlowStyle();
}

SurfaceArrhenius::SurfaceArrhenius()
    : m_b(0.0)
    , m_E(0.0)
    , m_A(0.0)
    , m_acov(0.0)
    , m_ecov(0.0)
    , m_mcov(0.0)
{
}

SurfaceArrhenius::SurfaceArrhenius(double A, double b, double Ta)
    : m_b(b)
    , m_E(Ta)
    , m_A(A)
    , m_acov(0.0)
    , m_ecov(0.0)
    , m_mcov(0.0)
{
}

void SurfaceArrhenius::addCoverageDependence(size_t k, doublereal a,
                               doublereal m, doublereal e)
{
    m_sp.push_back(k);
    m_ac.push_back(a);
    m_ec.push_back(e);
    if (m != 0.0) {
        m_msp.push_back(k);
        m_mc.push_back(m);
    }
}

PlogRate::PlogRate()
    : logP_(-1000)
    , logP1_(1000)
    , logP2_(-1000)
    , rDeltaP_(-1.0)
{
}

PlogRate::PlogRate(const std::multimap<double, Arrhenius3>& rates)
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
    auto rate_units = units.product();
    if (!node.hasKey("rate-constants")) {
        PlogRate::setRateParameters(std::vector<AnyMap> (), node.units(), rate_units);
        return;
    }

    setRateParameters(node.at("rate-constants").asVector<AnyMap>(),
                      node.units(), rate_units);
}

void PlogRate::setRateParameters(const std::vector<AnyMap>& rates,
                                 const UnitSystem& units, const Units& rate_units)
{
    std::multimap<double, Arrhenius3> multi_rates;
    if (rates.size()) {
        for (const auto& rate : rates) {
            multi_rates.insert({rate.convert("P", "Pa"),
                Arrhenius3(AnyValue(rate), units, rate_units)});
        }
    } else {
        // ensure that reaction rate can be evaluated (but returns NaN)
        multi_rates.insert({1.e-7, Arrhenius3(NAN, NAN, NAN)});
        multi_rates.insert({1.e14, Arrhenius3(NAN, NAN, NAN)});
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
    std::multimap<double, Arrhenius3> rates2;
    for (const auto& item : rates) {
        rates2.emplace(item.first, item.second);
    }
    setRates(rates2);
}

void PlogRate::setRates(const std::multimap<double, Arrhenius3>& rates)
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
            if (k < 0) {
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

std::multimap<double, Arrhenius3> PlogRate::getRates() const
{
    std::multimap<double, Arrhenius3> rateMap;
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

ChebyshevRate::ChebyshevRate(double Tmin, double Tmax, double Pmin, double Pmax,
                             const Array2D& coeffs) : ChebyshevRate()
{
    setLimits(Tmin, Tmax, Pmin, Pmax);
    setData(coeffs);
}

void ChebyshevRate::setParameters(const AnyMap& node, const UnitStack& units)
{
    m_rate_units = units.product();
    const UnitSystem& unit_system = node.units();
    Array2D coeffs;
    if (node.hasKey("data")) {
        const auto& T_range = node["temperature-range"].asVector<AnyValue>(2);
        const auto& P_range = node["pressure-range"].asVector<AnyValue>(2);
        auto& vcoeffs = node["data"].asVector<vector_fp>();
        coeffs = Array2D(vcoeffs.size(), vcoeffs[0].size());
        for (size_t i = 0; i < coeffs.nRows(); i++) {
            if (vcoeffs[i].size() != vcoeffs[0].size()) {
                throw InputFileError("ChebyshevRate::setParameters", node["data"],
                    "Inconsistent number of coefficients in row {} of matrix", i + 1);
            }
            for (size_t j = 0; j < coeffs.nColumns(); j++) {
                coeffs(i, j) = vcoeffs[i][j];
            }
        }
        if (m_rate_units.factor()) {
            coeffs(0, 0) += std::log10(unit_system.convertTo(1.0, m_rate_units));
        }
        setLimits(
            unit_system.convert(T_range[0], "K"),
            unit_system.convert(T_range[1], "K"),
            unit_system.convert(P_range[0], "Pa"),
            unit_system.convert(P_range[1], "Pa")
        );
    } else {
        // ensure that reaction rate can be evaluated (but returns NaN)
        coeffs = Array2D(1, 1);
        coeffs(0, 0) = NAN;
        setLimits(290., 3000., 1.e-7, 1.e14);
    }

    setData(coeffs);
}

void ChebyshevRate::setup(double Tmin, double Tmax, double Pmin, double Pmax,
                          const Array2D& coeffs)
{
    warn_deprecated("ChebyshevRate::setup", "Deprecated in Cantera 2.6; "
        "replaceable with setLimits() and setData().");
    setLimits(Tmin, Tmax, Pmin, Pmax);
    setData(coeffs);
}

void ChebyshevRate::setLimits(double Tmin, double Tmax, double Pmin, double Pmax)
{
    double logPmin = std::log10(Pmin);
    double logPmax = std::log10(Pmax);
    double TminInv = 1.0 / Tmin;
    double TmaxInv = 1.0 / Tmax;

    TrNum_ = - TminInv - TmaxInv;
    TrDen_ = 1.0 / (TmaxInv - TminInv);
    PrNum_ = - logPmin - logPmax;
    PrDen_ = 1.0 / (logPmax - logPmin);

    Tmin_ = Tmin;
    Tmax_ = Tmax;
    Pmin_ = Pmin;
    Pmax_ = Pmax;
}

void ChebyshevRate::setData(const Array2D& coeffs)
{
    m_coeffs = coeffs;
    dotProd_.resize(coeffs.nRows());

    // convert to row major for legacy output
    // note: chebCoeffs_ is not used internally (@todo: remove after Cantera 2.6)
    size_t rows = m_coeffs.nRows();
    size_t cols = m_coeffs.nColumns();
    chebCoeffs_.resize(rows * cols);
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            chebCoeffs_[cols * i + j] = m_coeffs(i, j);
        }
    }
}

void ChebyshevRate::getParameters(AnyMap& rateNode) const
{
    rateNode["type"] = type();
    if (!m_coeffs.data().size() || std::isnan(m_coeffs(0, 0))) {
        // object not fully set up
        return;
    }
    rateNode["temperature-range"].setQuantity({Tmin(), Tmax()}, "K");
    rateNode["pressure-range"].setQuantity({Pmin(), Pmax()}, "Pa");
    size_t nT = m_coeffs.nRows();
    size_t nP = m_coeffs.nColumns();
    std::vector<vector_fp> coeffs2d(nT, vector_fp(nP));
    for (size_t i = 0; i < nT; i++) {
        for (size_t j = 0; j < nP; j++) {
            coeffs2d[i][j] = m_coeffs(i, j);
        }
    }
    // Unit conversions must take place later, after the destination unit system
    // is known. A lambda function is used here to override the default behavior
    Units rate_units2 = m_rate_units;
    auto converter = [rate_units2](AnyValue& coeffs, const UnitSystem& units) {
        if (rate_units2.factor() != 0.0) {
            coeffs.asVector<vector_fp>()[0][0] += \
                std::log10(units.convertFrom(1.0, rate_units2));
        } else if (units.getDelta(UnitSystem()).size()) {
            throw CanteraError("ChebyshevRate::getParameters lambda",
                "Cannot convert rate constant with unknown dimensions to a "
                "non-default unit system");
        }
    };
    AnyValue coeffs;
    coeffs = std::move(coeffs2d);
    rateNode["data"].setQuantity(coeffs, converter);
}

}
