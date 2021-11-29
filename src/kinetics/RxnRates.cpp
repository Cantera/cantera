//! @file RxnRates.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/RxnRates.h"
#include "cantera/base/AnyMap.h"

namespace Cantera
{
Arrhenius2::Arrhenius2()
    : ArrheniusBase()
    , m_logA(-1.0E300)
{
    m_b = 0.0;
    m_A = 0.0;
}

Arrhenius2::Arrhenius2(doublereal A, doublereal b, doublereal E)
    : ArrheniusBase(A, b, E * GasConstant)
{
    if (m_A  <= 0.0) {
        m_logA = -1.0E300;
    } else {
        m_logA = std::log(m_A);
    }
}

Arrhenius2::Arrhenius2(const AnyValue& rate,
                       const UnitSystem& units, const Units& rate_units)
{
    setRateParameters(rate, units, rate_units);
}

void Arrhenius2::setRateParameters(const AnyValue& rate,
                                   const UnitSystem& units, const Units& rate_units)
{
    UnitStack units_stack(rate_units);
    ArrheniusBase::setRateParameters(rate, units, units_stack);
    if (m_A <= 0.0) {
        m_logA = -1.0E300;
    } else {
        m_logA = std::log(m_A);
    }
}

BlowersMasel2::BlowersMasel2()
    : m_logA(-1.0E300)
    , m_b(0.0)
    , m_A(0.0)
    , m_w(0.0)
    , m_E0(0.0)
{
}

BlowersMasel2::BlowersMasel2(double A, double b, double E0, double w)
    : m_b(b)
    , m_A(A)
    , m_w(w)
    , m_E0(E0)
{
    if (m_A  <= 0.0) {
        m_logA = -1.0E300;
    } else {
        m_logA = std::log(m_A);
    }
}

void BlowersMasel2::getParameters(AnyMap& rateNode, const Units& rate_units) const
{
    if (rate_units.factor() != 0.0) {
        rateNode["A"].setQuantity(preExponentialFactor(), rate_units);
    } else {
        rateNode["A"] = preExponentialFactor();
        // This can't be converted to a different unit system because the dimensions of
        // the rate constant were not set. Can occur if the reaction was created outside
        // the context of a Kinetics object and never added to a Kinetics object.
        rateNode["__unconvertible__"] = true;
    }

    rateNode["b"] = temperatureExponent();
    rateNode["Ea0"].setQuantity(activationEnergy_R0(), "K", true);
    rateNode["w"].setQuantity(bondEnergy(), "K", true);
    rateNode.setFlowStyle();
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

PlogRate::PlogRate(const std::multimap<double, Arrhenius>& rates)
    : PlogRate()
{
    setRates(rates);
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
    std::multimap<double, Arrhenius> multi_rates;
    if (rates.size()) {
        for (const auto& rate : rates) {
            multi_rates.insert({rate.convert("P", "Pa"),
                Arrhenius(AnyValue(rate), units, rate_units)});
        }
    } else {
        // ensure that reaction rate can be evaluated (but returns NaN)
        multi_rates.insert({1.e-7, Arrhenius(NAN, NAN, NAN)});
        multi_rates.insert({1.e14, Arrhenius(NAN, NAN, NAN)});
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
            r.second.getParameters(rateNode_);
        } else {
            r.second.getParameters(rateNode_, rate_units);
        }
        rateList.push_back(std::move(rateNode_));
    }
    rateNode["rate-constants"] = std::move(rateList);
}

void PlogRate::setup(const std::multimap<double, Arrhenius>& rates)
{
    warn_deprecated("PlogRate::setup", "Deprecated in Cantera 2.6; "
        "renamed to setRates.");
    setRates(rates);
}

void PlogRate::setRates(const std::multimap<double, Arrhenius>& rates)
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
                k += rates_.at(p).updateRC(log(T[i]), 1.0/T[i]);
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

std::vector<std::pair<double, Arrhenius> > PlogRate::rates() const
{
    warn_deprecated("PlogRate::rates", "Behavior to change after Cantera 2.6; "
        "see getRates() for new behavior.");
    auto rateMap = getRates();
    return std::vector<std::pair<double, Arrhenius>>(rateMap.begin(), rateMap.end());
}

std::multimap<double, Arrhenius> PlogRate::getRates() const
{
    std::multimap<double, Arrhenius> rateMap;
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

ChebyshevRate3::ChebyshevRate3(double Tmin, double Tmax, double Pmin, double Pmax,
                               const Array2D& coeffs) : ChebyshevRate3()
{
    setLimits(Tmin, Tmax, Pmin, Pmax);
    setData(coeffs);
}

void ChebyshevRate3::setParameters(const AnyMap& node, const UnitStack& units)
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
                throw InputFileError("ChebyshevRate3::setParameters", node["data"],
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

void ChebyshevRate3::setup(double Tmin, double Tmax, double Pmin, double Pmax,
                      const Array2D& coeffs)
{
    warn_deprecated("ChebyshevRate3::setup", "Deprecated in Cantera 2.6; "
        "replaceable with setLimits() and setData().");
    setLimits(Tmin, Tmax, Pmin, Pmax);
    setData(coeffs);
}

void ChebyshevRate3::setLimits(double Tmin, double Tmax, double Pmin, double Pmax)
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

void ChebyshevRate3::setData(const Array2D& coeffs)
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

void ChebyshevRate3::getParameters(AnyMap& rateNode) const
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
            throw CanteraError("ChebyshevRate3::getParameters lambda",
                "Cannot convert rate constant with unknown dimensions to a "
                "non-default unit system");
        }
    };
    AnyValue coeffs;
    coeffs = std::move(coeffs2d);
    rateNode["data"].setQuantity(coeffs, converter);
}
BMSurfaceArrhenius::BMSurfaceArrhenius()
    : m_b(0.0)
    , m_A(0.0)
    , m_E0(0.0)
    , m_w(0.0)
    , m_acov(0.0)
    , m_ecov(0.0)
    , m_mcov(0.0)
{
}

BMSurfaceArrhenius::BMSurfaceArrhenius(double A, double b, double Ta, double w)
    : m_b(b)
    , m_A(A)
    , m_E0(Ta)
    , m_w(w)
    , m_acov(0.0)
    , m_ecov(0.0)
    , m_mcov(0.0)
{
}

void BMSurfaceArrhenius::addCoverageDependence(size_t k, double a,
                               double m, double e)
{
    m_sp.push_back(k);
    m_ac.push_back(a);
    m_ec.push_back(e);
    if (m != 0.0) {
        m_msp.push_back(k);
        m_mc.push_back(m);
    }
}

}
