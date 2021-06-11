//! @file RxnRates.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/RxnRates.h"
#include "cantera/base/Array.h"
#include "cantera/base/AnyMap.h"
#include "cantera/base/global.h"

namespace Cantera
{
Arrhenius::Arrhenius()
    : m_logA(-1.0E300)
    , m_b(0.0)
    , m_E(0.0)
    , m_A(0.0)
{
}

Arrhenius::Arrhenius(doublereal A, doublereal b, doublereal E)
    : m_b(b)
    , m_E(E)
    , m_A(A)
{
    if (m_A  <= 0.0) {
        m_logA = -1.0E300;
    } else {
        m_logA = std::log(m_A);
    }
}

Arrhenius::Arrhenius(const AnyValue& rate,
                     const UnitSystem& units, const Units& rate_units)
{
    setParameters(rate, units, rate_units);
}

void Arrhenius::setParameters(const AnyValue& rate,
                              const UnitSystem& units, const Units& rate_units)
{
    if (rate.empty()) {
        m_A = SNAN;
        m_b = SNAN;
        m_E = SNAN;
    } else if (rate.is<AnyMap>()) {
        auto& rate_map = rate.as<AnyMap>();
        m_A = units.convert(rate_map["A"], rate_units);
        m_b = rate_map["b"].asDouble();
        m_E = units.convertActivationEnergy(rate_map["Ea"], "K");
    } else {
        auto& rate_vec = rate.asVector<AnyValue>(3);
        m_A = units.convert(rate_vec[0], rate_units);
        m_b = rate_vec[1].asDouble();
        m_E = units.convertActivationEnergy(rate_vec[2], "K");
    }

    if (m_A <= 0.0) {
        m_logA = -1.0E300;
    } else {
        m_logA = std::log(m_A);
    }
}

void Arrhenius::getParameters(AnyMap& rateNode, const Units& rate_units) const
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
    rateNode["Ea"].setQuantity(activationEnergy_R(), "K", true);
    rateNode.setFlowStyle();
}

BlowersMasel::BlowersMasel()
    : m_logA(-1.0E300)
    , m_b(0.0)
    , m_A(0.0)
    , m_w(0.0)
    , m_E0(0.0)
{
}

BlowersMasel::BlowersMasel(double A, double b, double E0, double w)
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

void BlowersMasel::getParameters(AnyMap& rateNode, const Units& rate_units) const
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

Plog::Plog()
    : logP_(-1000)
    , logP1_(1000)
    , logP2_(-1000)
    , rDeltaP_(-1.0)
{
}

Plog::Plog(const std::multimap<double, Arrhenius>& rates)
    : Plog()
{
    setup(rates);
}

void Plog::setParameters(const std::vector<AnyMap>& rates,
                         const UnitSystem& units, const Units& rate_units)
{
    std::multimap<double, Arrhenius> multi_rates;
    if (rates.size()) {
        for (const auto& rate : rates) {
            multi_rates.insert({rate.convert("P", "Pa"),
                Arrhenius(AnyValue(rate), units, rate_units)});
        }
    } else {
        // ensure that reaction rate can be evaluated (but returns SNAN)
        multi_rates.insert({1.e-7, Arrhenius(SNAN, SNAN, SNAN)});
        multi_rates.insert({1.e14, Arrhenius(SNAN, SNAN, SNAN)});
    }
    setup(multi_rates);
}

void Plog::getParameters(AnyMap& rateNode, const Units& rate_units) const
{
    std::vector<AnyMap> rateList;
    for (const auto& r : rates()) {
        AnyMap rateNode_;
        rateNode_["P"].setQuantity(r.first, "Pa");
        r.second.getParameters(rateNode_, rate_units);
        rateList.push_back(std::move(rateNode_));
    }
    rateNode["rate-constants"] = std::move(rateList);
}

void Plog::setup(const std::multimap<double, Arrhenius>& rates)
{
    size_t j = 0;
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

void Plog::validate(const std::string& equation)
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
                format_to(err_reactions,
                          "\nInvalid rate coefficient for reaction '{}'\n"
                          "at P = {:.5g}, T = {:.1f}\n",
                          equation, std::exp(iter->first), T[i]);
            }
        }
    }
    if (err_reactions.size()) {
        throw CanteraError("Plog::validate", to_string(err_reactions));
    }
}

std::vector<std::pair<double, Arrhenius> > Plog::rates() const
{
    std::vector<std::pair<double, Arrhenius> > R;
    // initial preincrement to skip rate for P --> 0
    for (auto iter = ++pressures_.begin();
         iter->first < 1000; // skip rates for (P --> infinity)
         ++iter) {
        for (size_t i = iter->second.first;
             i < iter->second.second;
             i++) {
            R.emplace_back(std::exp(iter->first), rates_[i]);
        }
    }
    return R;
}

Chebyshev::Chebyshev(double Tmin, double Tmax, double Pmin, double Pmax,
                     const Array2D& coeffs)
    : Tmin_(Tmin)
    , Tmax_(Tmax)
    , Pmin_(Pmin)
    , Pmax_(Pmax)
    , nP_(coeffs.nColumns())
    , nT_(coeffs.nRows())
    , chebCoeffs_(coeffs.nColumns() * coeffs.nRows(), 0.0)
    , dotProd_(coeffs.nRows())
{
    setup(Tmin, Tmax, Pmin, Pmax, coeffs);
}

void Chebyshev::setParameters(const AnyMap& node,
                              const UnitSystem& units, const Units& rate_units)
{
    Array2D coeffs;
    if (!node.empty()) {
        const auto& T_range = node["temperature-range"].asVector<AnyValue>(2);
        const auto& P_range = node["pressure-range"].asVector<AnyValue>(2);
        auto& vcoeffs = node["data"].asVector<vector_fp>();
        coeffs = Array2D(vcoeffs.size(), vcoeffs[0].size());
        for (size_t i = 0; i < coeffs.nRows(); i++) {
            if (vcoeffs[i].size() != vcoeffs[0].size()) {
                throw InputFileError("Chebyshev::setParameters", node["data"],
                    "Inconsistent number of coefficients in row {} of matrix", i + 1);
            }
            for (size_t j = 0; j < coeffs.nColumns(); j++) {
                coeffs(i, j) = vcoeffs[i][j];
            }
        }
        coeffs(0, 0) += std::log10(units.convertTo(1.0, rate_units));

        Tmin_ = units.convert(T_range[0], "K");
        Tmax_ = units.convert(T_range[1], "K");
        Pmin_ = units.convert(P_range[0], "Pa");
        Pmax_ = units.convert(P_range[1], "Pa");
        nP_ = coeffs.nColumns();
        nT_ = coeffs.nRows();
    } else {
        // ensure that reaction rate can be evaluated (but returns SNAN)
        coeffs = Array2D(1, 1);
        coeffs(0, 0) = SNAN;
        Tmin_ = 290.;
        Tmax_ = 3000.;
        Pmin_ = 1.e-7;
        Pmax_ = 1.e14;
        nP_ = 1;
        nT_ = 1;
    }

    chebCoeffs_.resize(nP_ * nT_);
    dotProd_.resize(nT_);
    setup(Tmin_, Tmax_, Pmin_, Pmax_, coeffs);
}

void Chebyshev::setup(double Tmin, double Tmax, double Pmin, double Pmax,
                      const Array2D& coeffs)
{
    double logPmin = std::log10(Pmin);
    double logPmax = std::log10(Pmax);
    double TminInv = 1.0 / Tmin;
    double TmaxInv = 1.0 / Tmax;

    TrNum_ = - TminInv - TmaxInv;
    TrDen_ = 1.0 / (TmaxInv - TminInv);
    PrNum_ = - logPmin - logPmax;
    PrDen_ = 1.0 / (logPmax - logPmin);

    for (size_t t = 0; t < nT_; t++) {
        for (size_t p = 0; p < nP_; p++) {
            chebCoeffs_[nP_*t + p] = coeffs(t,p);
        }
    }
}

void Chebyshev::getParameters(AnyMap& rateNode, const Units& rate_units) const
{
    rateNode["temperature-range"].setQuantity({Tmin(), Tmax()}, "K");
    rateNode["pressure-range"].setQuantity({Pmin(), Pmax()}, "Pa");
    const auto& coeffs1d = coeffs();
    size_t nT = nTemperature();
    size_t nP = nPressure();
    std::vector<vector_fp> coeffs2d(nT, vector_fp(nP));
    for (size_t i = 0; i < nT; i++) {
        for (size_t j = 0; j < nP; j++) {
            coeffs2d[i][j] = coeffs1d[nP*i + j];
        }
    }
    // Unit conversions must take place later, after the destination unit system
    // is known. A lambda function is used here to override the default behavior
    Units rate_units2 = rate_units;
    auto converter = [rate_units2](AnyValue& coeffs, const UnitSystem& units) {
        if (rate_units2.factor() != 0.0) {
            coeffs.asVector<vector_fp>()[0][0] += std::log10(units.convertFrom(1.0, rate_units2));
        } else if (units.getDelta(UnitSystem()).size()) {
            throw CanteraError("ChebyshevReaction::getParameters lambda",
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

ChebyshevRate::ChebyshevRate()
    : Chebyshev()
{
    warn_deprecated("ChebyshevRate::ChebyshevRate",
                    "Renamed to Chebyshev. Behavior will change after Cantera 2.6. "
                    "For future behavior, refer to ChebyshevRate3");
}

ChebyshevRate::ChebyshevRate(double Tmin, double Tmax, double Pmin, double Pmax,
                             const Array2D& coeffs)
    : Chebyshev(Tmin, Tmax, Pmin, Pmax, coeffs)
{
    warn_deprecated("ChebyshevRate::ChebyshevRate",
                    "Renamed to Chebyshev. Behavior will change after Cantera 2.6. "
                    "For future behavior, refer to ChebyshevRate3");}
}
