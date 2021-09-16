//! @file RxnRates.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/RxnRates.h"
#include "cantera/base/AnyMap.h"

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
        m_A = NAN;
        m_b = NAN;
        m_E = NAN;
    } else if (rate.is<AnyMap>()) {
        auto& rate_map = rate.as<AnyMap>();
        if (rate_units.factor() == 0 && rate_map["A"].is<std::string>()) {
            // A zero rate units factor is used as a sentinel to detect
            // stand-alone reaction rate objects
            throw InputFileError("Arrhenius::setParameters", rate_map,
                "Specification of units is not supported for pre-exponential factor "
                "when\ncreating a standalone 'ReactionRate' object.");
        }
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
    double A = preExponentialFactor();
    if (std::isnan(A)) {
        // Return empty/unmodified AnyMap
        return;
    } else if (rate_units.factor() != 0.0) {
        rateNode["A"].setQuantity(A, rate_units);
    } else {
        rateNode["A"] = A;
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
    setRates(rates);
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
        // ensure that reaction rate can be evaluated (but returns NaN)
        multi_rates.insert({1.e-7, Arrhenius(NAN, NAN, NAN)});
        multi_rates.insert({1.e14, Arrhenius(NAN, NAN, NAN)});
    }
    setRates(multi_rates);
}

void Plog::getParameters(AnyMap& rateNode, const Units& rate_units) const
{
    std::vector<AnyMap> rateList;
    double A = rates_[1].preExponentialFactor();
    if (std::isnan(A)) {
        // Return empty/unmodified AnyMap
        return;
    }
    for (const auto& r : getRates()) {
        AnyMap rateNode_;
        rateNode_["P"].setQuantity(r.first, "Pa");
        r.second.getParameters(rateNode_, rate_units);
        rateList.push_back(std::move(rateNode_));
    }
    rateNode["rate-constants"] = std::move(rateList);
}

void Plog::setup(const std::multimap<double, Arrhenius>& rates)
{
    warn_deprecated("Plog::setup", "Deprecated in Cantera 2.6; "
        "renamed to setRates.");
    setRates(rates);
}

void Plog::setRates(const std::multimap<double, Arrhenius>& rates)
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
    warn_deprecated("Plog::rates", "Behavior to change after Cantera 2.6; "
        "see getRates() for new behavior.");
    auto rateMap = getRates();
    return std::vector<std::pair<double, Arrhenius>>(rateMap.begin(), rateMap.end());
}

std::multimap<double, Arrhenius> Plog::getRates() const
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

Chebyshev::Chebyshev(double Tmin, double Tmax, double Pmin, double Pmax,
                     const Array2D& coeffs)
{
    setLimits(Tmin, Tmax, Pmin, Pmax);
    setData(coeffs);
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
        if (rate_units.factor()) {
            coeffs(0, 0) += std::log10(units.convertTo(1.0, rate_units));
        }
        setLimits(
            units.convert(T_range[0], "K"), units.convert(T_range[1], "K"),
            units.convert(P_range[0], "Pa"), units.convert(P_range[1], "Pa"));
    } else {
        // ensure that reaction rate can be evaluated (but returns NaN)
        coeffs = Array2D(1, 1);
        coeffs(0, 0) = NAN;
        setLimits(290., 3000., 1.e-7, 1.e14);
    }

    setData(coeffs);
}

void Chebyshev::setup(double Tmin, double Tmax, double Pmin, double Pmax,
                      const Array2D& coeffs)
{
    warn_deprecated("Chebyshev::setup", "Deprecated in Cantera 2.6; "
        "replaceable with setLimits() and setData().");
    setLimits(Tmin, Tmax, Pmin, Pmax);
    setData(coeffs);
}

void Chebyshev::setLimits(double Tmin, double Tmax, double Pmin, double Pmax)
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

void Chebyshev::setData(const Array2D& coeffs)
{
    m_coeffs = coeffs;
    dotProd_.resize(coeffs.nRows());

    // convert to row major for legacy output
    // note: chebCoeffs_ is not used internally (@TODO: remove after Cantera 2.6)
    size_t rows = m_coeffs.nRows();
    size_t cols = m_coeffs.nColumns();
    chebCoeffs_.resize(rows * cols);
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            chebCoeffs_[cols * i + j] = m_coeffs(i, j);
        }
    }
}

void Chebyshev::getParameters(AnyMap& rateNode, const Units& rate_units) const
{
    if (std::isnan(m_coeffs(0, 0))) {
        // Return empty/unmodified AnyMap
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
    Units rate_units2 = rate_units;
    auto converter = [rate_units2](AnyValue& coeffs, const UnitSystem& units) {
        if (rate_units2.factor() != 0.0) {
            coeffs.asVector<vector_fp>()[0][0] += std::log10(units.convertFrom(1.0, rate_units2));
        } else if (units.getDelta(UnitSystem()).size()) {
            throw CanteraError("Chebyshev::getParameters lambda",
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
