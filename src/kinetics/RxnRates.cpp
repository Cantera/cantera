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
