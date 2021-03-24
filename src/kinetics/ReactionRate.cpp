//! @file ReactionRate.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/ReactionRate.h"
#include "cantera/numerics/Func1.h"
#include "cantera/base/AnyMap.h"

namespace Cantera
{

ArrheniusRate::ArrheniusRate()
    : Arrhenius()
    , allow_negative_pre_exponential_factor(false)
{
}

ArrheniusRate::ArrheniusRate(double A, double b, double E)
    : Arrhenius(A, b, E / GasConstant)
    , allow_negative_pre_exponential_factor(false)
{
}

ArrheniusRate::ArrheniusRate(const AnyMap& node, const Units& rate_units) {
    setParameters(node, rate_units);
}

bool ArrheniusRate::setParameters(const AnyMap& node, const Units& rate_units) {
    allow_negative_pre_exponential_factor = node.getBool("negative-A", false);
    if (!node.hasKey("rate-constant")) {
        return false;
    }
    Arrhenius::setParameters(node["rate-constant"], node.units(), rate_units);
    return true;
}

void ArrheniusRate::getParameters(AnyMap& rateNode,
                                  const Units& rate_units) const {
    Arrhenius::getParameters(rateNode, rate_units);
}

void ArrheniusRate::validate(const std::string& equation) {
    if (!allow_negative_pre_exponential_factor && preExponentialFactor() < 0) {
        throw CanteraError("ArrheniusRate::validate",
            "Undeclared negative pre-exponential factor found in reaction '"
            + equation + "'");
    }
}

PlogRate::PlogRate()
    : Plog() {
}

PlogRate::PlogRate(const std::multimap<double, Arrhenius>& rates)
    : Plog(rates) {
}

PlogRate::PlogRate(const AnyMap& node, const Units& rate_units)
    : Plog() {
    setParameters(node, rate_units);
}

bool PlogRate::setParameters(const AnyMap& node, const Units& rate_units) {
    if (!node.hasKey("rate-constants")) {
        return false;
    }
    Plog::setParameters(node.at("rate-constants").asVector<AnyMap>(),
                        node.units(), rate_units);
    return true;
}

void PlogRate::getParameters(AnyMap& rateNode,
                             const Units& rate_units) const {
    Plog::getParameters(rateNode, rate_units);
}

ChebyshevRate3::ChebyshevRate3()
    : ChebyshevRate() {
}

ChebyshevRate3::ChebyshevRate3(double Tmin, double Tmax, double Pmin, double Pmax,
                               const Array2D& coeffs)
    : ChebyshevRate(Tmin, Tmax, Pmin, Pmax, coeffs) {
}

ChebyshevRate3::ChebyshevRate3(const AnyMap& node, const Units& rate_units)
    : ChebyshevRate() {
    setParameters(node, rate_units);
}

bool ChebyshevRate3::setParameters(const AnyMap& node, const Units& rate_units) {
    if (!node.hasKey("rate-constants")) {
        return false;
    }
    ChebyshevRate::setParameters(node.at("rate-constants").asVector<AnyMap>(),
                                 node.units(), rate_units);
    return true;
}

void ChebyshevRate3::getParameters(AnyMap& rateNode,
                                   const Units& rate_units) const {
    throw CanteraError("ChebyshevRate3::getParameters",
        "@todo");
}

void ChebyshevRate3::validate(const std::string& equation) {
}

CustomFunc1Rate::CustomFunc1Rate() : m_ratefunc(0) {}

void CustomFunc1Rate::setRateFunction(shared_ptr<Func1> f) {
    m_ratefunc = f;
}

double CustomFunc1Rate::eval(const CustomFunc1Data& shared_data,
                             double concm) const {
    if (m_ratefunc) {
        return m_ratefunc->eval(shared_data.m_temperature);
    }
    throw CanteraError("CustomFunc1Rate::eval",
                       "Custom rate function is not initialized.");
}

}
