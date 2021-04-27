//! @file ReactionRate.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/ReactionRate.h"
#include "cantera/numerics/Func1.h"
#include "cantera/base/AnyMap.h"
#include "cantera/base/Units.h"

namespace Cantera
{

AnyMap ReactionRateBase::parameters(const Units& rate_units) const
{
    AnyMap out;
    getParameters(out, rate_units);
    return out;
}

AnyMap ReactionRateBase::parameters() const
{
    AnyMap out;
    getParameters(out, units);
    return out;
}

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

ArrheniusRate::ArrheniusRate(const Arrhenius& arr, bool allow_negative_A)
    : Arrhenius(arr.preExponentialFactor(),
                arr.temperatureExponent(),
                arr.activationEnergy_R())
    , allow_negative_pre_exponential_factor(allow_negative_A)
{
}

ArrheniusRate::ArrheniusRate(const AnyMap& node, const Units& rate_units) {
    setParameters(node, rate_units);
}

void ArrheniusRate::setParameters(const AnyMap& node, const Units& rate_units) {
    units = rate_units;
    allow_negative_pre_exponential_factor = node.getBool("negative-A", false);
    if (!node.hasKey("rate-constant")) {
        return;
    }
    Arrhenius::setParameters(node["rate-constant"], node.units(), rate_units);
}

void ArrheniusRate::getParameters(AnyMap& rateNode,
                                  const Units& rate_units) const {
    if (allow_negative_pre_exponential_factor) {
        rateNode["negative-A"] = true;
    }
    AnyMap node;
    Arrhenius::getParameters(node, rate_units);
    rateNode["rate-constant"] = std::move(node);
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

void PlogRate::setParameters(const AnyMap& node, const Units& rate_units) {
    units = rate_units;
    if (!node.hasKey("rate-constants")) {
        // ensure that Plog has defined state and produces zero reaction rate
        AnyMap rate = AnyMap::fromYamlString(
            "rate-constants:\n"
            "- {P: 1e-7, A: 0., b: 0., Ea: 0.}\n"
            "- {P: 1e7, A: 0., b: 0., Ea: 0.}");
        Plog::setParameters(rate.at("rate-constants").asVector<AnyMap>(),
                            node.units(), rate_units);
        return;
    }
    Plog::setParameters(node.at("rate-constants").asVector<AnyMap>(),
                        node.units(), rate_units);
}

void PlogRate::getParameters(AnyMap& rateNode,
                             const Units& rate_units) const {
    Plog::getParameters(rateNode, rate_units);
}

ChebyshevRate3::ChebyshevRate3()
    : Chebyshev() {
}

ChebyshevRate3::ChebyshevRate3(double Tmin, double Tmax, double Pmin, double Pmax,
                               const Array2D& coeffs)
    : Chebyshev(Tmin, Tmax, Pmin, Pmax, coeffs) {
}

ChebyshevRate3::ChebyshevRate3(const AnyMap& node, const Units& rate_units)
    : Chebyshev() {
    setParameters(node, rate_units);
}

void ChebyshevRate3::setParameters(const AnyMap& node, const Units& rate_units) {
    units = rate_units;
    if (!node.hasKey("data")) {
        // ensure that Chebyshev has defined state and produces zero reaction rate
        AnyMap rate = AnyMap::fromYamlString(
            "temperature-range: [290, 3000]\n"
            "pressure-range: [1.e-7, 1.e7]\n"
            "data: [[-16.]]\n");
        Chebyshev::setParameters(rate, node.units(), rate_units);
        return;
    }
    Chebyshev::setParameters(node, node.units(), rate_units);
}

void ChebyshevRate3::getParameters(AnyMap& rateNode,
                                   const Units& rate_units) const {
    Chebyshev::getParameters(rateNode, rate_units);
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
