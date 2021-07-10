//! @file ReactionRate.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/ReactionRate.h"
#include "cantera/kinetics/MultiRate.h"
#include "cantera/base/Array.h"
#include "cantera/numerics/Func1.h"

namespace Cantera
{

void ReactionRateBase::setParameters(const AnyMap& node, const Units& rate_units)
{
    units = rate_units;
    input = node;
}

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

void ReactionRateBase::linkEvaluator(size_t index,
                                     const shared_ptr<MultiRateBase> evaluator)
{
    m_index = index;
    m_evaluator = evaluator;
}

void ReactionRateBase::releaseEvaluator()
{
    m_index = npos;
    m_evaluator.reset();
}

size_t ReactionRateBase::index()
{
    if (m_evaluator) {
        return m_index;
    }
    throw CanteraError("ReactionRateBase::index", "Not applicable, as reaction rate "
        "is not linked to Kinetics object with assoicated rate evaluator");
}

ArrheniusRate::ArrheniusRate()
    : Arrhenius(NAN, NAN, NAN)
    , allow_negative_pre_exponential_factor(false)
{
}

ArrheniusRate::ArrheniusRate(double A, double b, double E)
    : Arrhenius(A, b, E / GasConstant)
    , allow_negative_pre_exponential_factor(false)
{
}

ArrheniusRate::ArrheniusRate(const AnyMap& node, const Units& rate_units)
{
    setParameters(node, rate_units);
}

ArrheniusRate::ArrheniusRate(const AnyMap& node)
{
    setParameters(node, Units(1.));
}

ArrheniusRate::ArrheniusRate(const Arrhenius& arr, bool allow_negative_A)
    : Arrhenius(arr.preExponentialFactor(),
                arr.temperatureExponent(),
                arr.activationEnergy_R())
    , allow_negative_pre_exponential_factor(allow_negative_A)
{
}

void ArrheniusRate::setParameters(const AnyMap& node, const Units& rate_units)
{
    ReactionRateBase::setParameters(node, rate_units);
    allow_negative_pre_exponential_factor = node.getBool("negative-A", false);
    if (!node.hasKey("rate-constant")) {
        Arrhenius::setParameters(AnyValue(), node.units(), rate_units);
        return;
    }
    Arrhenius::setParameters(node["rate-constant"], node.units(), rate_units);
}

void ArrheniusRate::getParameters(AnyMap& rateNode,
                                  const Units& rate_units) const
{
    if (allow_negative_pre_exponential_factor) {
        rateNode["negative-A"] = true;
    }
    AnyMap node;
    Arrhenius::getParameters(node, rate_units);
    if (!node.empty()) {
        // Arrhenius object is configured
        rateNode["rate-constant"] = std::move(node);
    }
}

void ArrheniusRate::setPreExponentialFactor(double A)
{
    m_A = A;
    if (m_evaluator) {
        dynamic_cast<ArrheniusRate&>(
            m_evaluator->rate(m_index)).setPreExponentialFactor(A);
    }
}

void ArrheniusRate::setTemperatureExponent(double b)
{
    m_b = b;
    if (m_evaluator) {
        dynamic_cast<ArrheniusRate&>(
            m_evaluator->rate(m_index)).setTemperatureExponent(b);
    }
}

void ArrheniusRate::setActivationEnergy(double E)
{
    m_E = E / GasConstant;
    if (m_evaluator) {
        dynamic_cast<ArrheniusRate&>(
            m_evaluator->rate(m_index)).setActivationEnergy(E);
    }
}

void ArrheniusRate::validate(const std::string& equation)
{
    if (!allow_negative_pre_exponential_factor && preExponentialFactor() < 0) {
        throw CanteraError("ArrheniusRate::validate",
            "Undeclared negative pre-exponential factor found in reaction '"
            + equation + "'");
    }
}

PlogRate::PlogRate(const std::vector<std::pair<double, Arrhenius>>& rates)
    : Plog(rates)
{
}

PlogRate::PlogRate(const AnyMap& node, const Units& rate_units)
{
    setParameters(node, rate_units);
}

PlogRate::PlogRate(const AnyMap& node)
{
    setParameters(node, Units(1.));
}

void PlogRate::setParameters(const AnyMap& node, const Units& rate_units)
{
    // @TODO  implementation of Plog::setParameters should be transferred here
    //     when the Plog class is removed from RxnRates.h after Cantera 2.6
    ReactionRateBase::setParameters(node, rate_units);
    if (!node.hasKey("rate-constants")) {
        Plog::setParameters(std::vector<AnyMap> (), node.units(), rate_units);
        return;
    }
    Plog::setParameters(node.at("rate-constants").asVector<AnyMap>(),
                        node.units(), rate_units);
}

void PlogRate::getParameters(AnyMap& rateNode, const Units& rate_units) const
{
    // @TODO  implementation of Plog::getParameters should be transferred here
    //     when the Plog class is removed from RxnRates.h after Cantera 2.6
    Plog::getParameters(rateNode, rate_units);
}

void PlogRate::setRates(const std::vector<std::pair<double, Arrhenius>>& rates)
{
    Plog::setRates(rates);
    if (m_evaluator) {
        dynamic_cast<PlogRate&>(
            m_evaluator->rate(m_index)).setRates(rates);
    }
}

ChebyshevRate3::ChebyshevRate3(double Tmin, double Tmax, double Pmin, double Pmax,
                               const Array2D& coeffs)
    : Chebyshev(Tmin, Tmax, Pmin, Pmax, coeffs)
{
}

ChebyshevRate3::ChebyshevRate3(const AnyMap& node, const Units& rate_units)
{
    setParameters(node, rate_units);
}

ChebyshevRate3::ChebyshevRate3(const AnyMap& node)
{
    setParameters(node, Units(1.));
}

void ChebyshevRate3::setParameters(const AnyMap& node, const Units& rate_units)
{
    ReactionRateBase::setParameters(node, rate_units);
    if (!node.hasKey("data")) {
        Chebyshev::setParameters(AnyMap(), node.units(), rate_units);
        return;
    }
    // @TODO  implementation of Chebyshev::setParameters should be transferred here
    //     when the Chebyshev class is removed from RxnRates.h after Cantera 2.6
    Chebyshev::setParameters(node, node.units(), rate_units);
}

void ChebyshevRate3::getParameters(AnyMap& rateNode,
                                   const Units& rate_units) const
{
    // @TODO  implementation of Chebyshev::getParameters should be transferred here
    //     when the Chebyshev class is removed from RxnRates.h after Cantera 2.6
    Chebyshev::getParameters(rateNode, rate_units);
}

const Array2D& ChebyshevRate3::coeffs() const
{
    return ChebyshevRate3::coeffs2D();
}

void ChebyshevRate3::setCoeffs(const Array2D& coeffs)
{
    Chebyshev::setCoeffs(coeffs);
    if (m_evaluator) {
        dynamic_cast<ChebyshevRate3&>(
            m_evaluator->rate(m_index)).setCoeffs(coeffs);
    }
}

void ChebyshevRate3::validate(const std::string& equation)
{
}

CustomFunc1Rate::CustomFunc1Rate() : m_ratefunc(0) {}

void CustomFunc1Rate::setRateFunction(shared_ptr<Func1> f)
{
    m_ratefunc = f;
}

double CustomFunc1Rate::eval(const CustomFunc1Data& shared_data,
                             double concm) const
{
    if (m_ratefunc) {
        return m_ratefunc->eval(shared_data.m_temperature);
    }
    return NAN;
}

}
