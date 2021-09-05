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

void ReactionRateBase::setUnits(const Units& rate_units)
{
    units = rate_units;
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

unique_ptr<MultiRateBase> ArrheniusRate::newMultiRate() const
{
    return unique_ptr<MultiRateBase>(new MultiBulkRate<ArrheniusRate, ArrheniusData>);
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

void ArrheniusRate::validate(const std::string& equation)
{
    if (!allow_negative_pre_exponential_factor && preExponentialFactor() < 0) {
        throw CanteraError("ArrheniusRate::validate",
            "Undeclared negative pre-exponential factor found in reaction '"
            + equation + "'");
    }
}

PlogRate::PlogRate(const std::multimap<double, Arrhenius>& rates)
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

unique_ptr<MultiRateBase> PlogRate::newMultiRate() const
{
    return unique_ptr<MultiRateBase>(new MultiBulkRate<PlogRate, PlogData>);
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
    rateNode["type"] = type();
}

ChebyshevRate3::ChebyshevRate3(
    const std::pair<double, double> Trange,
    const std::pair<double, double> Prange,
    const Array2D& coeffs) : Chebyshev(Trange, Prange, coeffs)
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

unique_ptr<MultiRateBase> ChebyshevRate3::newMultiRate() const
{
    return unique_ptr<MultiRateBase>(
        new MultiBulkRate<ChebyshevRate3, ChebyshevData>);
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
    rateNode["type"] = type();
}

void ChebyshevRate3::validate(const std::string& equation)
{
}

CustomFunc1Rate::CustomFunc1Rate() : m_ratefunc(0) {}

unique_ptr<MultiRateBase> CustomFunc1Rate::newMultiRate() const
{
    return unique_ptr<MultiRateBase>(
        new MultiBulkRate<CustomFunc1Rate, CustomFunc1Data>);
}

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
