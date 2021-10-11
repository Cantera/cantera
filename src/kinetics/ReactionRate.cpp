//! @file ReactionRate.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/ReactionRate.h"
#include "cantera/kinetics/MultiRate.h"
#include "cantera/base/Array.h"
#include "cantera/numerics/Func1.h"

namespace Cantera
{

PlogRate::PlogRate(const std::multimap<double, Arrhenius>& rates)
    : Plog(rates)
{
}

PlogRate::PlogRate(const AnyMap& node, const UnitsVector& units)
{
    setParameters(node, units);
}

unique_ptr<MultiRateBase> PlogRate::newMultiRate() const
{
    return unique_ptr<MultiRateBase>(new MultiBulkRate<PlogRate, PlogData>);
}

void PlogRate::setParameters(const AnyMap& node, const UnitsVector& units)
{
    // @TODO  implementation of Plog::setParameters should be transferred here
    //     when the Plog class is removed from RxnRates.h after Cantera 2.6
    ReactionRateBase::setParameters(node, units);
    auto rate_units = Units::product(units);
    if (!node.hasKey("rate-constants")) {
        Plog::setParameters(std::vector<AnyMap> (), node.units(), rate_units);
        return;
    }

    Plog::setParameters(node.at("rate-constants").asVector<AnyMap>(),
                        node.units(), rate_units);
}

void PlogRate::getParameters(AnyMap& rateNode) const
{
    // @TODO  implementation of Plog::getParameters should be transferred here
    //     when the Plog class is removed from RxnRates.h after Cantera 2.6
    Plog::getParameters(rateNode);
    rateNode["type"] = type();
}

ChebyshevRate3::ChebyshevRate3(double Tmin, double Tmax, double Pmin, double Pmax,
                               const Array2D& coeffs)
    : Chebyshev(Tmin, Tmax, Pmin, Pmax, coeffs)
{
}

ChebyshevRate3::ChebyshevRate3(const AnyMap& node, const UnitsVector& units)
{
    setParameters(node, units);
}

unique_ptr<MultiRateBase> ChebyshevRate3::newMultiRate() const
{
    return unique_ptr<MultiRateBase>(
        new MultiBulkRate<ChebyshevRate3, ChebyshevData>);
}

void ChebyshevRate3::setParameters(const AnyMap& node, const UnitsVector& units)
{
    ReactionRateBase::setParameters(node, units);
    m_rate_units = Units::product(units);
    if (!node.hasKey("data")) {
        Chebyshev::setParameters(AnyMap(), node.units(), m_rate_units);
        return;
    }
    // @TODO  implementation of Chebyshev::setParameters should be transferred here
    //     when the Chebyshev class is removed from RxnRates.h after Cantera 2.6
    Chebyshev::setParameters(node, node.units(), m_rate_units);
}

void ChebyshevRate3::getParameters(AnyMap& rateNode) const
{
    // @TODO  implementation of Chebyshev::getParameters should be transferred here
    //     when the Chebyshev class is removed from RxnRates.h after Cantera 2.6
    Chebyshev::getParameters(rateNode, m_rate_units);
    rateNode["type"] = type();
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

double CustomFunc1Rate::eval(const CustomFunc1Data& shared_data) const
{
    if (m_ratefunc) {
        return m_ratefunc->eval(shared_data.temperature);
    }
    return NAN;
}

}
