//! @file ReactionRate.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/ReactionRate.h"
#include "cantera/kinetics/MultiRate.h"
#include "cantera/base/Array.h"
#include "cantera/numerics/Func1.h"

namespace Cantera
{

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
