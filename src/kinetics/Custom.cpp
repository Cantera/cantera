//! @file Custom.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/Custom.h"
#include "cantera/numerics/Func1.h"

namespace Cantera
{

CustomFunc1Rate::CustomFunc1Rate()
    : m_ratefunc(0)
{
}

void CustomFunc1Rate::setRateFunction(shared_ptr<Func1> f)
{
    m_ratefunc = f;
}

double CustomFunc1Rate::evalFromStruct(const ArrheniusData& shared_data) const
{
    if (m_ratefunc) {
        return m_ratefunc->eval(shared_data.temperature);
    }
    return NAN;
}

void CustomFunc1Rate::getParameters(AnyMap& rateNode, const Units& rate_units) const
{
    throw NotImplementedError("CustomFunc1Rate::getParameters",
                              "Not implemented by '{}' object.", type());
}

}