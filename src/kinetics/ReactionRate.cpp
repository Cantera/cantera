//! @file ReactionRate.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/ReactionRate.h"
#include "cantera/numerics/Func1.h"
#include "cantera/thermo/ThermoPhase.h"

namespace Cantera
{

void CustomFunc1Data::update(const ThermoPhase& bulk) {
    m_temperature = bulk.temperature();
}

CustomFunc1Rate::CustomFunc1Rate() : m_ratefunc(0) {}

void CustomFunc1Rate::setRateFunction(shared_ptr<Func1> f) {
    m_ratefunc = f;
}

double CustomFunc1Rate::eval(const CustomFunc1Data& shared_data) const {
    if (m_ratefunc) {
        return m_ratefunc->eval(shared_data.m_temperature);
    }
    throw CanteraError("CustomFunc1Rate::eval",
                       "Custom rate function is not initialized.");
}

void ArrheniusData::update(const ThermoPhase& bulk) {
    m_temperature = bulk.temperature();
    m_logT = std::log(m_temperature);
    m_recipT = 1./m_temperature;
}

ArrheniusRate::ArrheniusRate(double A, double b, double E)
    : Arrhenius(A, b, E) {
}

ArrheniusRate::ArrheniusRate(const AnyMap& node, const Units& rc_units) {
    setup(node["rate-constant"], node.units(), rc_units);
}

}
