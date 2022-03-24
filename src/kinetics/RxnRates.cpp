//! @file RxnRates.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/RxnRates.h"
#include "cantera/base/AnyMap.h"

namespace Cantera
{
Arrhenius2::Arrhenius2()
    : ArrheniusRate()
{
    m_b = 0.0;
    m_A = 0.0;
    m_logA = -1.0E300;
}

Arrhenius2::Arrhenius2(doublereal A, doublereal b, doublereal E)
    : ArrheniusRate(A, b, E * GasConstant)
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

Arrhenius2::Arrhenius2(const ArrheniusRate& other)
    : ArrheniusRate(other.preExponentialFactor(),
                 other.temperatureExponent(),
                 other.activationEnergy())
{
}

void Arrhenius2::setRateParameters(const AnyValue& rate,
                                   const UnitSystem& units, const Units& rate_units)
{
    UnitStack units_stack(rate_units);
    ArrheniusRate::setRateParameters(rate, units, units_stack);
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

}
