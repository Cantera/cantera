//! @file ReactorSurface.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/ReactorSurface.h"
#include "cantera/zeroD/ReactorNet.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/kinetics/Kinetics.h"

namespace Cantera
{

ReactorSurface::ReactorSurface()
    : m_area(1.0)
    , m_thermo(nullptr)
    , m_kinetics(nullptr)
    , m_reactor(nullptr)
{
}

double ReactorSurface::area() const
{
    return m_area;
}

void ReactorSurface::setArea(double a)
{
    m_area = a;
}

void ReactorSurface::setKinetics(Kinetics* kin) {
    m_kinetics = kin;
    if (kin == nullptr) {
        m_thermo = nullptr;
        return;
    }

    size_t i = kin->surfacePhaseIndex();
    if (i == npos) {
        throw CanteraError("ReactorSurface::setKinetics",
            "Specified surface kinetics manager does not represent a surface "
            "kinetics mechanism.");
    }
    m_thermo = dynamic_cast<SurfPhase*>(&kin->thermo(i));
    m_cov.resize(m_thermo->nSpecies());
    m_thermo->getCoverages(m_cov.data());
}

void ReactorSurface::setReactor(ReactorBase* reactor)
{
    m_reactor = reactor;
}

void ReactorSurface::setCoverages(const double* cov)
{
    copy(cov, cov + m_cov.size(), m_cov.begin());
}

void ReactorSurface::setCoverages(const Composition& cov)
{
    m_thermo->setCoveragesByName(cov);
    m_thermo->getCoverages(m_cov.data());
}

void ReactorSurface::setCoverages(const std::string& cov)
{
    m_thermo->setCoveragesByName(cov);
    m_thermo->getCoverages(m_cov.data());
}

void ReactorSurface::getCoverages(double* cov) const
{
    copy(m_cov.begin(), m_cov.end(), cov);
}

void ReactorSurface::syncCoverages()
{
    m_thermo->setCoveragesNoNorm(m_cov.data());
}

void ReactorSurface::addSensitivityReaction(size_t i)
{
    if (i >= m_kinetics->nReactions()) {
        throw CanteraError("ReactorSurface::addSensitivityReaction",
                           "Reaction number out of range ({})", i);
    }
    size_t p = m_reactor->network().registerSensitivityParameter(
        m_kinetics->reactionString(i), 1.0, 1.0);
    m_params.emplace_back(
        SensitivityParameter{i, p, 1.0, SensParameterType::reaction});
}

void ReactorSurface::setSensitivityParameters(const double* params)
{
    for (auto& p : m_params) {
        p.value = m_kinetics->multiplier(p.local);
        m_kinetics->setMultiplier(p.local, p.value*params[p.global]);
    }
}

void ReactorSurface::resetSensitivityParameters()
{
    for (auto& p : m_params) {
        m_kinetics->setMultiplier(p.local, p.value);
    }
}

}
