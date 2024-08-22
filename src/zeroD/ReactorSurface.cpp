//! @file ReactorSurface.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/ReactorSurface.h"
#include "cantera/zeroD/ReactorNet.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/kinetics/InterfaceKinetics.h"
#include "cantera/kinetics/Reaction.h"
#include "cantera/base/Solution.h"

namespace Cantera
{

ReactorSurface::ReactorSurface(shared_ptr<Solution> sol, const string& name)
    : ReactorNode(name)
{
    if (!sol || !(sol->thermo())) {
        warn_deprecated("ReactorSurface::ReactorSurface",
            "Creation of empty reactor objects is deprecated in Cantera 3.1 and will "
            "raise\nexceptions thereafter; reactor contents should be provided in the "
            "constructor.");
        return;
    } else if (!std::dynamic_pointer_cast<SurfPhase>(sol->thermo())) {
        throw CanteraError("ReactorSurface::ReactorSurface",
            "Solution object must have a SurfPhase object as the thermo manager.");
    }
    m_solution = sol;

    if (!sol->kinetics() ) {
        throw CanteraError("ReactorSurface::ReactorSurface",
            "Solution object must have kinetics manager.");
    } else if (!std::dynamic_pointer_cast<InterfaceKinetics>(sol->kinetics())) {
        throw CanteraError("ReactorSurface::ReactorSurface",
            "Kinetics manager must be an InterfaceKinetics object.");
    }
    // todo: move all member variables to use shared pointers after Cantera 3.1
    m_kinetics = m_solution->kinetics().get();
    m_solution->thermo()->addSpeciesLock();

    m_thermo = dynamic_cast<SurfPhase*>(&m_kinetics->thermo(0));
    if (m_thermo == nullptr) {
        // This should never happen, as the Solution should already contain the correct
        // SurfPhase object
        throw CanteraError("ReactorSurface::ReactorSurface",
            "Specified kinetics manager does not represent a surface "
            "kinetics mechanism.");
    }
    m_cov.resize(m_thermo->nSpecies());
    m_thermo->getCoverages(m_cov.data());
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

    m_thermo = dynamic_cast<SurfPhase*>(&kin->thermo(0));
    if (m_thermo == nullptr) {
        throw CanteraError("ReactorSurface::setKinetics",
            "Specified kinetics manager does not represent a surface "
            "kinetics mechanism.");
    }
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

void ReactorSurface::setCoverages(const string& cov)
{
    m_thermo->setCoveragesByName(cov);
    m_thermo->getCoverages(m_cov.data());
}

void ReactorSurface::getCoverages(double* cov) const
{
    copy(m_cov.begin(), m_cov.end(), cov);
}

void ReactorSurface::syncState()
{
    m_thermo->setTemperature(m_reactor->temperature());
    m_thermo->setCoveragesNoNorm(m_cov.data());
}

void ReactorSurface::addSensitivityReaction(size_t i)
{
    if (i >= m_kinetics->nReactions()) {
        throw CanteraError("ReactorSurface::addSensitivityReaction",
                           "Reaction number out of range ({})", i);
    }
    size_t p = m_reactor->network().registerSensitivityParameter(
        m_kinetics->reaction(i)->equation(), 1.0, 1.0);
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
