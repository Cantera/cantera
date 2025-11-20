//! @file ReactorSurface.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/ReactorSurface.h"
#include "cantera/zeroD/ReactorNet.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/kinetics/InterfaceKinetics.h"
#include "cantera/kinetics/Reaction.h"

namespace Cantera
{

ReactorSurface::ReactorSurface(shared_ptr<Solution> soln,
                               const vector<shared_ptr<ReactorBase>>& reactors,
                               bool clone,
                               const string& name)
    : ReactorBase(name)
{
    if (reactors.empty()) {
        throw CanteraError("ReactorSurface::ReactorSurface",
            "Surface requires at least one adjacent reactor.");
    }
    vector<shared_ptr<Solution>> adjacent;
    for (auto R : reactors) {
        adjacent.push_back(R->phase());
        m_reactors.push_back(R.get());
        R->addSurface(this);
    }
    if (clone) {
        m_solution = soln->clone(adjacent, true, false);
    } else {
        m_solution = soln;
    }
    m_solution->thermo()->addSpeciesLock();
    m_surf = std::dynamic_pointer_cast<SurfPhase>(m_solution->thermo());
    if (!m_surf) {
        throw CanteraError("ReactorSurface::ReactorSurface",
            "Solution object must have a SurfPhase object as the thermo manager.");
    }

    if (!soln->kinetics() ) {
        throw CanteraError("ReactorSurface::ReactorSurface",
            "Solution object must have kinetics manager.");
    } else if (!std::dynamic_pointer_cast<InterfaceKinetics>(soln->kinetics())) {
        throw CanteraError("ReactorSurface::ReactorSurface",
            "Kinetics manager must be an InterfaceKinetics object.");
    }
    m_kinetics = m_solution->kinetics();
    m_thermo = m_surf.get();
    m_cov.resize(m_surf->nSpecies());
    m_surf->getCoverages(m_cov.data());
}

double ReactorSurface::area() const
{
    return m_area;
}

void ReactorSurface::setArea(double a)
{
    m_area = a;
}

void ReactorSurface::setCoverages(const double* cov)
{
    copy(cov, cov + m_cov.size(), m_cov.begin());
}

void ReactorSurface::setCoverages(const Composition& cov)
{
    m_surf->setCoveragesByName(cov);
    m_surf->getCoverages(m_cov.data());
}

void ReactorSurface::setCoverages(const string& cov)
{
    m_surf->setCoveragesByName(cov);
    m_surf->getCoverages(m_cov.data());
}

void ReactorSurface::getCoverages(double* cov) const
{
    copy(m_cov.begin(), m_cov.end(), cov);
}

void ReactorSurface::restoreState()
{
    m_surf->setTemperature(m_reactors[0]->temperature());
    m_surf->setCoveragesNoNorm(m_cov.data());
}

void ReactorSurface::syncState()
{
    warn_user("ReactorSurface::syncState", "Behavior changed in Cantera 3.2 for "
        "consistency with ReactorBase. To set SurfPhase state from ReactorSurface "
        "object, use restoreState().");
    m_surf->getCoverages(m_cov.data());
}

void ReactorSurface::addSensitivityReaction(size_t rxn)
{
    if (rxn >= m_kinetics->nReactions()) {
        throw CanteraError("ReactorSurface::addSensitivityReaction",
                           "Reaction number out of range ({})", rxn);
    }
    size_t p = m_reactors[0]->network().registerSensitivityParameter(
        m_kinetics->reaction(rxn)->equation(), 1.0, 1.0);
    m_sensParams.emplace_back(
        SensitivityParameter{rxn, p, 1.0, SensParameterType::reaction});
}

void ReactorSurface::setSensitivityParameters(const double* params)
{
    for (auto& p : m_sensParams) {
        p.value = m_kinetics->multiplier(p.local);
        m_kinetics->setMultiplier(p.local, p.value*params[p.global]);
    }
}

void ReactorSurface::resetSensitivityParameters()
{
    for (auto& p : m_sensParams) {
        m_kinetics->setMultiplier(p.local, p.value);
    }
}

}
