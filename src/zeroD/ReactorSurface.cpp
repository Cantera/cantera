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
    m_nv = m_surf->nSpecies();
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
    m_surf->setCoveragesNoNorm(cov);
}

void ReactorSurface::setCoverages(const Composition& cov)
{
    m_surf->setCoveragesByName(cov);
}

void ReactorSurface::setCoverages(const string& cov)
{
    m_surf->setCoveragesByName(cov);
}

void ReactorSurface::getCoverages(double* cov) const
{
    m_surf->getCoverages(cov);
}

void ReactorSurface::getState(double* y)
{
    m_surf->getCoverages(y);
}

void ReactorSurface::initialize(double t0)
{
    m_sdot.resize(m_kinetics->nTotalSpecies(), 0.0);
    // Sync the surface temperature and pressure to that of the first adjacent reactor
    m_thermo->setState_TP(m_reactors[0]->temperature(), m_reactors[0]->pressure());
}

void ReactorSurface::updateState(double* y)
{
    m_surf->setCoveragesNoNorm(y);
    m_thermo->setState_TP(m_reactors[0]->temperature(), m_reactors[0]->pressure());
    m_kinetics->getNetProductionRates(m_sdot.data());
}

void ReactorSurface::eval(double t, double* LHS, double* RHS)
{
    size_t nsp = m_surf->nSpecies();
    double rs0 = 1.0 / m_surf->siteDensity();
    double sum = 0.0;
    for (size_t k = 1; k < nsp; k++) {
        RHS[k] = m_sdot[k] * rs0 * m_surf->size(k);
        sum -= RHS[k];
    }
    RHS[0] = sum;
}

void ReactorSurface::applySensitivity(double* params)
{
    if (!params) {
        return;
    }
    for (auto& p : m_sensParams) {
        if (p.type == SensParameterType::reaction) {
            p.value = m_kinetics->multiplier(p.local);
            m_kinetics->setMultiplier(p.local, p.value*params[p.global]);
        } else if (p.type == SensParameterType::enthalpy) {
            m_thermo->modifyOneHf298SS(p.local, p.value + params[p.global]);
        }
    }
    m_thermo->invalidateCache();
    m_kinetics->invalidateCache();
}

void ReactorSurface::resetSensitivity(double* params)
{
    if (!params) {
        return;
    }
    for (auto& p : m_sensParams) {
        if (p.type == SensParameterType::reaction) {
            m_kinetics->setMultiplier(p.local, p.value);
        } else if (p.type == SensParameterType::enthalpy) {
            m_thermo->resetHf298(p.local);
        }
    }
    m_thermo->invalidateCache();
    m_kinetics->invalidateCache();
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

size_t ReactorSurface::componentIndex(const string& nm) const
{
    return m_surf->speciesIndex(nm);
}

string ReactorSurface::componentName(size_t k)
{
    return m_surf->speciesName(k);
}

}
