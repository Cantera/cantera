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
    // TODO: After Cantera 3.2, raise an exception of 'reactors' is empty
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
    setThermo(*m_solution->thermo());
    if (!std::dynamic_pointer_cast<SurfPhase>(soln->thermo())) {
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
    // todo: move all member variables to use shared pointers after Cantera 3.2
    m_kinetics = m_solution->kinetics().get();
    m_thermo = m_solution->thermo().get();
    m_surf = dynamic_cast<SurfPhase*>(m_thermo);
    m_cov.resize(m_surf->nSpecies());
    m_surf->getCoverages(m_cov.data());
}

ReactorSurface::ReactorSurface(shared_ptr<Solution> sol, const string& name)
    : ReactorBase(sol, false, name)
{
    if (!std::dynamic_pointer_cast<SurfPhase>(sol->thermo())) {
        throw CanteraError("ReactorSurface::ReactorSurface",
            "Solution object must have a SurfPhase object as the thermo manager.");
    }

    if (!sol->kinetics() ) {
        throw CanteraError("ReactorSurface::ReactorSurface",
            "Solution object must have kinetics manager.");
    } else if (!std::dynamic_pointer_cast<InterfaceKinetics>(sol->kinetics())) {
        throw CanteraError("ReactorSurface::ReactorSurface",
            "Kinetics manager must be an InterfaceKinetics object.");
    }
    // todo: move all member variables to use shared pointers after Cantera 3.2
    m_kinetics = m_solution->kinetics().get();
    m_thermo = m_solution->thermo().get();
    m_surf = dynamic_cast<SurfPhase*>(m_thermo);
    m_cov.resize(m_surf->nSpecies());
    m_surf->getCoverages(m_cov.data());
}

ReactorSurface::ReactorSurface(shared_ptr<Solution> sol, bool clone, const string& name)
    : ReactorSurface(sol, name)
{
    // TODO: Remove after Cantera 3.2 when removing ReactorSurface from ReactorFactory
    if (clone) {
        throw CanteraError("ReactorSurface::ReactorSurface", "Bad constructor arguments."
            " When clone=true, list of adjacent reactors must be provided");
    }
}

double ReactorSurface::area() const
{
    return m_area;
}

void ReactorSurface::setArea(double a)
{
    m_area = a;
}

void ReactorSurface::setKinetics(Kinetics* kin)
{
    warn_deprecated("ReactorSurface::setKinetics",
                    "To be removed after Cantera 3.2.");
    m_kinetics = kin;
    if (kin == nullptr) {
        m_surf = nullptr;
        return;
    }

    m_surf = dynamic_cast<SurfPhase*>(&kin->thermo(0));
    if (m_surf == nullptr) {
        throw CanteraError("ReactorSurface::setKinetics",
            "Specified kinetics manager does not represent a surface "
            "kinetics mechanism.");
    }
    m_cov.resize(m_surf->nSpecies());
    m_surf->getCoverages(m_cov.data());
}

void ReactorSurface::setKinetics(Kinetics& kin)
{
    setKinetics(&kin);
}

void ReactorSurface::setReactor(ReactorBase* reactor)
{
    if (std::find(m_reactors.begin(), m_reactors.end(), reactor) != m_reactors.end()) {
        // Ignore case where reactor has already been added in the opposite direction
        return;
    }
    warn_deprecated("ReactorSurface::setReactor", "To be removed after Cantera 3.2. "
        "Superseded by constructor taking a list of adjacent reactors.");
    m_reactors.resize(1);
    m_reactors[0] = reactor;
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
    if (m_reactors.empty()) {
        // TODO: Remove this check after Cantera 3.2, since it will no longer be
        // possible to create the ReactorSurface without specifying the adjacent
        // reactors in the constructor.
        throw CanteraError("ReactorSurface::syncState",
                           "Surface is not installed on any Reactor");
    }
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
