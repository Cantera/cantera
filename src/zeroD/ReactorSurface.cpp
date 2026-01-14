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
        m_reactors.push_back(R);
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
    m_kinetics = std::dynamic_pointer_cast<InterfaceKinetics>(m_solution->kinetics());
    m_thermo = m_surf.get();
    m_nsp = m_nv = m_surf->nSpecies();
    m_sdot.resize(m_kinetics->nTotalSpecies(), 0.0);
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
    // Sync the surface temperature and pressure to that of the first adjacent reactor
    m_thermo->setState_TP(m_reactors[0]->temperature(), m_reactors[0]->pressure());
}

vector<size_t> ReactorSurface::initializeSteady()
{
    return {0}; // sum of coverages constraint
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

void ReactorSurface::evalSteady(double t, double* LHS, double* RHS)
{
    eval(t, LHS, RHS);
    vector<double> cov(m_nsp);
    m_surf->getCoverages(cov.data());
    double sum = 0.0;
    for (size_t k = 0; k < m_nsp; k++) {
        sum += cov[k];
    }
    RHS[0] = 1.0 - sum;
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

double ReactorSurface::upperBound(size_t k) const
{
    return 1.0;
}

double ReactorSurface::lowerBound(size_t k) const
{
    return -Tiny;
}

void ReactorSurface::resetBadValues(double* y)
{
    for (size_t k = 0; k < m_nsp; k++) {
        y[k] = std::max(y[k], 0.0);
    }
    m_surf->setCoverages(y);
    m_surf->getCoverages(y);
}

// ------ MoleReactorSurface methods ------

void MoleReactorSurface::initialize(double t0)
{
    ReactorSurface::initialize(t0);
    m_cov_tmp.resize(m_nsp);
    m_f_energy.resize(m_kinetics->nTotalSpecies(), 0.0);
    m_f_species.resize(m_kinetics->nTotalSpecies(), 0.0);
    m_kin2net.resize(m_kinetics->nTotalSpecies(), -1);
    m_kin2reactor.resize(m_kinetics->nTotalSpecies(), nullptr);

    for (size_t k = 0; k < m_nsp; k++) {
        m_kin2net[k] = static_cast<int>(m_offset + k);
    }
    for (auto R : m_reactors) {
        size_t nsp = R->phase()->thermo()->nSpecies();
        size_t k0 = m_kinetics->speciesOffset(*R->phase()->thermo());
        int offset = static_cast<int>(R->offset() + R->speciesOffset());
        for (size_t k = 0; k < nsp; k++) {
            m_kin2net[k + k0] = offset + k;
            m_kin2reactor[k + k0] = R.get();
        }
    }
}

void MoleReactorSurface::getState(double* y)
{
    m_surf->getCoverages(y);
    double totalSites = m_surf->siteDensity() * m_area;
    for (size_t k = 0; k < m_nsp; k++) {
        y[k] *= totalSites / m_surf->size(k);
    }
}

void MoleReactorSurface::updateState(double* y)
{
    std::copy(y, y + m_nsp, m_cov_tmp.data());
    double totalSites = m_surf->siteDensity() * m_area;
    for (size_t k = 0; k < m_nsp; k++) {
        m_cov_tmp[k] *= m_surf->size(k) / totalSites;
    }
    m_surf->setCoveragesNoNorm(m_cov_tmp.data());
    m_thermo->setState_TP(m_reactors[0]->temperature(), m_reactors[0]->pressure());
    m_kinetics->getNetProductionRates(m_sdot.data());
}

void MoleReactorSurface::eval(double t, double* LHS, double* RHS)
{
    for (size_t k = 0; k < m_nsp; k++) {
        RHS[k] = m_sdot[k] * m_area / m_surf->size(k);
    }
}

void MoleReactorSurface::evalSteady(double t, double* LHS, double* RHS)
{
    eval(t, LHS, RHS);
    double cov_sum = 0.0;
    for (size_t k = 0; k < m_nsp; k++) {
        cov_sum += m_cov_tmp[k];
    }
    RHS[0] = 1.0 - cov_sum;
}

double MoleReactorSurface::upperBound(size_t k) const
{
    return BigNumber;
}

double MoleReactorSurface::lowerBound(size_t k) const
{
    return -Tiny;
}

void MoleReactorSurface::resetBadValues(double* y)
{
    for (size_t k = 0; k < m_nsp; k++) {
        y[k] = std::max(y[k], 0.0);
    }
}

void MoleReactorSurface::getJacobianElements(vector<Eigen::Triplet<double>>& trips)
{
    auto sdot_ddC = m_kinetics->netProductionRates_ddCi();
    for (auto R : m_reactors) {
        double f_species;
        size_t nsp_R = R->phase()->thermo()->nSpecies();
        size_t k0 = m_kinetics->speciesOffset(*R->phase()->thermo());
        R->getJacobianScalingFactors(f_species, &m_f_energy[k0]);
        std::fill(&m_f_species[k0], &m_f_species[k0 + nsp_R], m_area * f_species);
    }
    std::fill(&m_f_species[0], &m_f_species[m_nsp], 1.0); // surface species
    for (int k = 0; k < sdot_ddC.outerSize(); k++) {
        int col = m_kin2net[k];
        if (col == -1) {
            continue;
        }
        for (Eigen::SparseMatrix<double>::InnerIterator it(sdot_ddC, k); it; ++it) {
            int row = m_kin2net[it.row()];
            if (row == -1 || it.value() == 0.0) {
                continue;
            }
            ReactorBase* R = m_kin2reactor[it.row()];
            trips.emplace_back(row, col, it.value() * m_f_species[k]);
            if (m_f_energy[it.row()] != 0.0) {
                trips.emplace_back(R->offset(), col,
                                   it.value() * m_f_energy[it.row()] * m_f_species[k]);
            }
        }
    }
}

// ------ FlowReactorSurface methods ------

FlowReactorSurface::FlowReactorSurface(shared_ptr<Solution> soln,
                                       const vector<shared_ptr<ReactorBase>>& reactors,
                                       bool clone,
                                       const string& name)
    : ReactorSurface(soln, reactors, clone, name)
{
    m_area = -1.0; // default to perimeter of cylindrical reactor
}

double FlowReactorSurface::area() const {
    if (m_area > 0) {
        return m_area;
    }

    // Assuming a cylindrical cross section, P = 2 * pi * r and A = pi * r^2, so
    // P(A) = 2 * sqrt(pi * A)
    return 2.0 * sqrt(Pi * m_reactors[0]->area());
}

void FlowReactorSurface::evalDae(double t, double* y, double* ydot, double* residual)
{
    size_t nsp = m_surf->nSpecies();
    double sum = y[0];
    for (size_t k = 1; k < nsp; k++) {
        residual[k] = m_sdot[k];
        sum += y[k];
    }
    residual[0] = sum - 1.0;
}

void FlowReactorSurface::getStateDae(double* y, double* ydot)
{
    // Advance the surface to steady state to get consistent initial coverages
    m_kinetics->advanceCoverages(100.0, m_ss_rtol, m_ss_atol, 0, m_max_ss_steps,
                                 m_max_ss_error_fails);
    getCoverages(y);
    // Update the values in m_sdot that are needed by the adjacent FlowReactor
    updateState(y);
}

void FlowReactorSurface::getConstraints(double* constraints)
{
    // The species coverages are algebraic constraints
    std::fill(constraints, constraints + m_nsp, 1.0);
}

}
