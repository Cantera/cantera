//! @file IdealGasMoleReactor.cpp A constant volume zero-dimensional
//! reactor with moles as the state

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/IdealGasMoleReactor.h"
#include "cantera/zeroD/FlowDevice.h"
#include "cantera/zeroD/ReactorNet.h"
#include "cantera/zeroD/ReactorSurface.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/base/utilities.h"
#include "cantera/numerics/eigen_dense.h"
#include <limits>

namespace Cantera
{

void IdealGasMoleReactor::getState(span<double> y)
{
    // get mass for calculations
    m_mass = m_thermo->density() * m_vol;

    // set the first component to the temperature
    y[0] = m_thermo->temperature();

    // set the second component to the volume
    y[1] = m_vol;

    // get moles of species in remaining state
    getMoles(y.subspan(m_sidx));
}

size_t IdealGasMoleReactor::componentIndex(const string& nm) const
{
    if (nm == "temperature") {
        return 0;
    }
    if (nm == "volume") {
        return 1;
    }
    try {
        return m_thermo->speciesIndex(nm) + m_sidx;
    } catch (const CanteraError&) {
        throw CanteraError("IdealGasMoleReactor::componentIndex",
            "Component '{}' not found", nm);
    }
}

string IdealGasMoleReactor::componentName(size_t k)
{
    if (k == 0) {
        return "temperature";
    } else {
        return MoleReactor::componentName(k);
    }
}

void IdealGasMoleReactor::initialize(double t0)
{
    MoleReactor::initialize(t0);
    m_uk.resize(m_nsp, 0.0);
}

double IdealGasMoleReactor::upperBound(size_t k) const {
    if (k == 0) {
        //@todo: Revise pending resolution of https://github.com/Cantera/enhancements/issues/229
        return 1.5 * m_thermo->maxTemp();
    } else {
        return MoleReactor::upperBound(k);
    }
}

double IdealGasMoleReactor::lowerBound(size_t k) const {
    if (k == 0) {
        //@todo: Revise pending resolution of https://github.com/Cantera/enhancements/issues/229
        return 0.5 * m_thermo->minTemp();
    } else {
        return MoleReactor::lowerBound(k);
    }
}

vector<size_t> IdealGasMoleReactor::initializeSteady()
{
    m_initialVolume = m_vol;
    m_initialTemperature = m_thermo->temperature();
    if (energyEnabled()) {
        return {1}; // volume
    } else {
        return {0, 1}; // temperature and volume
    }
}

void IdealGasMoleReactor::updateState(span<const double> y)
{
    // the components of y are: [0] the temperature, [1] the volume, [2...K+1) are the
    // moles of each species, and [K+1...] are the moles of surface
    // species on each wall.
    setMassFromMoles(y.subspan(m_sidx));
    m_vol = y[1];
    // set state
    m_thermo->setMolesNoTruncate(y.subspan(m_sidx, m_nsp));
    m_thermo->setState_TD(y[0], m_mass / m_vol);
    m_thermo->getPartialMolarIntEnergies_TV(m_uk);
    m_TotalCv = m_mass * m_thermo->cv_mass();
    updateConnected(true);
}

void IdealGasMoleReactor::eval(double time, span<double> LHS, span<double> RHS)
{
    double& mcvdTdt = RHS[0]; // m * c_v * dT/dt
    auto dndt = RHS.subspan(m_sidx); // kmol per s

    evalWalls(time);
    updateSurfaceProductionRates();
    auto imw = m_thermo->inverseMolecularWeights();

    if (m_chem) {
        m_kin->getNetProductionRates(m_wdot); // "omega dot"
    }

    // external heat transfer and compression work
    mcvdTdt += - (m_pressure + m_thermo->internalPressure()) * m_vdot + m_Qdot;

    if (m_energy) {
        mcvdTdt += m_thermo->intrinsicHeating() * m_vol;
    }

    for (size_t n = 0; n < m_nsp; n++) {
        // heat release from gas phase and surface reactions
        mcvdTdt -= m_wdot[n] * m_uk[n] * m_vol;
        mcvdTdt -= m_sdot[n] * m_uk[n];
        // production in gas phase and from surfaces
        dndt[n] = (m_wdot[n] * m_vol + m_sdot[n]);
    }

    // add terms for outlets
    for (auto outlet : m_outlet) {
        for (size_t n = 0; n < m_nsp; n++) {
            // flow of species into system and dilution by other species
            dndt[n] -= outlet->outletSpeciesMassFlowRate(n) * imw[n];
        }
        double mdot = outlet->massFlowRate();
        mcvdTdt -= mdot * m_pressure * m_vol / m_mass; // flow work
    }

    // add terms for inlets
    for (auto inlet : m_inlet) {
        double mdot = inlet->massFlowRate();
        mcvdTdt += inlet->enthalpy_mass() * mdot;
        for (size_t n = 0; n < m_nsp; n++) {
            double mdot_spec = inlet->outletSpeciesMassFlowRate(n);
            // flow of species into system and dilution by other species
            dndt[n] += inlet->outletSpeciesMassFlowRate(n) * imw[n];
            mcvdTdt -= m_uk[n] * imw[n] * mdot_spec;
        }
    }

    RHS[1] = m_vdot;
    if (m_energy) {
        LHS[0] = m_TotalCv;
    } else {
        RHS[0] = 0;
    }
}

void IdealGasMoleReactor::evalSteady(double t, span<double> LHS, span<double> RHS)
{
    eval(t, LHS, RHS);
    if (!energyEnabled()) {
        RHS[0] = m_thermo->temperature() - m_initialTemperature;
    }
    RHS[1] = m_vol - m_initialVolume;
}

void IdealGasMoleReactor::getJacobianElements(vector<Eigen::Triplet<double>>& trips)
{
    // dnk_dnj represents d(dot(n_k)) / d (n_j) but is first assigned as
    // d (dot(omega)) / d c_j, it is later transformed appropriately.
    Eigen::SparseMatrix<double> dnk_dnj = m_kin->netProductionRates_ddCi();

    // add species to species derivatives  elements to the jacobian
    // calculate ROP derivatives, excluding the terms -n_i / (V * N) dc_i/dn_j
    // as it substantially reduces matrix sparsity
    for (int k = 0; k < dnk_dnj.outerSize(); k++) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(dnk_dnj, k); it; ++it) {
            trips.emplace_back(static_cast<int>(it.row() + m_offset + m_sidx),
                static_cast<int>(it.col() + m_offset + m_sidx), it.value());
        }
    }

    // Temperature Derivatives
    if (m_energy) {
        // getting perturbed state for finite difference
        double deltaTemp = m_thermo->temperature()
            * std::sqrt(std::numeric_limits<double>::epsilon());
        // finite difference temperature derivatives
        vector<double> lhsPerturbed(m_nv, 1.0), lhsCurrent(m_nv, 1.0);
        vector<double> rhsPerturbed(m_nv, 0.0), rhsCurrent(m_nv, 0.0);
        vector<double> yCurrent(m_nv);
        getState(yCurrent);
        vector<double> yPerturbed = yCurrent;
        // perturb temperature
        yPerturbed[0] += deltaTemp;
        // getting perturbed state
        updateState(yPerturbed);
        double time = (m_net != nullptr) ? m_net->time() : 0.0;
        eval(time, lhsPerturbed, rhsPerturbed);
        // reset and get original state
        updateState(yCurrent);
        eval(time, lhsCurrent, rhsCurrent);
        // d ydot_j/dT
        for (size_t j = 0; j < m_nv; j++) {
            double ydotPerturbed = rhsPerturbed[j] / lhsPerturbed[j];
            double ydotCurrent = rhsCurrent[j] / lhsCurrent[j];
            trips.emplace_back(static_cast<int>(j + m_offset), m_offset,
                               (ydotPerturbed - ydotCurrent) / deltaTemp);
        }
        // d T_dot/dnj
        Eigen::VectorXd netProductionRates = Eigen::VectorXd::Zero(m_nsp);
        Eigen::VectorXd internal_energy = Eigen::VectorXd::Zero(m_nsp);
        Eigen::VectorXd specificHeat = Eigen::VectorXd::Zero(m_nsp);
        // getting species data
        m_thermo->getPartialMolarIntEnergies(asSpan(internal_energy));
        m_kin->getNetProductionRates(asSpan(netProductionRates));
        m_thermo->getPartialMolarCv_TV(asSpan(specificHeat));
        for (size_t i = 0; i < m_nsp; i++) {
            netProductionRates[i] *= m_vol;
        }
        // scale net production rates by  volume to get molar rate
        double qdot = internal_energy.dot(netProductionRates);
        double denom = 1 / (m_TotalCv * m_TotalCv);
        Eigen::VectorXd uk_dnkdnj_sums = dnk_dnj.transpose() * internal_energy;
        // add derivatives to jacobian
        for (size_t j = 0; j < m_nsp; j++) {
            trips.emplace_back(m_offset, static_cast<int>(j + m_offset + m_sidx),
                (specificHeat[j] * qdot - m_TotalCv * uk_dnkdnj_sums[j]) * denom);
        }
    }
}

void IdealGasMoleReactor::getJacobianScalingFactors(
    double& f_species, span<double> f_energy)
{
    f_species = 1.0 / m_vol;
    for (size_t k = 0; k < m_nsp; k++) {
        f_energy[k] = - m_uk[k] / m_TotalCv;
    }
}

}
