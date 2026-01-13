//! @file IdealGasConstPressureMoleReactor.cpp A constant pressure
//! zero-dimensional reactor with moles as the state

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/IdealGasConstPressureMoleReactor.h"
#include "cantera/zeroD/FlowDevice.h"
#include "cantera/zeroD/ReactorNet.h"
#include "cantera/zeroD/ReactorSurface.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/base/utilities.h"
#include <limits>

namespace Cantera
{

void IdealGasConstPressureMoleReactor::getState(double* y)
{
    // get mass for calculations
    m_mass = m_thermo->density() * m_vol;
    // set the first component to the temperature
    y[0] = m_thermo->temperature();
    // get moles of species in remaining state
    getMoles(y + m_sidx);
}

void IdealGasConstPressureMoleReactor::initialize(double t0)
{
    if (m_thermo->type() != "ideal-gas" && m_thermo->type() != "plasma") {
        throw CanteraError("IdealGasConstPressureMoleReactor::initialize",
                           "Incompatible phase type '{}' provided", m_thermo->type());
    }
    ConstPressureMoleReactor::initialize(t0);
    m_hk.resize(m_nsp, 0.0);
}

void IdealGasConstPressureMoleReactor::updateState(double* y)
{
    // the components of y are: [0] the temperature, [1...K+1) are the
    // moles of each species, and [K+1...] are the moles of surface
    // species on each wall.
    setMassFromMoles(y + m_sidx);
    m_thermo->setMolesNoTruncate(span<const double>(y + m_sidx, m_nsp));
    m_thermo->setState_TP(y[0], m_pressure);
    m_vol = m_mass / m_thermo->density();
    m_thermo->getPartialMolarEnthalpies(m_hk.data());
    m_TotalCp = m_mass * m_thermo->cp_mass();
    updateConnected(false);
}

void IdealGasConstPressureMoleReactor::eval(double time, double* LHS, double* RHS)
{
    double& mcpdTdt = RHS[0]; // m * c_p * dT/dt
    double* dndt = RHS + m_sidx; // kmol per s

    evalWalls(time);
    updateSurfaceProductionRates();
    auto imw = m_thermo->inverseMolecularWeights();

    if (m_chem) {
        m_kin->getNetProductionRates(&m_wdot[0]); // "omega dot"
    }

    // external heat transfer
    mcpdTdt += m_Qdot;

    if (m_energy) {
        mcpdTdt += m_thermo->intrinsicHeating() * m_vol;
    }

    for (size_t n = 0; n < m_nsp; n++) {
        // heat release from gas phase and surface reactions
        mcpdTdt -= m_wdot[n] * m_hk[n] * m_vol;
        mcpdTdt -= m_sdot[n] * m_hk[n];
        // production in gas phase and from surfaces
        dndt[n] = m_wdot[n] * m_vol + m_sdot[n];
    }

    // add terms for outlets
    for (auto outlet : m_outlet) {
        for (size_t n = 0; n < m_nsp; n++) {
            // flow of species out of system
            dndt[n] -= outlet->outletSpeciesMassFlowRate(n) * imw[n];
        }
    }

    // add terms for inlets
    for (auto inlet : m_inlet) {
        double mdot = inlet->massFlowRate();
        mcpdTdt += inlet->enthalpy_mass() * mdot;
        for (size_t n = 0; n < m_nsp; n++) {
            double mdot_spec = inlet->outletSpeciesMassFlowRate(n);
            // flow of species into system and dilution by other species
            dndt[n] += inlet->outletSpeciesMassFlowRate(n) * imw[n];
            mcpdTdt -= m_hk[n] * imw[n] * mdot_spec;
        }
    }

    if (m_energy) {
        LHS[0] = m_TotalCp;
    } else {
        RHS[0] = 0.0;
    }
}

void IdealGasConstPressureMoleReactor::getJacobianElements(
    vector<Eigen::Triplet<double>>& trips)
{
    // dnk_dnj represents d(dot(n_k)) / d (n_j) but is first assigned as
    // d (dot(omega)) / d c_j, it is later transformed appropriately.
    Eigen::SparseMatrix<double> dnk_dnj = m_kin->netProductionRates_ddCi();

    // add species to species derivatives  elements to the jacobian
    // calculate ROP derivatives, excluding the terms -n_i / (V * N) dc_i/dn_j
    // as it substantially reduces matrix sparsity
    size_t offset = static_cast<int>(m_offset + m_sidx);
    // double molarVol = m_thermo->molarVolume();
    for (int k = 0; k < dnk_dnj.outerSize(); k++) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(dnk_dnj, k); it; ++it) {
            trips.emplace_back(it.row() + offset, it.col() + offset, it.value());
        }
    }

    // Temperature Derivatives
    if (m_energy) {
        // getting perturbed state for finite difference
        double deltaTemp = m_thermo->temperature()
            * std::sqrt(std::numeric_limits<double>::epsilon());
        // get current state
        vector<double> yCurrent(m_nv);
        getState(yCurrent.data());
        // finite difference temperature derivatives
        vector<double> lhsPerturbed(m_nv, 1.0), lhsCurrent(m_nv, 1.0);
        vector<double> rhsPerturbed(m_nv, 0.0), rhsCurrent(m_nv, 0.0);
        vector<double> yPerturbed = yCurrent;
        // perturb temperature
        yPerturbed[0] += deltaTemp;
        // getting perturbed state
        updateState(yPerturbed.data());
        double time = (m_net != nullptr) ? m_net->time() : 0.0;
        eval(time, lhsPerturbed.data(), rhsPerturbed.data());
        // reset and get original state
        updateState(yCurrent.data());
        eval(time, lhsCurrent.data(), rhsCurrent.data());
        // d ydot_j/dT
        for (size_t j = 0; j < m_nv; j++) {
            double ydotPerturbed = rhsPerturbed[j] / lhsPerturbed[j];
            double ydotCurrent = rhsCurrent[j] / lhsCurrent[j];
            trips.emplace_back(static_cast<int>(j + m_offset), static_cast<int>(m_offset),
                                     (ydotPerturbed - ydotCurrent) / deltaTemp);
        }
        // d T_dot/dnj
        // allocate vectors for whole system
        Eigen::VectorXd netProductionRates = Eigen::VectorXd::Zero(m_nsp);
        Eigen::VectorXd enthalpy = Eigen::VectorXd::Zero(m_nsp);
        Eigen::VectorXd specificHeat = Eigen::VectorXd::Zero(m_nsp);
        // gas phase
        m_thermo->getPartialMolarCp(specificHeat.data());
        m_thermo->getPartialMolarEnthalpies(enthalpy.data());
        m_kin->getNetProductionRates(netProductionRates.data());
        // scale production rates by the volume for gas species
        for (size_t i = 0; i < m_nsp; i++) {
            netProductionRates[i] *= m_vol;
        }
        double qdot = enthalpy.dot(netProductionRates);
        double denom = 1 / (m_TotalCp * m_TotalCp);
        Eigen::VectorXd hk_dnkdnj_sums = dnk_dnj.transpose() * enthalpy;
        // Add derivatives to jac by spanning columns
        for (size_t j = 0; j < m_nsp; j++) {
            trips.emplace_back(m_offset, static_cast<int>(j + m_offset + m_sidx),
                (specificHeat[j] * qdot - m_TotalCp * hk_dnkdnj_sums[j]) * denom);
        }
    }
}

void IdealGasConstPressureMoleReactor::getJacobianScalingFactors(
    double& f_species, double* f_energy)
{
    f_species = 1.0 / m_vol;
    for (size_t k = 0; k < m_nsp; k++) {
        f_energy[k] = - m_hk[k] / m_TotalCp;
    }
}

size_t IdealGasConstPressureMoleReactor::componentIndex(const string& nm) const
{
    if (nm == "temperature") {
        return 0;
    }
    try {
        return m_thermo->speciesIndex(nm) + m_sidx;
    } catch (const CanteraError&) {
        throw CanteraError("IdealGasConstPressureReactor::componentIndex",
            "Component '{}' not found", nm);
    }
}

string IdealGasConstPressureMoleReactor::componentName(size_t k) {
    if (k == 0) {
        return "temperature";
    } else if (k >= m_sidx && k < neq()) {
        return m_thermo->speciesName(k - m_sidx);
    }
    throw IndexError("IdealGasConstPressureMoleReactor::componentName",
        "components", k, m_nv);
}

double IdealGasConstPressureMoleReactor::upperBound(size_t k) const {
    if (k == 0) {
        //@todo: Revise pending resolution of https://github.com/Cantera/enhancements/issues/229
        return 1.5 * m_thermo->maxTemp();
    } else {
        return BigNumber; // moles of a bulk or surface species
    }
}

double IdealGasConstPressureMoleReactor::lowerBound(size_t k) const {
    if (k == 0) {
        //@todo: Revise pending resolution of https://github.com/Cantera/enhancements/issues/229
        return 0.5 * m_thermo->minTemp();
    } else {
        return ConstPressureMoleReactor::lowerBound(k);
    }
}

}
