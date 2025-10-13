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
#include <limits>

namespace Cantera
{

void IdealGasMoleReactor::getState(double* y)
{
    if (m_thermo == 0) {
        throw CanteraError("IdealGasMoleReactor::getState",
                           "Error: reactor is empty.");
    }
    m_thermo->restoreState(m_state);

    // get mass for calculations
    m_mass = m_thermo->density() * m_vol;

    // set the first component to the temperature
    y[0] = m_thermo->temperature();

    // set the second component to the volume
    y[1] = m_vol;

    // get moles of species in remaining state
    getMoles(y + m_sidx);
    // set the remaining components to the surface species moles on
    // the walls
    getSurfaceInitialConditions(y + m_nsp + m_sidx);
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
        return speciesIndex(nm) + m_sidx;
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
    if (m_thermo->type() != "ideal-gas") {
        throw CanteraError("IdealGasMoleReactor::initialize",
                           "Incompatible phase type '{}' provided", m_thermo->type());
    }
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

vector<size_t> IdealGasMoleReactor::steadyConstraints() const
{
    if (nSurfs() != 0) {
        throw CanteraError("IdealGasMoleReactor::steadyConstraints",
            "Steady state solver cannot currently be used with IdealGasMoleReactor"
            " when reactor surfaces are present.\n"
            "See https://github.com/Cantera/enhancements/issues/234");
    }
    if (energyEnabled()) {
        return {1}; // volume
    } else {
        return {0, 1}; // temperature and volume
    }
}

void IdealGasMoleReactor::updateState(double* y)
{
    // the components of y are: [0] the temperature, [1] the volume, [2...K+1) are the
    // moles of each species, and [K+1...] are the moles of surface
    // species on each wall.
    setMassFromMoles(y + m_sidx);
    m_vol = y[1];
    // set state
    m_thermo->setMolesNoTruncate(y + m_sidx);
    m_thermo->setState_TD(y[0], m_mass / m_vol);
    updateConnected(true);
    updateSurfaceState(y + m_nsp + m_sidx);
}

void IdealGasMoleReactor::eval(double time, double* LHS, double* RHS)
{
    double& mcvdTdt = RHS[0]; // m * c_v * dT/dt
    double* dndt = RHS + m_sidx; // kmol per s

    evalWalls(time);

    m_thermo->restoreState(m_state);

    m_thermo->getPartialMolarIntEnergies(&m_uk[0]);
    const vector<double>& imw = m_thermo->inverseMolecularWeights();

    if (m_chem) {
        m_kin->getNetProductionRates(&m_wdot[0]); // "omega dot"
    }

    // evaluate surfaces
    evalSurfaces(LHS + m_nsp + m_sidx, RHS + m_nsp + m_sidx, m_sdot.data());

    // external heat transfer and compression work
    mcvdTdt += - m_pressure * m_vdot + m_Qdot;

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
        LHS[0] = m_mass * m_thermo->cv_mass();
    } else {
        RHS[0] = 0;
    }
}

Eigen::SparseMatrix<double> IdealGasMoleReactor::jacobian()
{
    if (m_nv == 0) {
        throw CanteraError("IdealGasMoleReactor::jacobian",
                           "Reactor must be initialized first.");
    }
    // clear former jacobian elements
    m_jac_trips.clear();
    // dnk_dnj represents d(dot(n_k)) / d (n_j) but is first assigned as
    // d (dot(omega)) / d c_j, it is later transformed appropriately.
    Eigen::SparseMatrix<double> dnk_dnj = m_kin->netProductionRates_ddCi();
    // species size that accounts for surface species
    size_t ssize = m_nv - m_sidx;
    // map derivatives from the surface chemistry jacobian
    // to the reactor jacobian
    if (!m_surfaces.empty()) {
        vector<Eigen::Triplet<double>> species_trips;
        for (int k = 0; k < dnk_dnj.outerSize(); k++) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(dnk_dnj, k); it; ++it) {
                species_trips.emplace_back(static_cast<int>(it.row()),
                                           static_cast<int>(it.col()), it.value());
            }
        }
        addSurfaceJacobian(species_trips);
        dnk_dnj.resize(ssize, ssize);
        dnk_dnj.setFromTriplets(species_trips.begin(), species_trips.end());
    }
    // add species to species derivatives  elements to the jacobian
    // calculate ROP derivatives, excluding the terms -n_i / (V * N) dc_i/dn_j
    // as it substantially reduces matrix sparsity
    for (int k = 0; k < dnk_dnj.outerSize(); k++) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(dnk_dnj, k); it; ++it) {
            m_jac_trips.emplace_back(static_cast<int>(it.row() + m_sidx),
                static_cast<int>(it.col() + m_sidx), it.value());
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
        getState(yCurrent.data());
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
            m_jac_trips.emplace_back(static_cast<int>(j), 0,
                                     (ydotPerturbed - ydotCurrent) / deltaTemp);
        }
        // d T_dot/dnj
        Eigen::VectorXd netProductionRates = Eigen::VectorXd::Zero(ssize);
        Eigen::VectorXd internal_energy = Eigen::VectorXd::Zero(ssize);
        Eigen::VectorXd specificHeat = Eigen::VectorXd::Zero(ssize);
        // getting species data
        m_thermo->getPartialMolarIntEnergies(internal_energy.data());
        m_kin->getNetProductionRates(netProductionRates.data());
        m_thermo->getPartialMolarCp(specificHeat.data());
        // convert Cp to Cv for ideal gas as Cp - Cv = R
        for (size_t i = 0; i < m_nsp; i++) {
            specificHeat[i] -= GasConstant;
            netProductionRates[i] *= m_vol;
        }
        // scale net production rates by  volume to get molar rate
        double qdot = internal_energy.dot(netProductionRates);
        // find the sum of n_i and cp_i
        double NCv = 0.0;
        double* moles = yCurrent.data() + m_sidx;
        for (size_t i = 0; i < ssize; i++) {
            NCv += moles[i] * specificHeat[i];
        }
        // make denominator beforehand
        double denom = 1 / (NCv * NCv);
        Eigen::VectorXd uk_dnkdnj_sums = dnk_dnj.transpose() * internal_energy;
        // add derivatives to jacobian
        for (size_t j = 0; j < ssize; j++) {
            m_jac_trips.emplace_back(0, static_cast<int>(j + m_sidx),
                (specificHeat[j] * qdot - NCv * uk_dnkdnj_sums[j]) * denom);
        }
    }
    // convert triplets to sparse matrix
    Eigen::SparseMatrix<double> jac(m_nv, m_nv);
    jac.setFromTriplets(m_jac_trips.begin(), m_jac_trips.end());
    return jac;
}

}
