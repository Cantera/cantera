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

void IdealGasConstPressureMoleReactor::setThermo(ThermoPhase& thermo)
{
    if (thermo.type() != "ideal-gas") {
        throw CanteraError("IdealGasConstPressureMoleReactor::setThermo",
                           "Incompatible phase type provided");
    }
    ConstPressureMoleReactor::setThermo(thermo);
}

void IdealGasConstPressureMoleReactor::getState(double* y)
{
    if (m_thermo == 0) {
        throw CanteraError("IdealGasConstPressureMoleReactor::getState",
                           "Error: reactor is empty.");
    }
    m_thermo->restoreState(m_state);
    // get mass for calculations
    m_mass = m_thermo->density() * m_vol;
    // set the first component to the temperature
    y[0] = m_thermo->temperature();
    // get moles of species in remaining state
    getMoles(y + m_sidx);
    // set the remaining components to the surface species moles on the walls
    getSurfaceInitialConditions(y + m_nsp + m_sidx);
}

void IdealGasConstPressureMoleReactor::initialize(double t0)
{
    ConstPressureMoleReactor::initialize(t0);
    m_hk.resize(m_nsp, 0.0);
}

void IdealGasConstPressureMoleReactor::updateState(double* y)
{
    // the components of y are: [0] the temperature, [1...K+1) are the
    // moles of each species, and [K+1...] are the moles of surface
    // species on each wall.
    setMassFromMoles(y + m_sidx);
    m_thermo->setMolesNoTruncate(y + m_sidx);
    m_thermo->setState_TP(y[0], m_pressure);
    m_vol = m_mass / m_thermo->density();
    updateConnected(false);
    updateSurfaceState(y + m_nsp + m_sidx);
}

void IdealGasConstPressureMoleReactor::eval(double time, double* LHS, double* RHS)
{
    double& mcpdTdt = RHS[0]; // m * c_p * dT/dt
    double* dndt = RHS + m_sidx; // kmol per s

    evalWalls(time);

    m_thermo->restoreState(m_state);

    m_thermo->getPartialMolarEnthalpies(&m_hk[0]);
    const vector<double>& imw = m_thermo->inverseMolecularWeights();

    if (m_chem) {
        m_kin->getNetProductionRates(&m_wdot[0]); // "omega dot"
    }

    // evaluate reactor surfaces
    evalSurfaces(LHS + m_nsp + m_sidx, RHS + m_nsp + m_sidx, m_sdot.data());

    // external heat transfer
    mcpdTdt += m_Qdot;

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
            // flow of species into system and dilution by other species
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
        LHS[0] = m_mass * m_thermo->cp_mass();
    } else {
        RHS[0] = 0.0;
    }
}

Eigen::SparseMatrix<double> IdealGasConstPressureMoleReactor::jacobian()
{
    if (m_nv == 0) {
        throw CanteraError("IdealGasConstPressureMoleReactor::jacobian",
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
        vector<Eigen::Triplet<double>> species_trips(dnk_dnj.nonZeros());
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
    // get net production rates
    Eigen::VectorXd netProductionRates = Eigen::VectorXd::Zero(ssize);
    // gas phase net production rates
    m_kin->getNetProductionRates(netProductionRates.data());
    // surface phase net production rates mapped to reactor gas phase
    for (auto &S: m_surfaces) {
        auto curr_kin = S->kinetics();
        vector<double> prod_rates(curr_kin->nTotalSpecies());
        curr_kin->getNetProductionRates(prod_rates.data());
        for (size_t i = 0; i < curr_kin->nTotalSpecies(); i++) {
            size_t row = speciesIndex(curr_kin->kineticsSpeciesName(i));
            if (row != npos) {
                netProductionRates[row] += prod_rates[i];
            }
        }
    }
    double molarVol = m_thermo->molarVolume();
    // add species to species derivatives  elements to the jacobian
    // calculate ROP derivatives, excluding the terms -n_i / (V * N) dc_i/dn_j
    // as it substantially reduces matrix sparsity
    for (int k = 0; k < dnk_dnj.outerSize(); k++) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(dnk_dnj, k); it; ++it) {
            // gas phase species need the addition of  V / N * omega_dot
            if (static_cast<size_t>(it.row()) < m_nsp) {
                it.valueRef() = it.value() + netProductionRates[it.row()] * molarVol;
            }
            m_jac_trips.emplace_back(static_cast<int>(it.row() + m_sidx),
                static_cast<int>(it.col() + m_sidx), it.value());
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
            m_jac_trips.emplace_back(static_cast<int>(j), 0,
                                     (ydotPerturbed - ydotCurrent) / deltaTemp);
        }
        // d T_dot/dnj
        // allocate vectors for whole system
        Eigen::VectorXd enthalpy = Eigen::VectorXd::Zero(ssize);
        Eigen::VectorXd specificHeat = Eigen::VectorXd::Zero(ssize);
        // gas phase
        m_thermo->getPartialMolarCp(specificHeat.data());
        m_thermo->getPartialMolarEnthalpies(enthalpy.data());
        // scale production rates by the volume for gas species
        for (size_t i = 0; i < m_nsp; i++) {
            netProductionRates[i] *= m_vol;
        }
        double qdot = enthalpy.dot(netProductionRates);
        // find denominator ahead of time
        double NCp = 0.0;
        double* moles = yCurrent.data() + m_sidx;
        for (size_t i = 0; i < ssize; i++) {
            NCp += moles[i] * specificHeat[i];
        }
        double denom = 1 / (NCp * NCp);
        Eigen::VectorXd hk_dnkdnj_sums = dnk_dnj.transpose() * enthalpy;
        // Add derivatives to jac by spanning columns
        for (size_t j = 0; j < ssize; j++) {
            m_jac_trips.emplace_back(0, static_cast<int>(j + m_sidx),
                (specificHeat[j] * qdot - NCp * hk_dnkdnj_sums[j]) * denom);
        }
    }
    // convert triplets to sparse matrix
    Eigen::SparseMatrix<double> jac(m_nv, m_nv);
    jac.setFromTriplets(m_jac_trips.begin(), m_jac_trips.end());
    return jac;
}

size_t IdealGasConstPressureMoleReactor::componentIndex(const string& nm) const
{
    size_t k = speciesIndex(nm);
    if (k != npos) {
        return k + m_sidx;
    } else if (nm == "temperature") {
        return 0;
    } else {
        return npos;
    }
}

string IdealGasConstPressureMoleReactor::componentName(size_t k) {
    if (k == 0) {
        return "temperature";
    } else if (k >= m_sidx && k < neq()) {
        k -= m_sidx;
        if (k < m_thermo->nSpecies()) {
            return m_thermo->speciesName(k);
        } else {
            k -= m_thermo->nSpecies();
        }
        for (auto& S : m_surfaces) {
            ThermoPhase* th = S->thermo();
            if (k < th->nSpecies()) {
                return th->speciesName(k);
            } else {
                k -= th->nSpecies();
            }
        }
    }
    throw CanteraError("IdealGasConstPressureMoleReactor::componentName",
                       "Index is out of bounds.");
}

}
