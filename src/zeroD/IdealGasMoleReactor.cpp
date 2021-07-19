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

using namespace std;

namespace Cantera
{

void IdealGasMoleReactor::setThermoMgr(ThermoPhase& thermo)
{
    if (thermo.type() != "IdealGas") {
        throw CanteraError("IdealGasMoleReactor::setThermoMgr",
                           "Incompatible phase type provided");
    }
    MoleReactor::setThermoMgr(thermo);
}

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

    // use inverse molecular weights
    const double* Y = m_thermo->massFractions();
    const vector_fp& imw = m_thermo->inverseMolecularWeights();
    double* ys = y + m_sidx;
    for (size_t i = 0; i < m_nsp; i++) {
        ys[i] = m_mass * imw[i] * Y[i];
    }
    // set the remaining components to the surface species moles on
    // the walls
    getSurfaceInitialConditions(y + m_nsp + m_sidx);
}

size_t IdealGasMoleReactor::componentIndex(const string& nm) const
{
    size_t k = speciesIndex(nm);
    if (k != npos) {
        return k + m_sidx;
    } else if (nm == "temperature") {
        return 0;
    } else if (nm == "volume") {
        return 1;
    } else {
        return npos;
    }
}

void IdealGasMoleReactor::initialize(double t0)
{
    MoleReactor::initialize(t0);
    m_uk.resize(m_nsp, 0.0);
}

void IdealGasMoleReactor::updateState(double* y)
{
    // the components of y are: [0] the temperature, [1] the volume, [2...K+1) are the
    // moles of each species, and [K+1...] are the moles of surface
    // species on each wall.
    const vector_fp& mw = m_thermo->molecularWeights();
    // calculate mass from moles
    m_mass = 0;
    for (size_t i = 0; i < m_nv - m_sidx; i++) {
        m_mass += y[i + m_sidx] * mw[i];
    }
    m_vol = y[1];
    // set state
    m_thermo->setMolesNoTruncate(y + m_sidx);
    m_thermo->setState_TR(y[0], m_mass / m_vol);
    updateConnected(true);
    updateSurfaceState(y + m_nsp + m_sidx);
}

void IdealGasMoleReactor::eval(double time, double* LHS, double* RHS)
{
    double& mcvdTdt = RHS[0]; // m * c_v * dT/dt
    double* dydt = RHS + m_sidx; // kmol per s

    evalWalls(time);

    m_thermo->restoreState(m_state);

    m_thermo->getPartialMolarIntEnergies(&m_uk[0]);
    const vector_fp& imw = m_thermo->inverseMolecularWeights();

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
        dydt[n] = (m_wdot[n] * m_vol + m_sdot[n]);
    }

    // add terms for outlets
    for (auto outlet : m_outlet) {
        for (size_t n = 0; n < m_nsp; n++) {
            // flow of species into system and dilution by other species
            dydt[n] -= outlet->outletSpeciesMassFlowRate(n) * imw[n];
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
            dydt[n] += inlet->outletSpeciesMassFlowRate(n) * imw[n];
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

Eigen::SparseMatrix<double> IdealGasMoleReactor::jacobian(double t, double* y)
{
    // clear former jacobian elements
    m_jac_trips.clear();
    // Determine Species Derivatives
    // get ROP derivatives
    double scalingFactor = m_vol/accumulate(y + m_sidx, y + m_nv, 0.0);
    Eigen::SparseMatrix<double> speciesDervs = scalingFactor * m_kin->netProductionRates_ddX();
    // add to preconditioner
    for (int k=0; k<speciesDervs.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(speciesDervs, k); it; ++it) {
            m_jac_trips.emplace_back(it.row() + m_sidx, it.col() + m_sidx, it.value());
        }
    }
    // Temperature Derivatives
    if (m_energy) {
        // getting perturbed state for finite difference
        double deltaTemp = y[0] * std::sqrt(std::numeric_limits<double>::epsilon());
        // finite difference temperature derivatives
        vector_fp yNext(m_nv);
        vector_fp ydotNext(m_nv);
        vector_fp yCurrent(m_nv);
        vector_fp ydotCurrent(m_nv);
        // copy y to current and next
        copy(y, y + m_nv, yCurrent.begin());
        copy(y, y + m_nv, yNext.begin());
        // perturb temperature
        yNext[0] += deltaTemp;
        // getting perturbed state
        updateState(yNext.data());
        eval(t, yNext.data(), ydotNext.data());
        // reset and get original state
        updateState(yCurrent.data());
        eval(t, yCurrent.data(), ydotCurrent.data());
        // d T_dot/dT
        m_jac_trips.emplace_back(0, 0, (ydotNext[0] - ydotCurrent[0]) / deltaTemp);
        // d omega_dot_j/dT
        for (size_t j = m_sidx; j < m_nv; j++) {
            m_jac_trips.emplace_back(j, 0, (ydotNext[j] - ydotCurrent[j]) / deltaTemp);
        }
        // find derivatives d T_dot/dNj
        vector_fp specificHeat(m_nsp);
        vector_fp netProductionRates(m_nsp);
        vector_fp internal_energy(m_nsp);
        // getting species data
        m_thermo->getPartialMolarCp(specificHeat.data());
        m_thermo->getPartialMolarIntEnergies(internal_energy.data());
        m_kin->getNetProductionRates(netProductionRates.data());
        // scale net production rates by  volume to get molar rate
        scale(netProductionRates.begin(), netProductionRates.end(), netProductionRates.begin(), m_vol);
        // convert Cp to Cv for ideal gas as Cp - Cv = R
        for (size_t i = 0; i < specificHeat.size(); i++) {
            specificHeat[i] -= GasConstant;
        }
        // finding a sum inside the derivative
        double uknkSum = 0;
        double NtotalCv = accumulate(y + m_sidx, y + m_nv, 0.0) * m_thermo->cv_mole();
        for (size_t i = 0; i < m_nsp; i++) {
            uknkSum += internal_energy[i] * netProductionRates[i];
        }
        // finding derivatives
        // spans columns
        for (size_t j = 0; j < m_nsp; j++) {
            double ukdnkdnjSum = 0;
            // spans rows
            for (size_t k = 0; k < m_nsp; k++) {
                ukdnkdnjSum += internal_energy[k] * speciesDervs.coeff(k, j);
            }
            // set appropriate column of preconditioner
            m_jac_trips.emplace_back(0, j + m_sidx, (-ukdnkdnjSum * NtotalCv +
            specificHeat[j] *  uknkSum) / (NtotalCv * NtotalCv));
        }
    }
    Eigen::SparseMatrix<double> jac (m_nv, m_nv);
    jac.setFromTriplets(m_jac_trips.begin(), m_jac_trips.end());
    return jac;
}

}
