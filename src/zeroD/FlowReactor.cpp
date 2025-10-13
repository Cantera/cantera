//! @file FlowReactor.cpp A steady-state, ideal-gas, adiabatic,
//!                       constant-area (cylindrical), frictionless plug flow reactor

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/FlowReactor.h"
#include "cantera/base/global.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/zeroD/ReactorSurface.h"
#include "cantera/zeroD/ReactorNet.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/numerics/DenseMatrix.h"
#include "cantera/kinetics/InterfaceKinetics.h"

namespace Cantera
{

void FlowReactor::getStateDae(double* y, double* ydot)
{
    if (m_thermo == nullptr) {
        throw CanteraError("FlowReactor::getStateDae", "Error: reactor is empty.");
    }
    m_thermo->restoreState(m_state);
    m_thermo->getMassFractions(y+m_offset_Y);
    const vector<double>& mw = m_thermo->molecularWeights();

    // set the first component to the initial density
    y[0] = m_rho;

    if (m_u < 0) {
        throw CanteraError("FlowReactor::getStateDae",
            "Set mass flow rate before initializing reactor");
    }

    // set the second component to the initial speed
    y[1] = m_u;

    // set the third component to the initial pressure
    y[2] = m_P;

    // set the fourth component to the initial temperature
    y[3] = m_T;

    if (m_chem) {
        m_kin->getNetProductionRates(m_wdot.data()); // "omega dot"
    }

    // need to advance the reactor surface to steady state to get the initial
    // coverages
    for (auto m_surf : m_surfaces) {
        m_surf->restoreState();
        auto kin = static_cast<InterfaceKinetics*>(m_surf->kinetics());
        kin->advanceCoverages(100.0, m_ss_rtol, m_ss_atol, 0, m_max_ss_steps,
                              m_max_ss_error_fails);
        auto& surf = dynamic_cast<SurfPhase&>(kin->thermo(0));
        vector<double> cov(surf.nSpecies());
        surf.getCoverages(cov.data());
        m_surf->setCoverages(cov.data());
    }

    // set the initial coverages
    getSurfaceInitialConditions(y + m_offset_Y + m_nsp);

    // reset ydot vector
    std::fill(ydot, ydot + m_nv, 0.0);
    // calculate initial d(coverage)/dt values
    evalSurfaces(ydot + m_nsp + 4, m_sdot.data());

    // next, we must solve for the initial derivative values
    // a . ydot = b
    // ydot -> [rho', u', P', T', Yk']
    // a -> coefficients of [rho', u', P', T', Yk'] in the governing eq's, given
    //      the initial state
    // b -> rhs constant of each conservation equation
    //
    // note that the species coverages are included in the algebraic constraints,
    // hence are not included here
    DenseMatrix a(m_offset_Y + m_nsp, m_offset_Y + m_nsp);

    // first row is the ideal gas equation, differentiated with respect to z
    a(0, 0) = - GasConstant * m_T / m_thermo->meanMolecularWeight();
    a(0, 2) = 1;
    a(0, 3) = - m_rho * GasConstant / m_thermo->meanMolecularWeight();
    for (size_t i = 0; i < m_nsp; ++i) {
        a(0, m_offset_Y + i) = - m_rho * m_T / mw[i];
    }

    // initialize the second row from mass conservation equation (Kee 16.48)
    a(1, 0) = m_u; // rho' * u
    a(1, 1) = m_rho; // u' * rho

    // initialize the third row from momentum conservation equation (Kee 16.62),
    // rho * u * u' + P' = -u * (P / A_c) * sum(sdot_k * W_k)
    a(2, 1) = m_rho * m_u; // rho * u * u'
    a(2, 2) = 1; // 1 * P'

    // initialize the fourth row from conservation of energy (Kee 16.58), adiabatic
    double cp_mass = m_thermo->cp_mass();
    a(3, 3) = m_rho * m_u * cp_mass; // rho * u * cp * T'

    // initialize the next rows from the mass-fraction equations (Kee 16.51)
    for (size_t i = 0; i < m_nsp; ++i) {
        a(m_offset_Y + i, m_offset_Y + i) = m_rho * m_u; // rho * u * Yk'
    }

    // now set the RHS vector

    // get (perim / Ac) * sum(sk' * wk), used multiple places
    double h_sk_wk = 0;
    double hydraulic = surfaceAreaToVolumeRatio();
    for (size_t i = 0; i < m_nsp; ++i) {
        h_sk_wk += hydraulic * m_sdot[i] * mw[i];
    }

    // RHS of ideal gas eq. is zero
    ydot[0] = 0;

    // mass continuity, Kee 16.48
    // (perim / Ac) * sum(sk' * wk)
    ydot[1] = h_sk_wk;

    // momentum conservation, Kee 16.62, no friction
    // - u * (perim / Ac) * sum(sk' * wk)
    ydot[2] = -m_u * h_sk_wk;

    if (m_energy) {
        // energy conservation, Kee 16.58, adiabatic
        // -sum(wk' * Wk * hk) - (perim / Ac) * sum(sk' * Wk * hk)
        // Note: the partial molar enthalpies are in molar units, while Kee uses mass
        //       units, where:
        //              h_mass = h_mole / Wk
        //       hence:
        //              h_mass * Wk = h_mol
        m_thermo->getPartialMolarEnthalpies(m_hk.data());
        for (size_t i = 0; i < m_nsp; ++i) {
            ydot[3] -= m_hk[i] * (m_wdot[i] + hydraulic * m_sdot[i]);
        }
    } else {
        ydot[3] = 0;
    }

    // mass-fraction equations Kee 16.51
    // - Yk * (perim / Ac) * sum(sk' * Wk) + wk' * Wk + (perim / Ac) * sk' * Wk
    for (size_t i = 0; i < m_nsp; ++i) {
        ydot[m_offset_Y + i] = -y[m_offset_Y + i] * h_sk_wk +
            mw[i] * (m_wdot[i] + hydraulic * m_sdot[i]);
    }

    // and solve
    solve(a, ydot, 1, 0);
}

void FlowReactor::initialize(double t0)
{
    Reactor::initialize(t0);
    m_thermo->restoreState(m_state);
    // initialize state
    m_T = m_thermo->temperature();
    m_rho = m_thermo->density();
    m_P = m_thermo->pressure();
    m_T = m_thermo->temperature();
    // resize temporary arrays
    m_wdot.resize(m_nsp);
    m_hk.resize(m_nsp);
    // set number of variables to the number of non-species equations
    // i.e., density, velocity, pressure and temperature
    // plus the number of species in the gas phase
    m_nv = m_offset_Y + m_nsp;
    if (m_surfaces.size()) {
        size_t n_surf_species = 0;
        size_t n_total_species = 0;
        for (auto const &m_surf : m_surfaces) {
            // finally, add the number of species the surface for the site-fraction
            // equations
            n_surf_species += m_surf->thermo()->nSpecies();
            n_total_species += m_surf->kinetics()->nTotalSpecies();
        }
        m_nv += n_surf_species;
        m_sdot_temp.resize(n_total_species);
    }
}

void FlowReactor::syncState()
{
    Reactor::syncState();
    m_rho = m_thermo->density();
    m_P = m_thermo->pressure();
    m_T = m_thermo->temperature();
}

void FlowReactor::updateState(double* y)
{
    // Set the mass fractions and density of the mixture.
    m_rho = y[0];
    m_u = y[1];
    m_P = y[2];
    m_T = y[3];
    double* mss = y + m_offset_Y;
    m_thermo->setMassFractions_NoNorm(mss);
    m_thermo->setState_TP(m_T, m_P);

    // update surface
    updateSurfaceState(y + m_nsp + m_offset_Y);

    m_thermo->saveState(m_state);
}

void FlowReactor::setMassFlowRate(double mdot)
{
    m_rho = m_thermo->density();
    m_u = mdot/(m_rho * m_area);
}

void FlowReactor::setArea(double area) {
    double mdot = m_rho * m_u * m_area;
    m_area = area;
    setMassFlowRate(mdot);
}

double FlowReactor::surfaceAreaToVolumeRatio() const {
    if (m_sa_to_vol > 0) {
        return m_sa_to_vol;
    }

    // assuming a cylinder, volume = Pi * r^2 * L, and perimeter = 2 * Pi * r * L
    // where L is the length, and r is the radius of the reactor
    // hence, perimeter / area = 2 * Pi * r * L / (Pi * L * r^2) = 2 / r
    return 2.0 / sqrt(m_area / Pi);
}

void FlowReactor::updateSurfaceState(double* y)
{
    size_t loc = 0;
    for (auto& S : m_surfaces) {
        // set the coverages without normalization to avoid over-constraining
        // the system.
        // note: the ReactorSurface class doesn't normalize when calling setCoverages
        S->setCoverages(y+loc);
        S->restoreState();
        loc += S->thermo()->nSpecies();
    }
}

void FlowReactor::evalDae(double time, double* y, double* ydot, double* residual)
{
    m_thermo->restoreState(m_state);

    evalSurfaces(ydot + m_nsp + 4, m_sdot.data());
    const vector<double>& mw = m_thermo->molecularWeights();
    double sk_wk = 0;
    for (size_t i = 0; i < m_nsp; ++i) {
        sk_wk = m_sdot[i] * mw[i];
    }
    m_thermo->getPartialMolarEnthalpies(m_hk.data());
    // get net production
    if (m_chem) {
        m_kin->getNetProductionRates(m_wdot.data());
    }

    // set dphi/dz variables
    double drhodz = ydot[0];
    double dudz = ydot[1];
    double dPdz = ydot[2];
    double dTdz = ydot[3];

    // use equation of state for density residual
    residual[0] = m_rho - m_thermo->density();

    //! use mass continuity for velocity residual
    //! Kee.'s Chemically Reacting Flow, Eq. 16.48
    double hydraulic = surfaceAreaToVolumeRatio();
    residual[1] = m_u * drhodz + m_rho * dudz - sk_wk * hydraulic;

    //! Use conservation of momentum for pressure residual
    //! Kee.'s Chemically Reacting Flow, Eq. 16.62
    //! [NOTE: friction is neglected]
    residual[2] = m_rho * m_u * dudz + m_u * hydraulic * sk_wk + dPdz;

    //! use conservation of energy for temperature residual
    //! Kee.'s Chemically Reacting Flow, Eq. 16.58
    //! [NOTE: adiabatic]
    //  Also, as above:
    //              h_mass = h_mole / Wk
    //       hence:
    //              h_mass * Wk = h_mol
    if (m_energy) {
        double cp_mass = m_thermo->cp_mass();
        residual[3] = m_rho * m_u * cp_mass * dTdz;
        for (size_t i = 0; i < m_nsp; ++i) {
            residual[3] += m_hk[i] * (m_wdot[i] + hydraulic * m_sdot[i]);
        }
    } else {
        residual[3] = dTdz;
    }

    //! species conservation equations
    //! Kee.'s Chemically Reacting Flow, Eq. 16.51
    for (size_t i = 0; i < m_nsp; ++i) {
        residual[i + m_offset_Y] = m_rho * m_u * ydot[i + m_offset_Y] +
            y[i + m_offset_Y] * hydraulic * sk_wk -
            mw[i] * (m_wdot[i] + hydraulic * m_sdot[i]);
    }

    // surface algebraic constraints
    if (m_surfaces.size()) {
        size_t loc = m_offset_Y + m_nsp; // offset into residual vector
        for (auto &m_surf : m_surfaces) {
            Kinetics* kin = m_surf->kinetics();
            size_t nk = m_surf->thermo()->nSpecies();
            kin->getNetProductionRates(m_sdot_temp.data());
            double sum = y[loc];
            for (size_t i = 1; i < nk; ++i) {
                //! net surface production rate residuals
                //! Kee.'s Chemically Reacting Flow, Eq. 16.63
                residual[loc + i] = m_sdot_temp[i];
                //! next, evaluate the algebraic constraint to explicitly conserve
                //! surface site fractions.
                //! Kee.'s Chemically Reacting Flow, Eq. 16.64
                sum += y[loc + i];
            }
            // set residual of the first species sum - 1 to ensure
            // surface site conservation.
            // Note: the first species is selected to be consistent with
            // Reactor::evalSurfaces
            residual[loc] = sum - 1.0;
            loc += nk;
        }
    }
}

void FlowReactor::getConstraints(double* constraints) {
    // mark all variables differential equations unless otherwise specified
    std::fill(constraints, constraints + m_nv, 1.0);
    // the species coverages are algebraic constraints
    std::fill(constraints + m_offset_Y + m_nsp, constraints + m_nv, 0.0);
}

size_t FlowReactor::componentIndex(const string& nm) const
{
    if (nm == "density") {
        return 0;
    }
    if (nm == "speed") {
        return 1;
    }
    if (nm == "pressure") {
        return 2;
    }
    if (nm == "temperature") {
        return 3;
    }
    try {
        return speciesIndex(nm) + m_offset_Y;
    } catch (const CanteraError&) {
        throw CanteraError("FlowReactor::componentIndex",
            "Unknown component '{}'", nm);
    }
}

string FlowReactor::componentName(size_t k)
{
    if (k == 0) {
        return "density";
    } else if (k == 1) {
        return "speed";
    } else if (k == 2) {
        return "pressure";
    } else if (k == 3) {
        return "temperature";
    } else if (k >= 4 && k < neq()) {
        k -= 4;
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
    throw CanteraError("FlowReactor::componentName", "Index {} is out of bounds.", k);
}

}
