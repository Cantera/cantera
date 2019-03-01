//! @file FlowReactor.cpp A steady-state, ideal-gas, adiabatic,
//!                       constant-area (cylindrical), frictionless plug flow reactor

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/FlowReactor.h"
#include "cantera/base/global.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/zeroD/ReactorSurface.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/numerics/DenseMatrix.h"
#include "cantera/kinetics/InterfaceKinetics.h"

using namespace std;

namespace Cantera
{

FlowReactor::FlowReactor(double area, double sa_to_vol,
                         double ss_atol, double ss_rtol,
                         int max_ss_steps,
                         int max_ss_error_fails) :
    m_u(0.0),
    m_T(0.0),
    m_area(area),
    m_sa_to_vol(sa_to_vol),
    m_sdot_temp(0),
    m_hk(0),
    m_ss_rtol(ss_rtol),
    m_ss_atol(ss_atol),
    m_max_ss_steps(max_ss_steps),
    m_max_ss_error_fails(max_ss_error_fails)
{
}

void FlowReactor::getState(double* y, double* ydot)
{
    if (m_thermo == 0) {
        throw CanteraError("FlowReactor::getState",
                           "Error: reactor is empty.");
    }
    m_thermo->restoreState(m_state);
    m_thermo->getMassFractions(y+m_non_spec_eq);
    const vector_fp& mw = m_thermo->molecularWeights();

    // set the first component to the initial density
    y[0] = m_rho;

    // set the second component to the initial speed
    y[1] = m_u;

    // set the fourth component to the initial pressure
    y[2] = m_P;

    // set the third component to the initial temperature
    y[3] = m_T;

    if (m_chem) {
        m_kin->getNetProductionRates(&m_wdot[0]); // "omega dot"
    }

    // need to advance the reactor surface to steady state to get the initial
    // coverages
    for (auto &m_surf : m_surfaces)
    {
        auto *kin = static_cast<InterfaceKinetics*>(m_surf->kinetics());
        kin->advanceCoverages(100.0, m_ss_rtol, m_ss_atol, 0, m_max_ss_steps,
                              m_max_ss_error_fails);
        auto *surf = static_cast<SurfPhase*>(&kin->thermo(kin->reactionPhaseIndex()));
        vector_fp cov(surf->nSpecies());
        surf->getCoverages(&cov[0]);
        m_surf->setCoverages(&cov[0]);
        m_surf->syncCoverages();
    }

    // set the initial coverages
    getSurfaceInitialConditions(y + m_non_spec_eq + m_nsp);

    // reset ydot vector
    std::fill(ydot, ydot + m_nv, 0.0);
    // calculate initial d(coverage)/dt values
    evalSurfaces(0, ydot + m_non_spec_eq + m_nsp);

    // next, we must solve for the initial derivative values
    // a . ydot = b
    // ydot -> [rho', u', P', T', Yk']
    // a -> coefficients of [rho', u', P', T', Yk'] in the governing eq's, given
    //      the initial state
    // b -> rhs constant of each conservation equation
    //
    // note that the species coverages are included in the algebraic constraints,
    // hence are not included ehre
    DenseMatrix a;
    a.resize(m_non_spec_eq + m_nsp, m_non_spec_eq + m_nsp, 0.0);

    // first row is the ideal gas equation
    a(0, 0) = - GasConstant * m_T / m_thermo->meanMolecularWeight();
    a(0, 2) = 1;
    a(0, 3) = - m_rho * GasConstant / m_thermo->meanMolecularWeight();
    for (size_t i = 0; i < m_nsp; ++i)
    {
        a(0, m_non_spec_eq + i) = - m_rho * m_T / mw[i];
    }


    // initialize the second row from mass conservation equation
    // Kee 16.48
    a(1, 0) = m_u; // rho' * u
    a(1, 1) = m_rho; // u' * rho

    // initialize the third row from momentum conservation equation
    // Kee 16.62
    // -> rho * u * u' + P' = -u * (P / A_c) * sum(sdot_k * W_k)
    a(2, 1) = m_rho * m_u; // rho * u * u'
    a(2, 2) = 1; // 1 * P'

    // initialize the fourth row from conservation of energy
    // Kee 16.58, adiabatic
    const doublereal cp_mass = m_thermo->cp_mass();
    a(3, 3) = m_rho * m_u * cp_mass; // rho * u * cp * T'

    // initialize the next rows from the mass-fraction equations
    // Kee 16.51
    for (size_t i = 0; i < m_nsp; ++i)
    {
        a(m_non_spec_eq + i, m_non_spec_eq + i) = m_rho * m_u; // rho * u * Yk'
    }

    // now set the RHS vector

    // get (perim / Ac) * sum(sk' * wk), used multiple places
    double h_sk_wk = 0;
    const doublereal hydraulic = surfaceAreaToVolumeRatio();
    for (size_t i = 0; i < m_nsp; ++i)
    {
        h_sk_wk += hydraulic * m_sdot[i] * mw[i];
    }

    // RHS of ideal gas eq. is zero
    ydot[0] = 0;

    // mass continutity, Kee 16.48
    // (perim / Ac) * sum(sk' * wk)
    ydot[1] = h_sk_wk;

    // momentum conservation, Kee 16.62, no friction
    // - u * (perim / Ac) * sum(sk' * wk)
    ydot[2] = -m_u * h_sk_wk;

    // energy conservation, Kee 16.58, adiabatic
    // -sum(wk' * Wk * hk) - (perim / Ac) * sum(sk' * Wk * hk)
    // Note: the partial molar enthalpies are in molar units, while Kee uses mass
    //       units, where:
    //              h_mass = h_mole / Wk
    //       hence:
    //              h_mass * Wk = h_mol
    m_thermo->getPartialMolarEnthalpies(&m_hk[0]);
    for (size_t i = 0; i < m_nsp; ++i)
    {
        ydot[3] -= m_hk[i] * (m_wdot[i] + hydraulic * m_sdot[i]);
    }

    // mass-fraction equations Kee 16.51
    // - Yk * (perim / Ac) * sum(sk' * Wk) + wk' * Wk + (perim / Ac) * sk' * Wk
    for (size_t i = 0; i < m_nsp; ++i)
    {
        ydot[m_non_spec_eq + i] = -y[m_non_spec_eq + i] * h_sk_wk +
            mw[i] * (m_wdot[i] + hydraulic * m_sdot[i]);
    }

    // and solve
    solve(a, ydot, 1, 0);
}


void FlowReactor::initialize(doublereal t0)
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
    m_nv = m_non_spec_eq + m_nsp;
    if (m_surfaces.size())
    {
        size_t n_surf_species = 0;
        size_t n_total_species = 0;
        for (auto const &m_surf : m_surfaces)
        {
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
    ReactorBase::syncState();
    m_rho = m_thermo->density();
    m_P = m_thermo->pressure();
    m_T = m_thermo->temperature();
}

void FlowReactor::updateState(doublereal* y)
{
    // Set the mass fractions and density of the mixture.
    m_rho = y[0];
    m_u = y[1];
    m_P = y[2];
    m_T = y[3];
    doublereal* mss = y + m_non_spec_eq;
    m_thermo->setMassFractions_NoNorm(mss);
    m_thermo->setState_TP(m_T, m_P);

    // update surface
    updateSurfaceState(y + m_nsp + m_non_spec_eq);

    m_thermo->saveState(m_state);
}

void FlowReactor::setMassFlowRate(double mdot)
{
    m_rho0 = m_thermo->density();
    m_u = mdot/m_rho0;
    m_u0 = m_u;
    m_T = m_thermo->temperature();
    m_P0 = m_thermo->pressure() + m_rho0*m_u*m_u;
    m_h0 = m_thermo->enthalpy_mass() + 0.5*m_u*m_u;
}

void FlowReactor::updateSurfaceState(double* y)
{
    size_t loc = 0;
    for (auto& S : m_surfaces) {
        // set the coverages without normalization to avoid over-constraining
        // the system.
        // note: the ReactorSurface class doesn't normalize when calling setCoverages
        S->setCoverages(y+loc);
        S->syncCoverages();
        loc += S->thermo()->nSpecies();
    }
}

void FlowReactor::eval(double time, double* y,
                          double* ydot, double* params,
                          double* residual)
{
    m_thermo->restoreState(m_state);

    evalSurfaces(t, &m_sdot_temp[0]);
    const vector_fp& mw = m_thermo->molecularWeights();
    double sk_wk = 0;
    for (size_t i = 0; i < m_nsp; ++i)
    {
        sk_wk = m_sdot[i] * mw[i];
    }
    m_thermo->getPartialMolarEnthalpies(&m_hk[0]);
    // get net production
    if (m_chem) {
        m_kin->getNetProductionRates(&m_wdot[0]);
    }

    // set dphi/dz variables
    doublereal m_drhodz = ydot[0];
    doublereal m_dudz = ydot[1];
    doublereal m_dPdz = ydot[2];
    doublereal m_dTdz = ydot[3];

    // use equation of state for density residual
    const double R_specific = GasConstant / m_thermo->meanMolecularWeight(); // specific gas constant
    // P' - rho' * (R/M) * T - rho * (R/M) * T'
    residual[0] = m_dPdz - m_drhodz * R_specific * m_T - m_rho * R_specific * m_dTdz;
    // next, subtract off rho * T * sum(Yi' / Wi)
    for (size_t i = 0; i < m_nsp; ++i)
    {
        residual[0] -= m_rho * m_T * ydot[m_non_spec_eq + i] / mw[i];
    }

    //! use mass continuity for velocity residual
    //! Kee.'s Chemically Reacting Flow, Eq. 16.48
    const doublereal hydraulic = surfaceAreaToVolumeRatio();
    residual[1] = m_u * m_drhodz + m_rho * m_dudz - sk_wk * hydraulic;

    //! Use conservation of momentum for pressure residual
    //! Kee.'s Chemically Reacting Flow, Eq. 16.62
    //! [NOTE: friction is neglected]
    residual[2] = m_rho * m_u * m_dudz + m_u * hydraulic * sk_wk + m_dPdz;

    //! use conservation of energy for temperature residual
    //! Kee.'s Chemically Reacting Flow, Eq. 16.58
    //! [NOTE: adiabatic]
    //  Also, as above:
    //              h_mass = h_mole / Wk
    //       hence:
    //              h_mass * Wk = h_mol
    const doublereal cp_mass = m_thermo->cp_mass();
    residual[3] = m_rho * m_u * cp_mass * m_dTdz;
    for (size_t i = 0; i < m_nsp; ++i)
    {
        residual[3] += m_hk[i] * (m_wdot[i] + hydraulic * m_sdot[i]);
    }

    //! species conservation equations
    //! Kee.'s Chemically Reacting Flow, Eq. 16.51
    for (size_t i = 0; i < m_nsp; ++i)
    {
        residual[i + m_non_spec_eq] = m_rho * m_u * ydot[i + m_non_spec_eq] +
            y[i + m_non_spec_eq] * hydraulic * sk_wk -
            mw[i] * (m_wdot[i] + hydraulic * m_sdot[i]);
    }

    // surface algebraic constraints
    if (m_surfaces.size())
    {
        size_t loc = m_non_spec_eq + m_nsp; // offset into residual vector
        for (auto &m_surf : m_surfaces)
        {
            Kinetics* kin = m_surf->kinetics();
            SurfPhase* surf = m_surf->thermo();
            size_t nk = surf->nSpecies();
            kin->getNetProductionRates(&m_sdot_temp[0]);
            size_t ns = kin->surfacePhaseIndex();
            size_t surfloc = kin->kineticsSpeciesIndex(0,ns);
            doublereal sum = y[loc];
            for (size_t i = 1; i < nk; ++i)
            {
                //! net surface production rate residuals
                //! Kee.'s Chemically Reacting Flow, Eq. 16.63
                residual[loc + i] = m_sdot_temp[surfloc + i];
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

size_t FlowReactor::componentIndex(const string& nm) const
{
    // check for a gas species name
    size_t k = m_thermo->speciesIndex(nm);
    if (k != npos) {
        return k + m_non_spec_eq;
    } else if (nm == "rho" || nm == "density") {
        return 0;
    } else if (nm == "U" || nm == "velocity") {
        return 1;
    } else if (nm == "P" || nm == "pressure") {
        return 2;
    } else if (nm == "T" || nm == "temperature") {
        return 3;
    } else {
        return npos;
    }
}

}
