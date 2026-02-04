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

FlowReactor::FlowReactor(shared_ptr<Solution> sol, const string& name)
    : FlowReactor(sol, true, name)
{
}

FlowReactor::FlowReactor(shared_ptr<Solution> sol, bool clone, const string& name)
    : IdealGasReactor(sol, clone, name)
{
    m_nv = 4 + m_nsp; // rho, u, P, T, and species mass fractions
    m_rho = m_thermo->density();
    // resize temporary arrays
    m_wdot.resize(m_nsp);
    m_hk.resize(m_nsp);
}

void FlowReactor::getStateDae(double* y, double* ydot)
{
    m_thermo->getMassFractions(span<double>(y + m_offset_Y, m_nsp));
    auto mw = m_thermo->molecularWeights();

    // set the first component to the initial density
    y[0] = m_thermo->density();

    if (m_u < 0) {
        throw CanteraError("FlowReactor::getStateDae",
            "Set mass flow rate before initializing reactor");
    }

    // set the second component to the initial speed
    y[1] = m_u;

    // set the third component to the initial pressure
    y[2] = m_thermo->pressure();

    // set the fourth component to the initial temperature
    y[3] = m_thermo->temperature();

    if (m_chem) {
        m_kin->getNetProductionRates(m_wdot.data()); // "omega dot"
    }

    // set the initial coverages
    updateSurfaceProductionRates();

    // reset ydot vector
    std::fill(ydot, ydot + m_nv, 0.0);

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
    a(0, 0) = - GasConstant * m_thermo->temperature() / m_thermo->meanMolecularWeight();
    a(0, 2) = 1;
    a(0, 3) = - m_rho * GasConstant / m_thermo->meanMolecularWeight();
    for (size_t i = 0; i < m_nsp; ++i) {
        a(0, m_offset_Y + i) = - m_rho * m_thermo->temperature() / mw[i];
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
    for (size_t i = 0; i < m_nsp; ++i) {
        h_sk_wk += m_sdot[i] * mw[i] / m_area;
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
        m_thermo->getPartialMolarEnthalpies(m_hk);
        for (size_t i = 0; i < m_nsp; ++i) {
            ydot[3] -= m_hk[i] * (m_wdot[i] + m_sdot[i] / m_area);
        }
    } else {
        ydot[3] = 0;
    }

    // mass-fraction equations Kee 16.51
    // - Yk * (perim / Ac) * sum(sk' * Wk) + wk' * Wk + (perim / Ac) * sk' * Wk
    for (size_t i = 0; i < m_nsp; ++i) {
        ydot[m_offset_Y + i] = -y[m_offset_Y + i] * h_sk_wk +
            mw[i] * (m_wdot[i] + m_sdot[i] / m_area);
    }

    // and solve
    solve(a, ydot, 1, 0);
}

void FlowReactor::updateState(double* y)
{
    // Set the mass fractions and density of the mixture.
    m_rho = y[0];
    m_u = y[1];
    m_thermo->setMassFractions_NoNorm(span<const double>(y + m_offset_Y, m_nsp));
    m_thermo->setState_TP(y[3], y[2]);
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

void FlowReactor::evalDae(double time, double* y, double* ydot, double* residual)
{
    updateSurfaceProductionRates();
    auto mw = m_thermo->molecularWeights();
    double sk_wk = 0;
    for (size_t i = 0; i < m_nsp; ++i) {
        sk_wk += m_sdot[i] * mw[i] / m_area;
    }
    m_thermo->getPartialMolarEnthalpies(m_hk);
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
    residual[1] = m_u * drhodz + m_rho * dudz - sk_wk;

    //! Use conservation of momentum for pressure residual
    //! Kee.'s Chemically Reacting Flow, Eq. 16.62
    //! [NOTE: friction is neglected]
    residual[2] = m_rho * m_u * dudz + m_u * sk_wk + dPdz;

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
            residual[3] += m_hk[i] * (m_wdot[i] + m_sdot[i] / m_area);
        }
    } else {
        residual[3] = dTdz;
    }

    //! species conservation equations
    //! Kee.'s Chemically Reacting Flow, Eq. 16.51
    double dSumYdz = 0;
    for (size_t i = 0; i < m_nsp; ++i) {
        residual[i + m_offset_Y] = m_rho * m_u * ydot[i + m_offset_Y] +
            y[i + m_offset_Y] * sk_wk -
            mw[i] * (m_wdot[i] + m_sdot[i] / m_area);
            dSumYdz += ydot[i + m_offset_Y];
    }
    // Spread d/dz(sum(Y)) = 0 constraint across all species equations. `scale` is
    // defined to make the size of the error in sum(Y) comparable to the overall rtol.
    // TODO: find a better way to access a value of rtol.
    double rtol = 1e-7;
    double scale = 0.1 * m_rho * m_u / rtol;
    for (size_t i = 0; i < m_nsp; ++i) {
        residual[i + m_offset_Y] += scale * std::max(0.0, y[i + m_offset_Y]) * dSumYdz;
    }
}

void FlowReactor::getConstraints(double* constraints) {
    // mark all variables differential equations unless otherwise specified
    std::fill(constraints, constraints + m_nv, 1.0);
}

size_t FlowReactor::componentIndex(const string& nm) const
{
    if (nm == "density") {
        return 0;
    } else if (nm == "speed") {
        return 1;
    } else if (nm == "pressure") {
        return 2;
    } else if (nm == "temperature") {
        return 3;
    }
    try {
        return m_thermo->speciesIndex(nm) + m_offset_Y;
    } catch (const CanteraError&) {
        throw CanteraError("FlowReactor::componentIndex",
            "Component '{}' not found", nm);
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
        return m_thermo->speciesName(k - 4);
    }
    throw IndexError("FlowReactor::componentName", "component", k, m_nv);
}

}
