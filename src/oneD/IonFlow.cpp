//! @file IonFlow.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/oneD/IonFlow.h"
#include "cantera/oneD/StFlow.h"
#include "cantera/base/ctml.h"
#include "cantera/transport/TransportBase.h"
#include "cantera/numerics/funcs.h"
#include "cantera/oneD/Domain1D.h"


using namespace std;

namespace Cantera
{

IonFlow::IonFlow(IdealGasPhase* ph, size_t nsp, size_t points) :
    FreeFlame(ph, nsp, points),
    m_do_electric(false),
    m_solnPhase(1)
{
    // make a local copy of species charge
    for (size_t k = 0; k < m_nsp; k++) {
        m_speciesCharge.push_back(m_thermo->charge(k));
    }

    // Find indices for charge of species
    for (size_t k = 0; k < m_nsp; k++){
        if (m_speciesCharge[k] != 0){
            m_kCharge.push_back(k);
        } else {
            m_kNeutral.push_back(k);
        }
    }

    // Find the index of electron 
    if (m_thermo->speciesIndex("E") < m_nsp ) {
        m_kElectron = m_thermo->speciesIndex("E");
    }

    // mass fraction bounds (strict bound for ions)
    for (size_t k : m_kCharge) {
        setBounds(c_offset_Y+k, -1.0e-20, 1e-5);
    }
    setBounds(c_offset_P, -1.0e20, 1.0e20);
    m_refiner->setActive(c_offset_P, false);

    m_mobi.resize(m_nsp*m_points);
    m_do_poisson.resize(m_points,false);
    m_do_velocity.resize(m_points,true);
}

void IonFlow::resize(size_t components, size_t points){
    StFlow::resize(components, points);
    m_mobi.resize(m_nsp*m_points);
    m_do_species.resize(m_nsp,true);
    m_do_poisson.resize(m_points,false);
    m_do_velocity.resize(m_points,true);
    m_fixedMassFrac.resize(m_points*m_nsp);
    m_fixedElecPoten.resize(m_points,0.0);
    m_fixedVelocity.resize(m_points);
}

void IonFlow::updateTransport(doublereal* x, size_t j0, size_t j1)
{
    StFlow::updateTransport(x,j0,j1);
    for (size_t j = j0; j < j1; j++) {
        setGasAtMidpoint(x,j);
        m_trans->getMobilities(&m_mobi[j*m_nsp]);
        m_mobi[m_kElectron+m_nsp*j] = 0.4;
        m_diff[m_kElectron+m_nsp*j] = 0.4*(Boltzmann * T(x,j)) / ElectronCharge;
    }
}
void IonFlow::updateDiffFluxes(const doublereal* x, size_t j0, size_t j1)
{
    if (m_solnPhase == 1) {
        phaseOneDiffFluxes(x,j0,j1);
    } else if (m_solnPhase == 2) {
        phaseTwoDiffFluxes(x,j0,j1);
    } else {
        phaseThreeDiffFluxes(x,j0,j1);
    }
}

void IonFlow::phaseOneDiffFluxes(const doublereal* x, size_t j0, size_t j1)
{
    for (size_t j = j0; j < j1; j++) {
        double wtm = m_wtm[j];
        double rho = density(j);
        double dz = z(j+1) - z(j);
        double sum = 0.0;
        for (size_t k : m_kNeutral) {
            m_flux(k,j) = m_wt[k]*(rho*m_diff[k+m_nsp*j]/wtm);
            m_flux(k,j) *= (X(x,k,j) - X(x,k,j+1))/dz;
            sum -= m_flux(k,j);
        }

        // correction flux to insure that \sum_k Y_k V_k = 0.
        for (size_t k : m_kNeutral) {  
            m_flux(k,j) += sum*Y(x,k,j);
        }

        // flux for ions
        for (size_t k : m_kCharge) {
            m_flux(k,j) = 0;
        }
    }
}

void IonFlow::phaseTwoDiffFluxes(const doublereal* x, size_t j0, size_t j1)
{
    for (size_t j = j0; j < j1; j++) {
        double wtm = m_wtm[j];
        double rho = density(j);
        double dz = z(j+1) - z(j);
        // mixture-average diffusion
        double sum_flux = 0.0;
        for (size_t k = 0; k < m_nsp; k++) {       
            m_flux(k,j) = m_wt[k]*(rho*m_diff[k+m_nsp*j]/wtm);
            m_flux(k,j) *= (X(x,k,j) - X(x,k,j+1))/dz;
            sum_flux -= m_flux(k,j);
        }

        // ambipolar diffusion
        double sum_chargeFlux = 0.0;
        double sum = 0.0;
        for (size_t k : m_kCharge) {
            double Xav = 0.5 * (X(x,k,j+1) + X(x,k,j));
            int q_k = m_speciesCharge[k];
            sum_chargeFlux += m_speciesCharge[k] / m_wt[k] * m_flux(k,j);
            sum += m_mobi[k+m_nsp*j] * Xav * q_k * q_k;
        }
        double drift;
        double sum_drift = 0.0;
        for (size_t k : m_kCharge) {
            double Xav = 0.5 * (X(x,k,j+1) + X(x,k,j));
            int q_k = m_speciesCharge[k];
            drift = q_k * q_k * m_mobi[k+m_nsp*j] * Xav / sum;
            drift *= -sum_chargeFlux * m_wt[k] / q_k;
            m_flux(k,j) += drift;
            sum_drift -= drift;
        }

        // correction flux
        for (size_t k = 0; k < m_nsp; k++) {       
            m_flux(k,j) += Y(x,k,j) * sum_flux;
        }
    }
}

void IonFlow::phaseThreeDiffFluxes(const doublereal* x, size_t j0, size_t j1)
{
    for (size_t j = j0; j < j1; j++) {
        double wtm = m_wtm[j];
        double rho = density(j);
        double dz = z(j+1) - z(j);
        // mixture-average diffusion
        double sum = 0.0;
        for (size_t k = 0; k < m_nsp; k++) {       
            m_flux(k,j) = m_wt[k]*(rho*m_diff[k+m_nsp*j]/wtm);
            m_flux(k,j) *= (X(x,k,j) - X(x,k,j+1))/dz;
            sum -= m_flux(k,j); 
        }

        // correction flux
        for (size_t k = 0; k < m_nsp; k++) {       
            m_flux(k,j) += Y(x,k,j) * sum;
        }

        // ambipolar diffusion
        double drift;
        double E_ambi = E(x,j);
        sum = 0.0;
        for (size_t k : m_kCharge) {
            double Yav = 0.5 * (Y(x,k,j) + Y(x,k,j+1));
            drift = rho * Yav * E_ambi;
            drift *= m_speciesCharge[k] * m_mobi[k+m_nsp*j]; 
            m_flux(k,j) += drift;
            sum -= drift;
        }

        // correction drift
        for (size_t k : m_kCharge) {
            m_flux(k,j) += Y(x,k,j) * sum;
        }
    }   
}

void IonFlow::enableElectric(bool withElectric)
{
    m_do_electric = withElectric;
}

void IonFlow::setSolvingPhase(const size_t phase)
{
    m_solnPhase = phase;
}

void IonFlow::eval(size_t jg, doublereal* xg,
                  doublereal* rg, integer* diagg, doublereal rdt)
{
    StFlow::eval(jg, xg, rg, diagg, rdt);
    // start of local part of global arrays
    doublereal* x = xg + loc();
    doublereal* rsd = rg + loc();
    integer* diag = diagg + loc();

    size_t jmin, jmax;
    if (jg == npos) { // evaluate all points
        jmin = 0;
        jmax = m_points - 1;
    } else { // evaluate points for Jacobian
        size_t jpt = (jg == 0) ? 0 : jg - firstPoint();
        jmin = std::max<size_t>(jpt, 1) - 1;
        jmax = std::min(jpt+1,m_points-1);
    }
    // the boundary points are not applied
    for (size_t j = jmin; j <= jmax; j++) {
        if (j == 0) {
            rsd[index(c_offset_P, j)] = -phi(x,j);
            diag[index(c_offset_P, j)] = 0;
            for ( size_t k : m_kCharge) {
                rsd[index(c_offset_Y + k, j)] = Y(x,k,j); 
                diag[index(c_offset_Y + k, j)] = 0;
            }
        } else if (j == m_points - 1) {
            rsd[index(c_offset_P, j)] = -phi(x,j);
            diag[index(c_offset_P, j)] = 0;
            for ( size_t k : m_kCharge) {
                rsd[index(c_offset_Y + k, j)] = Y(x,k,j); 
                diag[index(c_offset_Y + k, j)] = 0;
            }
        } else {
            if (!m_do_velocity[j]) {
                rsd[index(c_offset_U, j)] = u(x,j) - u_fixed(j);
                diag[index(c_offset_U, j)] = 0;
            }
            for (size_t k = 0; k < m_nsp; k++) {
                if (!m_do_species[k]) {
                rsd[index(c_offset_Y + k, j)] = Y(x,k,j) - Y_fixed(k,j);
                rsd[index(c_offset_Y + k, j)] -= rdt*(Y(x,k,j) - Y_prev(k,j));
                diag[index(c_offset_Y + k, j)] = 1;            
                }
            }
        }
    }
    
    // convinent method due to interference
    for (size_t j = jmin; j <= jmax; j++) {
        if (j == 0) {
            rsd[index(c_offset_P, j)] = -phi(x,j);
            diag[index(c_offset_P, j)] = 0;
        } else if (j == m_points - 1) {
            rsd[index(c_offset_P, j)] = -phi(x,j);
            diag[index(c_offset_P, j)] = 0;
        } else {
            if (m_do_poisson[j]) {
                evalPoisson(j,x,rsd,diag,rdt);
            } else {
                rsd[index(c_offset_P, j)] = phi(x,j) - phi_fixed(j);
                diag[index(c_offset_P, j)] = 0;
            }
        }
    }
}

void IonFlow::evalPoisson(size_t j, doublereal* x, doublereal* rsd, integer* diag, doublereal rdt)
{
    //-----------------------------------------------
    //    Poisson's equation
    //
    //    dE/dz = e/eps_0 * sum(q_k*n_k)
    //
    //    E = -dV/dz 
    //-----------------------------------------------
    doublereal chargeDensity = 0.0;
    for (size_t k : m_kCharge) {
        chargeDensity += m_speciesCharge[k] * ElectronCharge * ND(x,k,j);
    }
    rsd[index(c_offset_P, j)] = dEdz(x,j) - chargeDensity / epsilon_0;
    diag[index(c_offset_P, j)] = 0;
}

void IonFlow::solveSpeciesEqn(size_t k)
{
    bool changed = false;
    if (k == npos) {
        for (size_t i = 0; i < m_nsp; i++) {
            if (!m_do_energy[i]) {
                changed = true;
            }
            m_do_species[i] = true;
        }
    } else {
        if (!m_do_species[k]) {
            changed = true;
        }
        m_do_species[k] = true;
    }
    if (changed) {
        needJacUpdate();
    }
}

void IonFlow::fixSpeciesMassFrac(size_t k)
{
    bool changed = false;
    if (k == npos) {
        for (size_t i = 0; i < m_nsp; i++) {
            if (m_do_species[i]) {
                changed = true;
            }
            m_do_species[i] = false;
        }
    } else {
        if (m_do_species[k]) {
            changed = true;
        }
        m_do_species[k] = false;
    }
    if (changed) {
        needJacUpdate();
    }
}

void IonFlow::solvePoissonEqn(size_t j)
{
    bool changed = false;
    if (j == npos) {
        for (size_t i = 0; i < m_points; i++) {
            if (!m_do_poisson[i]) {
                changed = true;
            }
            m_do_poisson[i] = true;
        }
    } else {
        if (!m_do_poisson[j]) {
            changed = true;
        }
        m_do_poisson[j] = true;
    }
    m_refiner->setActive(0, true);
    m_refiner->setActive(1, true);
    m_refiner->setActive(2, true);
    if (changed) {
        needJacUpdate();
    }
}

void IonFlow::fixElectricPotential(size_t j)
{
    bool changed = false;
    if (j == npos) {
        for (size_t i = 0; i < m_points; i++) {
            if (m_do_poisson[i]) {
                changed = true;
            }
            m_do_poisson[i] = false;
        }
    } else {
        if (m_do_poisson[j]) {
            changed = true;
        }
        m_do_poisson[j] = false;
    }
    m_refiner->setActive(0, false);
    m_refiner->setActive(1, false);
    m_refiner->setActive(2, false);
    if (changed) {
        needJacUpdate();
    }
}

void IonFlow::solveVelocity(size_t j)
{
    bool changed = false;
    if (j == npos) {
        for (size_t i = 0; i < m_points; i++) {
            if (!m_do_velocity[i]) {
                changed = true;
            }
            m_do_velocity[i] = true;
        }
    } else {
        if (!m_do_velocity[j]) {
            changed = true;
        }
        m_do_velocity[j] = true;
    }
    m_refiner->setActive(0, true);
    m_refiner->setActive(1, true);
    m_refiner->setActive(2, true);
    if (changed) {
        needJacUpdate();
    }
}

void IonFlow::fixVelocity(size_t j)
{
    bool changed = false;
    if (j == npos) {
        for (size_t i = 0; i < m_points; i++) {
            if (m_do_velocity[i]) {
                changed = true;
            }
            m_do_velocity[i] = false;
        }
    } else {
        if (m_do_velocity[j]) {
            changed = true;
        }
        m_do_velocity[j] = false;
    }
    m_refiner->setActive(0, false);
    m_refiner->setActive(1, false);
    m_refiner->setActive(2, false);
    if (changed) {
        needJacUpdate();
    }
}

void IonFlow::_finalize(const doublereal* x)
{
    FreeFlame::_finalize(x);

    for (size_t k = 0; k < m_nsp; k++) {
        bool y = m_do_species[k];
        if (!y) {
            for (size_t j = 0; j < m_points; j++) {
                m_fixedMassFrac[m_points*k+j] = Y(x,k,j);
            }
        }
    }

    // This method is still not tested
    // not sure why you want to return to original state
    // if not doing on point zero
    bool p = m_do_poisson[0];
    for (size_t j = 0; j < m_points; j++) {
        if (!p) {
            m_fixedElecPoten[j] = phi(x, j);
        }
    }
    if (p) {
        solvePoissonEqn();
    }

    bool v = m_do_velocity[0];
    for (size_t j = 0; j < m_points; j++) {
        if (!v) {
            m_fixedVelocity[j] = u(x,j);
        }
    }
    if (v) {
        solveVelocity();
    }
}

}