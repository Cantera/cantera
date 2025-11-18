//! @file IonFlow.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/oneD/IonFlow.h"
#include "cantera/oneD/Flow1D.h"
#include "cantera/oneD/refine.h"
#include "cantera/transport/Transport.h"
#include "cantera/numerics/funcs.h"
#include "cantera/numerics/polyfit.h"
#include "cantera/base/utilities.h"
#include "cantera/base/global.h"

namespace Cantera
{

void IonFlow::_init(ThermoPhase* ph, size_t nsp, size_t points)
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
    if (m_thermo->speciesIndex("E", false) != npos ) {
        m_kElectron = m_thermo->speciesIndex("E", true);
    }

    // no bound for electric potential
    setBounds(c_offset_E, -1.0e20, 1.0e20);

    // Set tighter negative species limit on charged species to avoid
    // instabilities. Tolerance on electrons is even tighter to account for the
    // low "molecular" weight.
    for (size_t k : m_kCharge) {
        setBounds(c_offset_Y + k, -1e-10, 1.0);
    }
    setBounds(c_offset_Y + m_kElectron, -1e-14, 1.0);

    m_refiner->setActive(c_offset_E, false);
    m_mobility.resize(m_nsp*m_points);
}

IonFlow::IonFlow(shared_ptr<Solution> phase, const string& id, size_t points)
    : Flow1D(phase, id, points)
{
    _init(phase->thermo().get(), phase->thermo()->nSpecies(), points);
    m_solution = phase;
    m_solution->thermo()->addSpeciesLock();
    m_id = id;
    m_kin = m_solution->kinetics().get();
    m_trans = m_solution->transport().get();
    if (m_trans->transportModel() == "none") {
        throw CanteraError("IonFlow::IonFlow",
            "An appropriate transport model\nshould be set when instantiating the "
            "Solution ('gas') object.");
    }
    m_solution->registerChangedCallback(this, [this]() {
        _setKinetics(m_solution->kinetics());
        _setTransport(m_solution->transport());
    });
}

string IonFlow::domainType() const {
    if (m_isFree) {
        return "free-ion-flow";
    }
    if (m_usesLambda) {
        return "axisymmetric-ion-flow";
    }
    return "unstrained-ion-flow";
}

void IonFlow::resize(size_t components, size_t points){
    Flow1D::resize(components, points);
    m_mobility.resize(m_nsp*m_points);
}

bool IonFlow::componentActive(size_t n) const
{
    if (n == c_offset_E) {
        return true;
    } else {
        return Flow1D::componentActive(n);
    }
}

void IonFlow::updateTransport(double* x, size_t j0, size_t j1)
{
    Flow1D::updateTransport(x,j0,j1);
    for (size_t j = j0; j < j1; j++) {
        setGasAtMidpoint(x,j);
        m_trans->getMobilities(&m_mobility[j*m_nsp]);
        if (m_import_electron_transport) {
            size_t k = m_kElectron;
            double tlog = log(m_thermo->temperature());
            m_mobility[k+m_nsp*j] = poly5(tlog, m_mobi_e_fix.data());
            double rho = m_thermo->density();
            double wtm = m_thermo->meanMolecularWeight();
            m_diff[k+m_nsp*j] = m_wt[k]*rho*poly5(tlog, m_diff_e_fix.data())/wtm;
        }
    }
}

void IonFlow::updateDiffFluxes(const double* x, size_t j0, size_t j1)
{
    if (m_do_electric_field) {
        electricFieldMethod(x,j0,j1);
    } else {
        frozenIonMethod(x,j0,j1);
    }
}

void IonFlow::frozenIonMethod(const double* x, size_t j0, size_t j1)
{
    for (size_t j = j0; j < j1; j++) {
        double dz = z(j+1) - z(j);
        double sum = 0.0;
        for (size_t k : m_kNeutral) {
            m_flux(k,j) = m_diff[k+m_nsp*j];
            m_flux(k,j) *= (X(x,k,j) - X(x,k,j+1))/dz;
            sum -= m_flux(k,j);
        }

        // correction flux to insure that \sum_k Y_k V_k = 0.
        for (size_t k : m_kNeutral) {
            m_flux(k,j) += sum*Y(x,k,j);
        }

        // flux for ions
        // Set flux to zero to prevent some fast charged species (such electrons)
        // to run away
        for (size_t k : m_kCharge) {
            m_flux(k,j) = 0;
        }
    }
}

void IonFlow::electricFieldMethod(const double* x, size_t j0, size_t j1)
{
    for (size_t j = j0; j < j1; j++) {
        double rho = density(j);
        double dz = z(j+1) - z(j);

        // mixture-average diffusion
        for (size_t k = 0; k < m_nsp; k++) {
            m_flux(k,j) = m_diff[k+m_nsp*j];
            m_flux(k,j) *= (X(x,k,j) - X(x,k,j+1))/dz;
        }

        // ambipolar diffusion
        double E_ambi = E(x,j);
        for (size_t k : m_kCharge) {
            double Yav = 0.5 * (Y(x,k,j) + Y(x,k,j+1));
            double drift = rho * Yav * E_ambi
                           * m_speciesCharge[k] * m_mobility[k+m_nsp*j];
            m_flux(k,j) += drift;
        }

        // correction flux
        double sum_flux = 0.0;
        for (size_t k = 0; k < m_nsp; k++) {
            sum_flux -= m_flux(k,j); // total net flux
        }
        double sum_ion = 0.0;
        for (size_t k : m_kCharge) {
            sum_ion += Y(x,k,j);
        }
        // The portion of correction for ions is taken off
        for (size_t k : m_kNeutral) {
            m_flux(k,j) += Y(x,k,j) / (1-sum_ion) * sum_flux;
        }
    }
}

//! Evaluate the electric field equation residual
void IonFlow::evalElectricField(double* x, double* rsd, int* diag,
                                double rdt, size_t jmin, size_t jmax)
{
    Flow1D::evalElectricField(x, rsd, diag, rdt, jmin, jmax);
    if (!m_do_electric_field) {
        return;
    }

    if (jmin == 0) { // left boundary
        rsd[index(c_offset_E, jmin)] = E(x,jmin);
    }

    if (jmax == m_points - 1) { // right boundary
        rsd[index(c_offset_E, jmax)] = dEdz(x,jmax) - rho_e(x,jmax) / epsilon_0;
    }

    // j0 and j1 are constrained to only interior points
    size_t j0 = std::max<size_t>(jmin, 1);
    size_t j1 = std::min(jmax, m_points - 2);
    for (size_t j = j0; j <= j1; j++) {
        rsd[index(c_offset_E, j)] = dEdz(x,j) - rho_e(x,j) / epsilon_0;
        diag[index(c_offset_E, j)] = 0;
    }
}

void IonFlow::evalSpecies(double* x, double* rsd, int* diag,
                          double rdt, size_t jmin, size_t jmax)
{
    Flow1D::evalSpecies(x, rsd, diag, rdt, jmin, jmax);
    if (!m_do_electric_field) {
        return;
    }

    if (jmin == 0) { // left boundary
        // enforcing the flux for charged species is difficult
        // since charged species are also affected by electric
        // force, so Neumann boundary condition is used.
        for (size_t k : m_kCharge) {
            rsd[index(c_offset_Y + k, jmin)] = Y(x,k,jmin) - Y(x,k,jmin + 1);
        }
    }
}

void IonFlow::solveElectricField()
{
    if (!m_do_electric_field) {
        needJacUpdate();
    }
    m_refiner->setActive(c_offset_U, true);
    m_refiner->setActive(c_offset_V, true);
    m_refiner->setActive(c_offset_T, true);
    m_refiner->setActive(c_offset_E, true);
    m_do_electric_field = true;
}

void IonFlow::fixElectricField()
{
    if (m_do_electric_field) {
        needJacUpdate();
    }
    m_refiner->setActive(c_offset_U, false);
    m_refiner->setActive(c_offset_V, false);
    m_refiner->setActive(c_offset_T, false);
    m_refiner->setActive(c_offset_E, false);
    m_do_electric_field = false;
}

void IonFlow::setElectronTransport(vector<double>& tfix, vector<double>& diff_e,
                                   vector<double>& mobi_e)
{
    m_import_electron_transport = true;
    size_t degree = 5;
    size_t n = tfix.size();
    vector<double> tlog;
    for (size_t i = 0; i < n; i++) {
        tlog.push_back(log(tfix[i]));
    }
    vector<double> w(n, -1.0);
    m_diff_e_fix.resize(degree + 1);
    m_mobi_e_fix.resize(degree + 1);
    polyfit(n, degree, tlog.data(), diff_e.data(), w.data(), m_diff_e_fix.data());
    polyfit(n, degree, tlog.data(), mobi_e.data(), w.data(), m_mobi_e_fix.data());
}

}
