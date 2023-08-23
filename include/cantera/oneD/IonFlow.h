//! @file IonFlow.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_IONFLOW_H
#define CT_IONFLOW_H

#include "cantera/oneD/StFlow.h"

namespace Cantera
{
/**
 * This class models the ion transportation in a flame. There are three
 * stages of the simulation.
 *
 * The first stage turns off the diffusion of ions due to the fast
 * diffusion rate of electron without internal electric forces (ambi-
 * polar diffusion effect).
 *
 * The second stage evaluates drift flux from electric field calculated from
 * Poisson's equation, which is solved together with other equations. Poisson's
 * equation is coupled because the total charge densities depends on the species'
 * concentration. See Pedersen and Brown @cite pedersen1993 for details.
 *
 * @ingroup flowGroup
 */
class IonFlow : public StFlow
{
public:
    IonFlow(ThermoPhase* ph = 0, size_t nsp = 1, size_t points = 1);

    //! Create a new flow domain.
    //! @param sol  Solution object used to evaluate all thermodynamic, kinetic, and
    //!     transport properties
    //! @param id  name of flow domain
    //! @param points  initial number of grid points
    IonFlow(shared_ptr<Solution> sol, const string& id="", size_t points = 1);

    string type() const override;

    size_t getSolvingStage() const override {
        return m_stage;
    }
    void setSolvingStage(const size_t stage) override;

    void resize(size_t components, size_t points) override;
    bool componentActive(size_t n) const override;

    void _finalize(const double* x) override;

    void solveElectricField(size_t j=npos) override;
    void fixElectricField(size_t j=npos) override;
    bool doElectricField(size_t j) const override {
        return m_do_electric_field[j];
    }

    /**
     * Sometimes it is desired to carry out the simulation using a specified
     * electron transport profile, rather than assuming it as a constant (0.4).
     * See Bisetti and El Morsli @cite bisetti2012.
     * If in the future the class GasTransport is improved, this method may
     * be discarded. This method specifies this profile.
    */
    void setElectronTransport(vector<double>& tfix,
                              vector<double>& diff_e,
                              vector<double>& mobi_e);

protected:
    /*!
     * This function overloads the original function. The residual function
     * of electric field is added.
     */
    void evalResidual(double* x, double* rsd, int* diag,
                      double rdt, size_t jmin, size_t jmax) override;
    void updateTransport(double* x, size_t j0, size_t j1) override;
    void updateDiffFluxes(const double* x, size_t j0, size_t j1) override;
    //! Solving phase one: the fluxes of charged species are turned off
    void frozenIonMethod(const double* x, size_t j0, size_t j1);
    //! Solving phase two: the electric field equation is added coupled
    //! by the electrical drift
    void electricFieldMethod(const double* x, size_t j0, size_t j1);
    //! flag for solving electric field or not
    vector<bool> m_do_electric_field;

    //! flag for importing transport of electron
    bool m_import_electron_transport = false;

    //! electrical properties
    vector<double> m_speciesCharge;

    //! index of species with charges
    vector<size_t> m_kCharge;

    //! index of neutral species
    vector<size_t> m_kNeutral;

    //! coefficients of polynomial fitting of fixed electron transport profile
    vector<double> m_mobi_e_fix;
    vector<double> m_diff_e_fix;

    //! mobility
    vector<double> m_mobility;

    //! solving stage
    size_t m_stage = 1;

    //! index of electron
    size_t m_kElectron = npos;

    //! electric field
    double E(const double* x, size_t j) const {
        return x[index(c_offset_E, j)];
    }

    double dEdz(const double* x, size_t j) const {
        return (E(x,j)-E(x,j-1))/(z(j)-z(j-1));
    }

    //! number density
    double ND(const double* x, size_t k, size_t j) const {
        return Avogadro * m_rho[j] * Y(x,k,j) / m_wt[k];
    }

    //! total charge density
    double rho_e(double* x, size_t j) const {
        double chargeDensity = 0.0;
        for (size_t k : m_kCharge) {
            chargeDensity += m_speciesCharge[k] * ElectronCharge * ND(x,k,j);
        }
        return chargeDensity;
    }
};

}

#endif
