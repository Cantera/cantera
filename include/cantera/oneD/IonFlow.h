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
 * concentration.
 * Reference:
 * Pederson, Timothy, and R. C. Brown.
 * "Simulation of electric field effects in premixed methane flames."
 * Combustion and Flames 94.4(1993): 433-448.
 * @ingroup onedim
 */
class IonFlow : public StFlow
{
public:
    IonFlow(IdealGasPhase* ph = 0, size_t nsp = 1, size_t points = 1);
    //! set the solving stage
    virtual void setSolvingStage(const size_t phase);

    virtual void resize(size_t components, size_t points);

    virtual void _finalize(const double* x);
    //! set to solve electric field on a point
    void solveElectricField(size_t j=npos);
    //! set to fix voltage on a point
    void fixElectricField(size_t j=npos);
    bool doElectricField(size_t j) {
        return m_do_electric_field[j];
    }

    /**
     * Sometimes it is desired to carry out the simulation using a specified
     * electron transport profile, rather than assuming it as a constant (0.4).
     * Reference:
     * Bisetti, Fabrizio, and Mbark El Morsli.
     * "Calculation and analysis of the mobility and diffusion coefficient
     * of thermal electrons in methane/air premixed flames."
     * Combustion and flame 159.12 (2012): 3518-3521.
     * If in the future the class GasTranport is improved, this method may
     * be discard. This method specifies this profile.
    */
    void setElectronTransport(vector_fp& tfix,
                              vector_fp& diff_e,
                              vector_fp& mobi_e);

protected:
    /*!
     * This function overloads the original function. The residual function
     * of electric field is added.
     */
    virtual void evalResidual(double* x, double* rsd, int* diag,
                              double rdt, size_t jmin, size_t jmax);
    virtual void updateTransport(double* x, size_t j0, size_t j1);
    virtual void updateDiffFluxes(const double* x, size_t j0, size_t j1);
    //! Solving phase one: the fluxes of charged species are turned off
    virtual void frozenIonMethod(const double* x, size_t j0, size_t j1);
    //! Solving phase two: the electric field equation is added coupled
    //! by the electrical drift
    virtual void electricFieldMethod(const double* x, size_t j0, size_t j1);
    //! flag for solving electric field or not
    std::vector<bool> m_do_electric_field;

    //! flag for importing transport of electron
    bool m_import_electron_transport;

    //! electrical properties
    vector_fp m_speciesCharge;

    //! index of species with charges
    std::vector<size_t> m_kCharge;

    //! index of neutral species
    std::vector<size_t> m_kNeutral;

    //! coefficients of polynomial fitting of fixed electron transport profile
    vector_fp m_mobi_e_fix;
    vector_fp m_diff_e_fix;

    //! mobility
    vector_fp m_mobility;

    //! solving stage
    size_t m_stage;

    //! index of electron
    size_t m_kElectron;

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
