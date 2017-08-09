//! @file IonFlow.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

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
 * The second stage uses charge neutrality model, which assume zero charge
 * flux throughout the domain, to calculate drift flux. The drift flux is
 * added to the total flux of ions.
 * Reference:
 * Prager, J., U. Riedel, and J. Warnatz.
 * "Modeling ion chemistry and charged species diffusion in lean
 * methane–oxygen flames."
 * Proceedings of the Combustion Institute 31.1 (2007): 1129-1137.
 *
 * The third stage evaluates drift flux from electric field calculated from
 * Poisson's equation, which is solved together with other equations. Poisson's
 * equation is coupled because the total charge densities depends on the species'
 * concentration.
 * Reference:
 * Pederson, Timothy, and R. C. Brown.
 * "Simulation of electric field effects in premixed methane flames."
 * Combustion and Flames 94.4(1993): 433-448.
 * @ingroup onedim
 */
class IonFlow : public FreeFlame
{
public:
    IonFlow(IdealGasPhase* ph = 0, size_t nsp = 1, size_t points = 1);
    //! set the solving stage
    virtual void setSolvingStage(const size_t phase);
    //! set electric voltage at inlet and outlet
    virtual void setElectricPotential(const double v1, const double v2);

    virtual void resize(size_t components, size_t points);

    virtual void _finalize(const double* x);
    //! set to solve Poisson's equation on a point
    void solvePoissonEqn(size_t j=npos);
    //! set to fix voltage on a point
    void fixElectricPotential(size_t j=npos);
    bool doPoisson(size_t j) {
        return m_do_poisson[j];
    }
    //! set to solve velocity on a point
    void solveVelocity(size_t j=npos);
    //! set to fix velocity on a point
    void fixVelocity(size_t j=npos);
    bool doVelocity(size_t j) {
        return m_do_velocity[j];
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
    void setElectronTransport(vector_fp& zfixed, vector_fp& diff_e_fixed,
                              vector_fp& mobi_e_fixed);

protected:
    virtual void evalResidual(double* x, double* rsd, int* diag,
                              double rdt, size_t jmin, size_t jmax);
    virtual void updateTransport(double* x, size_t j0, size_t j1);
    virtual void updateDiffFluxes(const double* x, size_t j0, size_t j1);
    //! Solving phase one: the fluxes of charged species are turned off
    virtual void frozenIonMethod(const double* x, size_t j0, size_t j1);
    //! Solving phase two: the Prager's ambipolar-diffusion model is used
    virtual void chargeNeutralityModel(const double* x, size_t j0, size_t j1);
    //! Solving phase three: the Poisson's equation is added coupled by the electrical drift
    virtual void poissonEqnMethod(const double* x, size_t j0, size_t j1);
    //! flag for solving poisson's equation or not
    std::vector<bool> m_do_poisson;
    //! flag for solving the velocity or not
    std::vector<bool> m_do_velocity;

    //! flag for importing transport of electron
    bool m_import_electron_transport;

    //! flag for overwrite transport of electron or not
    bool m_overwrite_eTransport;

    //! electrical properties
    vector_int m_speciesCharge;

    //! index of species with charges
    std::vector<size_t> m_kCharge;

    //! index of neutral species
    std::vector<size_t> m_kNeutral;

    //! fixed transport profile of electron
    vector_fp m_elecMobility;
    vector_fp m_elecDiffCoeff;

    //! mobility
    vector_fp m_mobility;

    //! solving stage
    int m_stage;

    //! The voltage
    double m_inletVoltage;
    double m_outletVoltage;

    //! index of electron
    size_t m_kElectron;

    //! fixed electric potential value
    vector_fp m_fixedElecPoten;

    //! fixed velocity value
    vector_fp m_fixedVelocity;

    //! fixed electron transport values
    vector_fp m_ztfix;
    vector_fp m_diff_e_fix;
    vector_fp m_mobi_e_fix;

    //! The fixed electric potential value at point j
    double phi_fixed(size_t j) const {
        return m_fixedElecPoten[j];
    }

    //! The fixed velocity value at point j
    double u_fixed(size_t j) const {
        return m_fixedVelocity[j];
    }

    //! electric potential
    double phi(const double* x, size_t j) const {
        return x[index(c_offset_P, j)];
    }

    //! electric field
    double E(const double* x, size_t j) const {
        return -(phi(x,j+1)-phi(x,j))/(z(j+1)-z(j));
    }

    double dEdz(const double* x, size_t j) const {
        return 2*(E(x,j)-E(x,j-1))/(z(j+1)-z(j-1));
    }

    //! number density
    double ND(const double* x, size_t k, size_t j) const {
        return Avogadro * m_rho[j] * Y(x,k,j) / m_wt[k];
    }

    //! total number density
    double ND_t(size_t j) const {
        return Avogadro * m_rho[j] / m_wtm[j];
    }
};

}

#endif
