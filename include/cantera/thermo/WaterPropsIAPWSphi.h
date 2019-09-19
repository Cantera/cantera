/**
 * @file WaterPropsIAPWSphi.h
 * Header for Lowest level of the classes which support a real water model
 * (see class \link Cantera::WaterPropsIAPWS WaterPropsIAPWS\endlink and class
 * \link Cantera::WaterPropsIAPWSphi WaterPropsIAPWSphi\endlink).
 *
 * This class calculates dimensionless quantities.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef WATERPROPSIAPWSPHI_H
#define WATERPROPSIAPWSPHI_H

#include "cantera/base/config.h"

namespace Cantera
{

//! Low level class for the real description of water.
/*!
 * The reference is W. Wagner, A. Pruss, "The IAPWS Formulation 1995 for the
 * Thermodynamic Properties of Ordinary Water Substance for General and
 * Scientific Use," J. Phys. Chem. Ref. Dat, 31, 387, 2002.
 *
 * Units Note: This class works with reduced units exclusively.
 */
class WaterPropsIAPWSphi
{
public:
    //! Base constructor
    WaterPropsIAPWSphi();

    //! Calculate the Phi function, which is the base function
    /*!
     * The phi function is basically the Helmholtz free energy Eqn. (6.4) All
     * internal polynomials are recalculated.
     *
     * @param tau     Dimensionless temperature = T_c/T
     * @param delta   Dimensionless density =  delta = rho / Rho_c
     */
    doublereal phi(doublereal tau, doublereal delta);

    //! Calculate derivative of phi wrt delta
    /*!
     * @param tau     Dimensionless temperature = T_c/T
     * @param delta   Dimensionless density =  delta = rho / Rho_c
     */
    doublereal phi_d(doublereal tau, doublereal delta);

    //! 2nd derivative of phi wrt delta
    /*!
     * @param tau     Dimensionless temperature = T_c/T
     * @param delta   Dimensionless density =  delta = rho / Rho_c
     */
    doublereal phi_dd(doublereal tau, doublereal delta);

    //! First derivative of phi wrt tau
    /*!
     * @param tau     Dimensionless temperature = T_c/T
     * @param delta   Dimensionless density =  delta = rho / Rho_c
     */
    doublereal phi_t(doublereal tau, doublereal delta);

    //! Second derivative of phi wrt tau
    /*!
     * @param tau     Dimensionless temperature = T_c/T
     * @param delta   Dimensionless density =  delta = rho / Rho_c
     */
    doublereal phi_tt(doublereal tau, doublereal delta);

    //! Calculate the dimensionless pressure at tau and delta;
    /*!
     *       pM/(rhoRT) = delta * phi_d() = 1.0 + delta phiR_d()
     *
     * @param tau     Dimensionless temperature = T_c/T
     * @param delta   Dimensionless density =  delta = rho / Rho_c
     *
     * note: this is done so much, we have a separate routine.
     */
    doublereal pressureM_rhoRT(doublereal tau, doublereal delta);

    //! Dimensionless derivative of p wrt rho at constant T
    /*!
     *  dp/drho * 1/RT = (2. delta phi_d() + delta**2 phi_dd())
     *                   (1.0 + 2. delta phiR_d() + delta**2 phiR_dd())
     *
     * @param tau     Dimensionless temperature = T_c/T
     * @param delta   Dimensionless density =  delta = rho / Rho_c
     */
    doublereal dimdpdrho(doublereal tau, doublereal delta);

    //! Dimensionless derivative of p wrt T at constant rho
    /*!
     *  dp/dT * M/(Rho R) = (1.0 + delta phiR_d()
     *                   -  tau delta (phiR_dt())
     *
     * @param tau     Dimensionless temperature = T_c/T
     * @param delta   Dimensionless density =  delta = rho / Rho_c
     */
    doublereal dimdpdT(doublereal tau, doublereal delta);

    /**
     * This function computes the reduced density, given the reduced pressure
     * and the reduced temperature, tau. It takes an initial guess, deltaGuess.
     * DeltaGuess is important as this is a multivalued function below the
     * critical point.
     *
     * @param p_red       Value of the dimensionless pressure
     * @param tau         Dimensionless temperature = T_c/T
     * @param deltaGuess  Initial guess for the dimensionless density
     *
     * @returns the dimensionless density.
     */
    doublereal dfind(doublereal p_red, doublereal tau, doublereal deltaGuess);

    //! Calculate the dimensionless Gibbs free energy
    doublereal gibbs_RT() const;

    //! Calculate the dimensionless enthalpy, h/RT
    doublereal enthalpy_RT() const;

    //! Calculate the dimensionless entropy, s/R
    doublereal entropy_R() const;

    //! Calculate the dimensionless internal energy, u/RT
    doublereal intEnergy_RT() const;

    //! Calculate the dimensionless constant volume heat capacity, Cv/R
    doublereal cv_R() const;

    //! Calculate the dimensionless constant pressure heat capacity, Cv/R
    doublereal cp_R() const;

    //! Calculates internal polynomials in tau and delta.
    /*!
     * This routine is used to store the internal state of tau and delta
     * for later use by the other routines in the class.
     *
     * @param tau     Dimensionless temperature = T_c/T
     * @param delta   Dimensionless density =  delta = rho / Rho_c
     */
    void tdpolycalc(doublereal tau, doublereal delta);

    /*!
     * Calculate Equation 6.6 for phiR, the residual part of the
     * dimensionless Helmholtz free energy.
     */
    doublereal phiR() const;

protected:
    //! Calculate Equation 6.5 for phi0, the ideal gas part of the
    //! dimensionless Helmholtz free energy.
    doublereal phi0() const;
    //! Calculate d_phiR_d(delta), the first derivative of phiR wrt delta
    doublereal phiR_d() const;
    //! Calculate d_phi0_d(delta), the first derivative of phi0 wrt delta
    doublereal phi0_d() const;
    //! Calculate d2_phiR_dd(delta), the second derivative of phiR wrt delta
    doublereal phiR_dd() const;
    //! Calculate d2_phi0_dd(delta), the second derivative of phi0 wrt delta
    doublereal phi0_dd() const;
    //! Calculate d_phi0/d(tau)
    doublereal phi0_t() const;
    //! Calculate Equation 6.6 for dphiRdtau, the derivative residual part of
    //! the dimensionless Helmholtz free energy wrt temperature
    doublereal phiR_t() const;
    //! Calculate Equation 6.6 for dphiRdtau, the second derivative residual
    //! part of the dimensionless Helmholtz free energy wrt temperature
    doublereal phiR_tt() const;
    //! Calculate d2_phi0/dtau2
    doublereal phi0_tt() const;
    //! Calculate the mixed derivative d2_phiR/(dtau ddelta)
    doublereal phiR_dt() const;
    //! Calculate the mixed derivative d2_phi0/(dtau ddelta)
    doublereal phi0_dt() const;

    //! Value of internally calculated polynomials of powers of TAU
    doublereal TAUp[52];

    //! Value of internally calculated polynomials of powers of delta
    doublereal DELTAp[16];

    //! Last tau that was used to calculate polynomials
    doublereal TAUsave;

    //! sqrt of TAU
    doublereal TAUsqrt;

    //! Last delta that was used to calculate polynomials
    doublereal DELTAsave;
};

} // namespace Cantera
#endif
