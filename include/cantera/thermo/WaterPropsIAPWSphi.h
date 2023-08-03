/**
 * @file WaterPropsIAPWSphi.h
 * Header for Lowest level of the classes which support a real water model
 * (see class @link Cantera::WaterPropsIAPWS WaterPropsIAPWS@endlink and class
 * @link Cantera::WaterPropsIAPWSphi WaterPropsIAPWSphi@endlink).
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
 * This is a helper class for WaterSSTP and PDSS_Water and does not constitute
 * a complete implementation of a thermo phase by itself (see @ref thermoprops
 * and classes @link Cantera::WaterSSTP WaterSSTP@endlink and
 * @link Cantera::PDSS_Water PDSS_Water@endlink).
 *
 * The reference is Wagner and Pruss @cite wagner2002.
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
    double phi(double tau, double delta);

    //! Calculate derivative of phi wrt delta
    /*!
     * @param tau     Dimensionless temperature = T_c/T
     * @param delta   Dimensionless density =  delta = rho / Rho_c
     */
    double phi_d(double tau, double delta);

    //! 2nd derivative of phi wrt delta
    /*!
     * @param tau     Dimensionless temperature = T_c/T
     * @param delta   Dimensionless density =  delta = rho / Rho_c
     */
    double phi_dd(double tau, double delta);

    //! First derivative of phi wrt tau
    /*!
     * @param tau     Dimensionless temperature = T_c/T
     * @param delta   Dimensionless density =  delta = rho / Rho_c
     */
    double phi_t(double tau, double delta);

    //! Second derivative of phi wrt tau
    /*!
     * @param tau     Dimensionless temperature = T_c/T
     * @param delta   Dimensionless density =  delta = rho / Rho_c
     */
    double phi_tt(double tau, double delta);

    //! Calculate the dimensionless pressure at tau and delta;
    /*!
     *       pM/(rhoRT) = delta * phi_d() = 1.0 + delta phiR_d()
     *
     * @param tau     Dimensionless temperature = T_c/T
     * @param delta   Dimensionless density =  delta = rho / Rho_c
     *
     * note: this is done so much, we have a separate routine.
     */
    double pressureM_rhoRT(double tau, double delta);

    //! Dimensionless derivative of p wrt rho at constant T
    /*!
     *  dp/drho * 1/RT = (2. delta phi_d() + delta**2 phi_dd())
     *                   (1.0 + 2. delta phiR_d() + delta**2 phiR_dd())
     *
     * @param tau     Dimensionless temperature = T_c/T
     * @param delta   Dimensionless density =  delta = rho / Rho_c
     */
    double dimdpdrho(double tau, double delta);

    //! Dimensionless derivative of p wrt T at constant rho
    /*!
     *  dp/dT * M/(Rho R) = (1.0 + delta phiR_d()
     *                   -  tau delta (phiR_dt())
     *
     * @param tau     Dimensionless temperature = T_c/T
     * @param delta   Dimensionless density =  delta = rho / Rho_c
     */
    double dimdpdT(double tau, double delta);

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
    double dfind(double p_red, double tau, double deltaGuess);

    //! Calculate the dimensionless Gibbs free energy
    double gibbs_RT() const;

    //! Calculate the dimensionless enthalpy, h/RT
    double enthalpy_RT() const;

    //! Calculate the dimensionless entropy, s/R
    double entropy_R() const;

    //! Calculate the dimensionless internal energy, u/RT
    double intEnergy_RT() const;

    //! Calculate the dimensionless constant volume heat capacity, Cv/R
    double cv_R() const;

    //! Calculate the dimensionless constant pressure heat capacity, Cv/R
    double cp_R() const;

    //! Calculates internal polynomials in tau and delta.
    /*!
     * This routine is used to store the internal state of tau and delta
     * for later use by the other routines in the class.
     *
     * @param tau     Dimensionless temperature = T_c/T
     * @param delta   Dimensionless density =  delta = rho / Rho_c
     */
    void tdpolycalc(double tau, double delta);

    /**
     * Calculate Equation 6.6 for phiR, the residual part of the
     * dimensionless Helmholtz free energy.
     */
    double phiR() const;

protected:
    //! Calculate Equation 6.5 for phi0, the ideal gas part of the
    //! dimensionless Helmholtz free energy.
    double phi0() const;
    //! Calculate d_phiR_d(delta), the first derivative of phiR wrt delta
    double phiR_d() const;
    //! Calculate d_phi0_d(delta), the first derivative of phi0 wrt delta
    double phi0_d() const;
    //! Calculate d2_phiR_dd(delta), the second derivative of phiR wrt delta
    double phiR_dd() const;
    //! Calculate d2_phi0_dd(delta), the second derivative of phi0 wrt delta
    double phi0_dd() const;
    //! Calculate d_phi0/d(tau)
    double phi0_t() const;
    //! Calculate Equation 6.6 for dphiRdtau, the derivative residual part of
    //! the dimensionless Helmholtz free energy wrt temperature
    double phiR_t() const;
    //! Calculate Equation 6.6 for dphiRdtau, the second derivative residual
    //! part of the dimensionless Helmholtz free energy wrt temperature
    double phiR_tt() const;
    //! Calculate d2_phi0/dtau2
    double phi0_tt() const;
    //! Calculate the mixed derivative d2_phiR/(dtau ddelta)
    double phiR_dt() const;
    //! Calculate the mixed derivative d2_phi0/(dtau ddelta)
    double phi0_dt() const;

    //! Value of internally calculated polynomials of powers of TAU
    double TAUp[52];

    //! Value of internally calculated polynomials of powers of delta
    double DELTAp[16];

    //! Last tau that was used to calculate polynomials
    double TAUsave = -1.0;

    //! sqrt of TAU
    double TAUsqrt = -1.0;

    //! Last delta that was used to calculate polynomials
    double DELTAsave = -1.0;
};

} // namespace Cantera
#endif
