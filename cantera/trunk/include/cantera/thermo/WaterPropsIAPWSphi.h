/**
 * @file WaterPropsIAPWSphi.h
 *  Header for Lowest level of the classes which support a real water model
 *  (see class \link Cantera::WaterPropsIAPWS WaterPropsIAPWS\endlink and  class \link WaterPropsIAPWSphi WaterPropsIAPWSphi\endlink).
 *
 *   This class calculates dimensionless quantities.
 */
/*
 * Copyright (2006) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
#ifndef WATERPROPSIAPWSPHI_H
#define WATERPROPSIAPWSPHI_H

#include "cantera/base/config.h"

/*!
 * the WaterPropsIAPSWSphi class support low level calls for
 * the real description of water.
 *
 *  The reference is W. Wagner, A. Prub, "The IAPWS Formulation 1995 for the Thermodynamic
 *  Properties of Ordinary Water Substance for General and Scientific Use,"
 *  J. Phys. Chem. Ref. Dat, 31, 387, 2002.
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
     * The phi function is basically the helmholtz free energy
     * Eqn. (6.4)
     * All internal polynomials are recalculated.
     *
     * @param tau     Dimensionless temperature = T_c/T
     * @param delta   Dimensionless density =  delta = rho / Rho_c
     */
    doublereal phi(doublereal tau, doublereal delta);

    //! Delta derivative of phi
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

    //! Internal check # 1
    void   check1();

    //! Internal check # 2
    void   check2();

    //! Calculate the dimensionless pressure at tau and delta;
    /*!
     *
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
     * This program computes the reduced density, given the reduced pressure
     * and the reduced temperature, tau. It takes an initial guess, deltaGuess.
     * DeltaGuess is important as this is a multivalued function below the
     * critical point.
     *
     * @param p_red       Value of the dimensionless pressure
     * @param tau         Dimensionless temperature = T_c/T
     * @param deltaGuess Initial guess for the dimensionless density
     *
     * @return
     *   Returns the dimensionless density.
     */
    doublereal dfind(doublereal p_red, doublereal tau, doublereal deltaGuess);

    /**
     * Calculate the dimensionless gibbs free energy
     */
    doublereal gibbs_RT() const;

    /**
     * Calculate the dimensionless enthalpy, h/RT
     */
    doublereal enthalpy_RT() const;

    /**
     * Calculate the dimensionless entropy, s/R
     */
    doublereal entropy_R() const;

    /**
     * Calculate the dimensionless internal energy, u/RT
     */
    doublereal intEnergy_RT() const;

    /**
     * Calculate the dimensionless constant volume heat capacity, Cv/R
     */
    doublereal cv_R() const;

    /**
     * Calculate the dimensionless constant pressure heat capacity, Cv/R
     */
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

    //! Return the value of phiR(), res
    doublereal phiR() const;

private:

    //! nau calculation
    doublereal phi0() const;
    //! calculation of d_phiR/d_d
    doublereal phiR_d() const;
    //! calculation of d_nau/d_d
    doublereal phi0_d() const;
    //! calculation of d2_res/d_dd
    doublereal phiR_dd() const;
    //! calculation of d2_nau/d_dd
    doublereal phi0_dd() const;
    //! calculation of d_nau/d_t
    doublereal phi0_t() const;
    //! calculation of d_res/d_t
    doublereal phiR_t() const;
    //! calculation of d2_res/d_tt
    doublereal phiR_tt() const;
    //! calculation of d2_nau/d_tt
    doublereal phi0_tt() const;
    //! calculation of d2_res/d_dt
    doublereal phiR_dt() const;
    //! calculation of d2_nau/d_dt
    doublereal phi0_dt() const;

    /**
     * intCheck() calculates all of the functions at a one point and
     * prints out the result. It's used for conducting the internal
     * check.
     *
     * @param tau     Dimensionless temperature = T_c/T
     * @param delta   Dimensionless density =  delta = rho / Rho_c
     */
    void intCheck(doublereal tau, doublereal delta);

private:

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
#endif
