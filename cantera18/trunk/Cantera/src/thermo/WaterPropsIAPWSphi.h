/**
 * @file WaterPropsIAPWSphi.h
 *  Header for Lowest level of the classes which support a real water model
 *  (see class #WaterPropsIAPWS and class  #WaterPropsIAPWSphi).
 *
 *   This class calculates dimensionless quantitites.
 */
/*
 * Copywrite (2006) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
/*
 * $Id: WaterPropsIAPWSphi.h,v 1.8 2008/09/12 21:51:04 hkmoffa Exp $
 */

#ifndef WATERPROPSIAPWSPHI_H
#define WATERPROPSIAPWSPHI_H

#include "config.h"

/*!
 * the WaterPropsIAPSWSphi class support low level calls for
 * the real description of water.
 *
 *  The reference is W. Wagner, A. Prub, "The IAPWS Formulation 1995 for the Themodynamic
 *  Properties of Ordinary Water Substance for General and Scientific Use,"
 *  J. Phys. Chem. Ref. Dat, 31, 387, 2002.
 * 
 * Units Note: This class works with reduced units exclusively.
 */
class WaterPropsIAPWSphi {

public:

  //! Base constructor
  WaterPropsIAPWSphi();

  //! Calculate the Phi function, which is the base function
  /*!
   * The phi functino  is basically the helmholtz free energy
   * Eqn. (6.4)
   * All internal polynomials are recalculated.
   *
   * @param tau     Dimensionless temperature = T_c/T
   * @param delta   Dimensionless density =  delta = rho / Rho_c
   */
  double phi(double tau, double delta);

  //! Delta derivative of phi
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
   * note: this is done so much, we have a seperate routine.
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
  double dfind(double p_red, double tau, double deltaGuess);

  /**
   * Calculate the dimensionless gibbs free energy
   */
  double gibbs_RT() const;

  /**
   * Calculate the dimensionless enthalpy, h/RT
   */
  double enthalpy_RT() const;
    
  /**
   * Calculate the dimensionless entropy, s/R
   */
  double entropy_R() const;

  /**
   * Calculate the dimensionless internal energy, u/RT
   */
  double intEnergy_RT() const;

  /**
   * Calculate the dimensionless constant volume heat capacity, Cv/R
   */
  double cv_R() const;

  /**
   * Calculate the dimensionless constant pressure heat capacity, Cv/R
   */
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

  //! Return the value of phiR(), res
  double phiR() const;

private:

  //! nau calculation
  double phi0() const;
  //! calculation of d_phiR/d_d
  double phiR_d() const;
  //! calculation of d_nau/d_d
  double phi0_d() const;
  //! calculation of d2_res/d_dd
  double phiR_dd() const;
  //! calculation of d2_nau/d_dd
  double phi0_dd() const;
  //! calculation of d_nau/d_t
  double phi0_t() const;
  //! calculation of d_res/d_t
  double phiR_t() const;
  //! calculation of d2_res/d_tt
  double phiR_tt() const;
  //! calculation of d2_nau/d_tt
  double phi0_tt() const;
  //! calculation of d2_res/d_dt
  double phiR_dt() const;
  //! calculation of d2_nau/d_dt
  double phi0_dt() const;

  /**
   * intCheck() calculates all of the functions at a one point and
   * prints out the result. It's used for conducting the internal
   * check.
   *
   * @param tau     Dimensionless temperature = T_c/T
   * @param delta   Dimensionless density =  delta = rho / Rho_c
   */
  void intCheck(double tau, double delta);

private:

  //! Value of internally calculated polynomials of powers of TAU
  double TAUp[52];

 //! Value of internally calculated polynomials of powers of delta
  double DELTAp[16];

  //! Last tau that was used to calculate polynomials
  double TAUsave;

  //! sqrt of TAU
  double TAUsqrt;

  //! Last delta that was used to calculate polynomials
  double DELTAsave;
};
#endif
