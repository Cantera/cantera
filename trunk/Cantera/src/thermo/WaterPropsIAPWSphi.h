/**
 * @file WaterPropsIAPWSphi.h
 */
/*
 * Copywrite (2006) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
/*
 * $Id: WaterPropsIAPWSphi.h,v 1.1 2006/07/04 00:01:53 hkmoffa Exp $
 */

#ifndef WATERPROPSIAPWSPHI_H
#define WATERPROPSIAPWSPHI_H

/*
 * Units Note: This class works with reduced units exclusively.
 */

class WaterPropsIAPWSphi {

public:
  WaterPropsIAPWSphi();

  /*
   * Calculate the base phi's, recalculating the internal polynomials
   */
  double phi(double tau, double delta);
  double phi_d(double tau, double delta);  
  double phi_dd(double tau, double delta);
  double phi_t(double tau, double delta);
  double phi_tt(double tau, double delta);
  double phi_dt(double tau, double delta);

  void   check1();
  void   check2();

  /**
   * Calculate the dimensionless pressure, pred: 
   *       pred = pressure M / (rho RT)
   */
  double pressure_rhoRT(double tau, double delta);

  /**
   * This program computes the reduced density, given the reduced pressure
   * and the reduced temperature, tau. It takes an initial guess, deltaGuess.
   * DeltaGuess is important as this is a multivalued function below the
   * critical point.
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

  /**
   * Calculates internal polynomials in tau and delta. This
   * routine is used to store the internal state of tau and delta
   * for later use by the other routines in the class.
   */
  void tdpolycalc(double tau, double delta);

  double phiR() const;
private:
  double phi0() const;
  double phiR_d() const;
  double phi0_d() const;
  double phiR_dd() const;
  double phi0_dd() const;
  double phi0_t() const;
  double phiR_t() const;
  double phiR_tt() const;
  double phi0_tt() const;
  double phiR_dt() const;
  double phi0_dt() const;
  void intCheck(double tau, double delta);

protected:
  double TAUp[52];
  double DELTAp[16];
  double TAUsave;
  double TAUsqrt;
  double DELTAsave;
};
#endif
