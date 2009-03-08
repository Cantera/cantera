/**
 * @file WaterPropsIAPWS.h
 *
 */
/*
 * Copywrite (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
/*
 * $Id: WaterPropsIAPWS.h,v 1.1 2006/07/04 00:01:53 hkmoffa Exp $
 */

#ifndef WATERPROPSIAPWS_H
#define WATERPROPSIAPWS_H

#include "WaterPropsIAPWSphi.h"

/*
 * These constants are defined and used in the interphase
 * to describe desired phases.
 */
#define WATER_GAS       0
#define WATER_LIQUID    1
#define WATER_SUPERCRIT 2

/**
 * Class for calculating the properties of water. 
 *
 *
 * Note, the base thermodynamic state for this class is the one
 * used in the steam tables, i.e., the liquid at the triple point
 * for water has the following properties:
 *
 *   u(273.16, rho)    = 0.0
 *   s(273.16, rho)    = 0.0 
 *   psat(273.16)      = 611.655 Pascal
 *   rho(273.16, psat) = 999.793 kg m-3
 *
 */
class WaterPropsIAPWS {
public:
  WaterPropsIAPWS();
  WaterPropsIAPWS(const WaterPropsIAPWS &b);
  WaterPropsIAPWS & operator=(const WaterPropsIAPWS &b);
  ~WaterPropsIAPWS();


  void setState(double temperature, double rho);
  
  /**
   * Calculate the Helmholtz free energy in mks units of
   * J kmol-1 K-1.
   */
  double helmholtzFE(double temperature, double rho);
  double helmholtzFE() const;
  
  /**
   * Calculate the Gibbs free energy in mks units of
   * J kmol-1 K-1.
   */
  double Gibbs(double temperature, double rho);
  double Gibbs() const;

  /**
   * Calculate the enthalpy in mks units of
   * J kmol-1 
   */
  double enthalpy(double temperature, double rho);
  double enthalpy() const;
 
  /**
   * Calculate the internal energy in mks units of
   * J kmol-1 
   */
  double intEnergy(double temperature, double rho);
  double intEnergy() const;

  /**
   * Calculate the entropy in mks units of 
   * J kmol-1 K-1
   */
  double entropy(double temperature, double rho);
  double entropy() const;
   
  /**
   * Calculate the constant volume heat capacity
   * in mks units of J kmol-1 K-1
   */
  double cv(double temperature, double rho);
  double cv() const;

  /**
   * Calculate the constant pressure heat capacity
   * in mks units of J kmol-1 K-1
   */
  double cp(double temperature, double rho);
  double cp() const;

  double molarVolume(double temperature, double rho);
  double molarVolume() const;

  /**
   * Calculate the pressure (Pascals), given the temperature and density
   *  Temperature: kelvin
   *  rho: density in kg m-3
   */
  double pressure(double temperature, double rho);
  double pressure() const;

  /*
   * Calculates the density given the temperature and the pressure,
   * and a guess at the density. Note, below T_c, this is a
   * multivalued function. 
   *
   * parameters:
   *    temperature: Kelvin
   *    pressure   : Pressure in Pascals (Newton/m**2)
   *    phase      : guessed phase of water 
   *               : -1: no guessed phase
   *    rhoguess   : guessed density of the water
   *               : -1.0 no guessed density
   */
  double density(double temperature, double pressure, 
		 int phase = -1, double rhoguess = -1.0);
  double density() const;

  /**
   * This function returns an estimated value for the saturation
   * pressure. It does this via a polynomial fit of the vapor pressure
   * curve. 
   * units = (Pascals)
   */
  double psat_est(double temperature);

  /**
   * Returns the coefficient of thermal expansion as a function
   * of temperature and pressure.
   *           alpha = d (ln V) / dT at constant P.
   *
   * 
   */
  double coeffThermExp(double temperature, double pressure);

  /**
   * Returns the coefficient of isothermal compressibility as a function
   * of temperature and pressure.
   *     kappa = - d (ln V) / dP at constant T.
   *
   *  units - 1/Pascal
   */
  double isothermalCompressibility(double temperature, double pressure);

  /**
   * Utility routine in the calculation of the saturation pressure
   */
  void corr(double temperature, double pressure, double &densLiq, 
	      double &densGas, double &delGRT);
  void corr1(double temperature, double pressure, double &densLiq, 
	       double &densGas, double &pcorr);

  /**
   * This function returns the saturation pressure given the
   * temperature as an input parameter.
   * units = Pascal
   */
  double psat(double temperature);

  double Tcrit() { return 647.096;}
  double Pcrit() { return 22.064E6;}

  double Rhocrit() { return 322.;}


private:
  /**
   * Calculate the dimensionless temp and rho and store internally.
   */
  void calcDim(double temperature, double rho);

  /*
   * Dimensionless versions of thermo functions. Note these are 
   * private, because R value is specific to the class. We only
   * show the dimensional functions in the interface.
   */
  double helmholtzFE_RT() const;
  double Gibbs_RT() const;
  double enthalpy_RT() const;
  double intEnergy_RT() const;
  double entropy_R() const; 
  double cv_R() const;
  double cp_R() const;
  double pressure_rhoRT() const;

protected:
  WaterPropsIAPWSphi *m_phi;
    
  double tau;
  double delta;
  int iState;
};
#endif

