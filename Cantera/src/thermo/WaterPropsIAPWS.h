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
 * $Id$
 */

#ifndef WATERPROPSIAPWS_H
#define WATERPROPSIAPWS_H

#include "WaterPropsIAPWSphi.h"

/**
 *  @name Names for the phase regions
 *
 *  These constants are defined and used in the interface
 *  to describe desired phases.
 */ 
//@{
#define WATER_GAS       0
#define WATER_LIQUID    1
#define WATER_SUPERCRIT 2
//@}

/*!
 * Class for calculating the properties of water. 
 *
 *
 * Note, the base thermodynamic state for this class is the one
 * used in the steam tables, i.e., the liquid at the triple point
 * for water has the following properties:
 *
 * -  u(273.16, rho)    = 0.0
 * -  s(273.16, rho)    = 0.0 
 * -  psat(273.16)      = 611.655 Pascal
 * -  rho(273.16, psat) = 999.793 kg m-3
 *
 */
class WaterPropsIAPWS {
public:

  //! Base constructor
  WaterPropsIAPWS();

  //! Copy constructor
  WaterPropsIAPWS(const WaterPropsIAPWS &);

  //! assignment constructor
  WaterPropsIAPWS & operator=(const WaterPropsIAPWS &);

  //! destructor
  ~WaterPropsIAPWS();

  //! Set the internal state of the object wrt temperature and density
  /*!
   * @param temperature   temperature (kelvin)
   * @param rho           density  (kg m-3)
   */
  void setState(double temperature, double rho);
  
  
  //! Calculate the Helmholtz free energy in mks units of J kmol-1 K-1.
  /*!
   * @param temperature   temperature (kelvin)
   * @param rho           density  (kg m-3)
   */
  double helmholtzFE(double temperature, double rho);

  //! Calculate the Helmholtz free energy in mks units of J kmol-1 K-1, 
  //! using the last temperature and density
  double helmholtzFE() const;
  
  
  //! Calculate the Gibbs free energy in mks units of J kmol-1 K-1.
  /*!
   * @param temperature   temperature (kelvin)
   * @param rho           density  (kg m-3)
   */
  double Gibbs(double temperature, double rho);

  //! Calculate the Gibbs free energy in mks units of J kmol-1 K-1.
  //! using the last temperature and density
  double Gibbs() const;

  
  //!  Calculate the enthalpy in mks units of  J kmol-1 
  /*!
   * @param temperature   temperature (kelvin)
   * @param rho           density  (kg m-3)
   */
  double enthalpy(double temperature, double rho);

  //!  Calculate the enthalpy in mks units of  J kmol-1 
  //!  using the last temperature and density
  double enthalpy() const;
 
  //! Calculate the internal energy in mks units of J kmol-1 
  /*!
   * @param temperature   temperature (kelvin)
   * @param rho           density  (kg m-3)
   */
  double intEnergy(double temperature, double rho);

  //! Calculate the internal energy in mks units of J kmol-1 
  //! at the last internal energy
  double intEnergy() const;
  
  //! Calculate the entropy in mks units of  J kmol-1 K-1
  /*!
   * @param temperature   temperature (kelvin)
   * @param rho           density  (kg m-3)
   */
  double entropy(double temperature, double rho);

  //! Calculate the entropy in mks units of  J kmol-1 K-1
  //! at the last temperature and density
  double entropy() const;
   
  
  //! Calculate the constant volume heat capacity in mks units of J kmol-1 K-1
  /*!
   * @param temperature   temperature (kelvin)
   * @param rho           density  (kg m-3)
   */
  double cv(double temperature, double rho);

  //! Calculate the constant volume heat capacity in mks units of J kmol-1 K-1
  //! at the last temperature and density
  double cv() const;

  //! Calculate the constant pressure heat capacity in mks units of J kmol-1 K-1
  /*!
   * @param temperature   temperature (kelvin)
   * @param rho           density  (kg m-3)
   */
  double cp(double temperature, double rho);

  //! Calculate the constant pressure heat capacity in mks units of J kmol-1 K-1
  //! at the last temperature and density
  double cp() const;

  //! Calculate the molar volume (kmol m-3)
  /*!
   * @param temperature   temperature (kelvin)
   * @param rho           density  (kg m-3)
   */
  double molarVolume(double temperature, double rho);

  //! Calculate the molar volume (kmol m-3) 
  //! at the last temperature and density
  double molarVolume() const;
 
  //! Calculate the pressure (Pascals), given the temperature and density
  /*!
   *  @param temperature  input temperature kelvin
   *  @param rho          density in kg m-3
   *
   *  @return
   *    returns the pressure (Pascal)
   */
  double pressure(double temperature, double rho);

  //! Calculates the pressure (Pascals), given the current value of the
  //! temperature and density.
  /*!
   * The density is an independent variable in the underlying equation of state
   */
  double pressure() const;

  //! Calculates the density given the temperature and the pressure,
  //! and a guess at the density. 
  /*!
   *  Note, below T_c, this is a multivalued function. 
   *
   *
   *  @param  temperature: Kelvin
   *  @param  pressure   : Pressure in Pascals (Newton/m**2)
   *  @param  phase      : guessed phase of water 
   *                     : -1: no guessed phase
   *  @param rhoguess    : guessed density of the water
   *                     : -1.0 no guessed density
   *  @return
   *     Returns the density
   */
  double density(double temperature, double pressure, 
		 int phase = -1, double rhoguess = -1.0);

  //! Returns the density (kg m-3)
  /*!
   * The density is an independent variable in the underlying equation of state
   */
  double density() const;

  
  //! This function returns an estimated value for the saturation pressure. 
  /*!
   * It does this via a polynomial fit of the vapor pressure curve. 
   * units = (Pascals)
   *
   * @param temperature Input temperature (Kelvin)
   *
   * @return
   *   Returns the estimated saturation pressure
   */
  double psat_est(double temperature);

  //! Returns the coefficient of thermal expansion as a function of temperature and pressure.
  /*!
   *           alpha = d (ln V) / dT at constant P.
   *
   * @param temperature Input temperature (Kelvin)
   * @param pressure    Input pressure (Pa)
   * @return
   *    Returns the coefficient of thermal expansion
   */
  double coeffThermExp(double temperature, double pressure);

 
  //! Returns the coefficient of isothermal compressibility as a function
  //! of temperature and pressure.
  /*!
   *     kappa = - d (ln V) / dP at constant T.
   *
   *  units - 1/Pascal
   *
   * @param temperature Input temperature (Kelvin)
   * @param pressure    Input pressure (Pa)
   * @return 
   *    returns the isothermal compressibility
   */
  double isothermalCompressibility(double temperature, double pressure);

  
  //! Utility routine in the calculation of the saturation pressure
  /*!
   * @param temperature    temperature (kelvin)
   * @param pressure       pressure (Pascal)
   * @param densLiq        Output density of liquid
   * @param densGas        output Density of gas
   * @param delGRT         output delGRT
   */ 
  void corr(double temperature, double pressure, double &densLiq, 
	      double &densGas, double &delGRT);

  //! Utility routine in the calculation of the saturation pressure
  /*!
   * @param temperature    temperature (kelvin)
   * @param pressure       pressure (Pascal)
   * @param densLiq        Output density of liquid
   * @param densGas        output Density of gas
   * @param pcorr          output corrected pressure
   */ 
  void corr1(double temperature, double pressure, double &densLiq, 
	       double &densGas, double &pcorr);

 
  //!  This function returns the saturation pressure given the
  //! temperature as an input parameter.
  /*!
   * @param temperature   input temperature (kelvin)
   * @return 
   * Returns the saturation pressure
   *                units = Pascal
   */
  double psat(double temperature);

  //! Returns the critical temperature of water (Kelvin)
  /*!
   *  This is hard coded to the value 647.096 Kelvin
   */
  double Tcrit() { return 647.096;}

  //! Returns the critical pressure of water (22.064E6 Pa)
  /*!
   *  This is hard coded to the value of 22.064E6 pascals
   */
  double Pcrit() { return 22.064E6;}

  //! Return the critical density of water (kg m-3)
  /*!
   * This is equal to 322 kg m-3.
   */
  double Rhocrit() { return 322.;}

private:
  /**
   * Calculate the dimensionless temp and rho and store internally.
   *
   * @param temperature   input temperature (kelvin)
   *  @param rho          density in kg m-3
   */
  void calcDim(double temperature, double rho);

  /*
   * Dimensionless versions of thermo functions. Note these are 
   * private, because R value is specific to the class. We only
   * show the dimensional functions in the interface.
   */
  double helmholtzFE_RT() const;

  //! Returns the dimensionless gibbs free energy
  double Gibbs_RT() const;

  //! Returns the dimensionless enthalpy
  double enthalpy_RT() const;

  //! Returns the dimensionless internal energy
  double intEnergy_RT() const;

  //! Returns the dimensionless entropy
  double entropy_R() const; 

  //! Returns the dimensionless heat capacity at constant volume
  double cv_R() const;

  //! Returns the dimensionless heat capacity at constant pressure
  double cp_R() const;

  //! Return the current dimensionless pressure
  double pressure_rhoRT() const;

protected:

  //! pointer to the underlying object that does the calculations.
  WaterPropsIAPWSphi *m_phi;
    
  //! Dimensionless temperature
  /*!
   *   tau = T_C / T
   */
  double tau;

  //! Dimensionless density
  /*!
   *  delta = rho / rho_c
   */
  double delta;

  //! Current state of the system
  int iState;
};
#endif

