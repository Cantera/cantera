/**
 * @file WaterPropsIAPWS.h
 * Headers for a class for calculating the equation of state of water
 * from the IAPWS 1995 Formulation based on the steam tables thermodynamic
 * basis (See class \link WaterPropsIAPWS WaterPropsIAPWS\endlink).
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

//! Class for calculating the equation of state of water. 
/*!
 *
 *  The reference is W. Wagner, A. Prub, "The IAPWS Formulation 1995 for the Themodynamic
 *  Properties of Ordinary Water Substance for General and Scientific Use,"
 *  J. Phys. Chem. Ref. Dat, 31, 387, 2002.
 *
 * This class provides a very complicated polynomial for the specific helmholtz free
 * energy of water, as a function of temperature and density.
 *
 *    \f[
 *        \frac{M\hat{f}(\rho,T)}{R T} = \phi(\delta, \tau) =
 *                        \phi^o(\delta, \tau) +  \phi^r(\delta, \tau)
 *    \f]
 *
 *  where
 *
 *    \f[
 *         \delta = \rho / \rho_c \mbox{\qquad and \qquad} \tau = T_c / T 
 *    \f]
 *
 *  The following constants are assumed
 *
 *     \f[
 *         T_c = 647.096\mbox{\ K}
 *    \f]
 *    \f[
 *        \rho_c = 322 \mbox{\  kg\ m$^{-3}$}
 *     \f]
 *    \f[
 *        R/M = 0.46151805 \mbox{\ kJ\ kg$^{-1}$\ K$^{-1}$}
 *    \f]
 *
 *  The free energy is a unique single-valued function of the temperature and density
 *  over its entire range.
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
 *  Therefore, to use this class within %Cantera, offsets to u() and s() must be used
 *  to put the water class onto the same basis as other thermodynamic quantities.
 *  For example, in the WaterSSTP class, these offsets are calculated in the following way.
 *  The thermodynamic base state for water is set to the NIST basis here
 *  by specifying constants EW_Offset and SW_Offset. These offsets are
 *  calculated on the fly so that the following properties hold:
 *
 *   - Delta_Hfo_idealGas(298.15, 1bar) = -241.826 kJ/gmol
 *   - So_idealGas(298.15, 1bar)        = 188.835 J/gmolK
 *
 *   The offsets are calculated by actually computing the above quantities and then
 *   calculating the correction factor.
 *
 *  This class provides an interface to the #WaterPropsIAPWSphi class, which actually
 *  calculates the \f$ \phi^o(\delta, \tau)  \f$ and the  \f$ \phi^r(\delta, \tau) \f$
 *  polynomials in dimensionless form.
 *
 *  All thermodynamic results from this class are returned in dimensional form. This 
 *  is because the gas constant (and molecular weight) used within this class is allowed to be potentially
 *  different than that used elsewhere in %Cantera. Therefore, everything has to be
 *  in dimensional units. Note, however, the thermodynamic basis is set to that used
 *  in the steam tables. (u = s = 0 for liquid water at the triple point).
 *
 *  This class is not a %ThermoPhase. However, it does maintain an internal state of
 *  the object that is dependent on temperature and density. The internal state
 *  is characterized by an internally storred \f$ \tau\f$ and a \f$ \delta \f$ value,
 *  and an iState value, which indicates whether the point is a liquid, a gas,
 *  or a supercritical fluid.
 *  Along with that the  \f$ \tau\f$ and a \f$ \delta \f$ values are polynomials of
 *  \f$ \tau\f$ and a \f$ \delta \f$ that are kept by the #WaterPropsIAPWSphi class.
 *  Therefore, whenever  \f$ \tau\f$ or \f$ \delta \f$ is changed, the function setState()
 *  must be called in order for the internal state to be kept up to date. 
 *
 * The class is pretty straightfoward. However, one function deserves mention.
 * the #density() function calculates the density that is consistent with
 * a particular value of the temperature and pressure. It may therefore be 
 * multivalued or potentially there may be no answer from this function. It therefore
 * takes a phase guess and a density guess as optional parameters. If no guesses are
 * supplied to density(), a gas phase guess is assumed. This may or may not be what
 * is wanted. Therefore, density() should usually at leat be supplied with a phase
 * guess so that it may manufacture an appropriate density guess. 
 * #density() manufactures the initial density guess, nondimensionalizes everything,
 * and then calls #WaterPropsIAPWSphi::dfind(), which does the iterative calculation
 * to find the density condition that matches the desired input pressure.
 *
 * The phase guess defines are located in the .h file. they are
 *
 *   - WATER_GAS
 *   - WATER_LIQUID
 *   - WATER_SUPERCRIT
 *
 * @ingroup thermoprops
 *
 */
class WaterPropsIAPWS {
public:

  //! Base constructor
  WaterPropsIAPWS();

  //! Copy constructor
  /*!
   * @param right Object to be copied
   */
  WaterPropsIAPWS(const WaterPropsIAPWS &right);

  //! assignment constructor
  /*!
   * @param right Object to be copied
   */
  WaterPropsIAPWS & operator=(const WaterPropsIAPWS &right);

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
  //! and a guess at the density. Sets the internal state. 
  /*!
   *  Note, below T_c, this is a multivalued function. 
   *
   * The #density() function calculates the density that is consistent with
   * a particular value of the temperature and pressure. It may therefore be 
   * multivalued or potentially there may be no answer from this function. It therefore
   * takes a phase guess and a density guess as optional parameters. If no guesses are
   * supplied to density(), a gas phase guess is assumed. This may or may not be what
   * is wanted. Therefore, density() should usually at leat be supplied with a phase
   * guess so that it may manufacture an appropriate density guess. 
   * #density() manufactures the initial density guess, nondimensionalizes everything,
   * and then calls #WaterPropsIAPWSphi::dfind(), which does the iterative calculation
   * to find the density condition that matches the desired input pressure.
   *
   *  @param  temperature: Kelvin
   *  @param  pressure   : Pressure in Pascals (Newton/m**2)
   *  @param  phase      : guessed phase of water 
   *                     : -1: no guessed phase
   *  @param rhoguess    : guessed density of the water
   *                     : -1.0 no guessed density
   *  @return
   *     Returns the density. If an error is encountered in the calculation
   *     the value of -1.0 is returned.
   */
  double density(double temperature, double pressure, 
		 int phase = -1, double rhoguess = -1.0);

  //! Returns the density (kg m-3)
  /*!
   * The density is an independent variable in the underlying equation of state
   */
  double density() const;

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

private:

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
