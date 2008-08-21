/**
 *  @file WaterProps.h
 *   Header for a class used to house several approximation
 *   routines for properties of water.
 *  (see \ref thermoprops
 *   and class \link Cantera::WaterProps WaterProps\endlink).
 */
/*
 * Copywrite (2006) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
/*
 * $Id$
 */

#ifndef CT_WATERPROPS_H
#define CT_WATERPROPS_H


#include "ct_defs.h"
class WaterPropsIAPWS;
namespace Cantera {
  class WaterPDSS;
 
  /**
   * @defgroup relatedProps Electric Properties of Phases
   *
   *
   * These classes are used to compute the electrical and electrothermochemical properties of
   * phases of matter. The main property currently is the dielectric
   * constant, which is an important parameter for electolyte solutions.
   *
   * 
   * @ingroup phases
   */
  //@{


  //! The WaterProps class is used to 
  //! house several approximation routines for properties of water.
  /*!
   *  The class is also a wrapper around the WaterPropsIAPWS class
   *  which provides the calculations for the equation of
   *  state properties for water.
   *
   *  In particular, this class house routine for the calculation
   *  of the dielectric constant of water
   *
   * Most if not all of the member functions are static.
   */
  class WaterProps {

  public:

    //! Default constructor
    WaterProps();

    //! Constructor with pointer to Water PDSS object
    /*!
     * @param wptr Pointer to water standard state object
     */
    WaterProps(WaterPDSS *wptr);

    //! Copy Constructor
    /*!
     * @param b Object to be copied
     */
    WaterProps(const WaterProps &b);

    //! destructor
    virtual ~WaterProps();

    //! Assignment operator
    /*!
     * @param b Object to be copied
     */
    WaterProps& operator=(const WaterProps& b);

    
    //! Simple calculation of water density at atmospheric pressure.
    //! Valid up to boiling point.
    /*!
     * This formulation has no dependence on the pressure and shouldn't
     * be used where accuracy is needed.
     *
     * @param T temperature in kelvin
     * @param P Pressure in pascal
     * @param ifunc changes what's returned
     *
     * @return value returned depends on ifunc value:
     * ifunc = 0 Returns the density in kg/m^3
     * ifunc = 1 returns the derivative of the density wrt T.
     * ifunc = 2 returns the 2nd derivative of the density wrt T 
     * ifunc = 3 returns the derivative of the density wrt P.
     *
     * Verification:
     *   Agrees with the CRC values (6-10) for up to 4 sig digits.
     *
     * units = returns density in kg m-3.
     */
    static double density_T(double T, double P, int ifunc);

   
    //!     Bradley-Pitzer equation for the dielectric constant 
    //!     of water as a function of temperature and pressure.
    /*!
     *  Returns the dimensionless relative dielectric constant
     *  and its derivatives.
     * 
     *
     * Range of validity: 0 to 350C, 0 to 1 kbar pressure
     *
     * @param T temperature (kelvin)
     * @param P_pascal pressure in pascal
     * @param ifunc changes what's returned from the function
     *   - ifunc = 0 return value
     *   - ifunc = 1 return temperature derivative
     *   - ifunc = 2 return temperature second derivative
     *   - ifunc = 3 return pressure first derivative
     *   .
     *
     * @return Depends on the value of ifunc:
     *   - ifunc = 0 return value
     *   - ifunc = 1 return temperature derivative
     *   - ifunc = 2 return temperature second derivative
     *   - ifunc = 3 return pressure first derivative
     *   .
     *
     *  Validation:
     *   Numerical experiments indicate that this function agrees with
     *   the Archer and Wang data in the CRC p. 6-10 to all 4 significant
     *   digits shown (0 to 100C).
     * 
     *   value at 25C and 1 atm, relEps = 78.38
     * 
     */
    double relEpsilon(double T, double P_pascal,  int ifunc = 0);

    
    //! ADebye calculates the value of A_Debye as a function
    //! of temperature and pressure according to relations
    //! that take into account the temperature and pressure
    //! dependence of the water density and dieletric constant.
    /*!
     *  The A_Debye expression appears on the top of the
     *  ln actCoeff term in the general Debye-Huckel expression
     *  It depends on temperature and pressure. And, therefore,
     *  most be recalculated whenever T or P changes.
     *  The units returned by this expression are sqrt(kg/gmol).
     *     
     *
     *    \f[
     *      A_{Debye} = \frac{1}{8 \pi} \sqrt{\frac{2 N_{Avog} \rho_w}{1000}}
     *                        {\left(\frac{e^2}{\epsilon k_{boltz} T}\right)}^{\frac{3}{2}}
     *    \f]
     *
     *
     *   Nominal value at 25C and 1atm = 1.172576 sqrt(kg/gmol).
     *
     *                Based on:
     *                    epsilon/epsilon_0 = 78.54
     *                           (water at 25C)
     *                    epsilon_0 = 8.854187817E-12 C2 N-1 m-2
     *                    e = 1.60217653E-19 C
     *                    F = 9.6485309E7 C kmol-1
     *                    R = 8.314472E3 kg m2 s-2 kmol-1 K-1
     *                    T = 298.15 K
     *                    B_Debye = 3.28640E9 sqrt(kg/gmol)/m
     *                    Na = 6.0221415E26
     *
     *  @param T  Temperature (kelvin)
     *  @param P  pressure (pascal)
     *  @param ifunc Changes what's returned from the routine:
     *   - ifunc = 0 return value
     *   - ifunc = 1 return temperature derivative
     *   - ifunc = 2 return temperature second derivative
     *   - ifunc = 3 return pressure first derivative
     *   .
     *
     * @return Returns a single double whose meaning depends on ifunc:
     *   - ifunc = 0 return value
     *   - ifunc = 1 return temperature derivative
     *   - ifunc = 2 return temperature second derivative
     *   - ifunc = 3 return pressure first derivative
     *   .
     *
     *  Verification:
     *
     *    With the epsRelWater value from the Bradley-Pitzer relation,
     *    and the water density from the density_IAPWS() function,
     *    The A_Debye computed with this function agrees with
     *    the Pitzer table p. 99 to 4 significant digits at 25C.
     *    and 20C. (Aphi = ADebye/3)
     */
    double ADebye(double T, double P, int ifunc);


    //! Returns the saturation pressure given the temperature
    /*!
     * @param T temperature (kelvin)
     * @return returns the saturation pressure (pascal)
     */
    double satPressure(double T);
 

    //! Returns the density of water
    /*!
     * @param T Temperature (kelvin)
     * @param P pressure (pascal)
     */
    double density_IAPWS(double T, double P);

    //! returns the coefficient of thermal expansion
    /*!
     *  @param T Temperature (kelvin)
     *  @param P pressure (pascal)
     */
    double coeffThermalExp_IAPWS(double T, double P);

    //! Returns the isothermal compressibility of water
    /*!
     * @param T  temperature in kelvin
     * @param P  pressure in pascal
     */
    double isothermalCompressibility_IAPWS(double T, double P);

  protected:

    //! Pointer to the WaterPropsIAPWS object
    /*!
     *  this pointer points to the water object.
     */
    WaterPropsIAPWS *m_waterIAPWS;

    //! true if we own the WaterPropsIAPWS object
    bool m_own_sub;
  };

 //@}
}


#endif
