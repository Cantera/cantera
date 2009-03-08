/**
 *  @file WaterProps.h
 *
 */
/*
 * Copywrite (2006) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
/*
 * $Id: WaterProps.h,v 1.1 2006/07/04 00:01:53 hkmoffa Exp $
 */

#ifndef CT_WATERPROPS_H
#define CT_WATERPROPS_H


#include "ct_defs.h"
class WaterPropsIAPWS;
namespace Cantera {

  class WaterPDSS;
 

  /**
   * Definition of the WaterProps class. This class is used to 
   * house several approximation routines for properties of water
   * Most if not all of the member functions are static.
   */
  class WaterProps {

  public:

    WaterProps();

    WaterProps(WaterPDSS *wptr);

    WaterProps(const WaterProps &b);

    virtual ~WaterProps();

    WaterProps& operator=(const WaterProps&b);

    /*
     * Simple calculation of water density at atmospheric pressure.
     * Valid up to boiling point.
     *
     * ifunc = 0 Returns the density in kg/m^3
     * ifunc = 1 returns the derivative of the density wrt T.
     * ifunc = 2 returns the derivative of the density wrt P.
     * ifunc = 3 returns the 2nd derivative of the density wrt T 
     *
     * Note -> needs augmenting with a T,P implementation.
     *
     * Verification:
     *   Agrees with the CRC values (6-10) for up to 4 sig digits.
     *
     * units = returns density in kg m-3.
     */
    static double density_T(double T, double P, int ifunc);

    /**
     * Dielectric constant for water:
     *     Bradley-Pitzer equation for the dielectric constant 
     *     of water as a function of temperature and pressure.
     *
     *  ifunc = 0 value
     *  ifunc = 1 Temperature deriviative
     *  ifunc = 2 second temperature derivative
     *
     *  @param T temperature in Kelvin
     *  @param P Pressure in bar
     *
     * Range of validity 0 to 350C, 0 to 1 kbar pressure
     *
     * ifunc = 0 return value
     * ifunc = 1 return temperature derivative
     *
     *  Validation:
     *   Numerical experiments indicate that this function agrees with
     *   the Archer and Wang data in the CRC p. 6-10 to all 4 significant
     *   digits shown (0 to 100C).
     * 
     *   value at 25C, relEps = 78.38
     */
    static double relEpsilon(double T, double P_pascal,  int ifunc);

    /**
     * ADebye calculates the value of A_Debye as a function
     * of temperature and pressure according to relations
     * that take into account the temperature and pressure
     * dependence of the water density and dieletric constant.
     *
     * A_Debye -> this expression appears on the top of the
     *            ln actCoeff term in the general Debye-Huckel
     *            expression
     *            It depends on temperature. And, therefore,
     *            most be recalculated whenever T or P changes.
     *            
     *            A_Debye = (1/8Pi) sqrt(2Na dw/1000) 
     *                          (e e/(epsilon RT)^3/2
     *
     *            Units = sqrt(kg/gmol)
     *
     *            Nominal value = 1.172576 sqrt(kg/gmol)
     *                  based on:
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
     *  Verification:
     *    With the epsRelWater value from the BP relation,
     *    and the water density from the WaterDens function,
     *    The A_Debye computed with this function agrees with
     *    the Pitzer table p. 99 to 4 significant digits at 25C.
     *    and 20C. (Aphi = ADebye/3)
     */
    double ADebye(double T, double P, int ifunc);

    double satPressure(double T);
 

    double density_IAPWS(double T, double P);
    double coeffThermalExp_IAPWS(double T, double P);
    double isothermalCompressibility_IAPWS(double T, double P);

  protected:


    WaterPropsIAPWS *m_waterIAPWS;
    bool m_own_sub;
  };


}


#endif
