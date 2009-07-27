/**
 *  @file WaterProps.cpp
 */
/*
 * Copywrite (2006) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
/*
 * $Id: WaterProps.cpp,v 1.9 2008/12/17 17:04:47 hkmoffa Exp $
 */

//@{
#ifndef MAX
#define MAX(x,y)    (( (x) > (y) ) ? (x) : (y))
#endif
//@}

#include "WaterProps.h"
#include "ctml.h"
#include "PDSS_Water.h"
#include "WaterPropsIAPWS.h"

#include <cmath>

namespace Cantera {


  /*
   * default constructor -> object owns its own water evaluator
   */
  WaterProps::WaterProps():
    m_waterIAPWS(0),
    m_own_sub(false)
  {
    m_waterIAPWS = new WaterPropsIAPWS();
    m_own_sub = true;
  }

  /*
   * constructor -> object in slave mode, It doesn't own its
   *      own water evaluator.
   */
  WaterProps::WaterProps(PDSS_Water *wptr)  :
    m_waterIAPWS(0),
    m_own_sub(false)
  {
    if (wptr) {
      m_waterIAPWS = wptr->getWater();
      m_own_sub = false;
    } else {
      m_waterIAPWS = new WaterPropsIAPWS();
      m_own_sub = true;
    }
  }

  /**
   * Copy constructor
   */
  WaterProps::WaterProps(const WaterProps &b)  :
    m_waterIAPWS(0),
    m_own_sub(false)
  {
    *this = b;
  }

  /**
   * Destructor
   */
  WaterProps::~WaterProps() {
    if (m_own_sub) {
      delete m_waterIAPWS;
    }
  }

  /**
   * Assignment operator
   */
  WaterProps& WaterProps::operator=(const WaterProps&b) {
    if (&b == this) return *this;
   
    if (m_own_sub) {
      if (m_waterIAPWS) {
	delete m_waterIAPWS;
	m_waterIAPWS = 0;
      }
    }
    if (b.m_own_sub) {
      m_waterIAPWS = new WaterPropsIAPWS();
      m_own_sub = true;
    } else {
      m_waterIAPWS = b.m_waterIAPWS;
      m_own_sub = false;
    }

    return *this;
  }
  
  // Simple calculation of water density at atmospheric pressure.
  // Valid up to boiling point.
  /*
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
  double WaterProps::density_T(double T, double P, int ifunc) {
    double Tc = T - 273.15;
    const double U1 = 288.9414;
    const double U2 = 508929.2;
    const double U3 = 68.12963;
    const double U4 = -3.9863;
    
    double tmp1 = Tc + U1;
    double tmp4 = Tc + U4;
    double t4t4 = tmp4 * tmp4;
    double tmp3 = Tc + U3;
    double rho = 1000. * (1.0 - tmp1*t4t4/(U2 * tmp3));

    /*
     * Impose an ideal gas lower bound on rho. We need this
     * to ensure positivity of rho, even though it is
     * grossly unrepresentative.
     */
    double rhomin = P / (GasConstant * T);
    if (rho < rhomin) {
      rho = rhomin;
      if (ifunc == 1) {
	double drhodT = - rhomin / T;
	return drhodT;
      } else if (ifunc == 3) {
	double drhodP = rhomin / P;
	return drhodP;
      } else if (ifunc == 2) {
	double d2rhodT2 = 2.0 * rhomin / (T * T);
	return d2rhodT2;
      }
    }
    
    if (ifunc == 1) {
      double drhodT = 1000./U2 * (
				  - tmp4 * tmp4 / (tmp3)
				  - tmp1 * 2 * tmp4 / (tmp3)
				  + tmp1 * t4t4 / (tmp3*tmp3)
				  );
      return drhodT;
    } else if (ifunc == 3) {
      return 0.0;
    } else if (ifunc == 2) {
      double t3t3 = tmp3 * tmp3;
      double d2rhodT2 =  1000./U2 * 
	((-4.0*tmp4-2.0*tmp1)/tmp3 +
	 (2.0*t4t4 + 4.0*tmp1*tmp4)/t3t3 
	 - 2.0*tmp1 * t4t4/(t3t3*tmp3));
      return d2rhodT2;
    }
        
    return rho;
  }

  //  Bradley-Pitzer equation for the dielectric constant 
  //  of water as a function of temperature and pressure.
  /*!
   *  Returns the dimensionless relative dielectric constant
   *  and its derivatives.
   * 
   *  ifunc = 0 value
   *  ifunc = 1 Temperature deriviative
   *  ifunc = 2 second temperature derivative
   *  ifunc = 3 return pressure first derivative
   *
   * Range of validity 0 to 350C, 0 to 1 kbar pressure
   *
   * @param T temperature (kelvin)
   * @param P_pascal pressure in pascal
   * @param ifunc changes what's returned from the function
   *
   * @return Depends on the value of ifunc:
   * ifunc = 0 return value
   * ifunc = 1 return temperature derivative
   * ifunc = 2 return second temperature derivative
   * ifunc = 3 return pressure first derivative
   *
   *  Validation:
   *   Numerical experiments indicate that this function agrees with
   *   the Archer and Wang data in the CRC p. 6-10 to all 4 significant
   *   digits shown (0 to 100C).
   * 
   *   value at 25C, relEps = 78.38
   * 
   */
  double WaterProps::relEpsilon(double T, double P_pascal, 
				int ifunc) {
    const double U1 = 3.4279E2;
    const double U2 = -5.0866E-3;
    const double U3 = 9.4690E-7;
    const double U4 = -2.0525;
    const double U5 = 3.1159E3;
    const double U6 = -1.8289E2;
    const double U7 = -8.0325E3;
    const double U8 = 4.2142E6;
    const double U9 = 2.1417;
    double T2 = T * T;

    double eps1000 = U1 * exp(U2 * T + U3 * T2);
    double C = U4 + U5/(U6 + T);
    double B = U7 + U8/T + U9 * T;

    double Pbar = P_pascal * 1.0E-5;
    double tmpBpar = B + Pbar;
    double tmpB1000 = B + 1000.0;
    double ltmp =  log(tmpBpar/tmpB1000);
    double epsRel = eps1000 + C * ltmp;

    if (ifunc == 1 || ifunc == 2) {
      double tmpC = U6 + T;
      double dCdT = - U5/(tmpC * tmpC);

      double dBdT = - U8/(T * T) + U9;

      double deps1000dT = eps1000 * (U2 + 2.0 * U3 * T);

      double dltmpdT = (dBdT/tmpBpar - dBdT/tmpB1000);
      if (ifunc == 1) {
	double depsReldT = deps1000dT + dCdT * ltmp + C * dltmpdT;
	return depsReldT;
      }
      double T3     = T2 * T;
      double d2CdT2 = - 2.0 * dCdT / tmpC;
      double d2BdT2 =   2.0 * U8 / (T3);

      double d2ltmpdT2 = (d2BdT2*(1.0/tmpBpar - 1.0/tmpB1000) +
			  dBdT*dBdT*(1.0/(tmpB1000*tmpB1000) - 1.0/(tmpBpar*tmpBpar)));

      double d2eps1000dT2 =  (deps1000dT * (U2 + 2.0 * U3 * T) + eps1000  * (2.0 * U3));

      if (ifunc == 2) {
	double d2epsReldT2 = (d2eps1000dT2 + d2CdT2 * ltmp + 2.0 * dCdT * dltmpdT
			      + C * d2ltmpdT2);
	return d2epsReldT2;
      }
    }
    if (ifunc == 3) {
      double dltmpdP   = 1.0E-5 / tmpBpar; 
      double depsReldP = C * dltmpdP;
      return depsReldP;
    }

    return epsRel;
  }

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
   *            A_Debye = (1/(8 Pi)) sqrt(2 Na dw / 1000) 
   *                          (e e/(epsilon R T))^3/2
   *
   *            Units = sqrt(kg/gmol) ~ sqrt(1/I)
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
   * ifunc = 0 return value
   * ifunc = 1 return temperature derivative
   * ifunc = 2 return temperature second derivative
   * ifunc = 3 return pressure first derivative
   *
   *  Verification:
   *    With the epsRelWater value from the BP relation,
   *    and the water density from the WaterDens function,
   *    The A_Debye computed with this function agrees with
   *    the Pitzer table p. 99 to 4 significant digits at 25C.
   *    and 20C. (Aphi = ADebye/3)
   * 
   * (statically defined within the object)
   */
  double WaterProps::ADebye(double T, double P_input, int ifunc) {
    const double e =  1.60217653E-19;
    const double epsilon0 =  8.854187817E-12;
    const double R = 8.314472E3;
    double psat = satPressure(T);
    double P;
    if (psat > P_input) {
      //printf("ADebye WARNING: p_input < psat: %g %g\n",
      // P_input, psat);
      P = psat;
    } else {
      P = P_input;
    }
    double epsRelWater = relEpsilon(T, P, 0);
    //printf("releps calc = %g, compare to 78.38\n", epsRelWater);
    //double B_Debye = 3.28640E9;
    const double Na = 6.0221415E26;

    double epsilon = epsilon0 * epsRelWater;
    double dw = density_IAPWS(T, P);
    double tmp = sqrt( 2.0 * Na * dw / 1000.);
    double tmp2 = e * e * Na / (epsilon * R * T);
    double tmp3 = tmp2 * sqrt(tmp2);
    double A_Debye = tmp * tmp3 / (8.0 * Pi);


    /*
     *  dAdT = - 3/2 Ad/T + 1/2 Ad/dw d(dw)/dT - 3/2 Ad/eps d(eps)/dT
     *
     *  dAdT = - 3/2 Ad/T - 1/2 Ad/Vw d(Vw)/dT - 3/2 Ad/eps d(eps)/dT
     */
    if (ifunc == 1 || ifunc == 2) {
      double dAdT = - 1.5 * A_Debye / T;

      double depsRelWaterdT = relEpsilon(T, P, 1);
      dAdT -= A_Debye * (1.5 * depsRelWaterdT / epsRelWater);

      //int methodD = 1;
      //double ddwdT = density_T_new(T, P, 1);
      // double contrib1 = A_Debye * (0.5 * ddwdT / dw);
         
      /*
       * calculate d(lnV)/dT _constantP, i.e., the cte
       */
      double cte = coeffThermalExp_IAPWS(T, P);
      double contrib2 =  - A_Debye * (0.5 * cte);

      //dAdT += A_Debye * (0.5 * ddwdT / dw);
      dAdT += contrib2;

#ifdef DEBUG_HKM
      //printf("dAdT = %g, contrib1 = %g, contrib2 = %g\n", 
      //	 dAdT, contrib1, contrib2);
#endif
 
      if (ifunc == 1) {
	return dAdT;
      }

      if (ifunc == 2) {
	/*
	 * Get the second derivative of the dielectric constant wrt T
	 * -> we will take each of the terms in dAdT and differentiate
	 *    it again.
	 */
	double d2AdT2 = 1.5 / T * (A_Debye/T - dAdT);

	double d2epsRelWaterdT2 = relEpsilon(T, P, 2);

	//double dT = -0.01;
	//double TT = T + dT;
	//double depsRelWaterdTdel = relEpsilon(TT, P, 1);
	//double d2alt = (depsRelWaterdTdel- depsRelWaterdT ) / dT;
	//printf("diff %g %g\n",d2epsRelWaterdT2, d2alt); 
	// HKM -> checks out, i.e., they are the same.

	d2AdT2 += 1.5 * (- dAdT * depsRelWaterdT / epsRelWater 
			 - A_Debye / epsRelWater * 
			 (d2epsRelWaterdT2 - depsRelWaterdT * depsRelWaterdT / epsRelWater));
	    
	double deltaT = -0.1;
	double Tdel = T + deltaT;
	double cte_del =  coeffThermalExp_IAPWS(Tdel, P);
	double dctedT = (cte_del - cte) / Tdel;
	    
	   
	//double d2dwdT2 = density_T_new(T, P, 2);

	double contrib3 = 0.5 * ( -(dAdT * cte) -(A_Debye * dctedT));
	d2AdT2 += contrib3;

	return d2AdT2;
      }
    }
    /*
     *  A_Debye = (1/(8 Pi)) sqrt(2 Na dw / 1000) 
     *                          (e e/(epsilon R T))^3/2
     *
     *  dAdP =  + 1/2 Ad/dw d(dw)/dP - 3/2 Ad/eps d(eps)/dP
     *
     *  dAdP =  - 1/2 Ad/Vw d(Vw)/dP - 3/2 Ad/eps d(eps)/dP
     *
     *  dAdP =  + 1/2 Ad * kappa  - 3/2 Ad/eps d(eps)/dP
     *
     *  where kappa = - 1/Vw d(Vw)/dP_T (isothermal compressibility)
     */
    if (ifunc == 3) {
	  
      double dAdP = 0.0;
	  
      double depsRelWaterdP = relEpsilon(T, P, 3);
      dAdP -=  A_Debye * (1.5 * depsRelWaterdP / epsRelWater);
	  
      double kappa = isothermalCompressibility_IAPWS(T,P);

      //double ddwdP = density_T_new(T, P, 3);
      dAdP += A_Debye * (0.5 * kappa);

      return dAdP;
    }

    return A_Debye;
  }

  double WaterProps::satPressure(double T) {
    double pres = m_waterIAPWS->psat(T);
    return pres;
  }
 

  double WaterProps::density_IAPWS(double temp, double press) {

    double dens;
    dens = m_waterIAPWS->density(temp, press, WATER_LIQUID);
    return dens;
  }

  double WaterProps::coeffThermalExp_IAPWS(double temp, double press) {
    double dens = m_waterIAPWS->density(temp, press, WATER_LIQUID);
    if (dens < 0.0) {
      throw CanteraError("WaterProps::coeffThermalExp_IAPWS", 
			 "Unable to solve for density at T = " + fp2str(temp) + " and P = " + fp2str(press));
    }
    double cte = m_waterIAPWS->coeffThermExp();
    return cte;
  }

  double WaterProps::isothermalCompressibility_IAPWS(double temp, double press) {
    double dens = m_waterIAPWS->density(temp, press, WATER_LIQUID);
    if (dens < 0.0) {
      throw CanteraError("WaterProps::isothermalCompressibility_IAPWS", 
			 "Unable to solve for density at T = " + fp2str(temp) + " and P = " + fp2str(press));
    }
    double kappa = m_waterIAPWS->isothermalCompressibility();
    return kappa;
  }


}
