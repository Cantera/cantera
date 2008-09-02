/**
 * @file WaterPropsIAPWS.cpp
 * Definitions for a class for calculating the equation of state of water
 * from the IAPWS 1995 Formulation based on the steam tables thermodynamic
 * basis (See class \link WaterPropsIAPWS WaterPropsIAPWS\endlink).
 */
/*
 * Copywrite (2006) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
/*
 * $Id$
 */

#include "WaterPropsIAPWS.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
/*
 * Critical Point values of water in mks units
 */
//! Critical Temperature value (kelvin)
const double T_c = 647.096;  
//! Critical Pressure (Pascals)
static const double P_c = 22.064E6; 
//! Value of the Density at the critical point (kg  m-3)
const double Rho_c = 322.;    
//! Molecular Weight of water that is consistent with the paper (kg kmol-1)
static const double M_water = 18.015268; 

/*
 * Note, this is the Rgas value quoted in the paper. For consistency
 * we have to use that value and not the updated value
 *
 * The Ratio of R/M = 0.46151805 kJ kg-1 K-1 , which is Eqn. (6.3) in the paper.
 */
//static const double Rgas = 8.314472E3;   // Joules kmol-1 K-1
static const double Rgas = 8.314371E3;   // Joules kmol-1 K-1


WaterPropsIAPWS:: WaterPropsIAPWS() :
  m_phi(0),
  tau(-1.0),
  delta(-1.0),
  iState(-30000)
{
  m_phi = new WaterPropsIAPWSphi();
}

WaterPropsIAPWS::WaterPropsIAPWS(const WaterPropsIAPWS &b) :
  m_phi(0),
  tau(b.tau),
  delta(b.delta),
  iState(b.iState)
{
  m_phi = new WaterPropsIAPWSphi();
  m_phi->tdpolycalc(tau, delta);
}

WaterPropsIAPWS & WaterPropsIAPWS::operator=(const WaterPropsIAPWS &b) {
  if (this == &b) return *this;
  tau = b.tau;
  delta = b.delta;
  iState = b.iState;
  m_phi->tdpolycalc(tau, delta);
  return *this;
}

WaterPropsIAPWS::~WaterPropsIAPWS() {
  delete (m_phi);
  m_phi = 0;
}


void WaterPropsIAPWS::calcDim(double temperature, double rho) {
  tau = T_c / temperature;
  delta = rho / Rho_c;
  /*
   * Determine the internal state
   */
  if (temperature > T_c) {
      iState = WATER_SUPERCRIT;
  } else {
    if (delta < 1.0) {
      iState = WATER_GAS;
    } else {
      iState = WATER_LIQUID;
    }
  }
}

/*
 * Calculate the Helmholtz free energy in mks units of
 * J kmol-1 K-1.
 */
double  WaterPropsIAPWS::helmholtzFE(double temperature, double rho) {
  setState(temperature, rho);
  double retn = m_phi->phi(tau, delta);
  double RT = Rgas * temperature;
  return (retn * RT);
}
double  WaterPropsIAPWS::helmholtzFE() const{
  double retn = m_phi->phi(tau, delta);
  double temperature = T_c/tau;
  double RT = Rgas * temperature;
  return (retn * RT);
}

/*
 * Calculate the pressure (Pascals), given the temperature and density
 *  Temperature: kelvin
 *  rho: density in kg m-3 
 */
double  WaterPropsIAPWS::pressure(double temperature, double rho) {
  calcDim(temperature, rho);
  double retn = m_phi->pressureM_rhoRT(tau, delta);
  return (retn * rho * Rgas * temperature/M_water);
}

/*
 * Calculate the pressure (Pascals), using the 
 * current internally storred temperature and density
 *  Temperature: kelvin
 *  rho: density in kg m-3 
 */
double  WaterPropsIAPWS::pressure() const {
  double retn = m_phi->pressureM_rhoRT(tau, delta);
  double rho = delta * Rho_c;
  double temperature = T_c / tau;
  return (retn * rho * Rgas * temperature/M_water);
}

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
 *
 * If a problem is encountered, a negative 1 is returned.
 */
double WaterPropsIAPWS::
density(double temperature, double pressure, int phase, double rhoguess) {

  double deltaGuess = 0.0;
  if (rhoguess == -1.0) {
    if (phase != -1) {
      if (temperature > T_c) {
	rhoguess = pressure * M_water / (Rgas * temperature);
      } else {
	if (phase != WATER_LIQUID) {
	  rhoguess = pressure * M_water / (Rgas * temperature);
	} else {
	  /*
	   * Provide a guess about the liquid density that is 
	   * relatively high -> convergnce from above seems robust.
	   */
	  rhoguess = 1000.;
	}
      }
    } else {
      /*
       * Assume the Gas phase initial guess, if nothing is
       * specified to the routine
       */
      rhoguess = pressure * M_water / (Rgas * temperature);
    }

  }
  double p_red = pressure * M_water / (Rgas * temperature * Rho_c);
  deltaGuess = rhoguess / Rho_c;
  setState(temperature, rhoguess);
  double delta_retn = m_phi->dfind(p_red, tau, deltaGuess);
  double density_retn;
  if (delta_retn >0.0) {
    delta = delta_retn;
      
    /*
     * Dimensionalize the density before returning
     */
    density_retn = delta_retn * Rho_c;
    /*
     * Set the internal state -> this may be
     * a duplication. However, let's just be sure.
     */
    setState(temperature, density_retn);
    

  } else {
    density_retn = -1.0;
  }
  return density_retn;
}

double WaterPropsIAPWS::density() const {
  return (delta * Rho_c);
}

/*
 * psat_est provides a rough estimate of the saturation
 * pressure given the temperature. This is used as an initial
 * guess for refining the pressure.
 *
 * Input
 *   temperature (kelvin)
 *
 * return:
 *   psat (Pascals)
 */
double WaterPropsIAPWS::psat_est(double temperature) {
    
  static const double A[8] = {
    -7.8889166E0,
    2.5514255E0,
    -6.716169E0,
    33.2239495E0,
    -105.38479E0,
    174.35319E0,
    -148.39348E0,
    48.631602E0
  };
  double ps;
  if (temperature < 314.) {
    double pl = 6.3573118E0 - 8858.843E0 / temperature
      + 607.56335E0 * pow(temperature, -0.6);
    ps = 0.1 * exp(pl);
  } else {
    double v = temperature / 647.25;
    double w = fabs(1.0-v);
    double b = 0.0;
    for (int i = 0; i < 8; i++) {
      double z = i + 1;
      b += A[i] * pow(w, ((z+1.0)/2.0));
    }
    double q = b / v;
    ps = 22.093*exp(q);
  }
  /*
   * Original correlation was in cgs. Convert to mks
   */
  ps *= 1.0E6;
  return ps;
}

/*
 * Returns the coefficient of thermal expansion as a function
 * of temperature and pressure.
 *           alpha = d (ln V) / dT at constant P.
 *
 * Currently this function is calculated using a differencing scheme. 
 */
double WaterPropsIAPWS::coeffThermExp(double temperature, double pressure) {

  double deltaT = 0.01;
  double psat_at=0.0;
  double rhoguess = -1;
  int phase = -1;
  if (temperature > T_c) {
    rhoguess = pressure * M_water / (Rgas * temperature);
  } else {
    psat_at = psat(temperature);
    if (pressure >= psat_at) {
      phase = WATER_LIQUID;
      deltaT = -0.01;
    } else
      phase = WATER_GAS;
  }

  double dens_base = density(temperature, pressure, phase, rhoguess);

  if (dens_base == -1.0) {
    printf("problems\n");
    exit(-1);
  }
  double temp_del = temperature + deltaT;

  double dens_del = density(temp_del, pressure, phase, dens_base);

  double Vavg = 0.5 * (1./dens_del + 1./dens_base); 
  double retn = 1.0 / Vavg * (1./dens_del - 1.0/dens_base)/deltaT;
  return retn;
}

/*
 * Returns the coefficient of isothermal compressibility
 * of temperature and pressure.
 *          kappa = - d (ln V) / dP at constant T.
 *
 */
double WaterPropsIAPWS::isothermalCompressibility() const {
  double retn = m_phi->dimdpdrho(tau, delta);
  double temperature = T_c/tau;
  double dpdrho = retn * Rgas * temperature / M_water;
  double dens = delta * Rho_c;
  return (1.0 / (dens * dpdrho));
}

/*
 * Calculate the Gibbs free energy in mks units of
 * J kmol-1 K-1.
 */
double WaterPropsIAPWS::Gibbs(double temperature, double rho) {
  setState(temperature, rho);
  double gRT = m_phi->gibbs_RT();
  return (gRT * Rgas * temperature);
}

double WaterPropsIAPWS::Gibbs() const {
  double gRT = m_phi->gibbs_RT();
  double temperature = T_c/tau;
  return (gRT * Rgas * temperature);
}

/*
 * Calculate the Gibbs free energy in mks units of
 * J kmol-1 K-1.
 */
void  WaterPropsIAPWS::
corr(double temperature, double pressure, double &densLiq, 
     double &densGas, double &delGRT) {
    
  densLiq = density(temperature, pressure, WATER_LIQUID, densLiq);
  if (densLiq <= 0.0) {
    printf("error liq\n");
    exit(-1);
  }
  setState(temperature, densLiq);
  double gibbsLiqRT =  m_phi->gibbs_RT();

  densGas = density(temperature, pressure, WATER_GAS, densGas);
  if (densGas <= 0.0) {
    printf("error gas\n");
    exit(-1);
  }
  setState(temperature, densGas);
  double gibbsGasRT = m_phi->gibbs_RT();
    
  delGRT = gibbsLiqRT - gibbsGasRT;
}

void WaterPropsIAPWS::
corr1(double temperature, double pressure, double &densLiq, 
      double &densGas, double &pcorr) {
    
  densLiq = density(temperature, pressure, WATER_LIQUID, densLiq);
  setState(temperature, densLiq);
  double prL = m_phi->phiR();

  densGas = density(temperature, pressure, WATER_GAS, densGas);
  setState(temperature, densGas);
  double prG = m_phi->phiR();
    
  double rhs = (prL - prG) + log(densLiq/densGas);
  rhs /= (1.0/densGas - 1.0/densLiq);

  pcorr = rhs * Rgas * temperature / M_water;
}

/**
 * Calculate the saturation pressure given the temperature.
 * p : Pascals : Newtons/m**2
 */
static int method = 1;
double WaterPropsIAPWS::psat(double temperature) {
  double densLiq = -1.0, densGas = -1.0, delGRT = 0.0;
  double dp, pcorr;
  double p = psat_est(temperature);
  bool conv = false;
  for (int i = 0; i < 30; i++) {
    if (method == 1) {
      corr(temperature, p, densLiq, densGas, delGRT);
      double delV = M_water * (1.0/densLiq - 1.0/densGas);
      dp = - delGRT * Rgas * temperature / delV;
    } else {
      corr1(temperature, p, densLiq, densGas, pcorr);
      dp = pcorr - p;
    }
    p += dp;
      
    if ((method == 1) && delGRT < 1.0E-8) {
      conv = true;
      break;
    } else {
      if (fabs(dp/p) < 1.0E-9) {
	conv = true;
	break;
      }
    }
  }
  return p;
}

int WaterPropsIAPWS::phaseState() const {
  return iState;
}

/**
 * Sets the internal state of the object to the
 * specified temperature and density.
 */
void WaterPropsIAPWS::setState(double temperature, double rho) {
  calcDim(temperature, rho);
  m_phi->tdpolycalc(tau, delta);
}


/*
 * Calculate the enthalpy in mks units of
 * J kmol-1 K-1.
 */
double WaterPropsIAPWS::
enthalpy(double temperature, double rho) {
  setState(temperature, rho);
  double hRT =  m_phi->enthalpy_RT();
  return (hRT * Rgas * temperature);
}
double WaterPropsIAPWS::
enthalpy() const {
  double temperature = T_c/tau;
  double hRT =  m_phi->enthalpy_RT();
  return (hRT * Rgas * temperature);
}


/*
 * Calculate the internal Energy in mks units of
 * J kmol-1 K-1.
 */
double WaterPropsIAPWS::
intEnergy(double temperature, double rho) {
  setState(temperature, rho);
  double uRT = m_phi->intEnergy_RT();
  return (uRT * Rgas * temperature);
}
double WaterPropsIAPWS::
intEnergy() const{
  double temperature = T_c / tau;
  double uRT = m_phi->intEnergy_RT();
  return (uRT * Rgas * temperature);
}

/*
 * Calculate the enthalpy in mks units of
 * J kmol-1 K-1.
 */
double WaterPropsIAPWS::
entropy(double temperature, double rho) {
  setState(temperature, rho);
  double sR = m_phi->entropy_R();
  return (sR * Rgas);
}

/*
 * Calculate the enthalpy in mks units of
 * J kmol-1 K-1.
 */
double WaterPropsIAPWS::entropy() const {
  double sR = m_phi->entropy_R();
  return (sR * Rgas);
}

/*
 * Calculate heat capacity at constant volume
 * J kmol-1 K-1.
 */
double WaterPropsIAPWS::cv(double temperature, double rho) {
  setState(temperature, rho);
  double cvR = m_phi->cv_R();
  return (cvR * Rgas);
}

/*
 * Calculate heat capacity at constant pressure
 * J kmol-1 K-1.
 */
double WaterPropsIAPWS::
cp(double temperature, double rho) {
  setState(temperature, rho);
  double cpR = m_phi->cp_R();
  return (cpR * Rgas);
}

double  WaterPropsIAPWS::
cp() const {
  double cpR = m_phi->cp_R();
  return (cpR * Rgas);
}

double WaterPropsIAPWS::
molarVolume(double temperature, double rho) {
  setState(temperature, rho);
  return (M_water / rho);
}

double WaterPropsIAPWS::
molarVolume() const {
  double rho = delta * Rho_c;
  return (M_water / rho);
}
