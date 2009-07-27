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
 * $Id: WaterPropsIAPWS.cpp,v 1.17 2008/12/17 17:01:29 hkmoffa Exp $
 */

#include "WaterPropsIAPWS.h"
#include "ctexceptions.h"
#include "stringUtils.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
/*
 * Critical Point values of water in mks units
 */
//! Critical Temperature value (kelvin)
const doublereal T_c = 647.096;  
//! Critical Pressure (Pascals)
static const doublereal P_c = 22.064E6; 
//! Value of the Density at the critical point (kg  m-3)
const doublereal Rho_c = 322.;
//! Molecular Weight of water that is consistent with the paper (kg kmol-1)
static const doublereal M_water = 18.015268; 

/*
 * Note, this is the Rgas value quoted in the paper. For consistency
 * we have to use that value and not the updated value
 *
 * The Ratio of R/M = 0.46151805 kJ kg-1 K-1 , which is Eqn. (6.3) in the paper.
 */
//static const doublereal Rgas = 8.314472E3;   // Joules kmol-1 K-1
static const doublereal Rgas = 8.314371E3;   // Joules kmol-1 K-1
//@{
#ifndef MAX
# define MAX(x,y) (( (x) > (y) ) ? (x) : (y))
#endif

#ifndef MIN
# define MIN(x,y) (( (x) < (y) ) ? (x) : (y))
#endif
//@}

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


void WaterPropsIAPWS::calcDim(doublereal temperature, doublereal rho) {
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

doublereal  WaterPropsIAPWS::helmholtzFE() const {
  doublereal retn = m_phi->phi(tau, delta);
  doublereal temperature = T_c/tau;
  doublereal RT = Rgas * temperature;
  return (retn * RT);
}

/*
 * Calculate the pressure (Pascals), using the 
 * current internally storred temperature and density
 *  Temperature: kelvin
 *  rho: density in kg m-3 
 */
doublereal  WaterPropsIAPWS::pressure() const {
  doublereal retn = m_phi->pressureM_rhoRT(tau, delta);
  doublereal rho = delta * Rho_c;
  doublereal temperature = T_c / tau;
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
doublereal WaterPropsIAPWS::density(doublereal temperature, doublereal pressure,
				    int phase, doublereal rhoguess) {

  doublereal deltaGuess = 0.0;
  if (rhoguess == -1.0) {
    if (phase != -1) {
      if (temperature > T_c) {
	rhoguess = pressure * M_water / (Rgas * temperature);
      } else {
	if (phase == WATER_GAS || phase == WATER_SUPERCRIT) {
	  rhoguess = pressure * M_water / (Rgas * temperature);
	} else if (phase == WATER_LIQUID) {
	  /*
	   * Provide a guess about the liquid density that is 
	   * relatively high -> convergnce from above seems robust.
	   */
	  rhoguess = 1000.;
	} else if (phase == WATER_UNSTABLELIQUID || phase == WATER_UNSTABLEGAS) {
	  throw Cantera::CanteraError("WaterPropsIAPWS::density", 
				      "Unstable Branch finder is untested");
	} else {
	  throw Cantera::CanteraError("WaterPropsIAPWS::density", 
				      "unknown state: " + Cantera::int2str(phase));
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
  doublereal p_red = pressure * M_water / (Rgas * temperature * Rho_c);
  deltaGuess = rhoguess / Rho_c;
  setState_TR(temperature, rhoguess);
  doublereal delta_retn = m_phi->dfind(p_red, tau, deltaGuess);
  doublereal density_retn;
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
    setState_TR(temperature, density_retn);
    

  } else {
    density_retn = -1.0;
  }
  return density_retn;
}


doublereal WaterPropsIAPWS::density_const(doublereal pressure,
					  int phase, doublereal rhoguess) const {
  doublereal temperature = T_c / tau;
  doublereal deltaGuess = 0.0;
  doublereal deltaSave = delta;
  if (rhoguess == -1.0) {
    if (phase != -1) {
      if (temperature > T_c) {
	rhoguess = pressure * M_water / (Rgas * temperature);
      } else {
	if (phase == WATER_GAS || phase == WATER_SUPERCRIT) {
	  rhoguess = pressure * M_water / (Rgas * temperature);
	} else if (phase == WATER_LIQUID) {
	  /*
	   * Provide a guess about the liquid density that is 
	   * relatively high -> convergnce from above seems robust.
	   */
	  rhoguess = 1000.;
	} else if (phase == WATER_UNSTABLELIQUID || phase == WATER_UNSTABLEGAS) {
	  throw Cantera::CanteraError("WaterPropsIAPWS::density", 
				      "Unstable Branch finder is untested");
	} else {
	  throw Cantera::CanteraError("WaterPropsIAPWS::density", 
				      "unknown state: " + Cantera::int2str(phase));
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
  doublereal p_red = pressure * M_water / (Rgas * temperature * Rho_c);
  deltaGuess = rhoguess / Rho_c;

  delta = deltaGuess;
  m_phi->tdpolycalc(tau, delta);
  //  setState_TR(temperature, rhoguess);

  doublereal delta_retn = m_phi->dfind(p_red, tau, deltaGuess);
  doublereal density_retn;
  if (delta_retn > 0.0) {
    delta = delta_retn;
      
    /*
     * Dimensionalize the density before returning
     */
    density_retn = delta_retn * Rho_c;
  
  } else {
    density_retn = -1.0;
  }

  delta = deltaSave;
  m_phi->tdpolycalc(tau, delta);
  return density_retn;
}



doublereal WaterPropsIAPWS::density() const {
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
doublereal WaterPropsIAPWS::psat_est(doublereal temperature) const {
    
  static const doublereal A[8] = {
    -7.8889166E0,
    2.5514255E0,
    -6.716169E0,
    33.2239495E0,
    -105.38479E0,
    174.35319E0,
    -148.39348E0,
    48.631602E0
  };
  doublereal ps;
  if (temperature < 314.) {
    doublereal pl = 6.3573118E0 - 8858.843E0 / temperature
      + 607.56335E0 * pow(temperature, -0.6);
    ps = 0.1 * exp(pl);
  } else {
    doublereal v = temperature / 647.25;
    doublereal w = fabs(1.0-v);
    doublereal b = 0.0;
    for (int i = 0; i < 8; i++) {
      doublereal z = i + 1;
      b += A[i] * pow(w, ((z+1.0)/2.0));
    }
    doublereal q = b / v;
    ps = 22.093*exp(q);
  }
  /*
   * Original correlation was in cgs. Convert to mks
   */
  ps *= 1.0E6;
  return ps;
}

/*
 * Returns the coefficient of isothermal compressibility
 * of temperature and pressure.
 *          kappa = - d (ln V) / dP at constant T.
 */
doublereal WaterPropsIAPWS::isothermalCompressibility() const {
  doublereal dpdrho_val = dpdrho();
  doublereal dens = delta * Rho_c;
  return (1.0 / (dens * dpdrho_val));
}

doublereal WaterPropsIAPWS::dpdrho() const {
  doublereal retn = m_phi->dimdpdrho(tau, delta);
  doublereal temperature = T_c/tau;
  doublereal val = retn * Rgas * temperature / M_water;
  return val;
}

doublereal WaterPropsIAPWS:: coeffPresExp() const {
  doublereal retn = m_phi->dimdpdT(tau, delta);
  return (retn);
}

doublereal WaterPropsIAPWS:: coeffThermExp() const {
  doublereal kappa = isothermalCompressibility();
  doublereal beta = coeffPresExp();
  doublereal dens = delta * Rho_c;
  return (kappa * dens * Rgas * beta / M_water);
}

doublereal WaterPropsIAPWS::Gibbs() const {
  doublereal gRT = m_phi->gibbs_RT();
  doublereal temperature = T_c/tau;
  return (gRT * Rgas * temperature);
}

/*
 * Calculate the Gibbs free energy in mks units of
 * J kmol-1 K-1.
 */
void  WaterPropsIAPWS::
corr(doublereal temperature, doublereal pressure, doublereal &densLiq, 
     doublereal &densGas, doublereal &delGRT) {
    
  densLiq = density(temperature, pressure, WATER_LIQUID, densLiq);
  if (densLiq <= 0.0) {
    throw Cantera::CanteraError("WaterPropsIAPWS::corr",
				"Error occurred trying to find liquid density at (T,P) = "
				+ Cantera::fp2str(temperature) + "  " + Cantera::fp2str(pressure));
  }
  setState_TR(temperature, densLiq);
  doublereal gibbsLiqRT =  m_phi->gibbs_RT();

  densGas = density(temperature, pressure, WATER_GAS, densGas);
  if (densGas <= 0.0) {
    throw Cantera::CanteraError("WaterPropsIAPWS::corr",
				"Error occurred trying to find gas density at (T,P) = "
				+ Cantera::fp2str(temperature) + "  " + Cantera::fp2str(pressure));
  }
  setState_TR(temperature, densGas);
  doublereal gibbsGasRT = m_phi->gibbs_RT();
    
  delGRT = gibbsLiqRT - gibbsGasRT;
}

void WaterPropsIAPWS::
corr1(doublereal temperature, doublereal pressure, doublereal &densLiq, 
      doublereal &densGas, doublereal &pcorr) {
    
  densLiq = density(temperature, pressure, WATER_LIQUID, densLiq);
  if (densLiq <= 0.0) {
    throw Cantera::CanteraError("WaterPropsIAPWS::corr1",
				"Error occurred trying to find liquid density at (T,P) = "
				+ Cantera::fp2str(temperature) + "  " + Cantera::fp2str(pressure));
  }
  setState_TR(temperature, densLiq);
  doublereal prL = m_phi->phiR();

  densGas = density(temperature, pressure, WATER_GAS, densGas);
  if (densGas <= 0.0) {
    throw Cantera::CanteraError("WaterPropsIAPWS::corr1",
				"Error occurred trying to find gas density at (T,P) = "
				+ Cantera::fp2str(temperature) + "  " + Cantera::fp2str(pressure));
  }
  setState_TR(temperature, densGas);
  doublereal prG = m_phi->phiR();
    
  doublereal rhs = (prL - prG) + log(densLiq/densGas);
  rhs /= (1.0/densGas - 1.0/densLiq);

  pcorr = rhs * Rgas * temperature / M_water;
}

/**
 * Calculate the saturation pressure given the temperature.
 * p : Pascals : Newtons/m**2
 */
static int method = 1;

doublereal WaterPropsIAPWS::psat(doublereal temperature, int waterState) {
  doublereal densLiq = -1.0, densGas = -1.0, delGRT = 0.0;
  doublereal dp, pcorr;
  if (temperature >= T_c) {
    densGas = density(temperature, P_c, WATER_SUPERCRIT);
    setState_TR(temperature, densGas);
    return P_c;
  }
  doublereal p = psat_est(temperature);
  bool conv = false;
  for (int i = 0; i < 30; i++) {
    if (method == 1) {
      corr(temperature, p, densLiq, densGas, delGRT);
      doublereal delV = M_water * (1.0/densLiq - 1.0/densGas);
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
  // Put the fluid in the desired end condition
  if (waterState == WATER_LIQUID) {
    setState_TR(temperature, densLiq);
  } else if (waterState == WATER_GAS) {
    setState_TR(temperature, densGas);
  } else {
    throw Cantera::CanteraError("WaterPropsIAPWS::psat",
                                "unknown water state input: " + Cantera::int2str(waterState));
  }
  return p;
}

int WaterPropsIAPWS::phaseState(bool checkState) const {
  if (checkState) {
    if (tau <= 1.0) {
      iState = WATER_SUPERCRIT;
    } else {
      doublereal T = T_c / tau;
      doublereal rho = delta * Rho_c;
      //doublereal psatTable = psat_est(T);
      doublereal rhoMidAtm = 0.5 * (1.01E5 * M_water / (8314.472 * 373.15) + 1.0E3);
      doublereal rhoMid = Rho_c + (T - T_c) * (Rho_c - rhoMidAtm) / (T_c - 373.15);
      int iStateGuess = WATER_LIQUID;
      if (rho < rhoMid) {
	iStateGuess = WATER_GAS;
      }
      doublereal kappa = isothermalCompressibility();
      if (kappa >= 0.0) {
	iState = iStateGuess;
      } else {
	// When we are here we are between the spinodal curves
	doublereal rhoDel = rho * 1.000001;

	//setState_TR(T, rhoDel);
	doublereal deltaSave = delta;
	doublereal deltaDel = rhoDel / Rho_c;
	delta = deltaDel;
	m_phi->tdpolycalc(tau, deltaDel);

	doublereal kappaDel = isothermalCompressibility();
	doublereal d2rhodp2 = (rhoDel * kappaDel - rho * kappa) / (rhoDel - rho);
	if  (d2rhodp2 > 0.0) {
	  iState = WATER_UNSTABLELIQUID;
	} else {
	  iState = WATER_UNSTABLEGAS;
	}
	//setState_TR(T, rho);
	delta = deltaSave;

	m_phi->tdpolycalc(tau, delta);
      }
    }
  }
  return iState;
}


// Find the water spinodal density
doublereal WaterPropsIAPWS::densSpinodalWater() const {
  doublereal temperature = T_c/tau;
  doublereal delta_save = delta;
  // return the critical density if we are above or even just a little below
  // the critical temperature. We just don't want to worry about the critical
  // point at this juncture.
  if (temperature >= T_c - 0.001) {
    return Rho_c;
  }
  doublereal p = psat_est(temperature);
  doublereal rho_low = 0.0;
  doublereal rho_high = 1000;
  
  doublereal densSatLiq = density_const(p, WATER_LIQUID);
  doublereal dens_old = densSatLiq;
  delta = dens_old / Rho_c;
  m_phi->tdpolycalc(tau, delta);
  doublereal dpdrho_old = dpdrho();
  if (dpdrho_old > 0.0) {
    rho_high = MIN(dens_old, rho_high);
  } else {
    rho_low = MAX(rho_low, dens_old);
  }
  doublereal dens_new = densSatLiq* (1.0001);
  delta = dens_new / Rho_c;
  m_phi->tdpolycalc(tau, delta);
  doublereal dpdrho_new = dpdrho();
  if (dpdrho_new > 0.0) {
    rho_high = MIN(dens_new, rho_high);
  } else {
    rho_low = MAX(rho_low, dens_new);
  }
  bool conv = false;
 
  for (int it = 0; it < 50; it++) {
    doublereal slope = (dpdrho_new - dpdrho_old)/(dens_new - dens_old);
    if (slope >= 0.0) {
      slope = MAX(slope, dpdrho_new *5.0/ dens_new);
    } else {  
      slope = -dpdrho_new;
      //slope = MIN(slope, dpdrho_new *5.0 / dens_new);
      // shouldn't be here for liquid spinodal
    }
    doublereal delta_rho =  - dpdrho_new / slope;
    if (delta_rho > 0.0) {
      delta_rho = MIN(delta_rho, dens_new * 0.1);
    } else {
      delta_rho = MAX(delta_rho, - dens_new * 0.1);
    }
    doublereal dens_est =  dens_new + delta_rho;
    if (dens_est < rho_low) {
      dens_est = 0.5 * (rho_low + dens_new); 
    }
    if (dens_est > rho_high) {
      dens_est = 0.5 * (rho_high + dens_new); 
    }

    
    dens_old = dens_new;
    dpdrho_old = dpdrho_new;
    dens_new = dens_est;

    delta = dens_new / Rho_c;
    m_phi->tdpolycalc(tau, delta);
    dpdrho_new = dpdrho();
    if (dpdrho_new > 0.0) {
      rho_high = MIN(dens_new, rho_high);
    } else  if (dpdrho_new < 0.0) {
      rho_low = MAX(rho_low, dens_new);
    } else {
      conv = true;
      break;
    }
    
    if (fabs(dpdrho_new) < 1.0E-5) {
      conv = true;
      break;
    }
  }

  if (!conv) {
    throw Cantera::CanteraError(" WaterPropsIAPWS::densSpinodalWater()",
				" convergence failure");
  }
  // Restore the original delta
  delta = delta_save;
  m_phi->tdpolycalc(tau, delta);

  return dens_new;
}


// Find the steam spinodal density
doublereal WaterPropsIAPWS::densSpinodalSteam() const {
  doublereal temperature = T_c/tau;
  doublereal delta_save = delta;
  // return the critical density if we are above or even just a little below
  // the critical temperature. We just don't want to worry about the critical
  // point at this juncture.
  if (temperature >= T_c - 0.001) {
    return Rho_c;
  }
  doublereal p = psat_est(temperature);
  doublereal rho_low = 0.0;
  doublereal rho_high = 1000;
  
  doublereal densSatGas = density_const(p, WATER_GAS);
  doublereal dens_old = densSatGas;
  delta = dens_old / Rho_c;
  m_phi->tdpolycalc(tau, delta);
  doublereal dpdrho_old = dpdrho();
  if (dpdrho_old < 0.0) {
    rho_high = MIN(dens_old, rho_high);
  } else {
    rho_low = MAX(rho_low, dens_old);
  }
  doublereal dens_new = densSatGas * (0.99);
  delta = dens_new / Rho_c;
  m_phi->tdpolycalc(tau, delta);
  doublereal dpdrho_new = dpdrho();
  if (dpdrho_new < 0.0) {
    rho_high = MIN(dens_new, rho_high);
  } else {
    rho_low = MAX(rho_low, dens_new);
  }
  bool conv = false;
 
  for (int it = 0; it < 50; it++) {
    doublereal slope = (dpdrho_new - dpdrho_old)/(dens_new - dens_old);
    if (slope >= 0.0) {
      slope = dpdrho_new;
      //slope = MAX(slope, dpdrho_new *5.0/ dens_new);
      // shouldn't be here for gas spinodal
    } else {  
      //slope = -dpdrho_new;
      slope = MIN(slope, dpdrho_new *5.0 / dens_new);
    
    }
    doublereal delta_rho = - dpdrho_new / slope;
    if (delta_rho > 0.0) {
      delta_rho = MIN(delta_rho, dens_new * 0.1);
    } else {
      delta_rho = MAX(delta_rho, - dens_new * 0.1);
    }
    doublereal dens_est =  dens_new + delta_rho;
    if (dens_est < rho_low) {
      dens_est = 0.5 * (rho_low + dens_new); 
    }
    if (dens_est > rho_high) {
      dens_est = 0.5 * (rho_high + dens_new); 
    }

    
    dens_old = dens_new;
    dpdrho_old = dpdrho_new;
    dens_new = dens_est;

    delta = dens_new / Rho_c;
    m_phi->tdpolycalc(tau, delta);
    dpdrho_new = dpdrho();
    if (dpdrho_new < 0.0) {
      rho_high = MIN(dens_new, rho_high);
    } else  if (dpdrho_new > 0.0) {
      rho_low = MAX(rho_low, dens_new);
    } else {
      conv = true;
      break;
    }
    
    if (fabs(dpdrho_new) < 1.0E-5) {
      conv = true;
      break;
    }
  }

  if (!conv) {
    throw Cantera::CanteraError(" WaterPropsIAPWS::densSpinodalSteam()",
				" convergence failure");
  }
  // Restore the original delta
  delta = delta_save;
  m_phi->tdpolycalc(tau, delta);

  return dens_new;
}



/**
 * Sets the internal state of the object to the
 * specified temperature and density.
 */
void WaterPropsIAPWS::setState_TR(doublereal temperature, doublereal rho) {
  calcDim(temperature, rho);
  m_phi->tdpolycalc(tau, delta);
}


/*
 * Calculate the enthalpy in mks units of
 * J kmol-1 K-1.
 */
doublereal WaterPropsIAPWS::enthalpy() const {
  doublereal temperature = T_c/tau;
  doublereal hRT =  m_phi->enthalpy_RT();
  return (hRT * Rgas * temperature);
}


/*
 * Calculate the internal Energy in mks units of
 * J kmol-1 K-1.
 */
doublereal WaterPropsIAPWS::intEnergy() const {
  doublereal temperature = T_c / tau;
  doublereal uRT = m_phi->intEnergy_RT();
  return (uRT * Rgas * temperature);
}

/*
 * Calculate the enthalpy in mks units of356
 * J kmol-1 K-1.
 */
doublereal WaterPropsIAPWS::entropy() const {
  doublereal sR = m_phi->entropy_R();
  return (sR * Rgas);
}

/*
 * Calculate heat capacity at constant volume
 * J kmol-1 K-1.
 */
doublereal WaterPropsIAPWS::cv() const {
  doublereal cvR = m_phi->cv_R();
  return (cvR * Rgas);
}

doublereal  WaterPropsIAPWS::cp() const {
  doublereal cpR = m_phi->cp_R();
  return (cpR * Rgas);
}

doublereal WaterPropsIAPWS::molarVolume() const {
  doublereal rho = delta * Rho_c;
  return (M_water / rho);
}
