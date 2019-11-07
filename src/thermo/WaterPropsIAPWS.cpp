/**
 * @file WaterPropsIAPWS.cpp
 * Definitions for a class for calculating the equation of state of water
 * from the IAPWS 1995 Formulation based on the steam tables thermodynamic
 * basis (See class \link Cantera::WaterPropsIAPWS WaterPropsIAPWS\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/WaterPropsIAPWS.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/stringUtils.h"

namespace Cantera
{
// Critical Point values of water in mks units

//! Critical Temperature value (kelvin)
const doublereal T_c = 647.096;
//! Critical Pressure (Pascals)
static const doublereal P_c = 22.064E6;
//! Value of the Density at the critical point (kg m-3)
const doublereal Rho_c = 322.;
//! Molecular Weight of water that is consistent with the paper (kg kmol-1)
static const doublereal M_water = 18.015268;

//! Gas constant that is quoted in the paper
/*
 * Note, this is the Rgas value quoted in the paper. For consistency
 * we have to use that value and not the updated value
 *
 * The Ratio of R/M = 0.46151805 kJ kg-1 K-1 , which is Eqn. (6.3) in the paper.
 */
static const doublereal Rgas = 8.314371E3; // Joules kmol-1 K-1

// Base constructor
WaterPropsIAPWS::WaterPropsIAPWS() :
    tau(-1.0),
    delta(-1.0),
    iState(-30000)
{
}

void WaterPropsIAPWS::calcDim(doublereal temperature, doublereal rho)
{
    tau = T_c / temperature;
    delta = rho / Rho_c;

    // Determine the internal state
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

doublereal WaterPropsIAPWS::helmholtzFE() const
{
    doublereal retn = m_phi.phi(tau, delta);
    doublereal temperature = T_c/tau;
    doublereal RT = Rgas * temperature;
    return retn * RT;
}

doublereal WaterPropsIAPWS::pressure() const
{
    doublereal retn = m_phi.pressureM_rhoRT(tau, delta);
    doublereal rho = delta * Rho_c;
    doublereal temperature = T_c / tau;
    return retn * rho * Rgas * temperature/M_water;
}

doublereal WaterPropsIAPWS::density(doublereal temperature, doublereal pressure,
                                    int phase, doublereal rhoguess)
{
    doublereal deltaGuess = 0.0;
    if (rhoguess == -1.0) {
        if (phase != -1) {
            if (temperature > T_c) {
                rhoguess = pressure * M_water / (Rgas * temperature);
            } else {
                if (phase == WATER_GAS || phase == WATER_SUPERCRIT) {
                    rhoguess = pressure * M_water / (Rgas * temperature);
                } else if (phase == WATER_LIQUID) {
                    // Provide a guess about the liquid density that is
                    // relatively high -> convergence from above seems robust.
                    rhoguess = 1000.;
                } else if (phase == WATER_UNSTABLELIQUID || phase == WATER_UNSTABLEGAS) {
                    throw CanteraError("WaterPropsIAPWS::density",
                                       "Unstable Branch finder is untested");
                } else {
                    throw CanteraError("WaterPropsIAPWS::density",
                                       "unknown state: {}", phase);
                }
            }
        } else {
            // Assume the Gas phase initial guess, if nothing is specified to
            // the routine
            rhoguess = pressure * M_water / (Rgas * temperature);
        }
    }
    doublereal p_red = pressure * M_water / (Rgas * temperature * Rho_c);
    deltaGuess = rhoguess / Rho_c;
    setState_TR(temperature, rhoguess);
    doublereal delta_retn = m_phi.dfind(p_red, tau, deltaGuess);
    doublereal density_retn;
    if (delta_retn >0.0) {
        delta = delta_retn;

        // Dimensionalize the density before returning
        density_retn = delta_retn * Rho_c;

        // Set the internal state -> this may be a duplication. However, let's
        // just be sure.
        setState_TR(temperature, density_retn);
    } else {
        density_retn = -1.0;
    }
    return density_retn;
}

doublereal WaterPropsIAPWS::density_const(doublereal pressure,
        int phase, doublereal rhoguess) const
{
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
                    // Provide a guess about the liquid density that is
                    // relatively high -> convergence from above seems robust.
                    rhoguess = 1000.;
                } else if (phase == WATER_UNSTABLELIQUID || phase == WATER_UNSTABLEGAS) {
                    throw CanteraError("WaterPropsIAPWS::density_const",
                                       "Unstable Branch finder is untested");
                } else {
                    throw CanteraError("WaterPropsIAPWS::density_const",
                                       "unknown state: {}", phase);
                }
            }
        } else {
            // Assume the Gas phase initial guess, if nothing is specified to
            // the routine
            rhoguess = pressure * M_water / (Rgas * temperature);
        }
    }
    doublereal p_red = pressure * M_water / (Rgas * temperature * Rho_c);
    deltaGuess = rhoguess / Rho_c;

    delta = deltaGuess;
    m_phi.tdpolycalc(tau, delta);

    doublereal delta_retn = m_phi.dfind(p_red, tau, deltaGuess);
    doublereal density_retn;
    if (delta_retn > 0.0) {
        delta = delta_retn;

        // Dimensionalize the density before returning
        density_retn = delta_retn * Rho_c;

    } else {
        density_retn = -1.0;
    }

    delta = deltaSave;
    m_phi.tdpolycalc(tau, delta);
    return density_retn;
}

doublereal WaterPropsIAPWS::density() const
{
    return delta * Rho_c;
}

doublereal WaterPropsIAPWS::temperature() const
{
    return T_c / tau;
}

doublereal WaterPropsIAPWS::psat_est(doublereal temperature) const
{
    // Formula and constants from: "NBS/NRC Steam Tables: Thermodynamic and
    // Transport Properties and Computer Programs for Vapor and Liquid States of
    // Water in SI Units". L. Haar, J. S. Gallagher, G. S. Kell. Hemisphere
    // Publishing. 1984.
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

    // Original correlation was in cgs. Convert to mks
    ps *= 1.0E6;
    return ps;
}

doublereal WaterPropsIAPWS::isothermalCompressibility() const
{
    doublereal dpdrho_val = dpdrho();
    doublereal dens = delta * Rho_c;
    return 1.0 / (dens * dpdrho_val);
}

doublereal WaterPropsIAPWS::dpdrho() const
{
    doublereal retn = m_phi.dimdpdrho(tau, delta);
    doublereal temperature = T_c/tau;
    return retn * Rgas * temperature / M_water;
}

doublereal WaterPropsIAPWS::coeffPresExp() const
{
    return m_phi.dimdpdT(tau, delta);
}

doublereal WaterPropsIAPWS::coeffThermExp() const
{
    doublereal kappa = isothermalCompressibility();
    doublereal beta = coeffPresExp();
    doublereal dens = delta * Rho_c;
    return kappa * dens * Rgas * beta / M_water;
}

doublereal WaterPropsIAPWS::Gibbs() const
{
    doublereal gRT = m_phi.gibbs_RT();
    doublereal temperature = T_c/tau;
    return gRT * Rgas * temperature;
}

void WaterPropsIAPWS::corr(doublereal temperature, doublereal pressure,
    doublereal& densLiq, doublereal& densGas, doublereal& delGRT)
{
    densLiq = density(temperature, pressure, WATER_LIQUID, densLiq);
    if (densLiq <= 0.0) {
        throw CanteraError("WaterPropsIAPWS::corr",
            "Error occurred trying to find liquid density at (T,P) = {}  {}",
            temperature, pressure);
    }
    setState_TR(temperature, densLiq);
    doublereal gibbsLiqRT = m_phi.gibbs_RT();

    densGas = density(temperature, pressure, WATER_GAS, densGas);
    if (densGas <= 0.0) {
        throw CanteraError("WaterPropsIAPWS::corr",
            "Error occurred trying to find gas density at (T,P) = {}  {}",
            temperature, pressure);
    }
    setState_TR(temperature, densGas);
    doublereal gibbsGasRT = m_phi.gibbs_RT();

    delGRT = gibbsLiqRT - gibbsGasRT;
}

void WaterPropsIAPWS::corr1(doublereal temperature, doublereal pressure,
    doublereal& densLiq, doublereal& densGas, doublereal& pcorr)
{
    densLiq = density(temperature, pressure, WATER_LIQUID, densLiq);
    if (densLiq <= 0.0) {
        throw CanteraError("WaterPropsIAPWS::corr1",
            "Error occurred trying to find liquid density at (T,P) = {}  {}",
            temperature, pressure);
    }
    setState_TR(temperature, densLiq);
    doublereal prL = m_phi.phiR();

    densGas = density(temperature, pressure, WATER_GAS, densGas);
    if (densGas <= 0.0) {
        throw CanteraError("WaterPropsIAPWS::corr1",
            "Error occurred trying to find gas density at (T,P) = {}  {}",
            temperature, pressure);
    }
    setState_TR(temperature, densGas);
    doublereal prG = m_phi.phiR();
    doublereal rhs = (prL - prG) + log(densLiq/densGas);
    rhs /= (1.0/densGas - 1.0/densLiq);
    pcorr = rhs * Rgas * temperature / M_water;
}

doublereal WaterPropsIAPWS::psat(doublereal temperature, int waterState)
{
    static int method = 1;
    doublereal densLiq = -1.0, densGas = -1.0, delGRT = 0.0;
    doublereal dp, pcorr;
    if (temperature >= T_c) {
        densGas = density(temperature, P_c, WATER_SUPERCRIT);
        setState_TR(temperature, densGas);
        return P_c;
    }
    doublereal p = psat_est(temperature);
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
            break;
        } else {
            if (fabs(dp/p) < 1.0E-9) {
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
        throw CanteraError("WaterPropsIAPWS::psat",
                           "unknown water state input: {}", waterState);
    }
    return p;
}

int WaterPropsIAPWS::phaseState(bool checkState) const
{
    if (checkState) {
        if (tau <= 1.0) {
            iState = WATER_SUPERCRIT;
        } else {
            doublereal T = T_c / tau;
            doublereal rho = delta * Rho_c;
            doublereal rhoMidAtm = 0.5 * (OneAtm * M_water / (Rgas * 373.15) + 1.0E3);
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
                doublereal deltaSave = delta;
                doublereal deltaDel = rhoDel / Rho_c;
                delta = deltaDel;
                m_phi.tdpolycalc(tau, deltaDel);

                doublereal kappaDel = isothermalCompressibility();
                doublereal d2rhodp2 = (rhoDel * kappaDel - rho * kappa) / (rhoDel - rho);
                if (d2rhodp2 > 0.0) {
                    iState = WATER_UNSTABLELIQUID;
                } else {
                    iState = WATER_UNSTABLEGAS;
                }
                delta = deltaSave;
                m_phi.tdpolycalc(tau, delta);
            }
        }
    }
    return iState;
}

doublereal WaterPropsIAPWS::densSpinodalWater() const
{
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
    m_phi.tdpolycalc(tau, delta);
    doublereal dpdrho_old = dpdrho();
    if (dpdrho_old > 0.0) {
        rho_high = std::min(dens_old, rho_high);
    } else {
        rho_low = std::max(rho_low, dens_old);
    }
    doublereal dens_new = densSatLiq* (1.0001);
    delta = dens_new / Rho_c;
    m_phi.tdpolycalc(tau, delta);
    doublereal dpdrho_new = dpdrho();
    if (dpdrho_new > 0.0) {
        rho_high = std::min(dens_new, rho_high);
    } else {
        rho_low = std::max(rho_low, dens_new);
    }
    bool conv = false;

    for (int it = 0; it < 50; it++) {
        doublereal slope = (dpdrho_new - dpdrho_old)/(dens_new - dens_old);
        if (slope >= 0.0) {
            slope = std::max(slope, dpdrho_new *5.0/ dens_new);
        } else {
            slope = -dpdrho_new;
            // shouldn't be here for liquid spinodal
        }
        doublereal delta_rho = - dpdrho_new / slope;
        if (delta_rho > 0.0) {
            delta_rho = std::min(delta_rho, dens_new * 0.1);
        } else {
            delta_rho = std::max(delta_rho, - dens_new * 0.1);
        }
        doublereal dens_est = dens_new + delta_rho;
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
        m_phi.tdpolycalc(tau, delta);
        dpdrho_new = dpdrho();
        if (dpdrho_new > 0.0) {
            rho_high = std::min(dens_new, rho_high);
        } else if (dpdrho_new < 0.0) {
            rho_low = std::max(rho_low, dens_new);
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
        throw CanteraError("WaterPropsIAPWS::densSpinodalWater",
                           "convergence failure");
    }
    // Restore the original delta
    delta = delta_save;
    m_phi.tdpolycalc(tau, delta);
    return dens_new;
}

doublereal WaterPropsIAPWS::densSpinodalSteam() const
{
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
    m_phi.tdpolycalc(tau, delta);
    doublereal dpdrho_old = dpdrho();
    if (dpdrho_old < 0.0) {
        rho_high = std::min(dens_old, rho_high);
    } else {
        rho_low = std::max(rho_low, dens_old);
    }
    doublereal dens_new = densSatGas * (0.99);
    delta = dens_new / Rho_c;
    m_phi.tdpolycalc(tau, delta);
    doublereal dpdrho_new = dpdrho();
    if (dpdrho_new < 0.0) {
        rho_high = std::min(dens_new, rho_high);
    } else {
        rho_low = std::max(rho_low, dens_new);
    }
    bool conv = false;
    for (int it = 0; it < 50; it++) {
        doublereal slope = (dpdrho_new - dpdrho_old)/(dens_new - dens_old);
        if (slope >= 0.0) {
            slope = dpdrho_new;
            // shouldn't be here for gas spinodal
        } else {
            slope = std::min(slope, dpdrho_new *5.0 / dens_new);

        }
        doublereal delta_rho = - dpdrho_new / slope;
        if (delta_rho > 0.0) {
            delta_rho = std::min(delta_rho, dens_new * 0.1);
        } else {
            delta_rho = std::max(delta_rho, - dens_new * 0.1);
        }
        doublereal dens_est = dens_new + delta_rho;
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
        m_phi.tdpolycalc(tau, delta);
        dpdrho_new = dpdrho();
        if (dpdrho_new < 0.0) {
            rho_high = std::min(dens_new, rho_high);
        } else if (dpdrho_new > 0.0) {
            rho_low = std::max(rho_low, dens_new);
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
        throw CanteraError("WaterPropsIAPWS::densSpinodalSteam",
                           "convergence failure");
    }
    // Restore the original delta
    delta = delta_save;
    m_phi.tdpolycalc(tau, delta);
    return dens_new;
}

void WaterPropsIAPWS::setState_TR(doublereal temperature, doublereal rho)
{
    calcDim(temperature, rho);
    m_phi.tdpolycalc(tau, delta);
}

doublereal WaterPropsIAPWS::enthalpy() const
{
    doublereal temperature = T_c/tau;
    doublereal hRT = m_phi.enthalpy_RT();
    return hRT * Rgas * temperature;
}

doublereal WaterPropsIAPWS::intEnergy() const
{
    doublereal temperature = T_c / tau;
    doublereal uRT = m_phi.intEnergy_RT();
    return uRT * Rgas * temperature;
}

doublereal WaterPropsIAPWS::entropy() const
{
    doublereal sR = m_phi.entropy_R();
    return sR * Rgas;
}

doublereal WaterPropsIAPWS::cv() const
{
    doublereal cvR = m_phi.cv_R();
    return cvR * Rgas;
}

doublereal WaterPropsIAPWS::cp() const
{
    doublereal cpR = m_phi.cp_R();
    return cpR * Rgas;
}

doublereal WaterPropsIAPWS::molarVolume() const
{
    doublereal rho = delta * Rho_c;
    return M_water / rho;
}

}
