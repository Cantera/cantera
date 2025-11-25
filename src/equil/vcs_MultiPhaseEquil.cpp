/**
 *  @file vcs_MultiPhaseEquil.cpp
 *    Driver routine for the VCSnonideal equilibrium solver package
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/equil/vcs_MultiPhaseEquil.h"
#include "cantera/equil/vcs_VolPhase.h"
#include "cantera/base/stringUtils.h"
#include "cantera/thermo/speciesThermoTypes.h"
#include "cantera/thermo/IdealSolidSolnPhase.h"
#include "cantera/thermo/IdealMolalSoln.h"

namespace Cantera
{
vcs_MultiPhaseEquil::vcs_MultiPhaseEquil(MultiPhase* mix, int printLvl) :
    m_mix(mix),
    m_printLvl(printLvl),
    m_vsolve(mix, printLvl)
{
}

int vcs_MultiPhaseEquil::equilibrate_TV(int XY, double xtarget, int estimateEquil,
                                        int printLvl, double err,
                                        int maxsteps, int loglevel)
{
    double Vtarget = m_mix->volume();
    if ((XY != TV) && (XY != HV) && (XY != UV) && (XY != SV)) {
        throw CanteraError("vcs_MultiPhaseEquil::equilibrate_TV",
                           "Wrong XY flag: {}", XY);
    }
    int maxiter = 100;
    int iSuccess = 0;
    if (XY == TV) {
        m_mix->setTemperature(xtarget);
    }
    int strt = estimateEquil;
    double P1 = 0.0;
    double V1 = 0.0;
    double V2 = 0.0;
    double P2 = 0.0;
    double Tlow = 0.5 * m_mix->minTemp();
    double Thigh = 2.0 * m_mix->maxTemp();
    int printLvlSub = std::max(0, printLvl - 1);
    for (int n = 0; n < maxiter; n++) {
        double Pnow = m_mix->pressure();

        switch (XY) {
        case TV:
            iSuccess = equilibrate_TP(strt, printLvlSub, err, maxsteps, loglevel);
            break;
        case HV:
            iSuccess = equilibrate_HP(xtarget, HP, Tlow, Thigh, strt,
                                      printLvlSub, err, maxsteps, loglevel);
            break;
        case UV:
            iSuccess = equilibrate_HP(xtarget, UP, Tlow, Thigh, strt,
                                      printLvlSub, err, maxsteps, loglevel);
            break;
        case SV:
            iSuccess = equilibrate_SP(xtarget, Tlow, Thigh, strt,
                                      printLvlSub, err, maxsteps, loglevel);
            break;
        default:
            break;
        }
        strt = false;
        double Vnow = m_mix->volume();
        if (n == 0) {
            V2 = Vnow;
            P2 = Pnow;
        } else if (n == 1) {
            V1 = Vnow;
            P1 = Pnow;
        } else {
            P2 = P1;
            V2 = V1;
            P1 = Pnow;
            V1 = Vnow;
        }

        double Verr = fabs((Vtarget - Vnow)/Vtarget);
        if (Verr < err) {
            return iSuccess;
        }
        double Pnew;
        // find dV/dP
        if (n > 1) {
            double dVdP = (V2 - V1) / (P2 - P1);
            if (dVdP == 0.0) {
                throw CanteraError("vcs_MultiPhase::equilibrate_TV",
                                   "dVdP == 0.0");
            } else {
                Pnew = Pnow + (Vtarget - Vnow) / dVdP;
                if (Pnew < 0.2 * Pnow) {
                    Pnew = 0.2 * Pnow;
                }
                if (Pnew > 3.0 * Pnow) {
                    Pnew = 3.0 * Pnow;
                }
            }
        } else {
            m_mix->setPressure(Pnow*1.01);
            double dVdP = (m_mix->volume() - Vnow)/(0.01*Pnow);
            Pnew = Pnow + 0.5*(Vtarget - Vnow)/dVdP;
            if (Pnew < 0.5* Pnow) {
                Pnew = 0.5 * Pnow;
            }
            if (Pnew > 1.7 * Pnow) {
                Pnew = 1.7 * Pnow;
            }
        }
        m_mix->setPressure(Pnew);
    }
    throw CanteraError("vcs_MultiPhase::equilibrate_TV",
                       "No convergence for V");
}

int vcs_MultiPhaseEquil::equilibrate_HP(double Htarget, int XY, double Tlow,
    double Thigh, int estimateEquil, int printLvl, double err, int maxsteps,
    int loglevel)
{
    int maxiter = 100;
    int iSuccess;
    if (XY != HP && XY != UP) {
        throw CanteraError("vcs_MultiPhaseEquil::equilibrate_HP",
                           "Wrong XP", XY);
    }
    int strt = estimateEquil;

    // Lower bound on T. This will change as we progress in the calculation
    if (Tlow <= 0.0) {
        Tlow = 0.5 * m_mix->minTemp();
    }
    // Upper bound on T. This will change as we progress in the calculation
    if (Thigh <= 0.0 || Thigh > 1.0E6) {
        Thigh = 2.0 * m_mix->maxTemp();
    }

    double cpb = 1.0;
    double Hlow = Undef;
    double Hhigh = Undef;
    double Tnow = m_mix->temperature();
    int printLvlSub = std::max(printLvl - 1, 0);

    for (int n = 0; n < maxiter; n++) {
        // start with a loose error tolerance, but tighten it as we get
        // close to the final temperature
        try {
            Tnow = m_mix->temperature();
            iSuccess = equilibrate_TP(strt, printLvlSub, err, maxsteps, loglevel);
            strt = 0;
            double Hnow = (XY == UP) ? m_mix->IntEnergy() : m_mix->enthalpy();
            double pmoles[10];
            pmoles[0] = m_mix->phaseMoles(0);
            double Tmoles = pmoles[0];
            double HperMole = Hnow/Tmoles;
            if (printLvl > 0) {
                plogf("T = %g, Hnow = %g ,Tmoles = %g,  HperMole = %g\n",
                      Tnow, Hnow, Tmoles, HperMole);
            }

            // the equilibrium enthalpy monotonically increases with T;
            // if the current value is below the target, then we know the
            // current temperature is too low. Set the lower bounds.
            if (Hnow < Htarget) {
                if (Tnow > Tlow) {
                    Tlow = Tnow;
                    Hlow = Hnow;
                }
            } else {
                // the current enthalpy is greater than the target; therefore
                // the current temperature is too high. Set the high bounds.
                if (Tnow < Thigh) {
                    Thigh = Tnow;
                    Hhigh = Hnow;
                }
            }
            double dT;
            if (Hlow != Undef && Hhigh != Undef) {
                cpb = (Hhigh - Hlow)/(Thigh - Tlow);
                dT = (Htarget - Hnow)/cpb;
                double dTa = fabs(dT);
                double dTmax = 0.5*fabs(Thigh - Tlow);
                if (dTa > dTmax) {
                    dT *= dTmax/dTa;
                }
            } else {
                double Tnew = sqrt(Tlow*Thigh);
                dT = clip(Tnew - Tnow, -200.0, 200.0);
            }
            double acpb = std::max(fabs(cpb), 1.0E-6);
            double denom = std::max(fabs(Htarget), acpb);
            double Herr = Htarget - Hnow;
            double HConvErr = fabs((Herr)/denom);
            if (printLvl > 0) {
                plogf("   equilibrate_HP: It = %d, Tcurr  = %g Hcurr = %g, Htarget = %g\n",
                      n, Tnow, Hnow, Htarget);
                plogf("                   H rel error = %g, cp = %g, HConvErr = %g\n",
                      Herr, cpb, HConvErr);
            }

            if (HConvErr < err) { // || dTa < 1.0e-4) {
                if (printLvl > 0) {
                    plogf("   equilibrate_HP: CONVERGENCE: Hfinal  = %g Tfinal = %g, Its = %d \n",
                          Hnow, Tnow, n);
                    plogf("                   H rel error = %g, cp = %g, HConvErr = %g\n",
                          Herr, cpb, HConvErr);
                }
                return iSuccess;
            }
            double Tnew = Tnow + dT;
            if (Tnew < 0.0) {
                Tnew = 0.5*Tnow;
            }
            m_mix->setTemperature(Tnew);
        } catch (const CanteraError&) {
            if (!estimateEquil) {
                strt = -1;
            } else {
                double Tnew = 0.5*(Tnow + Thigh);
                if (fabs(Tnew - Tnow) < 1.0) {
                    Tnew = Tnow + 1.0;
                }
                m_mix->setTemperature(Tnew);
            }
        }
    }
    throw CanteraError("vcs_MultiPhaseEquil::equilibrate_HP",
                       "No convergence for T");
}

int vcs_MultiPhaseEquil::equilibrate_SP(double Starget, double Tlow, double Thigh,
    int estimateEquil, int printLvl, double err, int maxsteps, int loglevel)
{
    int maxiter = 100;
    int strt = estimateEquil;

    // Lower bound on T. This will change as we progress in the calculation
    if (Tlow <= 0.0) {
        Tlow = 0.5 * m_mix->minTemp();
    }
    // Upper bound on T. This will change as we progress in the calculation
    if (Thigh <= 0.0 || Thigh > 1.0E6) {
        Thigh = 2.0 * m_mix->maxTemp();
    }

    double cpb = 1.0, dT;
    double Slow = Undef;
    double Shigh = Undef;
    double Tnow = m_mix->temperature();
    Tlow = std::min(Tnow, Tlow);
    Thigh = std::max(Tnow, Thigh);
    int printLvlSub = std::max(printLvl - 1, 0);

    for (int n = 0; n < maxiter; n++) {
        // start with a loose error tolerance, but tighten it as we get
        // close to the final temperature
        try {
            Tnow = m_mix->temperature();
            int iSuccess = equilibrate_TP(strt, printLvlSub, err, maxsteps, loglevel);
            strt = 0;
            double Snow = m_mix->entropy();
            double pmoles[10];
            pmoles[0] = m_mix->phaseMoles(0);
            double Tmoles = pmoles[0];
            double SperMole = Snow/Tmoles;
            if (printLvl > 0) {
                plogf("T = %g, Snow = %g ,Tmoles = %g,  SperMole = %g\n",
                      Tnow, Snow, Tmoles, SperMole);
            }

            // the equilibrium entropy monotonically increases with T;
            // if the current value is below the target, then we know the
            // current temperature is too low. Set the lower bounds to the
            // current condition.
            if (Snow < Starget) {
                if (Tnow > Tlow) {
                    Tlow = Tnow;
                    Slow = Snow;
                } else {
                    if (Slow > Starget && Snow < Slow) {
                        Thigh = Tlow;
                        Shigh = Slow;
                        Tlow = Tnow;
                        Slow = Snow;
                    }
                }
            } else {
                // the current enthalpy is greater than the target; therefore
                // the current temperature is too high. Set the high bounds.
                if (Tnow < Thigh) {
                    Thigh = Tnow;
                    Shigh = Snow;
                }
            }
            if (Slow != Undef && Shigh != Undef) {
                cpb = (Shigh - Slow)/(Thigh - Tlow);
                dT = (Starget - Snow)/cpb;
                double Tnew = Tnow + dT;
                double dTa = fabs(dT);
                double dTmax = 0.5*fabs(Thigh - Tlow);
                if (Tnew > Thigh || Tnew < Tlow) {
                    dTmax = 1.5*fabs(Thigh - Tlow);
                }
                dTmax = std::min(dTmax, 300.);
                if (dTa > dTmax) {
                    dT *= dTmax/dTa;
                }
            } else {
                double Tnew = sqrt(Tlow*Thigh);
                dT = Tnew - Tnow;
            }

            double acpb = std::max(fabs(cpb), 1.0E-6);
            double denom = std::max(fabs(Starget), acpb);
            double Serr = Starget - Snow;
            double SConvErr = fabs((Serr)/denom);
            if (printLvl > 0) {
                plogf("   equilibrate_SP: It = %d, Tcurr  = %g Scurr = %g, Starget = %g\n",
                      n, Tnow, Snow, Starget);
                plogf("                   S rel error = %g, cp = %g, SConvErr = %g\n",
                      Serr, cpb, SConvErr);
            }

            if (SConvErr < err) { // || dTa < 1.0e-4) {
                if (printLvl > 0) {
                    plogf("   equilibrate_SP: CONVERGENCE: Sfinal  = %g Tfinal = %g, Its = %d \n",
                          Snow, Tnow, n);
                    plogf("                   S rel error = %g, cp = %g, HConvErr = %g\n",
                          Serr, cpb, SConvErr);
                }
                return iSuccess;
            }
            double Tnew = Tnow + dT;
            if (Tnew < 0.0) {
                Tnew = 0.5*Tnow;
            }
            m_mix->setTemperature(Tnew);
        } catch (const CanteraError&) {
            if (!estimateEquil) {
                strt = -1;
            } else {
                double Tnew = 0.5*(Tnow + Thigh);
                if (fabs(Tnew - Tnow) < 1.0) {
                    Tnew = Tnow + 1.0;
                }
                m_mix->setTemperature(Tnew);
            }
        }
    }
    throw CanteraError("vcs_MultiPhaseEquil::equilibrate_SP",
                       "No convergence for T");
}

int vcs_MultiPhaseEquil::equilibrate(int XY, int estimateEquil, int printLvl,
                                     double err, int maxsteps, int loglevel)
{
    double xtarget;
    if (XY == TP) {
        return equilibrate_TP(estimateEquil, printLvl, err, maxsteps, loglevel);
    } else if (XY == HP || XY == UP) {
        if (XY == HP) {
            xtarget = m_mix->enthalpy();
        } else {
            xtarget = m_mix->IntEnergy();
        }
        double Tlow = 0.5 * m_mix->minTemp();
        double Thigh = 2.0 * m_mix->maxTemp();
        return equilibrate_HP(xtarget, XY, Tlow, Thigh,
                              estimateEquil, printLvl, err, maxsteps, loglevel);
    } else if (XY == SP) {
        xtarget = m_mix->entropy();
        double Tlow = 0.5 * m_mix->minTemp();
        double Thigh = 2.0 * m_mix->maxTemp();
        return equilibrate_SP(xtarget, Tlow, Thigh,
                              estimateEquil, printLvl, err, maxsteps, loglevel);
    } else if (XY == TV) {
        xtarget = m_mix->temperature();
        return equilibrate_TV(XY, xtarget,
                              estimateEquil, printLvl, err, maxsteps, loglevel);
    } else if (XY == HV) {
        xtarget = m_mix->enthalpy();
        return equilibrate_TV(XY, xtarget,
                              estimateEquil, printLvl, err, maxsteps, loglevel);
    } else if (XY == UV) {
        xtarget = m_mix->IntEnergy();
        return equilibrate_TV(XY, xtarget,
                              estimateEquil, printLvl, err, maxsteps, loglevel);
    } else if (XY == SV) {
        xtarget = m_mix->entropy();
        return equilibrate_TV(XY, xtarget, estimateEquil,
                                  printLvl, err, maxsteps, loglevel);
    } else {
        throw CanteraError("vcs_MultiPhaseEquil::equilibrate",
                           "Unsupported Option");
    }
}

int vcs_MultiPhaseEquil::equilibrate_TP(int estimateEquil, int printLvl, double err,
                                        int maxsteps, int loglevel)
{
    int maxit = maxsteps;
    m_printLvl = printLvl;
    m_vsolve.m_printLvl = printLvl;
    m_vsolve.m_doEstimateEquil = estimateEquil;

    // Check obvious bounds on the temperature and pressure NOTE, we may want to
    // do more here with the real bounds given by the ThermoPhase objects.
    if (m_mix->temperature() <= 0.0) {
        throw CanteraError("vcs_MultiPhaseEquil::equilibrate_TP",
                           "Temperature less than zero on input");
    }
    if (m_mix->pressure() <= 0.0) {
        throw CanteraError("vcs_MultiPhaseEquil::equilibrate_TP",
                           "Pressure less than zero on input");
    }

    //! Call the thermo Program
    int ip1 = m_printLvl;
    int ipr = std::max(0, m_printLvl-1);
    if (m_printLvl >= 3) {
        ip1 = m_printLvl - 2;
    } else {
        ip1 = 0;
    }
    int iSuccess = m_vsolve.vcs(ipr, ip1, maxit);

    if (printLvl > 0) {
        vector<double> mu(m_mix->nSpecies());
        m_mix->getChemPotentials(mu.data());
        plogf("\n Results from vcs:\n");
        if (iSuccess != 0) {
            plogf("\nVCS FAILED TO CONVERGE!\n");
        }
        plogf("\n");
        plogf("Temperature = %g Kelvin\n", m_vsolve.m_temperature);
        plogf("Pressure    = %g Pa\n", m_vsolve.m_pressurePA);
        plogf("\n");
        plogf("----------------------------------------"
              "---------------------\n");
        plogf(" Name             Mole_Number(kmol)");
        plogf("  Mole_Fraction     Chem_Potential (J/kmol)\n");
        plogf("--------------------------------------------------"
              "-----------\n");
        for (size_t i = 0; i < m_mix->nSpecies(); i++) {
            plogf("%-12s", m_mix->speciesName(i));
            if (m_vsolve.m_speciesUnknownType[i] == VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
                plogf("  %15.3e %15.3e  ", 0.0, m_mix->moleFraction(i));
                plogf("%15.3e\n", mu[i]);
            } else {
                plogf("  %15.3e   %15.3e  ", m_mix->speciesMoles(i), m_mix->moleFraction(i));
                if (m_mix->speciesMoles(i) <= 0.0) {
                    size_t iph = m_vsolve.m_phaseID[i];
                    vcs_VolPhase* VPhase = m_vsolve.m_VolPhaseList[iph].get();
                    if (VPhase->nSpecies() > 1) {
                        plogf("     -1.000e+300\n");
                    } else {
                        plogf("%15.3e\n", mu[i]);
                    }
                } else {
                    plogf("%15.3e\n", mu[i]);
                }
            }
        }
        plogf("------------------------------------------"
              "-------------------\n");
    }
    return iSuccess;
}

}
