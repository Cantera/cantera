/**
 *  @file vcs_MultiPhaseEquil.cpp
 *    Driver routine for the VCSnonideal equilibrium solver package
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/equil/vcs_MultiPhaseEquil.h"
#include "cantera/equil/vcs_VolPhase.h"
#include "cantera/equil/vcs_species_thermo.h"
#include "cantera/base/clockWC.h"
#include "cantera/base/stringUtils.h"
#include "cantera/thermo/speciesThermoTypes.h"
#include "cantera/thermo/IdealSolidSolnPhase.h"
#include "cantera/thermo/IdealMolalSoln.h"

#include <cstdio>

using namespace std;

namespace Cantera
{
vcs_MultiPhaseEquil::vcs_MultiPhaseEquil() :
    m_vprob(0, 0, 0),
    m_mix(0),
    m_printLvl(0)
{
}

vcs_MultiPhaseEquil::vcs_MultiPhaseEquil(MultiPhase* mix, int printLvl) :
    m_vprob(mix->nSpecies(), mix->nElements(), mix->nPhases()),
    m_mix(0),
    m_printLvl(printLvl)
{
    m_mix = mix;
    m_vprob.m_printLvl = m_printLvl;

    // Work out the details of the VCS_VPROB construction and Transfer the
    // current problem to VCS_PROB object
    int res = vcs_Cantera_to_vprob(mix, &m_vprob);
    if (res != 0) {
        plogf("problems\n");
    }
}

int vcs_MultiPhaseEquil::equilibrate_TV(int XY, doublereal xtarget,
                                        int estimateEquil,
                                        int printLvl, doublereal err,
                                        int maxsteps, int loglevel)
{
    doublereal Vtarget = m_mix->volume();
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
    doublereal Tlow = 0.5 * m_mix->minTemp();
    doublereal Thigh = 2.0 * m_mix->maxTemp();
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

int vcs_MultiPhaseEquil::equilibrate_HP(doublereal Htarget,
                                        int XY, double Tlow, double Thigh,
                                        int estimateEquil,
                                        int printLvl, doublereal err,
                                        int maxsteps, int loglevel)
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

    doublereal cpb = 1.0;
    doublereal Hlow = Undef;
    doublereal Hhigh = Undef;
    doublereal Tnow = m_mix->temperature();
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
                plogf("T = %g, Hnow = %g ,Tmoles = %g,  HperMole = %g",
                      Tnow, Hnow, Tmoles, HperMole);
                plogendl();
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
        } catch (CanteraError err) {
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
    throw CanteraError("MultiPhase::equilibrate_HP",
                       "No convergence for T");
}

int vcs_MultiPhaseEquil::equilibrate_SP(doublereal Starget,
                                        double Tlow, double Thigh,
                                        int estimateEquil,
                                        int printLvl, doublereal err,
                                        int maxsteps, int loglevel)
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

    doublereal cpb = 1.0, dT;
    doublereal Slow = Undef;
    doublereal Shigh = Undef;
    doublereal Tnow = m_mix->temperature();
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
        } catch (CanteraError err) {
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
    throw CanteraError("MultiPhase::equilibrate_SP",
                       "No convergence for T");
}

int vcs_MultiPhaseEquil::equilibrate(int XY, int estimateEquil,
                                     int printLvl, doublereal err,
                                     int maxsteps, int loglevel)
{
    doublereal xtarget;
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
        throw CanteraError(" vcs_MultiPhaseEquil::equilibrate",
                           "Unsupported Option");
    }
}

int vcs_MultiPhaseEquil::equilibrate_TP(int estimateEquil,
                                        int printLvl, doublereal err,
                                        int maxsteps, int loglevel)
{
    int maxit = maxsteps;
    clockWC tickTock;
    m_printLvl = printLvl;
    m_vprob.m_printLvl = printLvl;

    // Extract the current state information from the MultiPhase object and
    // Transfer it to VCS_PROB object.
    int res = vcs_Cantera_update_vprob(m_mix, &m_vprob);
    if (res != 0) {
        plogf("problems\n");
    }

    // Set the estimation technique
    if (estimateEquil) {
        m_vprob.iest = estimateEquil;
    } else {
        m_vprob.iest = 0;
    }

    // Check obvious bounds on the temperature and pressure NOTE, we may want to
    // do more here with the real bounds given by the ThermoPhase objects.
    if (m_mix->temperature() <= 0.0) {
        throw CanteraError("vcs_MultiPhaseEquil::equilibrate",
                           "Temperature less than zero on input");
    }
    if (m_mix->pressure() <= 0.0) {
        throw CanteraError("vcs_MultiPhaseEquil::equilibrate",
                           "Pressure less than zero on input");
    }

    // Print out the problem specification from the point of
    // view of the vprob object.
    m_vprob.prob_report(m_printLvl);

    //! Call the thermo Program
    int ip1 = m_printLvl;
    int ipr = std::max(0, m_printLvl-1);
    if (m_printLvl >= 3) {
        ip1 = m_printLvl - 2;
    } else {
        ip1 = 0;
    }
    int iSuccess = m_vsolve.vcs(&m_vprob, 0, ipr, ip1, maxit);

    // Transfer the information back to the MultiPhase object. Note we don't
    // just call setMoles, because some multispecies solution phases may be
    // zeroed out, and that would cause a problem for that routine. Also, the
    // mole fractions of such zeroed out phases actually contain information
    // about likely reemergent states.
    m_mix->uploadMoleFractionsFromPhases();
    size_t kGlob = 0;
    for (size_t ip = 0; ip < m_vprob.NPhase; ip++) {
        double phaseMole = 0.0;
        ThermoPhase& tref = m_mix->phase(ip);
        for (size_t k = 0; k < tref.nSpecies(); k++, kGlob++) {
            phaseMole += m_vprob.w[kGlob];
        }
        m_mix->setPhaseMoles(ip, phaseMole);
    }

    double te = tickTock.secondsWC();
    if (printLvl > 0) {
        plogf("\n Results from vcs:\n");
        if (iSuccess != 0) {
            plogf("\nVCS FAILED TO CONVERGE!\n");
        }
        plogf("\n");
        plogf("Temperature = %g Kelvin\n", m_vprob.T);
        plogf("Pressure    = %g Pa\n", m_vprob.PresPA);
        plogf("\n");
        plogf("----------------------------------------"
              "---------------------\n");
        plogf(" Name             Mole_Number(kmol)");
        plogf("  Mole_Fraction     Chem_Potential (J/kmol)\n");
        plogf("--------------------------------------------------"
              "-----------\n");
        for (size_t i = 0; i < m_vprob.nspecies; i++) {
            plogf("%-12s", m_vprob.SpName[i]);
            if (m_vprob.SpeciesUnknownType[i] == VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
                plogf("  %15.3e %15.3e  ", 0.0, m_vprob.mf[i]);
                plogf("%15.3e\n", m_vprob.m_gibbsSpecies[i]);
            } else {
                plogf("  %15.3e   %15.3e  ", m_vprob.w[i], m_vprob.mf[i]);
                if (m_vprob.w[i] <= 0.0) {
                    size_t iph = m_vprob.PhaseID[i];
                    vcs_VolPhase* VPhase = m_vprob.VPhaseList[iph];
                    if (VPhase->nSpecies() > 1) {
                        plogf("     -1.000e+300\n");
                    } else {
                        plogf("%15.3e\n", m_vprob.m_gibbsSpecies[i]);
                    }
                } else {
                    plogf("%15.3e\n", m_vprob.m_gibbsSpecies[i]);
                }
            }
        }
        plogf("------------------------------------------"
              "-------------------\n");
        if (printLvl > 2 && m_vsolve.m_timing_print_lvl > 0) {
            plogf("Total time = %12.6e seconds\n", te);
        }
    }
    return iSuccess;
}

void vcs_MultiPhaseEquil::reportCSV(const std::string& reportFile)
{
    size_t nphase = m_vprob.NPhase;

    FILE* FP = fopen(reportFile.c_str(), "w");
    if (!FP) {
        throw CanteraError("vcs_MultiPhaseEquil::reportCSV",
                           "Failure to open file");
    }
    vector_fp& mf = m_vprob.mf;
    double* fe = &m_vprob.m_gibbsSpecies[0];
    vector_fp VolPM;
    vector_fp activity;
    vector_fp ac;
    vector_fp mu;
    vector_fp mu0;
    vector_fp molalities;

    double vol = 0.0;
    for (size_t iphase = 0; iphase < nphase; iphase++) {
        size_t istart = m_mix->speciesIndex(0, iphase);
        ThermoPhase& tref = m_mix->phase(iphase);
        size_t nSpecies = tref.nSpecies();
        VolPM.resize(nSpecies, 0.0);
        tref.getPartialMolarVolumes(&VolPM[0]);
        vcs_VolPhase* volP = m_vprob.VPhaseList[iphase];

        double TMolesPhase = volP->totalMoles();
        double VolPhaseVolumes = 0.0;
        for (size_t k = 0; k < nSpecies; k++) {
            VolPhaseVolumes += VolPM[k] * mf[istart + k];
        }
        VolPhaseVolumes *= TMolesPhase;
        vol += VolPhaseVolumes;
    }

    fprintf(FP,"--------------------- VCS_MULTIPHASE_EQUIL FINAL REPORT"
            " -----------------------------\n");
    fprintf(FP,"Temperature  = %11.5g kelvin\n", m_mix->temperature());
    fprintf(FP,"Pressure     = %11.5g Pascal\n", m_mix->pressure());
    fprintf(FP,"Total Volume = %11.5g m**3\n", vol);
    fprintf(FP,"Number Basis optimizations = %d\n", m_vprob.m_NumBasisOptimizations);
    fprintf(FP,"Number VCS iterations = %d\n", m_vprob.m_Iterations);

    for (size_t iphase = 0; iphase < nphase; iphase++) {
        size_t istart = m_mix->speciesIndex(0, iphase);
        ThermoPhase& tref = m_mix->phase(iphase);
        string phaseName = tref.name();
        vcs_VolPhase* volP = m_vprob.VPhaseList[iphase];
        double TMolesPhase = volP->totalMoles();
        size_t nSpecies = tref.nSpecies();
        activity.resize(nSpecies, 0.0);
        ac.resize(nSpecies, 0.0);
        mu0.resize(nSpecies, 0.0);
        mu.resize(nSpecies, 0.0);
        VolPM.resize(nSpecies, 0.0);
        molalities.resize(nSpecies, 0.0);
        int actConvention = tref.activityConvention();
        tref.getActivities(&activity[0]);
        tref.getActivityCoefficients(&ac[0]);
        tref.getStandardChemPotentials(&mu0[0]);
        tref.getPartialMolarVolumes(&VolPM[0]);
        tref.getChemPotentials(&mu[0]);
        double VolPhaseVolumes = 0.0;
        for (size_t k = 0; k < nSpecies; k++) {
            VolPhaseVolumes += VolPM[k] * mf[istart + k];
        }
        VolPhaseVolumes *= TMolesPhase;
        vol += VolPhaseVolumes;

        if (actConvention == 1) {
            MolalityVPSSTP* mTP = static_cast<MolalityVPSSTP*>(&tref);
            mTP->getMolalities(&molalities[0]);
            tref.getChemPotentials(&mu[0]);

            if (iphase == 0) {
                fprintf(FP,"        Name,      Phase,  PhaseMoles,  Mole_Fract, "
                        "Molalities,  ActCoeff,   Activity,"
                        "ChemPot_SS0,   ChemPot,   mole_num,       PMVol, Phase_Volume\n");

                fprintf(FP,"            ,           ,      (kmol),            , "
                        "          ,          ,           ,"
                        "   (J/kmol),  (J/kmol),     (kmol), (m**3/kmol),     (m**3)\n");
            }
            for (size_t k = 0; k < nSpecies; k++) {
                std::string sName = tref.speciesName(k);
                fprintf(FP,"%12s, %11s, %11.3e, %11.3e, %11.3e, %11.3e, %11.3e,"
                        "%11.3e, %11.3e, %11.3e, %11.3e, %11.3e\n",
                        sName.c_str(),
                        phaseName.c_str(), TMolesPhase,
                        mf[istart + k], molalities[k], ac[k], activity[k],
                        mu0[k]*1.0E-6, mu[k]*1.0E-6,
                        mf[istart + k] * TMolesPhase,
                        VolPM[k], VolPhaseVolumes);
            }
        } else {
            if (iphase == 0) {
                fprintf(FP,"        Name,       Phase,  PhaseMoles,  Mole_Fract,  "
                        "Molalities,   ActCoeff,    Activity,"
                        "  ChemPotSS0,     ChemPot,   mole_num,       PMVol, Phase_Volume\n");

                fprintf(FP,"            ,            ,      (kmol),            ,  "
                        "          ,           ,            ,"
                        "    (J/kmol),    (J/kmol),     (kmol), (m**3/kmol),       (m**3)\n");
            }
            for (size_t k = 0; k < nSpecies; k++) {
                molalities[k] = 0.0;
            }
            for (size_t k = 0; k < nSpecies; k++) {
                std::string sName = tref.speciesName(k);
                fprintf(FP,"%12s, %11s, %11.3e, %11.3e, %11.3e, %11.3e, %11.3e, "
                        "%11.3e, %11.3e,% 11.3e, %11.3e, %11.3e\n",
                        sName.c_str(),
                        phaseName.c_str(), TMolesPhase,
                        mf[istart + k], molalities[k], ac[k],
                        activity[k], mu0[k]*1.0E-6, mu[k]*1.0E-6,
                        mf[istart + k] * TMolesPhase,
                        VolPM[k], VolPhaseVolumes);
            }
        }

        // Check consistency: These should be equal
        tref.getChemPotentials(fe+istart);
        for (size_t k = 0; k < nSpecies; k++) {
            if (!vcs_doubleEqual(fe[istart+k], mu[k])) {
                fprintf(FP,"ERROR: incompatibility!\n");
                fclose(FP);
                throw CanteraError("vcs_MultiPhaseEquil::reportCSV", "incompatibility!");
            }
        }

    }
    fclose(FP);
}

// HKM -> Work on transferring the current value of the voltages into the
//        equilibrium problem.
int vcs_Cantera_to_vprob(MultiPhase* mphase, VCS_PROB* vprob)
{
    VCS_SPECIES_THERMO* ts_ptr = 0;

    // Calculate the total number of species and phases in the problem
    size_t totNumPhases = mphase->nPhases();
    size_t totNumSpecies = mphase->nSpecies();

    // Problem type has yet to be worked out.
    vprob->prob_type = 0;
    vprob->nspecies = totNumSpecies;
    vprob->ne = 0;
    vprob->NPhase = totNumPhases;
    // Set the initial estimate to a machine generated estimate for now
    // We will work out the details later.
    vprob->iest = -1;
    vprob->T = mphase->temperature();
    vprob->PresPA = mphase->pressure();
    vprob->Vol = mphase->volume();
    vprob->Title = "MultiPhase Object";

    int printLvl = vprob->m_printLvl;

    // Loop over the phases, transferring pertinent information
    int kT = 0;
    for (size_t iphase = 0; iphase < totNumPhases; iphase++) {
        // Get the ThermoPhase object - assume volume phase
        ThermoPhase* tPhase = &mphase->phase(iphase);
        size_t nelem = tPhase->nElements();

        // Query Cantera for the equation of state type of the current phase.
        std::string eos = tPhase->type();
        bool gasPhase = (eos == "IdealGas");

        // Find out the number of species in the phase
        size_t nSpPhase = tPhase->nSpecies();
        // Find out the name of the phase
        string phaseName = tPhase->name();

        // Call the basic vcs_VolPhase creation routine.
        // Properties set here:
        //    ->PhaseNum = phase number in the thermo problem
        //    ->GasPhase = Boolean indicating whether it is a gas phase
        //    ->NumSpecies = number of species in the phase
        //    ->TMolesInert = Inerts in the phase = 0.0 for cantera
        //    ->PhaseName  = Name of the phase
        vcs_VolPhase* VolPhase = vprob->VPhaseList[iphase];
        VolPhase->resize(iphase, nSpPhase, nelem, phaseName.c_str(), 0.0);
        VolPhase->m_gasPhase = gasPhase;

        // Tell the vcs_VolPhase pointer about cantera
        VolPhase->setPtrThermoPhase(tPhase);
        VolPhase->setTotalMoles(0.0);

        // Set the electric potential of the volume phase from the
        // ThermoPhase object's value.
        VolPhase->setElectricPotential(tPhase->electricPotential());

        // Query the ThermoPhase object to find out what convention
        // it uses for the specification of activity and Standard State.
        VolPhase->p_activityConvention = tPhase->activityConvention();

        // Assign the value of eqn of state. Handle conflicts here.
        if (eos == "IdealGas") {
            VolPhase->m_eqnState = VCS_EOS_IDEAL_GAS;
        } else if (eos == "ConstDensity") {
            VolPhase->m_eqnState = VCS_EOS_CONSTANT;
        } else if (eos == "StoichSubstance") {
            VolPhase->m_eqnState = VCS_EOS_STOICH_SUB;
        } else if (eos == "IdealSolidSoln") {
            VolPhase->m_eqnState = VCS_EOS_IDEAL_SOLN;
        } else if (eos == "Surf" || eos == "Edge") {
            throw CanteraError("VCSnonideal",
                               "Surface/edge phase not handled yet.");
        } else {
            if (printLvl > 1) {
                writelog("Unknown Cantera EOS to VCSnonideal: '{}'\n", eos);
            }
            VolPhase->m_eqnState = VCS_EOS_UNK_CANTERA;
        }

        // Transfer all of the element information from the ThermoPhase object
        // to the vcs_VolPhase object. Also decide whether we need a new charge
        // neutrality element in the phase to enforce a charge neutrality
        // constraint. We also decide whether this is a single species phase
        // with the voltage being the independent variable setting the chemical
        // potential of the electrons.
        VolPhase->transferElementsFM(tPhase);

        // Combine the element information in the vcs_VolPhase
        // object into the vprob object.
        vprob->addPhaseElements(VolPhase);
        VolPhase->setState_TP(vprob->T, vprob->PresPA);
        vector_fp muPhase(tPhase->nSpecies(),0.0);
        tPhase->getChemPotentials(&muPhase[0]);
        double tMoles = 0.0;

        // Loop through each species in the current phase
        for (size_t k = 0; k < nSpPhase; k++) {
            // Obtain the molecular weight of the species from the
            // ThermoPhase object
            vprob->WtSpecies[kT] = tPhase->molecularWeight(k);

            // Obtain the charges of the species from the ThermoPhase object
            vprob->Charge[kT] = tPhase->charge(k);

            // Set the phaseid of the species
            vprob->PhaseID[kT] = iphase;

            // Transfer the Species name
            string stmp = mphase->speciesName(kT);
            vprob->SpName[kT] = stmp;

            // Transfer the type of unknown
            vprob->SpeciesUnknownType[kT] = VolPhase->speciesUnknownType(k);
            if (vprob->SpeciesUnknownType[kT] == VCS_SPECIES_TYPE_MOLNUM) {
                // Set the initial number of kmoles of the species
                // and the mole fraction vector
                vprob->w[kT] = mphase->speciesMoles(kT);
                tMoles += vprob->w[kT];
                vprob->mf[kT] = mphase->moleFraction(kT);
             } else if (vprob->SpeciesUnknownType[kT] == VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
                vprob->w[kT] = tPhase->electricPotential();
                vprob->mf[kT] = mphase->moleFraction(kT);
             } else {
               throw CanteraError(" vcs_Cantera_to_vprob() ERROR",
                                  "Unknown species type: {}", vprob->SpeciesUnknownType[kT]);
             }

            // transfer chemical potential vector
            vprob->m_gibbsSpecies[kT] = muPhase[k];

            // Transfer the species information from the
            // volPhase structure to the VPROB structure
            // This includes:
            //      FormulaMatrix[][]
            //      VolPhase->IndSpecies[]
            vprob->addOnePhaseSpecies(VolPhase, k, kT);

            // Get a pointer to the thermo object
            ts_ptr = vprob->SpeciesThermo[kT];

            // Fill in the vcs_SpeciesProperty structure
            vcs_SpeciesProperties* sProp = VolPhase->speciesProperty(k);
            sProp->NumElements = vprob->ne;
            sProp->SpName = vprob->SpName[kT];
            sProp->SpeciesThermo = ts_ptr;
            sProp->WtSpecies = tPhase->molecularWeight(k);
            sProp->FormulaMatrixCol.resize(vprob->ne, 0.0);
            for (size_t e = 0; e < vprob->ne; e++) {
                sProp->FormulaMatrixCol[e] = vprob->FormulaMatrix(kT,e);
            }
            sProp->Charge = tPhase->charge(k);
            sProp->SurfaceSpecies = false;
            sProp->VolPM = 0.0;

            // Transfer the thermo specification of the species
            //              vprob->SpeciesThermo[]

            // Add lookback connectivity into the thermo object first
            ts_ptr->IndexPhase = iphase;
            ts_ptr->IndexSpeciesPhase = k;
            ts_ptr->OwningPhase = VolPhase;

            // get a reference to the Cantera species thermo.
            MultiSpeciesThermo& sp = tPhase->speciesThermo();

            int spType = sp.reportType(k);
            if (spType == SIMPLE) {
                double c[4];
                double minTemp, maxTemp, refPressure;
                sp.reportParams(k, spType, c, minTemp, maxTemp, refPressure);
                ts_ptr->SS0_Model = VCS_SS0_CONSTANT;
                ts_ptr->SS0_T0 = c[0];
                ts_ptr->SS0_H0 = c[1];
                ts_ptr->SS0_S0 = c[2];
                ts_ptr->SS0_Cp0 = c[3];
                if (gasPhase) {
                    ts_ptr->SSStar_Model = VCS_SSSTAR_IDEAL_GAS;
                    ts_ptr->SSStar_Vol_Model = VCS_SSVOL_IDEALGAS;
                } else {
                    ts_ptr->SSStar_Model = VCS_SSSTAR_CONSTANT;
                    ts_ptr->SSStar_Vol_Model = VCS_SSVOL_CONSTANT;
                }
            } else {
                if (vprob->m_printLvl > 2) {
                    plogf("vcs_Cantera_convert: Species Type %d not known \n",
                          spType);
                }
                ts_ptr->SS0_Model = VCS_SS0_NOTHANDLED;
                ts_ptr->SSStar_Model = VCS_SSSTAR_NOTHANDLED;
            }

            // Transfer the Volume Information -> NEEDS WORK
            if (gasPhase) {
                ts_ptr->SSStar_Vol_Model = VCS_SSVOL_IDEALGAS;
                ts_ptr->SSStar_Vol0 = 82.05 * 273.15 / 1.0;
            } else {
                vector_fp phaseTermCoeff(nSpPhase, 0.0);
                int nCoeff;
                tPhase->getParameters(nCoeff, &phaseTermCoeff[0]);
                ts_ptr->SSStar_Vol_Model = VCS_SSVOL_CONSTANT;
                ts_ptr->SSStar_Vol0 = phaseTermCoeff[k];
            }
            kT++;
        }

        // Now go back through the species in the phase and assign a valid mole
        // fraction to all phases, even if the initial estimate of the total
        // number of moles is zero.
        if (tMoles > 0.0) {
            for (size_t k = 0; k < nSpPhase; k++) {
                size_t kTa = VolPhase->spGlobalIndexVCS(k);
                vprob->mf[kTa] = vprob->w[kTa] / tMoles;
            }
        } else {
            // Perhaps, we could do a more sophisticated treatment below.
            // But, will start with this.
            for (size_t k = 0; k < nSpPhase; k++) {
                size_t kTa = VolPhase->spGlobalIndexVCS(k);
                vprob->mf[kTa]= 1.0 / (double) nSpPhase;
            }
        }

        VolPhase->setMolesFromVCS(VCS_STATECALC_OLD, &vprob->w[0]);

        // Now, calculate a sample naught Gibbs free energy calculation
        // at the specified temperature.
        for (size_t k = 0; k < nSpPhase; k++) {
            vcs_SpeciesProperties* sProp = VolPhase->speciesProperty(k);
            ts_ptr = sProp->SpeciesThermo;
            ts_ptr->SS0_feSave = VolPhase->G0_calc_one(k)/ GasConstant;
            ts_ptr->SS0_TSave = vprob->T;
        }
    }

    // Transfer initial element abundances to the vprob object.
    // We have to find the mapping index from one to the other
    vprob->gai.resize(vprob->ne, 0.0);
    vprob->set_gai();

    // Printout the species information: PhaseID's and mole nums
    if (vprob->m_printLvl > 1) {
        writeline('=', 80, true, true);
        writeline('=', 16, false);
        plogf(" Cantera_to_vprob: START OF PROBLEM STATEMENT ");
        writeline('=', 20);
        writeline('=', 80);
        plogf("             Phase IDs of species\n");
        plogf("            species     phaseID        phaseName   ");
        plogf(" Initial_Estimated_kMols\n");
        for (size_t i = 0; i < vprob->nspecies; i++) {
            size_t iphase = vprob->PhaseID[i];

            vcs_VolPhase* VolPhase = vprob->VPhaseList[iphase];
            plogf("%16s      %5d   %16s", vprob->SpName[i].c_str(), iphase,
                  VolPhase->PhaseName.c_str());
            if (vprob->SpeciesUnknownType[i] == VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
                plogf("     Volts = %-10.5g\n", vprob->w[i]);
            } else {
                plogf("             %-10.5g\n", vprob->w[i]);
            }
        }

        // Printout of the Phase structure information
        writeline('-', 80, true, true);
        plogf("             Information about phases\n");
        plogf("  PhaseName    PhaseNum SingSpec GasPhase EqnState NumSpec");
        plogf("  TMolesInert       Tmoles(kmol)\n");

        for (size_t iphase = 0; iphase < vprob->NPhase; iphase++) {
            vcs_VolPhase* VolPhase = vprob->VPhaseList[iphase];
            plogf("%16s %5d %5d %8d %16s %8d %16e ", VolPhase->PhaseName.c_str(),
                  VolPhase->VP_ID_, VolPhase->m_singleSpecies,
                  VolPhase->m_gasPhase, VolPhase->eos_name(),
                  VolPhase->nSpecies(), VolPhase->totalMolesInert());
            plogf("%16e\n", VolPhase->totalMoles());
        }

        writeline('=', 80, true, true);
        writeline('=', 16, false);
        plogf(" Cantera_to_vprob: END OF PROBLEM STATEMENT ");
        writeline('=', 20);
        writeline('=', 80);
        plogf("\n");
    }
    return VCS_SUCCESS;
}

int vcs_Cantera_update_vprob(MultiPhase* mphase, VCS_PROB* vprob)
{
    size_t totNumPhases = mphase->nPhases();
    size_t kT = 0;
    vector_fp tmpMoles;
    // Problem type has yet to be worked out.
    vprob->prob_type = 0;
    // Whether we have an estimate or not gets overwritten on
    // the call to the equilibrium solver.
    vprob->iest = -1;
    vprob->T = mphase->temperature();
    vprob->PresPA = mphase->pressure();
    vprob->Vol = mphase->volume();

    for (size_t iphase = 0; iphase < totNumPhases; iphase++) {
        ThermoPhase* tPhase = &mphase->phase(iphase);
        vcs_VolPhase* volPhase = vprob->VPhaseList[iphase];

        // Set the electric potential of the volume phase from the
        // ThermoPhase object's value.
        volPhase->setElectricPotential(tPhase->electricPotential());

        volPhase->setState_TP(vprob->T, vprob->PresPA);
        vector_fp muPhase(tPhase->nSpecies(),0.0);
        tPhase->getChemPotentials(&muPhase[0]);

        // Loop through each species in the current phase
        size_t nSpPhase = tPhase->nSpecies();
        tmpMoles.resize(nSpPhase);
        for (size_t k = 0; k < nSpPhase; k++) {
            tmpMoles[k] = mphase->speciesMoles(kT);
            vprob->w[kT] = mphase->speciesMoles(kT);
            vprob->mf[kT] = mphase->moleFraction(kT);

            // transfer chemical potential vector
            vprob->m_gibbsSpecies[kT] = muPhase[k];

            kT++;
        }
        if (volPhase->phiVarIndex() != npos) {
            size_t kphi = volPhase->phiVarIndex();
            size_t kglob = volPhase->spGlobalIndexVCS(kphi);
            vprob->w[kglob] = tPhase->electricPotential();
        }
        volPhase->setMolesFromVCS(VCS_STATECALC_OLD, &vprob->w[0]);
        if ((nSpPhase == 1) && (volPhase->phiVarIndex() == 0)) {
            volPhase->setExistence(VCS_PHASE_EXIST_ALWAYS);
        } else if (volPhase->totalMoles() > 0.0) {
            volPhase->setExistence(VCS_PHASE_EXIST_YES);
        } else {
            volPhase->setExistence(VCS_PHASE_EXIST_NO);
        }
    }

    // Transfer initial element abundances to the vprob object. Put them in the
    // front of the object. There may be more constraints than there are
    // elements. But, we know the element abundances are in the front of the
    // vector.
    vprob->set_gai();

    // Printout the species information: PhaseID's and mole nums
    if (vprob->m_printLvl > 1) {
        writeline('=', 80, true, true);
        writeline('=', 20, false);
        plogf(" Cantera_to_vprob: START OF PROBLEM STATEMENT ");
        writeline('=', 20);
        writeline('=', 80);
        plogf("\n");
        plogf("             Phase IDs of species\n");
        plogf("            species     phaseID        phaseName   ");
        plogf(" Initial_Estimated_kMols\n");
        for (size_t i = 0; i < vprob->nspecies; i++) {
            size_t iphase = vprob->PhaseID[i];
            vcs_VolPhase* VolPhase = vprob->VPhaseList[iphase];
            plogf("%16s      %5d   %16s", vprob->SpName[i].c_str(), iphase,
                  VolPhase->PhaseName.c_str());
            if (vprob->SpeciesUnknownType[i] == VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
                plogf("     Volts = %-10.5g\n", vprob->w[i]);
            } else {
                plogf("             %-10.5g\n", vprob->w[i]);
            }
        }

        // Printout of the Phase structure information
        writeline('-', 80, true, true);
        plogf("             Information about phases\n");
        plogf("  PhaseName    PhaseNum SingSpec GasPhase EqnState NumSpec");
        plogf("  TMolesInert       Tmoles(kmol)\n");

        for (size_t iphase = 0; iphase < vprob->NPhase; iphase++) {
            vcs_VolPhase* VolPhase = vprob->VPhaseList[iphase];
            plogf("%16s %5d %5d %8d %16s %8d %16e ", VolPhase->PhaseName.c_str(),
                  VolPhase->VP_ID_, VolPhase->m_singleSpecies,
                  VolPhase->m_gasPhase, VolPhase->eos_name(),
                  VolPhase->nSpecies(), VolPhase->totalMolesInert());
            plogf("%16e\n", VolPhase->totalMoles());
        }

        writeline('=', 80, true, true);
        writeline('=', 20, false);
        plogf(" Cantera_to_vprob: END OF PROBLEM STATEMENT ");
        writeline('=', 20);
        writeline('=', 80);
        plogf("\n");
    }
    return VCS_SUCCESS;
}

void vcs_MultiPhaseEquil::getStoichVector(size_t rxn, vector_fp& nu)
{
    warn_deprecated("vcs_MultiPhaseEquil::getStoichVector",
                    "Unused. To be removed after Cantera 2.3.");
    size_t nsp = m_vsolve.m_numSpeciesTot;
    nu.resize(nsp, 0.0);
    for (size_t i = 0; i < nsp; i++) {
        nu[i] = 0.0;
    }
    size_t nc = numComponents();
    const std::vector<size_t>& indSpecies = m_vsolve.m_speciesMapIndex;
    if (rxn > nsp - nc) {
        return;
    }
    size_t j = indSpecies[rxn + nc];
    nu[j] = 1.0;
    for (size_t kc = 0; kc < nc; kc++) {
        j = indSpecies[kc];
        nu[j] = m_vsolve.m_stoichCoeffRxnMatrix(kc,rxn);
    }
}

size_t vcs_MultiPhaseEquil::numComponents() const
{
    warn_deprecated("vcs_MultiPhaseEquil::numComponents",
                    "Unused. To be removed after Cantera 2.3.");
    return m_vsolve.m_numComponents;
}

size_t vcs_MultiPhaseEquil::numElemConstraints() const
{
    warn_deprecated("vcs_MultiPhaseEquil::numElemConstraints",
                    "Unused. To be removed after Cantera 2.3.");
    return m_vsolve.m_numElemConstraints;
}

size_t vcs_MultiPhaseEquil::component(size_t m) const
{
    warn_deprecated("vcs_MultiPhaseEquil::component",
                    "Unused. To be removed after Cantera 2.3.");
    size_t nc = numComponents();
    if (m < nc) {
        return m_vsolve.m_speciesMapIndex[m];
    } else {
        return npos;
    }
}

int vcs_MultiPhaseEquil::determine_PhaseStability(int iph, double& funcStab, int printLvl, int loglevel)
{
    warn_deprecated("vcs_MultiPhaseEquil::determine_PhaseStability",
                    "Broken and unused. To be removed after Cantera 2.3.");
    clockWC tickTock;
    m_printLvl = printLvl;
    m_vprob.m_printLvl = printLvl;

    // Extract the current state information from the MultiPhase object and
    // Transfer it to VCS_PROB object.
    int res = vcs_Cantera_update_vprob(m_mix, &m_vprob);
    if (res != 0) {
        plogf("problems\n");
    }

    // Check obvious bounds on the temperature and pressure
    // NOTE, we may want to do more here with the real bounds
    // given by the ThermoPhase objects.
    double T = m_mix->temperature();
    if (T <= 0.0) {
        throw CanteraError("vcs_MultiPhaseEquil::determine_PhaseStability",
                           "Temperature less than zero on input");
    }
    double pres = m_mix->pressure();
    if (pres <= 0.0) {
        throw CanteraError("vcs_MultiPhaseEquil::determine_PhaseStability",
                           "Pressure less than zero on input");
    }

    // Print out the problem specification from the point of
    // view of the vprob object.
    m_vprob.prob_report(m_printLvl);

    // Call the thermo Program
    int iStable = m_vsolve.vcs_PS(&m_vprob, iph, printLvl, funcStab);

    // Transfer the information back to the MultiPhase object. Note we don't
    // just call setMoles, because some multispecies solution phases may be
    // zeroed out, and that would cause a problem for that routine. Also, the
    // mole fractions of such zeroed out phases actually contain information
    // about likely reemergent states.
    m_mix->uploadMoleFractionsFromPhases();
    m_mix->getChemPotentials(m_vprob.m_gibbsSpecies.data());

    double te = tickTock.secondsWC();
    if (printLvl > 0) {
        plogf("\n Results from vcs_PS:\n");
        plogf("\n");
        plogf("Temperature = %g Kelvin\n", m_vprob.T);
        plogf("Pressure    = %g Pa\n", m_vprob.PresPA);
        std::string sss = m_mix->phaseName(iph);
        if (iStable) {
            plogf("Phase %d named %s is     stable, function value = %g > 0\n", iph, sss.c_str(), funcStab);
        } else {
            plogf("Phase %d named %s is not stable + function value = %g < 0\n", iph, sss.c_str(), funcStab);
        }
        plogf("\n");
        plogf("----------------------------------------"
              "---------------------\n");
        plogf(" Name             Mole_Number(kmol)");
        plogf("  Mole_Fraction     Chem_Potential (J/kmol)\n");
        plogf("-------------------------------------------------------------\n");
        for (size_t i = 0; i < m_vprob.nspecies; i++) {
            plogf("%-12s", m_vprob.SpName[i]);
            if (m_vprob.SpeciesUnknownType[i] == VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
                plogf("  %15.3e %15.3e  ", 0.0, m_vprob.mf[i]);
                plogf("%15.3e\n", m_vprob.m_gibbsSpecies[i]);
            } else {
                plogf("  %15.3e   %15.3e  ", m_vprob.w[i], m_vprob.mf[i]);
                if (m_vprob.w[i] <= 0.0) {
                    plogf("%15.3e\n", m_vprob.m_gibbsSpecies[i]);
                } else {
                    plogf("%15.3e\n", m_vprob.m_gibbsSpecies[i]);
                }
            }
        }
        plogf("------------------------------------------"
              "-------------------\n");
        if (printLvl > 2 && m_vsolve.m_timing_print_lvl > 0) {
            plogf("Total time = %12.6e seconds\n", te);
        }
    }
    return iStable;
}

}
