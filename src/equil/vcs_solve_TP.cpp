/**
 * @file vcs_solve_TP.cpp Implementation file that contains the
 *     main algorithm for finding an equilibrium
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/equil/vcs_solve.h"
#include "cantera/equil/vcs_species_thermo.h"
#include "cantera/equil/vcs_VolPhase.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/clockWC.h"
#include "cantera/base/stringUtils.h"
#include "cantera/numerics/DenseMatrix.h"

#include <cstdio>

using namespace std;
namespace {
enum stages {MAIN, EQUILIB_CHECK, ELEM_ABUND_CHECK,
             RECHECK_DELETED, RETURN_A, RETURN_B};
}

namespace Cantera
{

void VCS_SOLVE::checkDelta1(double* const dsLocal,
                            double* const delTPhMoles, size_t kspec)
{
    vector_fp dchange(m_numPhases, 0.0);
    for (size_t k = 0; k < kspec; k++) {
        if (m_speciesUnknownType[k] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
            size_t iph = m_phaseID[k];
            dchange[iph] += dsLocal[k];
        }
    }
    for (size_t iphase = 0; iphase < m_numPhases; iphase++) {
        double denom = max(m_totalMolNum, 1.0E-4);
        if (!vcs_doubleEqual(dchange[iphase]/denom, delTPhMoles[iphase]/denom)) {
            throw CanteraError("VCS_SOLVE::checkDelta1",
                               "we have found a problem");
        }
    }
}

int VCS_SOLVE::vcs_solve_TP(int print_lvl, int printDetails, int maxit)
{
    int stage = MAIN;
    bool allMinorZeroedSpecies = false;
    size_t it1 = 0;
    size_t iti;
    int rangeErrorFound = 0;
    bool giveUpOnElemAbund = false;
    int finalElemAbundAttempts = 0;
    bool uptodate_minors = true;
    int forceComponentCalc = 1;

    char ANOTE[128];
    // Set the debug print lvl to the same as the print lvl.
    m_debug_print_lvl = printDetails;
    if (printDetails > 0 && print_lvl == 0) {
        print_lvl = 1;
    }
    // Initialize and set up all counters
    vcs_counters_init(0);
    clockWC ticktock;

    // temporary space for usage in this routine and in subroutines
    m_sm.assign(m_nelem * m_nelem, 0.0);
    m_ss.assign(m_nelem, 0.0);
    m_sa.assign(m_nelem, 0.0);
    m_aw.assign(m_nsp, 0.0);
    m_wx.assign(m_nelem, 0.0);

    int solveFail = false;

    // Evaluate the elemental composition
    vcs_elab();

    // Printout the initial conditions for problem
    if (print_lvl != 0) {
        plogf("VCS CALCULATION METHOD\n\n ");
        plogf("MultiPhase Object\n");
        plogf("\n\n%5d SPECIES\n%5d ELEMENTS\n", m_nsp, m_nelem);
        plogf("%5d COMPONENTS\n", m_numComponents);
        plogf("%5d PHASES\n", m_numPhases);
        plogf(" PRESSURE%22.8g %3s\n", m_pressurePA, "Pa ");
        plogf(" TEMPERATURE%19.3f K\n", m_temperature);
        vcs_VolPhase* Vphase = m_VolPhaseList[0].get();
        if (Vphase->nSpecies() > 0) {
            plogf(" PHASE1 INERTS%17.3f\n", TPhInertMoles[0]);
        }
        if (m_numPhases > 1) {
            plogf(" PHASE2 INERTS%17.3f\n", TPhInertMoles[1]);
        }
        plogf("\n ELEMENTAL ABUNDANCES             CORRECT");
        plogf("          FROM ESTIMATE           Type\n\n");
        for (size_t i = 0; i < m_nelem; ++i) {
            writeline(' ', 26, false);
            plogf("%-2.2s", m_elementName[i]);
            plogf("%20.12E%20.12E     %3d\n", m_elemAbundancesGoal[i], m_elemAbundances[i],
                  m_elType[i]);
        }
        if (m_doEstimateEquil < 0) {
            plogf("\n MODIFIED LINEAR PROGRAMMING ESTIMATE OF EQUILIBRIUM - forced\n");
        } else if (m_doEstimateEquil > 0) {
            plogf("\n MODIFIED LINEAR PROGRAMMING ESTIMATE OF EQUILIBRIUM - where necessary\n");
        }
        if (m_doEstimateEquil == 0) {
            plogf("\n USER ESTIMATE OF EQUILIBRIUM\n");
        }
        plogf(" Stan. Chem. Pot. in J/kmol\n");
        plogf("\n SPECIES            FORMULA VECTOR   ");
        writeline(' ', 41, false);
        plogf("   STAN_CHEM_POT  EQUILIBRIUM_EST.  Species_Type\n\n");
        writeline(' ', 20, false);
        for (size_t i = 0; i < m_nelem; ++i) {
            plogf("%-4.4s    ", m_elementName[i]);
        }
        plogf("   PhaseID\n");
        double RT = GasConstant * m_temperature;
        for (size_t i = 0; i < m_nsp; ++i) {
            plogf(" %-18.18s", m_speciesName[i]);
            for (size_t j = 0; j < m_nelem; ++j) {
                plogf("% -7.3g ", m_formulaMatrix(i,j));
            }
            plogf("  %3d  ", m_phaseID[i]);
            writeline(' ', std::max(55-int(m_nelem)*8, 0), false);
            plogf("%12.5E  %12.5E", RT * m_SSfeSpecies[i], m_molNumSpecies_old[i]);
            if (m_speciesUnknownType[i] == VCS_SPECIES_TYPE_MOLNUM) {
                plogf("       Mol_Num\n");
            } else if (m_speciesUnknownType[i] == VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
                plogf("       Voltage\n");
            } else {
                plogf("       Unknown\n");
            }
        }
    }

    for (size_t i = 0; i < m_nsp; ++i) {
        if (m_molNumSpecies_old[i] < 0.0) {
            plogf("On Input species %-12s has a negative MF, setting it small\n",
                  m_speciesName[i]);
            size_t iph = m_phaseID[i];
            double tmp = m_tPhaseMoles_old[iph] * VCS_RELDELETE_SPECIES_CUTOFF * 10;
            tmp = std::max(tmp, VCS_DELETE_MINORSPECIES_CUTOFF*10.);
            m_molNumSpecies_old[i] = tmp;
        }
    }

    // Evaluate the total moles of species in the problem
    vcs_tmoles();

    // Evaluate all chemical potentials at the old mole numbers at the outset of
    // the calculation.
    vcs_setFlagsVolPhases(false, VCS_STATECALC_OLD);
    vcs_dfe(VCS_STATECALC_OLD, 0, 0, m_numSpeciesRdc);

    bool lec;
    while (true) {
        if (stage == MAIN) {
            // DETERMINE BASIS SPECIES, EVALUATE STOICHIOMETRY
            if (forceComponentCalc) {
                int retn = solve_tp_component_calc(allMinorZeroedSpecies);
                if (retn != VCS_SUCCESS) {
                    return retn;
                }
                it1 = 1;
                forceComponentCalc = 0;
                iti = 0;
            }
            // Check on too many iterations. If we have too many iterations,
            // Clean up and exit code even though we haven't converged.
            //     -> we have run out of iterations!
            if (m_VCount->Its > maxit) {
                return -1;
            }
            solve_tp_inner(iti, it1, uptodate_minors, allMinorZeroedSpecies,
                           forceComponentCalc, stage, printDetails > 0, ANOTE);
            lec = false;
        } else if (stage == EQUILIB_CHECK) {
            // EQUILIBRIUM CHECK FOR MAJOR SPECIES
            solve_tp_equilib_check(allMinorZeroedSpecies, uptodate_minors,
                                   giveUpOnElemAbund, solveFail, iti, it1,
                                   maxit, stage, lec);
        } else if (stage == ELEM_ABUND_CHECK) {
            // CORRECT ELEMENTAL ABUNDANCES
            solve_tp_elem_abund_check(iti, stage, lec, giveUpOnElemAbund,
                                      finalElemAbundAttempts, rangeErrorFound);
        } else if (stage == RECHECK_DELETED) {
            // RECHECK DELETED SPECIES
            //
            // We are here for two reasons. One is if we have achieved
            // convergence, but some species have been eliminated from the
            // problem because they were in multispecies phases and their mole
            // fractions drifted less than VCS_DELETE_SPECIES_CUTOFF. The other
            // reason why we are here is because all of the non-component
            // species in the problem have been eliminated for one reason or
            // another.
            size_t npb = vcs_recheck_deleted();

            // If we haven't found any species that needed adding we are done.
            if (npb <= 0) {
                stage = RETURN_B;
            } else {
                // If we have found something to add, recalculate everything for
                // minor species and go back to do a full iteration
                vcs_setFlagsVolPhases(false, VCS_STATECALC_OLD);
                vcs_dfe(VCS_STATECALC_OLD, 1, 0, m_numSpeciesRdc);
                vcs_deltag(0, false, VCS_STATECALC_OLD);
                iti = 0;
                stage = MAIN;
            }
        } else if (stage == RETURN_A) {
            // CLEANUP AND RETURN BLOCK
            size_t npb = vcs_recheck_deleted();

            // If we haven't found any species that needed adding we are done.
            if (npb > 0) {
                // If we have found something to add, recalculate everything for
                // minor species and go back to do a full iteration
                vcs_setFlagsVolPhases(false, VCS_STATECALC_OLD);
                vcs_dfe(VCS_STATECALC_OLD, 1, 0, m_numSpeciesRdc);
                vcs_deltag(0, false, VCS_STATECALC_OLD);
                iti = 0;
                stage = MAIN;
            } else {
                stage = RETURN_B;
            }
        } else if (stage == RETURN_B) {
            // Add back deleted species in non-zeroed phases. Estimate their
            // mole numbers.
            size_t npb = vcs_add_all_deleted();
            if (npb > 0) {
                iti = 0;
                if (m_debug_print_lvl >= 1) {
                    plogf("  --- add_all_deleted(): some rxns not converged. RETURNING TO LOOP!\n");
                }
                stage = MAIN;
            } else {
                break;
            }
        }
    }

    // Make sure the volume phase objects hold the same state and information as
    // the vcs object. This also update the Cantera objects with this
    // information.
    vcs_updateVP(VCS_STATECALC_OLD);

    // Evaluate the final mole fractions storing them in wt[]
    m_molNumSpecies_new.assign(m_molNumSpecies_new.size(), 0.0);
    for (size_t kspec = 0; kspec < m_nsp; ++kspec) {
        if (m_SSPhase[kspec]) {
            m_molNumSpecies_new[kspec] = 1.0;
        } else {
            size_t iph = m_phaseID[kspec];
            if (m_tPhaseMoles_old[iph] != 0.0) {
                m_molNumSpecies_new[kspec] = m_molNumSpecies_old[kspec] / m_tPhaseMoles_old[iph];
            } else {
                // For MultiSpecies phases that are zeroed out, return the mole
                // fraction vector from the VolPhase object. This contains the
                // mole fraction that would be true if the phase just pops into
                // existence.
                size_t i = m_speciesLocalPhaseIndex[kspec];
                m_molNumSpecies_new[kspec] = m_VolPhaseList[iph]->molefraction(i);
            }
        }
    }
    // Return an error code if a Range Space Error is thought to have occurred.
    if (rangeErrorFound) {
        solveFail = 1;
    }

    // Calculate counters
    double tsecond = ticktock.secondsWC();
    m_VCount->Time_vcs_TP = tsecond;
    m_VCount->T_Time_vcs_TP += m_VCount->Time_vcs_TP;
    m_VCount->T_Calls_vcs_TP++;
    m_VCount->T_Its += m_VCount->Its;
    m_VCount->T_Basis_Opts += m_VCount->Basis_Opts;
    m_VCount->T_Time_basopt += m_VCount->Time_basopt;

    // Return a Flag indicating whether convergence occurred
    return solveFail;
}

int VCS_SOLVE::solve_tp_component_calc(bool& allMinorZeroedSpecies)
{
    double test = -1.0e-10;
    bool usedZeroedSpecies;
    int retn = vcs_basopt(false, &m_aw[0], &m_sa[0], &m_sm[0], &m_ss[0],
                          test, &usedZeroedSpecies);
    if (retn != VCS_SUCCESS) {
        return retn;
    }

    // Update the phase objects with the contents of the soln vector
    vcs_updateVP(VCS_STATECALC_OLD);
    vcs_deltag(0, false, VCS_STATECALC_OLD);

    // EVALUATE INITIAL SPECIES STATUS VECTOR
    allMinorZeroedSpecies = vcs_evaluate_speciesType();

    // EVALUATE THE ELELEMT ABUNDANCE CHECK
    if (! vcs_elabcheck(0)) {
        debuglog("   --- Element Abundance check failed\n", m_debug_print_lvl >= 2);
        vcs_elcorr(&m_sm[0], &m_wx[0]);
        vcs_setFlagsVolPhases(false, VCS_STATECALC_OLD);
        vcs_dfe(VCS_STATECALC_OLD, 0, 0, m_numSpeciesRdc);
        // Update the phase objects with the contents of the soln vector
        vcs_updateVP(VCS_STATECALC_OLD);
        vcs_deltag(0, false, VCS_STATECALC_OLD);
    } else {
        debuglog("   --- Element Abundance check passed\n", m_debug_print_lvl >= 2);
    }
    return retn;
}

void VCS_SOLVE::solve_tp_inner(size_t& iti, size_t& it1,
                               bool& uptodate_minors,
                               bool& allMinorZeroedSpecies,
                               int& forceComponentCalc,
                               int& stage, bool printDetails, char* ANOTE)
{
    if (iti == 0) {
        // SET INITIAL VALUES FOR ITERATION
        // EVALUATE REACTION ADJUSTMENTS
        //
        // Evaluate the minor non-component species chemical potentials and
        // delta G for their formation reactions We have already evaluated the
        // major non-components
        if (!uptodate_minors) {
            vcs_setFlagsVolPhases(false, VCS_STATECALC_OLD);
            vcs_dfe(VCS_STATECALC_OLD, 1, 0, m_numSpeciesRdc);
            vcs_deltag(1, false, VCS_STATECALC_OLD);
        }
        uptodate_minors = true;
    } else {
        uptodate_minors = false;
    }

    if (printDetails) {
        plogf("\n");
        writeline('=', 110);
        plogf(" Iteration = %3d, Iterations since last evaluation of "
              "optimal basis = %3d",
              m_VCount->Its, it1 - 1);
        if (iti == 0) {
            plogf(" (all species)\n");
        } else {
            plogf(" (only major species)\n");
        }
    }

     // Calculate the total moles in each phase -> old solution
     //  -> Needed for numerical stability when phases disappear.
     //  -> the phase moles tend to drift off without this step.
    check_tmoles();
    vcs_tmoles();
    // COPY OLD into NEW and ZERO VECTORS
    // Copy the old solution into the new solution as an initial guess
    m_feSpecies_new = m_feSpecies_old;
    m_actCoeffSpecies_new = m_actCoeffSpecies_old;
    m_deltaGRxn_new = m_deltaGRxn_old;
    m_deltaGRxn_Deficient = m_deltaGRxn_old;
    m_tPhaseMoles_new = m_tPhaseMoles_old;

    // Zero out the entire vector of updates. We sometimes would query these
    // values below, and we want to be sure that no information is left from
    // previous iterations.
    m_deltaMolNumSpecies.assign(m_deltaMolNumSpecies.size(), 0.0);

    // DETERMINE IF DEAD PHASES POP INTO EXISTENCE
    //
    // First step is a major branch in the algorithm. We first determine if a
    // phase pops into existence.
    std::vector<size_t> phasePopPhaseIDs(0);
    size_t iphasePop = vcs_popPhaseID(phasePopPhaseIDs);
    if (iphasePop != npos) {
        int soldel = vcs_popPhaseRxnStepSizes(iphasePop);
        if (soldel == 3) {
            iphasePop = npos;
            if (m_debug_print_lvl >= 2) {
                plogf("   --- vcs_popPhaseRxnStepSizes() was called but stoich "
                      "prevented phase %d popping\n");
            }
        }
    }

    // DETERMINE THE REACTION STEP SIZES FOR MAIN STEP AND IF PHASES DIE
    //
    // Don't do this step if there is a phase pop
    size_t iphaseDelete = npos;
    size_t kspec;
    if (iphasePop == npos) {
        // Figure out the new reaction step sizes for the major species (do
        // minor species in the future too)
        kspec = npos;
        iphaseDelete = vcs_RxnStepSizes(forceComponentCalc, kspec);
    } else if (m_debug_print_lvl >= 2) {
        plogf("   --- vcs_RxnStepSizes not called because alternative"
              "phase creation delta was used instead\n");
    }
    size_t doPhaseDeleteKspec = npos;
    size_t doPhaseDeleteIph = npos;

    // Zero out the net change in moles of multispecies phases
    m_deltaPhaseMoles.assign(m_deltaPhaseMoles.size(), 0.0);

    // MAIN LOOP IN CALCULATION: LOOP OVER IRXN TO DETERMINE STEP SIZE
    //
    // Loop through all of the reactions, irxn, pertaining to the formation
    // reaction for species kspec in canonical form.
    //
    // At the end of this loop, we will have a new estimate for the mole numbers
    // for all species consistent with an extent of reaction for all
    // noncomponent species formation reactions. We will have also ensured that
    // all predicted non-component mole numbers are greater than zero.
    //
    //         Old_Solution               New_Solution             Description
    // -----------------------------------------------------------------------------
    //  m_molNumSpecies_old[kspec]   m_molNumSpecies_new[kspec]   Species Mole Numbers
    //                               m_deltaMolNumSpecies[kspec]  Delta in the Species Mole Numbers
    if (iphaseDelete != npos) {
        debuglog("   --- Main Loop Treatment -> Circumvented due to Phase Deletion\n", m_debug_print_lvl >= 2);
        for (size_t k = 0; k < m_nsp; k++) {
            m_molNumSpecies_new[k] = m_molNumSpecies_old[k] + m_deltaMolNumSpecies[k];
            size_t iph = m_phaseID[k];
            m_tPhaseMoles_new[iph] += m_deltaMolNumSpecies[k];
        }
        if (kspec >= m_numComponents) {
            if (m_SSPhase[kspec] == 1) {
                m_speciesStatus[kspec] = VCS_SPECIES_ZEROEDSS;
            } else {
                throw CanteraError("VCS_SOLVE::solve_tp_inner",
                                   "we shouldn't be here!");
            }
            ++m_numRxnMinorZeroed;
            allMinorZeroedSpecies = (m_numRxnMinorZeroed == m_numRxnRdc);
        }

        // Set the flags indicating the mole numbers in the vcs_VolPhase objects
        // are out of date.
        vcs_setFlagsVolPhases(false, VCS_STATECALC_NEW);

        // Calculate the new chemical potentials using the tentative solution
        // values. We only calculate a subset of these, because we have only
        // updated a subset of the W().
        vcs_dfe(VCS_STATECALC_NEW, 0, 0, m_nsp);

        // Evaluate DeltaG for all components if ITI=0, and for major components
        // only if ITI NE 0
        vcs_deltag(0, false, VCS_STATECALC_NEW);
    } else {
        if (m_debug_print_lvl >= 2) {
            plogf("   --- Main Loop Treatment of each non-component species ");
            if (iti == 0) {
                plogf("- Full Calculation:\n");
            } else {
                plogf("- Major Components Calculation:\n");
            }
            plogf("   --- Species     IC    ");
            plogf("  KMoles   Tent_KMoles  Rxn_Adj   |    Comment \n");
        }
        for (size_t irxn = 0; irxn < m_numRxnRdc; irxn++) {
            size_t kspec = m_indexRxnToSpecies[irxn];
            double* sc_irxn = m_stoichCoeffRxnMatrix.ptrColumn(irxn);
            size_t iph = m_phaseID[kspec];
            vcs_VolPhase* Vphase = m_VolPhaseList[iph].get();
            ANOTE[0] = '\0';
            double dx;

            if (iphasePop != npos) {
                if (iph == iphasePop) {
                    dx = m_deltaMolNumSpecies[kspec];
                    m_molNumSpecies_new[kspec] = m_molNumSpecies_old[kspec] + m_deltaMolNumSpecies[kspec];
                    sprintf(ANOTE, "Phase pop");
                } else {
                    dx = 0.0;
                    m_molNumSpecies_new[kspec] = m_molNumSpecies_old[kspec];
                }
            } else if (m_speciesStatus[kspec] == VCS_SPECIES_INTERFACIALVOLTAGE) {
                // VOLTAGE SPECIES
                bool soldel_ret;
                dx = vcs_minor_alt_calc(kspec, irxn, &soldel_ret, ANOTE);
                m_deltaMolNumSpecies[kspec] = dx;
            } else if (m_speciesStatus[kspec] < VCS_SPECIES_MINOR) {
                // ZEROED OUT SPECIES
                bool resurrect = (m_deltaMolNumSpecies[kspec] > 0.0);
                if (m_debug_print_lvl >= 3) {
                    plogf("   --- %s currently zeroed (SpStatus=%-2d):",
                          m_speciesName[kspec], m_speciesStatus[kspec]);
                    plogf("%3d DG = %11.4E WT = %11.4E W = %11.4E DS = %11.4E\n",
                          irxn, m_deltaGRxn_new[irxn], m_molNumSpecies_new[kspec],
                          m_molNumSpecies_old[kspec], m_deltaMolNumSpecies[kspec]);
                }
                if (m_deltaGRxn_new[irxn] >= 0.0 || !resurrect) {
                    m_molNumSpecies_new[kspec] = m_molNumSpecies_old[kspec];
                    m_deltaMolNumSpecies[kspec] = 0.0;
                    resurrect = false;
                    sprintf(ANOTE, "Species stays zeroed: DG = %11.4E", m_deltaGRxn_new[irxn]);
                    if (m_deltaGRxn_new[irxn] < 0.0) {
                        if (m_speciesStatus[kspec] == VCS_SPECIES_STOICHZERO) {
                            sprintf(ANOTE, "Species stays zeroed even though dg neg due to "
                                    "STOICH/PHASEPOP constraint: DG = %11.4E",
                                    m_deltaGRxn_new[irxn]);
                        } else {
                            sprintf(ANOTE, "Species stays zeroed even though dg neg: DG = %11.4E, ds zeroed",
                                    m_deltaGRxn_new[irxn]);
                        }
                    }
                } else {
                    for (size_t j = 0; j < m_nelem; ++j) {
                        int elType = m_elType[j];
                        if (elType == VCS_ELEM_TYPE_ABSPOS) {
                            double atomComp = m_formulaMatrix(kspec,j);
                            if (atomComp > 0.0) {
                                double maxPermissible = m_elemAbundancesGoal[j] / atomComp;
                                if (maxPermissible < VCS_DELETE_MINORSPECIES_CUTOFF) {
                                    sprintf(ANOTE, "Species stays zeroed even though dG "
                                            "neg, because of %s elemAbund",
                                            m_elementName[j].c_str());
                                    resurrect = false;
                                    break;
                                }
                            }
                        }
                    }
                }

                // Resurrect the species
                if (resurrect) {
                    if (Vphase->exists() == VCS_PHASE_EXIST_NO) {
                        if (m_debug_print_lvl >= 2) {
                            plogf("   --- Zeroed species changed to major: ");
                            plogf("%-12s\n", m_speciesName[kspec]);
                        }
                        m_speciesStatus[kspec] = VCS_SPECIES_MAJOR;
                        allMinorZeroedSpecies = false;
                    } else {
                        if (m_debug_print_lvl >= 2) {
                            plogf("   --- Zeroed species changed to minor: ");
                            plogf("%-12s\n", m_speciesName[kspec]);
                        }
                        m_speciesStatus[kspec] = VCS_SPECIES_MINOR;
                    }
                    if (m_deltaMolNumSpecies[kspec] > 0.0) {
                        dx = m_deltaMolNumSpecies[kspec] * 0.01;
                        m_molNumSpecies_new[kspec] = m_molNumSpecies_old[kspec] + dx;
                    } else {
                        m_molNumSpecies_new[kspec] = m_totalMolNum * VCS_DELETE_PHASE_CUTOFF * 10.;
                        dx = m_molNumSpecies_new[kspec] - m_molNumSpecies_old[kspec];
                    }
                    m_deltaMolNumSpecies[kspec] = dx;
                    sprintf(ANOTE, "Born:IC=-1 to IC=1:DG=%11.4E", m_deltaGRxn_new[irxn]);
                } else {
                    m_molNumSpecies_new[kspec] = m_molNumSpecies_old[kspec];
                    m_deltaMolNumSpecies[kspec] = 0.0;
                    dx = 0.0;
                }
            } else if (m_speciesStatus[kspec] == VCS_SPECIES_MINOR) {
                // MINOR SPECIES
                //
                // Unless ITI isn't equal to zero we zero out changes to minor
                // species.
                if (iti != 0) {
                    m_molNumSpecies_new[kspec] = m_molNumSpecies_old[kspec];
                    m_deltaMolNumSpecies[kspec] = 0.0;
                    dx = 0.0;
                    sprintf(ANOTE,"minor species not considered");
                    if (m_debug_print_lvl >= 2) {
                        plogf("   --- %-12s%3d% 11.4E %11.4E %11.4E | %s\n",
                              m_speciesName[kspec], m_speciesStatus[kspec], m_molNumSpecies_old[kspec], m_molNumSpecies_new[kspec],
                              m_deltaMolNumSpecies[kspec], ANOTE);
                    }
                    continue;
                }

                // Minor species alternative calculation
                //
                // This is based upon the following approximation:
                // The mole fraction changes due to these reactions don't affect
                // the mole numbers of the component species. Therefore the
                // following approximation is valid for an ideal solution
                //    0 = DG(I) + log(WT(I)/W(I))
                //    (DG contains the contribution from FF(I) + log(W(I)/TL) )
                // Thus,
                //     WT(I) = W(I) EXP(-DG(I))
                // If soldel is true on return, then we branch to the section
                // that deletes a species from the current set of active species.
                bool soldel_ret;
                dx = vcs_minor_alt_calc(kspec, irxn, &soldel_ret, ANOTE);
                m_deltaMolNumSpecies[kspec] = dx;
                m_molNumSpecies_new[kspec] = m_molNumSpecies_old[kspec] + dx;
                if (soldel_ret) {
                    // DELETE MINOR SPECIES LESS THAN  VCS_DELETE_SPECIES_CUTOFF
                    // MOLE NUMBER
                    if (m_debug_print_lvl >= 2) {
                        plogf("   --- Delete minor species in multispec phase: %-12s\n",
                              m_speciesName[kspec]);
                    }
                    m_deltaMolNumSpecies[kspec] = 0.0;

                    // Delete species, kspec. The alternate return is for the
                    // case where all species become deleted. Then, we need to
                    // branch to the code where we reevaluate the deletion of
                    // all species.
                    size_t lnospec = vcs_delete_species(kspec);
                    if (lnospec) {
                        stage = RECHECK_DELETED;
                        break;
                    }

                    // Go back to consider the next species in the list. Note,
                    // however, that the next species in the list is now in slot
                    // l. In deleting the previous species L, We have exchanged
                    // slot MR with slot l, and then have decremented MR.
                    // Therefore, we will decrement the species counter, here.
                    --irxn;
                    continue;
                }
            } else {
                // MAJOR SPECIES
                sprintf(ANOTE, "Normal Major Calc");

                // Check for superconvergence of the formation reaction. Do
                // nothing if it is superconverged. Skip to the end of the irxn
                // loop if it is superconverged.
                if (fabs(m_deltaGRxn_new[irxn]) <= m_tolmaj2) {
                    m_molNumSpecies_new[kspec] = m_molNumSpecies_old[kspec];
                    m_deltaMolNumSpecies[kspec] = 0.0;
                    dx = 0.0;
                    sprintf(ANOTE, "major species is converged");
                    if (m_debug_print_lvl >= 2) {
                        plogf("   --- %-12s %3d %11.4E %11.4E %11.4E | %s\n",
                              m_speciesName[kspec], m_speciesStatus[kspec], m_molNumSpecies_old[kspec], m_molNumSpecies_new[kspec],
                              m_deltaMolNumSpecies[kspec], ANOTE);
                    }
                    continue;
                }

                // Set the initial step size, dx, equal to the value produced by
                // the routine, vcs_RxnStepSize().
                //
                // Note the multiplication logic is to make sure that dg[]
                // didn't change sign due to w[] changing in the middle of the
                // iteration. (it can if a single species phase goes out of
                // existence).
                if ((m_deltaGRxn_new[irxn] * m_deltaMolNumSpecies[kspec]) <= 0.0) {
                    dx = m_deltaMolNumSpecies[kspec];
                } else {
                    dx = 0.0;
                    m_deltaMolNumSpecies[kspec] = 0.0;
                    sprintf(ANOTE, "dx set to 0, DG flipped sign due to "
                            "changed initial point");
                }

                //Form a tentative value of the new species moles
                m_molNumSpecies_new[kspec] = m_molNumSpecies_old[kspec] + dx;

                // Check for non-positive mole fraction of major species. If we
                // find one, we branch to a section below. Then, depending upon
                // the outcome, we branch to sections below, or we restart the
                // entire iteration.
                if (m_molNumSpecies_new[kspec] <= 0.0) {
                    sprintf(ANOTE, "initial nonpos kmoles= %11.3E",
                            m_molNumSpecies_new[kspec]);
                    // NON-POSITIVE MOLES OF MAJOR SPECIES
                    //
                    // We are here when a tentative value of a mole fraction
                    // created by a tentative value of M_DELTAMOLNUMSPECIES(*)
                    // is negative. We branch from here depending upon whether
                    // this species is in a single species phase or in a
                    // multispecies phase.
                    if (!m_SSPhase[kspec]) {
                        // Section for multispecies phases:
                        //   - Cut reaction adjustment for positive kmoles of
                        //     major species in multispecies phases. Decrease
                        //     its concentration by a factor of 10.
                        dx = -0.9 * m_molNumSpecies_old[kspec];
                        m_deltaMolNumSpecies[kspec] = dx;
                        m_molNumSpecies_new[kspec] = m_molNumSpecies_old[kspec] + dx;
                    } else {
                        // Section for single species phases:
                        //     Calculate a dx that will wipe out the
                        //     moles in the phase.
                        dx = -m_molNumSpecies_old[kspec];

                        // Calculate an update that doesn't create a negative
                        // mole number for a component species. Actually,
                        // restrict this a little more so that the component
                        // values can only be reduced by two 99%,
                        for (size_t j = 0; j < m_numComponents; ++j) {
                            if (sc_irxn[j] != 0.0) {
                                m_wx[j] = m_molNumSpecies_old[j] + sc_irxn[j] * dx;
                                if (m_wx[j] <= m_molNumSpecies_old[j] * 0.01 - 1.0E-150) {
                                    dx = std::max(dx, m_molNumSpecies_old[j] * -0.99 / sc_irxn[j]);
                                }
                            } else {
                                m_wx[j] = m_molNumSpecies_old[j];
                            }
                        }
                        m_molNumSpecies_new[kspec] = m_molNumSpecies_old[kspec] + dx;
                        if (m_molNumSpecies_new[kspec] > 0.0) {
                            m_deltaMolNumSpecies[kspec] = dx;
                            sprintf(ANOTE,
                                    "zeroing SS phase created a neg component species "
                                    "-> reducing step size instead");
                        } else {
                            // We are going to zero the single species phase.
                            // Set the existence flag
                            iph = m_phaseID[kspec];
                            Vphase = m_VolPhaseList[iph].get();
                            sprintf(ANOTE, "zeroing out SS phase: ");

                            // Change the base mole numbers for the iteration.
                            // We need to do this here, because we have decided
                            // to eliminate the phase in this special section
                            // outside the main loop.
                            m_molNumSpecies_new[kspec] = 0.0;
                            doPhaseDeleteIph = iph;

                            doPhaseDeleteKspec = kspec;
                            if (m_debug_print_lvl >= 2 && m_speciesStatus[kspec] >= 0) {
                                plogf("   --- SS species changed to zeroedss: ");
                                plogf("%-12s\n", m_speciesName[kspec]);
                            }
                            m_speciesStatus[kspec] = VCS_SPECIES_ZEROEDSS;
                            ++m_numRxnMinorZeroed;
                            allMinorZeroedSpecies = (m_numRxnMinorZeroed == m_numRxnRdc);

                            for (size_t kk = 0; kk < m_nsp; kk++) {
                                m_deltaMolNumSpecies[kk] = 0.0;
                                m_molNumSpecies_new[kk] = m_molNumSpecies_old[kk];
                            }
                            m_deltaMolNumSpecies[kspec] = dx;
                            m_molNumSpecies_new[kspec] = 0.0;

                            for (size_t k = 0; k < m_numComponents; ++k) {
                                m_deltaMolNumSpecies[k] = 0.0;
                            }
                            for (iph = 0; iph < m_numPhases; iph++) {
                                m_deltaPhaseMoles[iph] = 0.0;
                            }
                        }
                    }
                }
            } // End of Loop on ic[irxn] -> the type of species

            // CALCULATE KMOLE NUMBER CHANGE FOR THE COMPONENT BASIS
            if (dx != 0.0 && (m_speciesUnknownType[kspec] !=
                              VCS_SPECIES_TYPE_INTERFACIALVOLTAGE)) {

                // Change the amount of the component compounds according to the
                // reaction delta that we just computed. This should keep the
                // amount of material constant.
                AssertThrowMsg(fabs(m_deltaMolNumSpecies[kspec] -dx) <
                        1.0E-14*(fabs(m_deltaMolNumSpecies[kspec]) + fabs(dx) + 1.0E-32),
                        "VCS_SOLVE::solve_tp_inner",
                        "ds[kspec] = {}, dx = {}, kspec = {}\nwe have a problem!",
                        m_deltaMolNumSpecies[kspec], dx, kspec);
                for (size_t k = 0; k < m_numComponents; ++k) {
                    m_deltaMolNumSpecies[k] += sc_irxn[k] * dx;
                }

                // Calculate the tentative change in the total number of moles
                // in all of the phases
                for (iph = 0; iph < m_numPhases; iph++) {
                    m_deltaPhaseMoles[iph] += dx * m_deltaMolNumPhase(iph,irxn);
                }
            }

            checkDelta1(&m_deltaMolNumSpecies[0],
                        &m_deltaPhaseMoles[0], kspec+1);

            // Branch point for returning
            if (m_debug_print_lvl >= 2) {
                m_molNumSpecies_new[kspec] = m_molNumSpecies_old[kspec] + m_deltaMolNumSpecies[kspec];
                plogf("   --- %-12.12s%3d %11.4E %11.4E %11.4E | %s\n",
                      m_speciesName[kspec], m_speciesStatus[kspec],
                      m_molNumSpecies_old[kspec], m_molNumSpecies_new[kspec],
                      m_deltaMolNumSpecies[kspec], ANOTE);
            }

            if (doPhaseDeleteIph != npos) {
                if (m_debug_print_lvl >= 2) {
                    plogf("   --- %-12.12s Main Loop Special Case deleting phase with species:\n",
                          m_speciesName[doPhaseDeleteKspec]);
                }
                break;
            }
        } // END OF MAIN LOOP OVER FORMATION REACTIONS
        if (m_debug_print_lvl >= 2) {
            for (size_t k = 0; k < m_numComponents; k++) {
                plogf("   --- ");
                plogf("%-12.12s", m_speciesName[k]);
                plogf("  c %11.4E %11.4E %11.4E |\n",
                      m_molNumSpecies_old[k],
                      m_molNumSpecies_old[k]+m_deltaMolNumSpecies[k], m_deltaMolNumSpecies[k]);
            }
            plogf("   ");
            writeline('-', 80);
            plogf("   --- Finished Main Loop\n");
        }

        // LIMIT REDUCTION OF BASIS SPECIES TO 99%
        //
        // We have a tentative m_deltaMolNumSpecies[]. Now apply other criteria
        // to limit its magnitude.
        double par = 0.5;
        size_t ll;
        for (size_t k = 0; k < m_numComponents; ++k) {
            if (m_molNumSpecies_old[k] > 0.0) {
                double xx = -m_deltaMolNumSpecies[k] / m_molNumSpecies_old[k];
                if (par < xx) {
                    par = xx;
                    ll = k;
                }
            } else if (m_deltaMolNumSpecies[k] < 0.0) {
                // If we are here, we then do a step which violates element
                // conservation.
                size_t iph = m_phaseID[k];
                m_deltaPhaseMoles[iph] -= m_deltaMolNumSpecies[k];
                m_deltaMolNumSpecies[k] = 0.0;
            }
        }
        par = 1.0 / par;
        if (par <= 1.01 && par > 0.0) {
            // Reduce the size of the step by the multiplicative factor, par
            par *= 0.99;
            if (m_debug_print_lvl >= 2) {
                plogf("   --- Reduction in step size due to component %s going negative = %11.3E\n",
                    m_speciesName[ll], par);
            }
            for (size_t i = 0; i < m_nsp; ++i) {
                m_deltaMolNumSpecies[i] *= par;
            }
            for (size_t iph = 0; iph < m_numPhases; iph++) {
                m_deltaPhaseMoles[iph] *= par;
            }
        } else {
            par = 1.0;
        }
        checkDelta1(&m_deltaMolNumSpecies[0],
                    &m_deltaPhaseMoles[0], m_nsp);

        // Now adjust the wt[kspec]'s so that the reflect the decrease in the
        // overall length of m_deltaMolNumSpecies[kspec] just calculated. At the
        // end of this section wt[], m_deltaMolNumSpecies[], tPhMoles, and
        // tPhMoles1 should all be consistent with a new estimate of the state
        // of the system.
        for (size_t kspec = 0; kspec < m_nsp; ++kspec) {
            m_molNumSpecies_new[kspec] = m_molNumSpecies_old[kspec] + m_deltaMolNumSpecies[kspec];
            if (m_molNumSpecies_new[kspec] < 0.0 && (m_speciesUnknownType[kspec]
                    != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE)) {
                throw CanteraError("VCS_SOLVE::solve_tp_inner",
                    "vcs_solve_TP: ERROR on step change wt[{}:{}]: {} < 0.0",
                    kspec, m_speciesName[kspec], m_molNumSpecies_new[kspec]);
            }
        }

        // Calculate the tentative total mole numbers for each phase
        for (size_t iph = 0; iph < m_numPhases; iph++) {
            m_tPhaseMoles_new[iph] = m_tPhaseMoles_old[iph] + m_deltaPhaseMoles[iph];
        }

        // Set the flags indicating the mole numbers in the vcs_VolPhase objects
        // are out of date.
        vcs_setFlagsVolPhases(false, VCS_STATECALC_NEW);

        // Calculate the new chemical potentials using the tentative solution
        // values. We only calculate a subset of these, because we have only
        // updated a subset of the W().
        vcs_dfe(VCS_STATECALC_NEW, 0, 0, m_nsp);

        // Evaluate DeltaG for all components if ITI=0, and for major components
        // only if ITI NE 0
        vcs_deltag(0, false, VCS_STATECALC_NEW);

        // CONVERGENCE FORCER SECTION
        if (printDetails) {
            plogf("   --- Total Old       Dimensionless Gibbs Free Energy = %20.13E\n",
                  vcs_Total_Gibbs(&m_molNumSpecies_old[0], &m_feSpecies_old[0],
                                  &m_tPhaseMoles_old[0]));
            plogf("   --- Total tentative Dimensionless Gibbs Free Energy = %20.13E\n",
                  vcs_Total_Gibbs(&m_molNumSpecies_new[0], &m_feSpecies_new[0],
                                  &m_tPhaseMoles_new[0]));
        }

        bool forced = vcs_globStepDamp();

        // Print out the changes to the solution that FORCER produced
        if (printDetails && forced) {
            plogf(" -----------------------------------------------------\n");
            plogf("   --- FORCER SUBROUTINE changed the solution:\n");
            plogf("   --- SPECIES Status INIT MOLES TENT_MOLES");
            plogf("  FINAL KMOLES  INIT_DEL_G/RT  TENT_DEL_G/RT  FINAL_DELTA_G/RT\n");
            for (size_t i = 0; i < m_numComponents; ++i) {
                plogf("  --- %-12.12s", m_speciesName[i]);
                plogf("    %14.6E %14.6E %14.6E\n", m_molNumSpecies_old[i],
                      m_molNumSpecies_old[i] + m_deltaMolNumSpecies[i], m_molNumSpecies_new[i]);
            }
            for (size_t kspec = m_numComponents; kspec < m_numSpeciesRdc; ++kspec) {
                size_t irxn = kspec - m_numComponents;
                plogf("  --- %-12.12s", m_speciesName[kspec]);
                plogf(" %2d %14.6E%14.6E%14.6E%14.6E%14.6E%14.6E\n", m_speciesStatus[kspec],
                      m_molNumSpecies_old[kspec],
                      m_molNumSpecies_old[kspec]+m_deltaMolNumSpecies[kspec],
                      m_molNumSpecies_new[kspec], m_deltaGRxn_old[irxn],
                      m_deltaGRxn_tmp[irxn], m_deltaGRxn_new[irxn]);
            }
            writeline(' ', 26, false);
            plogf("Norms of Delta G():%14.6E%14.6E\n",
                  l2normdg(&m_deltaGRxn_old[0]),
                  l2normdg(&m_deltaGRxn_new[0]));
            plogf("   Total kmoles of gas    = %15.7E\n", m_tPhaseMoles_old[0]);
            if ((m_numPhases > 1) && (!m_VolPhaseList[1]->m_singleSpecies)) {
                plogf("   Total kmoles of liquid = %15.7E\n", m_tPhaseMoles_old[1]);
            } else {
                plogf("   Total kmoles of liquid = %15.7E\n", 0.0);
            }
            plogf("   Total New Dimensionless Gibbs Free Energy = %20.13E\n",
                  vcs_Total_Gibbs(&m_molNumSpecies_new[0], &m_feSpecies_new[0],
                                  &m_tPhaseMoles_new[0]));
            plogf(" -----------------------------------------------------\n");
        }

    }

    // ITERATION SUMMARY PRINTOUT SECTION
    if (printDetails) {
        plogf("   ");
        writeline('-', 103);
        plogf("   --- Summary of the Update ");
        if (iti == 0) {
            plogf(" (all species):");
        } else {
            plogf(" (only major species):");
        }
        plogf("\n");
        plogf("   ---      Species Status Initial_KMoles Final_KMoles Initial_Mu/RT");
        plogf("     Mu/RT     Init_Del_G/RT   Delta_G/RT\n");
        for (size_t i = 0; i < m_numComponents; ++i) {
            plogf("   ---   %-12.12s", m_speciesName[i]);
            plogf("    ");
            plogf("%14.6E%14.6E%14.6E%14.6E\n", m_molNumSpecies_old[i],
                  m_molNumSpecies_new[i], m_feSpecies_old[i], m_feSpecies_new[i]);
        }
        for (size_t i = m_numComponents; i < m_numSpeciesRdc; ++i) {
            size_t l1 = i - m_numComponents;
            plogf("   ---   %-12.12s", m_speciesName[i]);
            plogf(" %2d %14.6E%14.6E%14.6E%14.6E%14.6E%14.6E\n",
                  m_speciesStatus[i], m_molNumSpecies_old[i],
                  m_molNumSpecies_new[i], m_feSpecies_old[i], m_feSpecies_new[i],
                  m_deltaGRxn_old[l1], m_deltaGRxn_new[l1]);
        }
        for (size_t kspec = m_numSpeciesRdc; kspec < m_nsp; ++kspec) {
            size_t l1 = kspec - m_numComponents;
            plogf("   ---   %-12.12s", m_speciesName[kspec]);
            plogf(" %2d %14.6E%14.6E%14.6E%14.6E%14.6E%14.6E\n",
                  m_speciesStatus[kspec], m_molNumSpecies_old[kspec],
                  m_molNumSpecies_new[kspec], m_feSpecies_old[kspec], m_feSpecies_new[kspec],
                  m_deltaGRxn_old[l1], m_deltaGRxn_new[l1]);
        }
        plogf("   ---");
        writeline(' ', 56, false);
        plogf("Norms of Delta G():%14.6E%14.6E\n",
              l2normdg(&m_deltaGRxn_old[0]),
              l2normdg(&m_deltaGRxn_new[0]));

        plogf("   ---           Phase_Name    KMoles(after update)\n");
        plogf("   ---   ");
        writeline('-', 50);
        for (size_t iph = 0; iph < m_numPhases; iph++) {
            vcs_VolPhase* Vphase = m_VolPhaseList[iph].get();
            plogf("   ---   %18s = %15.7E\n", Vphase->PhaseName, m_tPhaseMoles_new[iph]);
        }
        plogf("   ");
        writeline('-', 103);
        plogf("   --- Total Old Dimensionless Gibbs Free Energy = %20.13E\n",
              vcs_Total_Gibbs(&m_molNumSpecies_old[0], &m_feSpecies_old[0],
                              &m_tPhaseMoles_old[0]));
        plogf("   --- Total New Dimensionless Gibbs Free Energy = %20.13E\n",
              vcs_Total_Gibbs(&m_molNumSpecies_new[0], &m_feSpecies_new[0],
                              &m_tPhaseMoles_new[0]));
        debuglog("   --- Troublesome solve\n", m_VCount->Its > 550);
    }

    // RESET VALUES AT END OF ITERATION
    // UPDATE MOLE NUMBERS
    //
    // If the solution wasn't changed in the forcer routine, then copy the
    // tentative mole numbers and Phase moles into the actual mole numbers and
    // phase moles. We will consider this current step to be completed.
    //
    // Accept the step. -> the tentative solution now becomes the real solution.
    // If FORCED is true, then we have already done this inside the FORCED loop.
    vcs_updateMolNumVolPhases(VCS_STATECALC_NEW);
    m_tPhaseMoles_old = m_tPhaseMoles_new;
    m_molNumSpecies_old = m_molNumSpecies_new;
    m_actCoeffSpecies_old = m_actCoeffSpecies_new;
    m_deltaGRxn_old = m_deltaGRxn_new;
    m_feSpecies_old = m_feSpecies_new;

    vcs_setFlagsVolPhases(true, VCS_STATECALC_OLD);

    // Increment the iteration counters
    ++m_VCount->Its;
    ++it1;
    if (m_debug_print_lvl >= 2) {
        plogf("   --- Increment counter increased, step is accepted: %4d\n",
              m_VCount->Its);
    }

    // HANDLE DELETION OF MULTISPECIES PHASES
    //
    // We delete multiphases, when the total moles in the multiphase is reduced
    // below a relative threshold. Set microscopic multispecies phases with
    // total relative number of moles less than VCS_DELETE_PHASE_CUTOFF to
    // absolute zero.
    bool justDeletedMultiPhase = false;
    for (size_t iph = 0; iph < m_numPhases; iph++) {
        if (!m_VolPhaseList[iph]->m_singleSpecies && m_tPhaseMoles_old[iph] != 0.0 &&
                    m_tPhaseMoles_old[iph]/m_totalMolNum <= VCS_DELETE_PHASE_CUTOFF) {
            if (m_debug_print_lvl >= 1) {
                plogf("   --- Setting microscopic phase %d to zero\n", iph);
            }
            justDeletedMultiPhase = true;
            vcs_delete_multiphase(iph);
        }
    }

    // If we have deleted a multispecies phase because the equilibrium moles
    // decreased, then we will update all the component basis calculation, and
    // therefore all of the thermo functions just to be safe.
    if (justDeletedMultiPhase) {
        bool usedZeroedSpecies;
        double test = -1.0e-10;
        int retn = vcs_basopt(false, &m_aw[0], &m_sa[0], &m_sm[0], &m_ss[0],
                              test, &usedZeroedSpecies);
        if (retn != VCS_SUCCESS) {
            throw CanteraError("VCS_SOLVE::solve_tp_inner",
                               "BASOPT returned with an error condition");
        }
        vcs_setFlagsVolPhases(false, VCS_STATECALC_OLD);
        vcs_dfe(VCS_STATECALC_OLD, 0, 0, m_numSpeciesRdc);
        vcs_deltag(0, true, VCS_STATECALC_OLD);
        iti = 0;
        return;
    }

    // CHECK FOR ELEMENT ABUNDANCE
    if (m_debug_print_lvl >= 2) {
        plogf("   --- Normal element abundance check");
    }
    vcs_elab();
    if (! vcs_elabcheck(0)) {
        debuglog(" - failed -> redoing element abundances.\n", m_debug_print_lvl >= 2);
        vcs_elcorr(&m_sm[0], &m_wx[0]);
        vcs_setFlagsVolPhases(false, VCS_STATECALC_OLD);
        vcs_dfe(VCS_STATECALC_OLD, 0, 0, m_numSpeciesRdc);
        vcs_deltag(0, true, VCS_STATECALC_OLD);
        uptodate_minors = true;
    } else {
        debuglog(" - passed\n", m_debug_print_lvl >= 2);
    }

    // CHECK FOR OPTIMUM BASIS
    for (size_t i = 0; i < m_numRxnRdc; ++i) {
        size_t k = m_indexRxnToSpecies[i];
        if (m_speciesUnknownType[k] == VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
            continue;
        }
        for (size_t j = 0; j < m_numComponents; ++j) {
            bool doSwap = false;
            if (m_SSPhase[j]) {
                doSwap = (m_molNumSpecies_old[k] * m_spSize[k]) >
                         (m_molNumSpecies_old[j] * m_spSize[j] * 1.01);
                if (!m_SSPhase[k] && doSwap) {
                    doSwap = m_molNumSpecies_old[k] > (m_molNumSpecies_old[j] * 1.01);
                }
            } else {
                if (m_SSPhase[k]) {
                    doSwap = (m_molNumSpecies_old[k] * m_spSize[k]) >
                             (m_molNumSpecies_old[j] * m_spSize[j] * 1.01);
                    if (!doSwap) {
                        doSwap = m_molNumSpecies_old[k] > (m_molNumSpecies_old[j] * 1.01);
                    }
                } else {
                    doSwap = (m_molNumSpecies_old[k] * m_spSize[k]) >
                             (m_molNumSpecies_old[j] * m_spSize[j] * 1.01);
                }
            }
            if (doSwap && m_stoichCoeffRxnMatrix(j,i) != 0.0) {
                if (m_debug_print_lvl >= 2) {
                    plogf("   --- Get a new basis because %s", m_speciesName[k]);
                    plogf(" is better than comp %s", m_speciesName[j]);
                    plogf(" and share nonzero stoic: %-9.1f\n",
                          m_stoichCoeffRxnMatrix(j,i));
                }
                forceComponentCalc = 1;
                return;
            }
        }
    }
    debuglog("   --- Check for an optimum basis passed\n", m_debug_print_lvl >= 2);
    stage = EQUILIB_CHECK;

    // RE-EVALUATE MAJOR-MINOR VECTOR IF NECESSARY
    //
    // Skip this section if we haven't done a full calculation. Go right to the
    // check equilibrium section
    if (iti != 0) {
        return;
    }
    if (m_debug_print_lvl >= 2) {
        plogf("   --- Reevaluate major-minor status of noncomponents:\n");
    }
    m_numRxnMinorZeroed = 0;
    for (size_t irxn = 0; irxn < m_numRxnRdc; irxn++) {
        size_t kspec = m_indexRxnToSpecies[irxn];
        int speciesType = vcs_species_type(kspec);
        if (speciesType < VCS_SPECIES_MINOR) {
            if (m_debug_print_lvl >= 2 && m_speciesStatus[kspec] >= VCS_SPECIES_MINOR) {
                plogf("   ---    major/minor species is now zeroed out: %s\n",
                      m_speciesName[kspec]);
            }
            ++m_numRxnMinorZeroed;
        } else if (speciesType == VCS_SPECIES_MINOR) {
            if (m_debug_print_lvl >= 2 && m_speciesStatus[kspec] != VCS_SPECIES_MINOR) {
                if (m_speciesStatus[kspec] == VCS_SPECIES_MAJOR) {
                    plogf("   ---   Noncomponent turned from major to minor: ");
                } else if (kspec < m_numComponents) {
                    plogf("   ---   Component turned into a minor species: ");
                } else {
                    plogf("   ---   Zeroed Species turned into a "
                          "minor species: ");
                }
                plogf("%s\n", m_speciesName[kspec]);
            }
            ++m_numRxnMinorZeroed;
        } else if (speciesType == VCS_SPECIES_MAJOR) {
            if (m_speciesStatus[kspec] != VCS_SPECIES_MAJOR) {
                if (m_debug_print_lvl >= 2) {
                    if (m_speciesStatus[kspec] == VCS_SPECIES_MINOR) {
                        plogf("   ---   Noncomponent turned from minor to major: ");
                    } else if (kspec < m_numComponents) {
                        plogf("   ---   Component turned into a major: ");
                    } else {
                        plogf("   ---   Noncomponent turned from zeroed to major: ");
                    }
                    plogf("%s\n", m_speciesName[kspec]);
                }
                m_speciesStatus[kspec] = VCS_SPECIES_MAJOR;
            }
        }
        m_speciesStatus[kspec] = speciesType;
    }

    // This logical variable indicates whether all current non-component species
    // are minor or nonexistent
    allMinorZeroedSpecies = (m_numRxnMinorZeroed == m_numRxnRdc);
}

void VCS_SOLVE::solve_tp_equilib_check(bool& allMinorZeroedSpecies,
                                       bool& uptodate_minors,
                                       bool& giveUpOnElemAbund, int& solveFail,
                                       size_t& iti, size_t& it1, int maxit,
                                       int& stage, bool& lec)
{
    if (! allMinorZeroedSpecies) {
        if (m_debug_print_lvl >= 2) {
            plogf("   --- Equilibrium check for major species: ");
        }
        for (size_t irxn = 0; irxn < m_numRxnRdc; ++irxn) {
            size_t kspec = irxn + m_numComponents;
            if (m_speciesStatus[kspec] == VCS_SPECIES_MAJOR && (fabs(m_deltaGRxn_new[irxn]) > m_tolmaj)) {
                if (m_VCount->Its >= maxit) {
                    solveFail = -1;

                    // Clean up and exit code even though we haven't converged.
                    // -> we have run out of iterations!
                    stage = RETURN_A;
                } else {
                    if (m_debug_print_lvl >= 2) {
                        plogf("%s failed\n", m_speciesName[m_indexRxnToSpecies[irxn]]);
                    }
                    // Convergence amongst major species has not been achieved.
                    // Go back and do another iteration with variable ITI
                    iti = ((it1/4) *4) - it1;
                    stage = MAIN;
                }
                return;
            }
        }
        debuglog(" MAJOR SPECIES CONVERGENCE achieved", m_debug_print_lvl >= 2);
    } else {
        debuglog(" MAJOR SPECIES CONVERGENCE achieved "
            "(because there are no major species)\n", m_debug_print_lvl >= 2);
    }
    // Convergence amongst major species has been achieved

    // EQUILIBRIUM CHECK FOR MINOR SPECIES
    if (m_numRxnMinorZeroed != 0) {
        // Calculate the chemical potential and reaction DeltaG for minor
        // species, if needed.
        if (iti != 0) {
            vcs_setFlagsVolPhases(false, VCS_STATECALC_OLD);
            vcs_dfe(VCS_STATECALC_OLD, 1, 0, m_numSpeciesRdc);
            vcs_deltag(1, false, VCS_STATECALC_OLD);
            uptodate_minors = true;
        }
        if (m_debug_print_lvl >= 2) {
            plogf("   --- Equilibrium check for minor species: ");
        }
        for (size_t irxn = 0; irxn < m_numRxnRdc; ++irxn) {
            size_t kspec = irxn + m_numComponents;
            if (m_speciesStatus[kspec] == VCS_SPECIES_MINOR && (fabs(m_deltaGRxn_new[irxn]) > m_tolmin)) {
                if (m_VCount->Its >= maxit) {
                    solveFail = -1;
                    // Clean up and exit code. -> Even though we have not
                    // converged, we have run out of iterations !
                    stage = RETURN_A;
                    return;
                }
                if (m_debug_print_lvl >= 2) {
                    plogf("%s failed\n", m_speciesName[m_indexRxnToSpecies[irxn]]);
                }

                // Set iti to zero to force a full calculation, and go back to
                // the main loop to do another iteration.
                iti = 0;
                stage = MAIN;
                return;
            }
        }
        if (m_debug_print_lvl >= 2) {
            plogf(" CONVERGENCE achieved\n");
        }
    }

    // FINAL ELEMENTAL ABUNDANCE CHECK
    // Recalculate the element abundance vector again
    vcs_updateVP(VCS_STATECALC_OLD);
    vcs_elab();

    // LEC is only true when we are near the end game
    if (lec) {
        if (!giveUpOnElemAbund) {
            if (m_debug_print_lvl >= 2) {
                plogf("   --- Check the Full Element Abundances: ");
            }

            // Final element abundance check: If we fail then we need to go back
            // and correct the element abundances, and then go do a major step
            if (! vcs_elabcheck(1)) {
                if (m_debug_print_lvl >= 2) {
                    if (! vcs_elabcheck(0)) {
                        plogf(" failed\n");
                    } else {
                        plogf(" passed for NC but failed for NE: RANGE ERROR\n");
                    }
                }
                // delete?
                stage = ELEM_ABUND_CHECK;
                return;
            }
            if (m_debug_print_lvl >= 2) {
                plogf(" passed\n");
            }
        }

        // If we have deleted a species then we need to recheck the the deleted
        // species, before exiting
        if (m_numSpeciesRdc != m_nsp) {
            stage = RECHECK_DELETED;
            return;
        }
        // Final checks are passed -> go check out
        stage = RETURN_A;
    }
    lec = true;
}

void VCS_SOLVE::solve_tp_elem_abund_check(size_t& iti, int& stage, bool& lec,
                                          bool& giveUpOnElemAbund,
                                          int& finalElemAbundAttempts,
                                          int& rangeErrorFound)
{

    // HKM - Put in an element abundance check. The element abundances were
    // being corrected even if they were perfectly OK to start with. This is
    // actually an expensive operation, so I took it out. Also vcs_dfe() doesn't
    // need to be called if no changes were made.
    rangeErrorFound = 0;
    if (! vcs_elabcheck(1)) {
        bool ncBefore = vcs_elabcheck(0);
        vcs_elcorr(&m_sm[0], &m_wx[0]);
        bool ncAfter = vcs_elabcheck(0);
        bool neAfter = vcs_elabcheck(1);

        // Go back to evaluate the total moles of gas and liquid.
        vcs_setFlagsVolPhases(false, VCS_STATECALC_OLD);
        vcs_dfe(VCS_STATECALC_OLD, 0, 0, m_numSpeciesRdc);
        vcs_deltag(0, false, VCS_STATECALC_OLD);
        if (!ncBefore) {
            if (ncAfter) {
                // We have breathed new life into the old problem. Now the
                // element abundances up to NC agree. Go back and restart the
                // main loop calculation, resetting the end conditions.
                lec = false;
                iti = 0;
                stage = MAIN;
            } else {
                // We are still hosed
                if (finalElemAbundAttempts >= 3) {
                    giveUpOnElemAbund = true;
                    stage = EQUILIB_CHECK;
                } else {
                    finalElemAbundAttempts++;
                    lec = false;
                    iti = 0;
                    stage = MAIN;
                }
            }
            return;
        } else if (ncAfter) {
            if (!neAfter) {
                // Probably an unrecoverable range error
                if (m_debug_print_lvl >= 2) {
                    plogf(" ---  vcs_solve_tp: RANGE SPACE ERROR ENCOUNTERED\n");
                    plogf(" ---  vcs_solve_tp: - Giving up on NE Element Abundance satisfaction \n");
                    plogf(" ---  vcs_solve_tp: - However, NC Element Abundance criteria is satisfied \n");
                    plogf(" ---  vcs_solve_tp: - Returning the calculated equilibrium condition\n");
                }
                rangeErrorFound = 1;
                giveUpOnElemAbund = true;
            }

            // Recovery of end element abundances -> go do equilibrium check
            // again and then check out.
            stage = EQUILIB_CHECK;
            return;
        }
    }
    // Calculate delta g's
    vcs_deltag(0, false, VCS_STATECALC_OLD);
    // Go back to equilibrium check as a prep to eventually checking out
    stage = EQUILIB_CHECK;
}

double VCS_SOLVE::vcs_minor_alt_calc(size_t kspec, size_t irxn, bool* do_delete,
                                     char* ANOTE) const
{
    double w_kspec = m_molNumSpecies_old[kspec];
    double dg_irxn = m_deltaGRxn_old[irxn];
    size_t iph = m_phaseID[kspec];

    *do_delete = false;
    if (m_speciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
        if (w_kspec <= 0.0) {
            w_kspec = VCS_DELETE_MINORSPECIES_CUTOFF;
        }
        dg_irxn = std::max(dg_irxn, -200.0);
        if (ANOTE) {
            sprintf(ANOTE,"minor species alternative calc");
        }
        if (dg_irxn >= 23.0) {
            double molNum_kspec_new = w_kspec * 1.0e-10;
            if (w_kspec < VCS_DELETE_MINORSPECIES_CUTOFF) {
                // delete the species from the current list of active species,
                // because its concentration has gotten too small.
                *do_delete = true;
                return - w_kspec;
            }
            return molNum_kspec_new - w_kspec;
        } else {
            if (fabs(dg_irxn) <= m_tolmin2) {
                return 0.0;
            }
        }

        // get the diagonal of the activity coefficient Jacobian
        double s = m_np_dLnActCoeffdMolNum(kspec,kspec) / m_tPhaseMoles_old[iph];

        // We fit it to a power law approximation of the activity coefficient
        //
        //     gamma = gamma_0 * ( x / x0)**a
        //
        // where a is forced to be a little bit greater than -1.
        // We do this so that the resulting expression is always nonnegative
        // We then solve the resulting calculation:
        //
        //     gamma * x  = gamma_0 * x0 exp (-deltaG/RT);
        double a = clip(w_kspec * s, -1.0+1e-8, 100.0);
        double tmp = clip(-dg_irxn / (1.0 + a), -200.0, 200.0);
        double wTrial = w_kspec * exp(tmp);
        double molNum_kspec_new = wTrial;

        if (wTrial > 100. * w_kspec) {
            double molNumMax = 0.0001 * m_tPhaseMoles_old[iph];
            if (molNumMax < 100. * w_kspec) {
                molNumMax = 100. * w_kspec;
            }
            if (wTrial > molNumMax) {
                molNum_kspec_new = molNumMax;
            } else {
                molNum_kspec_new = wTrial;
            }
        } else if (1.0E10 * wTrial < w_kspec) {
            molNum_kspec_new= 1.0E-10 * w_kspec;
        } else {
            molNum_kspec_new = wTrial;
        }

        if ((molNum_kspec_new) < VCS_DELETE_MINORSPECIES_CUTOFF) {
            // delete the species from the current list of active species,
            // because its concentration has gotten too small.
            *do_delete = true;
            return - w_kspec;
        }
        return molNum_kspec_new - w_kspec;

    } else {
        // Voltage calculation
        // Need to check the sign -> This is good for electrons
        double dx = m_deltaGRxn_old[irxn]/ m_Faraday_dim;
        if (ANOTE) {
            sprintf(ANOTE,"voltage species alternative calc");
        }
        return dx;
    }
}

int VCS_SOLVE::delta_species(const size_t kspec, double* const delta_ptr)
{
    size_t irxn = kspec - m_numComponents;
    int retn = 1;
    AssertThrowMsg(kspec >= m_numComponents, "VCS_SOLVE::delta_species",
        "delete_species() ERROR: called for a component {}", kspec);
    if (m_speciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
        // Attempt the given dx. If it doesn't work, try to see if a smaller one
        // would work,
        double dx = *delta_ptr;
        double* sc_irxn = m_stoichCoeffRxnMatrix.ptrColumn(irxn);
        for (size_t j = 0; j < m_numComponents; ++j) {
            if (m_molNumSpecies_old[j] > 0.0) {
                double tmp = sc_irxn[j] * dx;
                if (-tmp > m_molNumSpecies_old[j]) {
                    retn = 0;
                    dx = std::min(dx, - m_molNumSpecies_old[j] / sc_irxn[j]);
                }
            }

            // If the component has a zero concentration and is a reactant
            // in the formation reaction, then dx == 0.0, and we just return.
            if (m_molNumSpecies_old[j] <= 0.0 && sc_irxn[j] < 0.0) {
                *delta_ptr = 0.0;
                return 0;
            }
        }

        // ok, we found a positive dx. implement it.
        *delta_ptr = dx;
        m_molNumSpecies_old[kspec] += dx;
        size_t iph = m_phaseID[kspec];
        m_tPhaseMoles_old[iph] += dx;
        vcs_setFlagsVolPhase(iph, false, VCS_STATECALC_OLD);

        for (size_t j = 0; j < m_numComponents; ++j) {
            double tmp = sc_irxn[j] * dx;
            if (tmp != 0.0) {
                iph = m_phaseID[j];
                m_molNumSpecies_old[j] += tmp;
                m_tPhaseMoles_old[iph] += tmp;
                vcs_setFlagsVolPhase(iph, false, VCS_STATECALC_OLD);
                m_molNumSpecies_old[j] = std::max(m_molNumSpecies_old[j], 0.0);
            }
        }
    }
    return retn;
}

int VCS_SOLVE::vcs_zero_species(const size_t kspec)
{
    int retn = 1;

    // Calculate a delta that will eliminate the species.
    if (m_speciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
        double dx = -m_molNumSpecies_old[kspec];
        if (dx != 0.0) {
            retn = delta_species(kspec, &dx);
            if (!retn && m_debug_print_lvl >= 1) {
                plogf("vcs_zero_species: Couldn't zero the species %d, "
                      "did delta of %g. orig conc of %g\n",
                      kspec, dx, m_molNumSpecies_old[kspec] + dx);
            }
        }
    }
    return retn;
}

int VCS_SOLVE::vcs_delete_species(const size_t kspec)
{
    const size_t klast = m_numSpeciesRdc - 1;
    const size_t iph = m_phaseID[kspec];
    vcs_VolPhase* const Vphase = m_VolPhaseList[iph].get();
    const size_t irxn = kspec - m_numComponents;

    // Zero the concentration of the species.
    //     -> This zeroes w[kspec] and modifies m_tPhaseMoles_old[]
    const int retn = vcs_zero_species(kspec);
    if (!retn) {
        throw CanteraError("VCS_SOLVE::vcs_delete_species",
                           "Failed to delete a species!");
    }

    // Decrement the minor species counter if the current species is a minor
    // species
    if (m_speciesStatus[kspec] != VCS_SPECIES_MAJOR) {
        --m_numRxnMinorZeroed;
    }
    m_speciesStatus[kspec] = VCS_SPECIES_DELETED;
    m_deltaGRxn_new[irxn] = 0.0;
    m_deltaGRxn_old[irxn] = 0.0;
    m_feSpecies_new[kspec] = 0.0;
    m_feSpecies_old[kspec] = 0.0;
    m_molNumSpecies_new[kspec] = m_molNumSpecies_old[kspec];

    // Rearrange the data if the current species isn't the last active species.
    if (kspec != klast) {
        vcs_switch_pos(true, klast, kspec);
    }

    // Adjust the total moles in a phase downwards.
    Vphase->setMolesFromVCSCheck(VCS_STATECALC_OLD, &m_molNumSpecies_old[0],
                                 &m_tPhaseMoles_old[0]);

    // Adjust the current number of active species and reactions counters
    --m_numRxnRdc;
    --m_numSpeciesRdc;

    // Check to see whether we have just annihilated a multispecies phase. If it
    // is extinct, call the delete_multiphase() function.
    if (! m_SSPhase[klast] && Vphase->exists() != VCS_PHASE_EXIST_ALWAYS) {
        bool stillExists = false;
        for (size_t k = 0; k < m_numSpeciesRdc; k++) {
            if (m_speciesUnknownType[k] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE &&
                m_phaseID[k] == iph && m_molNumSpecies_old[k] > 0.0) {
                stillExists = true;
                break;
            }
        }
        if (!stillExists) {
            vcs_delete_multiphase(iph);
        }
    }

    // When the total number of noncomponent species is zero, we have to signal
    // the calling code
    return (m_numRxnRdc == 0);
}

void VCS_SOLVE::vcs_reinsert_deleted(size_t kspec)
{
    size_t iph = m_phaseID[kspec];
    if (m_debug_print_lvl >= 2) {
        plogf("   --- Add back a deleted species: %-12s\n", m_speciesName[kspec]);
    }

    // Set the species back to minor species status
    //  this adjusts m_molNumSpecies_old[] and m_tPhaseMoles_old[]
    // HKM -> make this a relative mole number!
    double dx = m_tPhaseMoles_old[iph] * VCS_RELDELETE_SPECIES_CUTOFF * 10.;
    delta_species(kspec, &dx);
    m_speciesStatus[kspec] = VCS_SPECIES_MINOR;

    if (m_SSPhase[kspec]) {
        m_speciesStatus[kspec] = VCS_SPECIES_MAJOR;
        --m_numRxnMinorZeroed;
    }

    vcs_VolPhase* Vphase = m_VolPhaseList[iph].get();
    Vphase->setMolesFromVCSCheck(VCS_STATECALC_OLD,
                                 &m_molNumSpecies_old[0],
                                 &m_tPhaseMoles_old[0]);

    // We may have popped a multispecies phase back into existence. If we did,
    // we have to check the other species in that phase. Take care of the
    // m_speciesStatus[] flag. The value of m_speciesStatus[] must change from
    // VCS_SPECIES_ZEROEDPHASE to VCS_SPECIES_ZEROEDMS for those other species.
    if (! m_SSPhase[kspec]) {
        if (Vphase->exists() == VCS_PHASE_EXIST_NO) {
            Vphase->setExistence(VCS_PHASE_EXIST_YES);
            for (size_t k = 0; k < m_nsp; k++) {
                if (m_phaseID[k] == iph && m_speciesStatus[k] != VCS_SPECIES_DELETED) {
                    m_speciesStatus[k] = VCS_SPECIES_MINOR;
                }
            }
        }
    } else {
        Vphase->setExistence(VCS_PHASE_EXIST_YES);
    }

    ++m_numRxnRdc;
    ++m_numSpeciesRdc;
    ++m_numRxnMinorZeroed;

    if (kspec != (m_numSpeciesRdc - 1)) {
        // Rearrange both the species and the non-component global data
        vcs_switch_pos(true, (m_numSpeciesRdc - 1), kspec);
    }
}

bool VCS_SOLVE::vcs_delete_multiphase(const size_t iph)
{
    vcs_VolPhase* Vphase = m_VolPhaseList[iph].get();
    bool successful = true;

    // set the phase existence flag to dead
    Vphase->setTotalMoles(0.0);
    if (m_debug_print_lvl >= 2) {
        plogf("   --- delete_multiphase %d, %s\n", iph, Vphase->PhaseName);
    }

    // Loop over all of the species in the phase.
    for (size_t kspec = m_numComponents; kspec < m_numSpeciesRdc; ++kspec) {
        if (m_phaseID[kspec] == iph) {
            if (m_speciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
                // calculate an extent of rxn, dx, that zeroes out the species.
                double dx = - m_molNumSpecies_old[kspec];
                double dxTent = dx;

                int retn = delta_species(kspec, &dxTent);
                if (retn != 1) {
                    successful = false;
                    if (m_debug_print_lvl >= 2) {
                        plogf("   --- delete_multiphase %d, %s ERROR problems deleting species %s\n",
                              iph, Vphase->PhaseName, m_speciesName[kspec]);
                        plogf("   ---     delta  attempted: %g  achieved: %g   "
                              "  Zeroing it manually\n", dx, dxTent);
                    }
                    m_molNumSpecies_old[kspec] = 0.0;
                    m_molNumSpecies_new[kspec] = 0.0;
                    m_deltaMolNumSpecies[kspec] = 0.0;
                    // recover the total phase moles.
                    vcs_tmoles();
                } else {
                    // Set the mole number of that species to zero.
                    m_molNumSpecies_old[kspec] = 0.0;
                    m_molNumSpecies_new[kspec] = 0.0;
                    m_deltaMolNumSpecies[kspec] = 0.0;
                }

                // Change the status flag of the species to that of an zeroed
                // phase
                m_speciesStatus[kspec] = VCS_SPECIES_ZEROEDMS;
            }
        }
    }

    double dxPerm = 0.0, dxPerm2 = 0.0;
    for (size_t kcomp = 0; kcomp < m_numComponents; ++kcomp) {
        if (m_phaseID[kcomp] == iph) {
            if (m_debug_print_lvl >= 2) {
                plogf("   --- delete_multiphase   One of the species is a component %d - %s with mole number %g\n",
                      kcomp, m_speciesName[kcomp], m_molNumSpecies_old[kcomp]);
            }
            if (m_molNumSpecies_old[kcomp] != 0.0) {
                for (size_t kspec = m_numComponents; kspec < m_numSpeciesRdc; ++kspec) {
                    size_t irxn = kspec - m_numComponents;
                    if (m_phaseID[kspec] != iph) {
                        if (m_stoichCoeffRxnMatrix(kcomp,irxn) != 0.0) {
                            double dxWant = -m_molNumSpecies_old[kcomp] / m_stoichCoeffRxnMatrix(kcomp,irxn);
                            if (dxWant + m_molNumSpecies_old[kspec]  < 0.0) {
                                dxPerm = -m_molNumSpecies_old[kspec];
                            }
                            for (size_t jcomp = 0; kcomp < m_numComponents; ++kcomp) {
                                if (jcomp != kcomp) {
                                    if (m_phaseID[jcomp] == iph) {
                                        dxPerm = 0.0;
                                    } else {
                                        double dj = dxWant * m_stoichCoeffRxnMatrix(jcomp,irxn);
                                        if (dj + m_molNumSpecies_old[kcomp]  < 0.0) {
                                            dxPerm2 = -m_molNumSpecies_old[kcomp] / m_stoichCoeffRxnMatrix(jcomp,irxn);
                                        }
                                        if (fabs(dxPerm2) < fabs(dxPerm)) {
                                            dxPerm = dxPerm2;
                                        }
                                    }
                                }
                            }
                        }
                        if (dxPerm != 0.0) {
                            delta_species(kspec, &dxPerm);
                        }
                    }
                }

            }
            if (m_molNumSpecies_old[kcomp] != 0.0) {
                if (m_debug_print_lvl >= 2) {
                    plogf("   --- delete_multiphase   One of the species is a component %d - %s still with mole number %g\n",
                          kcomp, m_speciesName[kcomp], m_molNumSpecies_old[kcomp]);
                    plogf("   ---                     zeroing it \n");
                }
                m_molNumSpecies_old[kcomp] = 0.0;
            }
            m_speciesStatus[kcomp] = VCS_SPECIES_ZEROEDMS;
        }
    }

    // Loop over all of the inactive species in the phase: Right now we
    // reinstate all species in a deleted multiphase. We may only want to
    // reinstate the "major ones" in the future. Note, species in phases with
    // zero mole numbers are still considered active. Whether the phase pops
    // back into existence or not is checked as part of the main iteration loop.
    for (size_t kspec = m_numSpeciesRdc; kspec < m_nsp; ++kspec) {
        if (m_phaseID[kspec] == iph) {
            m_molNumSpecies_old[kspec] = 0.0;
            m_molNumSpecies_new[kspec] = 0.0;
            m_deltaMolNumSpecies[kspec] = 0.0;
            m_speciesStatus[kspec] = VCS_SPECIES_ZEROEDMS;

            ++m_numRxnRdc;
            ++m_numSpeciesRdc;
            if (m_debug_print_lvl >= 2) {
                plogf("   ---    Make %s", m_speciesName[kspec]);
                plogf(" an active but zeroed species because its phase "
                      "was zeroed\n");
            }
            if (kspec != (m_numSpeciesRdc - 1)) {
                // Rearrange both the species and the non-component global data
                vcs_switch_pos(true, (m_numSpeciesRdc - 1), kspec);
            }
        }
    }

    // Zero out the total moles counters for the phase
    m_tPhaseMoles_old[iph] = 0.0;
    m_tPhaseMoles_new[iph] = 0.0;
    m_deltaPhaseMoles[iph] = 0.0;

    // Upload the state to the VP object
    Vphase->setTotalMoles(0.0);
    return successful;
}

int VCS_SOLVE::vcs_recheck_deleted()
{
    vector_fp& xtcutoff = m_TmpPhase;
    if (m_debug_print_lvl >= 2) {
        plogf("   --- Start rechecking deleted species in multispec phases\n");
    }
    if (m_numSpeciesRdc == m_nsp) {
        return 0;
    }

    // Use the standard chemical potentials for the chemical potentials of
    // deleted species. Then, calculate Delta G for for formation reactions.
    // Note: fe[] here includes everything except for the ln(x[i]) term
    for (size_t kspec = m_numSpeciesRdc; kspec < m_nsp; ++kspec) {
        size_t iph = m_phaseID[kspec];
        m_feSpecies_new[kspec] = (m_SSfeSpecies[kspec] + log(m_actCoeffSpecies_old[kspec])
                                  - m_lnMnaughtSpecies[kspec]
                                  + m_chargeSpecies[kspec] * m_Faraday_dim * m_phasePhi[iph]);
    }

    // Recalculate the DeltaG's of the formation reactions for the deleted
    // species in the mechanism
    vcs_deltag(0, true, VCS_STATECALC_NEW);

    for (size_t iph = 0; iph < m_numPhases; iph++) {
        if (m_tPhaseMoles_old[iph] > 0.0) {
            xtcutoff[iph] = log(1.0 / VCS_RELDELETE_SPECIES_CUTOFF);
        } else {
            xtcutoff[iph] = 0.0;
        }
    }

    // We are checking the equation:
    //
    //     sum_u = sum_j_comp [ sigma_i_j * u_j ]
    //           = u_i_O + log((AC_i * W_i)/m_tPhaseMoles_old)
    //
    // by first evaluating:
    //
    //     DG_i_O = u_i_O - sum_u.
    //
    // Then, if TL is zero, the phase pops into existence if DG_i_O < 0. Also,
    // if the phase exists, then we check to see if the species can have a mole
    // number larger than VCS_DELETE_SPECIES_CUTOFF (default value = 1.0E-32).
    //
    // HKM: This seems to be an inconsistency in the algorithm here that needs
    // correcting. The requirement above may bypass some multiphases which
    // should exist. The real requirement for the phase to exist is:
    //
    //     sum_i_in_phase [ exp(-DG_i_O) ] >= 1.0
    //
    // Thus, we need to amend the code. Also nonideal solutions will tend to
    // complicate matters severely also.
    int npb = 0;
    for (size_t irxn = m_numRxnRdc; irxn < m_numRxnTot; ++irxn) {
        size_t kspec = m_indexRxnToSpecies[irxn];
        size_t iph = m_phaseID[kspec];
        if (m_tPhaseMoles_old[iph] == 0.0) {
            if (m_deltaGRxn_new[irxn] < 0.0) {
                vcs_reinsert_deleted(kspec);
                npb++;
            } else {
                m_molNumSpecies_old[kspec] = 0.0;
            }
        } else if (m_tPhaseMoles_old[iph] > 0.0) {
            if (m_deltaGRxn_new[irxn] < xtcutoff[iph]) {
                vcs_reinsert_deleted(kspec);
                npb++;
            }
        }
    }
    return npb;
}

size_t VCS_SOLVE::vcs_add_all_deleted()
{
    if (m_numSpeciesRdc == m_nsp) {
        return 0;
    }

    // Use the standard chemical potentials for the chemical potentials of
    // deleted species. Then, calculate Delta G for for formation reactions. We
    // are relying here on a old saved value of m_actCoeffSpecies_old[kspec]
    // being sufficiently good. Note, we will recalculate everything at the end
    // of the routine.
    m_molNumSpecies_new = m_molNumSpecies_old;
    for (int cits = 0; cits < 3; cits++) {
        for (size_t kspec = m_numSpeciesRdc; kspec < m_nsp; ++kspec) {
            size_t iph = m_phaseID[kspec];
            vcs_VolPhase* Vphase = m_VolPhaseList[iph].get();
            if (m_molNumSpecies_new[kspec] == 0.0) {
                m_molNumSpecies_new[kspec] = VCS_DELETE_MINORSPECIES_CUTOFF * 1.0E-10;
            }
            if (!Vphase->m_singleSpecies) {
                Vphase->sendToVCS_ActCoeff(VCS_STATECALC_NEW, &m_actCoeffSpecies_new[0]);
            }
            m_feSpecies_new[kspec] = (m_SSfeSpecies[kspec] + log(m_actCoeffSpecies_new[kspec]) - m_lnMnaughtSpecies[kspec]
                                      + m_chargeSpecies[kspec] * m_Faraday_dim * m_phasePhi[iph]);
        }

        // Recalculate the DeltaG's of the formation reactions for the deleted
        // species in the mechanism
        vcs_deltag(0, true, VCS_STATECALC_NEW);
        for (size_t irxn = m_numRxnRdc; irxn < m_numRxnTot; ++irxn) {
            size_t kspec = m_indexRxnToSpecies[irxn];
            size_t iph = m_phaseID[kspec];
            if (m_tPhaseMoles_old[iph] > 0.0) {
                double maxDG = std::min(m_deltaGRxn_new[irxn], 690.0);
                double dx = m_tPhaseMoles_old[iph] * exp(- maxDG);
                m_molNumSpecies_new[kspec] = dx;
            }
        }
    }
    for (size_t irxn = m_numRxnRdc; irxn < m_numRxnTot; ++irxn) {
        size_t kspec = m_indexRxnToSpecies[irxn];
        size_t iph = m_phaseID[kspec];
        if (m_tPhaseMoles_old[iph] > 0.0) {
            double dx = m_molNumSpecies_new[kspec];
            size_t retn = delta_species(kspec, &dx);
            if (retn == 0) {
                if (m_debug_print_lvl) {
                    plogf("  --- add_deleted(): delta_species() failed for species %s (%d) with mol number %g\n",
                          m_speciesName[kspec], kspec, dx);
                }
                if (dx > 1.0E-50) {
                    dx = 1.0E-50;
                    retn = delta_species(kspec, &dx);
                    if (retn == 0 && m_debug_print_lvl) {
                        plogf("  --- add_deleted(): delta_species() failed for species %s (%d) with mol number %g\n",
                              m_speciesName[kspec], kspec, dx);
                    }
                }
            }
            if (m_debug_print_lvl >= 2) {
                if (retn != 0) {
                    plogf("  --- add_deleted():  species %s added back in with mol number %g\n",
                          m_speciesName[kspec], dx);
                } else {
                    plogf("  --- add_deleted():  species %s failed to be added  back in\n");
                }
            }
        }
    }
    vcs_setFlagsVolPhases(false, VCS_STATECALC_OLD);
    vcs_dfe(VCS_STATECALC_OLD, 0, 0, m_nsp);
    vcs_deltag(0, true, VCS_STATECALC_OLD);

    size_t retn = 0;
    for (size_t irxn = m_numRxnRdc; irxn < m_numRxnTot; ++irxn) {
        size_t kspec = m_indexRxnToSpecies[irxn];
        size_t iph = m_phaseID[kspec];
        if (m_tPhaseMoles_old[iph] > 0.0 && fabs(m_deltaGRxn_old[irxn]) > m_tolmin) {
            if (((m_molNumSpecies_old[kspec] * exp(-m_deltaGRxn_old[irxn])) >
                    VCS_DELETE_MINORSPECIES_CUTOFF) ||
                    (m_molNumSpecies_old[kspec] > VCS_DELETE_MINORSPECIES_CUTOFF)) {
                retn++;
                if (m_debug_print_lvl >= 2) {
                    plogf("  --- add_deleted():  species %s with mol number %g not converged: DG = %g\n",
                          m_speciesName[kspec], m_molNumSpecies_old[kspec],
                          m_deltaGRxn_old[irxn]);
                }
            }
        }
    }
    return retn;
}

bool VCS_SOLVE::vcs_globStepDamp()
{
    double* dptr = &m_deltaGRxn_new[0];

    // CALCULATE SLOPE AT END OF THE STEP
    double s2 = 0.0;
    for (size_t irxn = 0; irxn < m_numRxnRdc; ++irxn) {
        size_t kspec = irxn + m_numComponents;
        if (m_speciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
            s2 += dptr[irxn] * m_deltaMolNumSpecies[kspec];
        }
    }

    // CALCULATE ORIGINAL SLOPE
    double s1 = 0.0;
    dptr = &m_deltaGRxn_old[0];
    for (size_t irxn = 0; irxn < m_numRxnRdc; ++irxn) {
        size_t kspec = irxn + m_numComponents;
        if (m_speciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
            s1 += dptr[irxn] * m_deltaMolNumSpecies[kspec];
        }
    }

    if (m_debug_print_lvl >= 2) {
        plogf("   --- subroutine FORCE: Beginning Slope = %g\n", s1);
        plogf("   --- subroutine FORCE: End Slope       = %g\n", s2);
    }

    if (s1 > 0.0) {
        if (m_debug_print_lvl >= 2) {
            plogf("   --- subroutine FORCE produced no adjustments,");
            if (s1 < 1.0E-40) {
                plogf(" s1 positive but really small\n");
            } else {
                plogf(" failed s1 test\n");
            }
        }
        return false;
    }

    if (s2 <= 0.0) {
        debuglog("   --- subroutine FORCE produced no adjustments, s2 < 0\n", m_debug_print_lvl >= 2);
        return false;
    }

    // FIT PARABOLA
    double al = 1.0;
    if (fabs(s1 -s2) > 1.0E-200) {
        al = s1 / (s1 - s2);
    }
    if (al >= 0.95 || al < 0.0) {
        if (m_debug_print_lvl >= 2) {
            plogf("   --- subroutine FORCE produced no adjustments (al = %g)\n", al);
        }
        return false;
    }
    if (m_debug_print_lvl >= 2) {
        plogf("   --- subroutine FORCE produced a damping factor = %g\n", al);
    }

    // ADJUST MOLE NUMBERS, CHEM. POT
    if (m_debug_print_lvl >= 2) {
        m_deltaGRxn_tmp = m_deltaGRxn_new;
    }

    dptr = &m_molNumSpecies_new[0];
    for (size_t kspec = 0; kspec < m_numSpeciesRdc; ++kspec) {
        m_molNumSpecies_new[kspec] = m_molNumSpecies_old[kspec] +
                                     al * m_deltaMolNumSpecies[kspec];
    }
    for (size_t iph = 0; iph < m_numPhases; iph++) {
        m_tPhaseMoles_new[iph] = m_tPhaseMoles_old[iph] + al * m_deltaPhaseMoles[iph];
    }
    vcs_updateVP(VCS_STATECALC_NEW);

    if (m_debug_print_lvl >= 2) {
        plogf("   --- subroutine FORCE adjusted the mole "
              "numbers, AL = %10.3f\n", al);
    }

    // Because we changed the mole numbers, we need to calculate the chemical
    // potentials again. If a major-only step is being carried out, then we
    // don't need to update the minor noncomponents.
    vcs_setFlagsVolPhases(false, VCS_STATECALC_NEW);
    vcs_dfe(VCS_STATECALC_NEW, 0, 0, m_numSpeciesRdc);

    // Evaluate DeltaG for all components if ITI=0, and for major components
    // only if ITI NE 0
    vcs_deltag(0, false, VCS_STATECALC_NEW);

    dptr = &m_deltaGRxn_new[0];
    s2 = 0.0;
    for (size_t irxn = 0; irxn < m_numRxnRdc; ++irxn) {
        size_t kspec = irxn + m_numComponents;
        if (m_speciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
            s2 += dptr[irxn] * m_deltaMolNumSpecies[kspec];
        }
    }

    if (m_debug_print_lvl >= 2) {
        plogf("   --- subroutine FORCE: Adj End Slope   = %g\n", s2);
    }
    return true;
}

int VCS_SOLVE::vcs_basopt(const bool doJustComponents, double aw[], double sa[], double sm[],
                          double ss[], double test, bool* const usedZeroedSpecies)
{
    size_t k;
    size_t juse = npos;
    size_t jlose = npos;
    DenseMatrix C;
    clockWC tickTock;
    if (m_debug_print_lvl >= 2) {
        plogf("   ");
        for (size_t i=0; i<77; i++) {
            plogf("-");
        }
        plogf("\n");
        plogf("   --- Subroutine BASOPT called to ");
        if (doJustComponents) {
            plogf("calculate the number of components\n");
        } else {
            plogf("reevaluate the components\n");
        }
        if (m_debug_print_lvl >= 2) {
            plogf("\n");
            plogf("   ---       Formula Matrix used in BASOPT calculation\n");
            plogf("   ---      Active |   ");
            for (size_t j = 0; j < m_nelem; j++) {
                plogf("  %1d  ", m_elementActive[j]);
            }
            plogf("\n");
            plogf("   ---     Species | ");
            for (size_t j = 0; j < m_nelem; j++) {
                writelog(" {:>8.8s}", m_elementName[j]);
            }
            plogf("\n");
            for (k = 0; k < m_nsp; k++) {
                writelog("   --- {:>11.11s} | ", m_speciesName[k]);
                for (size_t j = 0; j < m_nelem; j++) {
                    plogf(" %8.2g", m_formulaMatrix(k,j));
                }
                plogf("\n");
            }
            writelogendl();
        }
    }

    // Calculate the maximum value of the number of components possible. It's
    // equal to the minimum of the number of elements and the number of total
    // species.
    size_t ncTrial = std::min(m_nelem, m_nsp);
    m_numComponents = ncTrial;
    *usedZeroedSpecies = false;
    vector_int ipiv(ncTrial);

    // Use a temporary work array for the mole numbers, aw[]
    std::copy(m_molNumSpecies_old.begin(),
              m_molNumSpecies_old.begin() + m_nsp, aw);

    // Take out the Voltage unknowns from consideration
    for (k = 0; k < m_nsp; k++) {
        if (m_speciesUnknownType[k] == VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
            aw[k] = test;
        }
    }

    size_t jr = 0;

    // Top of a loop of some sort based on the index JR. JR is the current
    // number of component species found.
    while (jr < ncTrial) {
        // Top of another loop point based on finding a linearly independent
        // species
        while (true) {
            // Search the remaining part of the mole fraction vector, AW, for
            // the largest remaining species. Return its identity in K. The
            // first search criteria is always the largest positive magnitude of
            // the mole number.
            k = vcs_basisOptMax(aw, jr, m_nsp);

            // The fun really starts when you have run out of species that have
            // a significant concentration. It becomes extremely important to
            // make a good choice of which species you want to pick to fill out
            // the basis. Basically, you don't want to use species with elements
            // abundances which aren't pegged to zero. This means that those
            // modes will never be allowed to grow. You want to have the best
            // chance that the component will grow positively.
            //
            // Suppose you start with CH4, N2, as the only species with nonzero
            // compositions. You have the following abundances:
            //
            //  Abundances:
            // ----------------
            //     C 2.0
            //     N 2.0
            //     H 4.0
            //     O 0.0
            //
            // For example, Make the following choice:
            //
            //   CH4   N2  O choose ->  OH
            // or
            //   CH4   N2  O choose ->  H2
            //
            // OH and H2 both fill out the basis. They will pass the algorithm.
            // However, choosing OH as the next species will create a situation
            // where H2 can not grow in concentration. This happened in
            // practice, btw. The reason is that the formation reaction for H2
            // will cause one of the component species to go negative.
            //
            // The basic idea here is to pick a simple species whose mole number
            // can grow according to the element compositions. Candidates are
            // still filtered according to their linear independence.
            //
            // Note, if there is electronic charge and the electron species, you
            // should probably pick the electron as a component, if it linearly
            // independent. The algorithm below will do this automagically.
            if ((aw[k] != test) && aw[k] < VCS_DELETE_MINORSPECIES_CUTOFF) {
                *usedZeroedSpecies = true;
                double maxConcPossKspec = 0.0;
                double maxConcPoss = 0.0;
                size_t kfound = npos;
                int minNonZeroes = 100000;
                int nonZeroesKspec = 0;
                for (size_t kspec = ncTrial; kspec < m_nsp; kspec++) {
                    if (aw[kspec] >= 0.0 && m_speciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
                        maxConcPossKspec = 1.0E10;
                        nonZeroesKspec = 0;
                        for (size_t j = 0; j < m_nelem; ++j) {
                            if (m_elementActive[j] && m_elType[j] == VCS_ELEM_TYPE_ABSPOS) {
                                double nu = m_formulaMatrix(kspec,j);
                                if (nu != 0.0) {
                                    nonZeroesKspec++;
                                    maxConcPossKspec = std::min(m_elemAbundancesGoal[j] / nu, maxConcPossKspec);
                                }
                            }
                        }
                        if ((maxConcPossKspec >= maxConcPoss) || (maxConcPossKspec > 1.0E-5)) {
                            if (nonZeroesKspec <= minNonZeroes) {
                                if (kfound == npos || nonZeroesKspec < minNonZeroes) {
                                    kfound = kspec;
                                } else {
                                    // ok we are sitting pretty equal here
                                    // decide on the raw ss Gibbs energy
                                    if (m_SSfeSpecies[kspec] <= m_SSfeSpecies[kfound]) {
                                        kfound = kspec;
                                    }
                                }
                            }
                            minNonZeroes = std::min(minNonZeroes, nonZeroesKspec);
                            maxConcPoss = std::max(maxConcPoss, maxConcPossKspec);
                        }
                    }
                }
                if (kfound == npos) {
                    double gmin = 0.0;
                    kfound = k;
                    for (size_t kspec = ncTrial; kspec < m_nsp; kspec++) {
                        if (aw[kspec] >= 0.0) {
                            size_t irxn = kspec - ncTrial;
                            if (m_deltaGRxn_new[irxn] < gmin) {
                                gmin = m_deltaGRxn_new[irxn];
                                kfound = kspec;
                            }
                        }
                    }
                }
                k = kfound;
            }

            if (aw[k] == test) {
                m_numComponents = jr;
                ncTrial = m_numComponents;
                size_t numPreDeleted = m_numRxnTot - m_numRxnRdc;
                if (numPreDeleted != (m_nsp - m_numSpeciesRdc)) {
                    throw CanteraError("VCS_SOLVE::vcs_basopt", "we shouldn't be here");
                }
                m_numRxnTot = m_nsp - ncTrial;
                m_numRxnRdc = m_numRxnTot - numPreDeleted;
                m_numSpeciesRdc = m_nsp - numPreDeleted;
                for (size_t i = 0; i < m_nsp; ++i) {
                    m_indexRxnToSpecies[i] = ncTrial + i;
                }
                if (m_debug_print_lvl >= 2) {
                    plogf("   ---   Total number of components found = %3d (ne = %d)\n ",
                          ncTrial, m_nelem);
                }
                goto L_END_LOOP;
            }

            // Assign a small negative number to the component that we have just
            // found, in order to take it out of further consideration.
            aw[k] = test;

            // CHECK LINEAR INDEPENDENCE WITH PREVIOUS SPECIES
            //
            // Modified Gram-Schmidt Method, p. 202 Dalquist
            // QR factorization of a matrix without row pivoting.
            size_t jl = jr;
            for (size_t j = 0; j < m_nelem; ++j) {
                sm[j + jr*m_nelem] = m_formulaMatrix(k,j);
            }
            if (jl > 0) {
                // Compute the coefficients of JA column of the the upper
                // triangular R matrix, SS(J) = R_J_JR this is slightly
                // different than Dalquist) R_JA_JA = 1
                for (size_t j = 0; j < jl; ++j) {
                    ss[j] = 0.0;
                    for (size_t i = 0; i < m_nelem; ++i) {
                        ss[j] += sm[i + jr*m_nelem] * sm[i + j*m_nelem];
                    }
                    ss[j] /= sa[j];
                }
                // Now make the new column, (*,JR), orthogonal to the previous
                // columns
                for (size_t j = 0; j < jl; ++j) {
                    for (size_t i = 0; i < m_nelem; ++i) {
                        sm[i + jr*m_nelem] -= ss[j] * sm[i + j*m_nelem];
                    }
                }
            }

            // Find the new length of the new column in Q. It will be used in
            // the denominator in future row calcs.
            sa[jr] = 0.0;
            for (size_t ml = 0; ml < m_nelem; ++ml) {
                sa[jr] += pow(sm[ml + jr*m_nelem], 2);
            }

            // IF NORM OF NEW ROW .LT. 1E-3 REJECT
            if (sa[jr] > 1.0e-6) {
              break;
            }
        }

        // REARRANGE THE DATA
        if (jr != k) {
            if (m_debug_print_lvl >= 2) {
                plogf("   ---   %-12.12s", m_speciesName[k]);
                if (m_speciesUnknownType[k] == VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
                   plogf("(Volts = %9.2g)", m_molNumSpecies_old[k]);
                } else {
                   plogf("(%9.2g)", m_molNumSpecies_old[k]);
                }
                plogf(" replaces %-12.12s", m_speciesName[jr]);
                if (m_speciesUnknownType[jr] == VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
                   plogf("(Volts = %9.2g)", m_molNumSpecies_old[jr]);
                } else {
                   plogf("(%9.2g)", m_molNumSpecies_old[jr]);
                }
                plogf(" as component %3d\n", jr);
            }
            vcs_switch_pos(false, jr, k);
            std::swap(aw[jr], aw[k]);
        } else if (m_debug_print_lvl >= 2) {
            plogf("   ---   %-12.12s", m_speciesName[k]);
            if (m_speciesUnknownType[k] == VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
                plogf("(Volts = %9.2g) remains            ", m_molNumSpecies_old[k]);
            } else {
                plogf("(%9.2g) remains            ", m_molNumSpecies_old[k]);
            }
            plogf("              as component %3d\n", jr);
        }
// entry point from up above
L_END_LOOP:
        ;
        // If we haven't found enough components, go back and find some more.
        jr++;
    }

    if (doJustComponents) {
        goto L_CLEANUP;
    }

    // EVALUATE THE STOICHIOMETRY
    //
    // Formulate the matrix problem for the stoichiometric
    // coefficients. CX + B = 0
    //
    // C will be an nc x nc matrix made up of the formula vectors for the
    // components. n RHS's will be solved for. Thus, B is an nc x n matrix.
    //
    // BIG PROBLEM 1/21/99:
    //
    // This algorithm makes the assumption that the first nc rows of the formula
    // matrix aren't rank deficient. However, this might not be the case. For
    // example, assume that the first element in m_formulaMatrix[] is argon.
    // Assume that no species in the matrix problem actually includes argon.
    // Then, the first row in sm[], below will be identically zero. bleh.
    //
    // What needs to be done is to perform a rearrangement of the ELEMENTS ->
    // i.e. rearrange, m_formulaMatrix, sp, and m_elemAbundancesGoal, such that
    // the first nc elements form in combination with the nc components create
    // an invertible sm[]. not a small project, but very doable.
    //
    // An alternative would be to turn the matrix problem below into an ne x nc
    // problem, and do QR elimination instead of Gauss-Jordan elimination. Note
    // the rearrangement of elements need only be done once in the problem. It's
    // actually very similar to the top of this program with ne being the
    // species and nc being the elements!!
    C.resize(ncTrial, ncTrial);
    for (size_t j = 0; j < ncTrial; ++j) {
        for (size_t i = 0; i < ncTrial; ++i) {
            C(i, j) = m_formulaMatrix(j,i);
        }
    }
    for (size_t i = 0; i < m_numRxnTot; ++i) {
        k = m_indexRxnToSpecies[i];
        for (size_t j = 0; j < ncTrial; ++j) {
            m_stoichCoeffRxnMatrix(j,i) = - m_formulaMatrix(k,j);
        }
    }
    // Solve the linear system to calculate the reaction matrix,
    // m_stoichCoeffRxnMatrix.
    solve(C, m_stoichCoeffRxnMatrix.ptrColumn(0), m_numRxnTot, m_nelem);

    // NOW, if we have interfacial voltage unknowns, what we did was just wrong
    // -> hopefully it didn't blow up. Redo the problem. Search for inactive E
    juse = npos;
    jlose = npos;
    for (size_t j = 0; j < m_nelem; j++) {
        if (!m_elementActive[j] && !strcmp(m_elementName[j].c_str(), "E")) {
            juse = j;
        }
    }
    for (size_t j = 0; j < m_nelem; j++) {
        if (m_elementActive[j] && !strncmp((m_elementName[j]).c_str(), "cn_", 3)) {
            jlose = j;
        }
    }
    for (k = 0; k < m_nsp; k++) {
        if (m_speciesUnknownType[k] == VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
            for (size_t j = 0; j < ncTrial; ++j) {
                for (size_t i = 0; i < ncTrial; ++i) {
                    if (i == jlose) {
                        C(i, j) = m_formulaMatrix(j,juse);
                    } else {
                        C(i, j) = m_formulaMatrix(j,i);
                    }
                }
            }
            for (size_t i = 0; i < m_numRxnTot; ++i) {
                k = m_indexRxnToSpecies[i];
                for (size_t j = 0; j < ncTrial; ++j) {
                    if (j == jlose) {
                        aw[j] = - m_formulaMatrix(k,juse);
                    } else {
                        aw[j] = - m_formulaMatrix(k,j);
                    }
                }
            }

            solve(C, aw, 1, m_nelem);
            size_t i = k - ncTrial;
            for (size_t j = 0; j < ncTrial; j++) {
                m_stoichCoeffRxnMatrix(j,i) = aw[j];
            }
        }
    }

    // Calculate the szTmp array for each formation reaction
    for (size_t i = 0; i < m_numRxnTot; i++) {
        double szTmp = 0.0;
        for (size_t j = 0; j < ncTrial; j++) {
            szTmp += fabs(m_stoichCoeffRxnMatrix(j,i));
        }
        m_scSize[i] = szTmp;
    }

    if (m_debug_print_lvl >= 2) {
        plogf("   ---                Components:");
        for (size_t j = 0; j < ncTrial; j++) {
            plogf("        %3d", j);
        }
        plogf("\n   ---          Components Moles:");
        for (size_t j = 0; j < ncTrial; j++) {
            if (m_speciesUnknownType[j] == VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
                plogf(" % -10.3E", 0.0);
            } else {
                plogf(" % -10.3E", m_molNumSpecies_old[j]);
            }
        }
        plogf("\n   ---   NonComponent|   Moles  |");
        for (size_t j = 0; j < ncTrial; j++) {
            plogf(" %10.10s", m_speciesName[j]);
        }
        plogf("\n");
        for (size_t i = 0; i < m_numRxnTot; i++) {
            plogf("   --- %3d ", m_indexRxnToSpecies[i]);
            plogf("%-10.10s", m_speciesName[m_indexRxnToSpecies[i]]);
            if (m_speciesUnknownType[m_indexRxnToSpecies[i]] == VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
                plogf("|% -10.3E|", 0.0);
            } else {
                plogf("|% -10.3E|", m_molNumSpecies_old[m_indexRxnToSpecies[i]]);
            }
            for (size_t j = 0; j < ncTrial; j++) {
                plogf("    %+7.3f", m_stoichCoeffRxnMatrix(j,i));
            }
            plogf("\n");
        }

        // Manual check on the satisfaction of the reaction matrix's ability to
        // conserve elements
        double sumMax = -1.0;
        size_t iMax = npos;
        size_t jMax = npos;
        for (size_t i = 0; i < m_numRxnTot; ++i) {
            k = m_indexRxnToSpecies[i];
            double sum;
            for (size_t j = 0; j < ncTrial; ++j) {
                if (j == jlose) {
                    sum = m_formulaMatrix(k,juse);
                    for (size_t n = 0; n < ncTrial; n++) {
                        double numElements = m_formulaMatrix(n,juse);
                        double coeff = m_stoichCoeffRxnMatrix(n,i);
                        sum += coeff * numElements;
                    }
                } else {
                    sum = m_formulaMatrix(k,j);
                    for (size_t n = 0; n < ncTrial; n++) {
                        double numElements = m_formulaMatrix(n,j);
                        double coeff = m_stoichCoeffRxnMatrix(n,i);
                        sum += coeff * numElements;
                    }
                }
                if (fabs(sum) > sumMax) {
                    sumMax = fabs(sum);
                    iMax = i;
                    jMax = j;
                    if (j == jlose) {
                        jMax = juse;
                    }
                }
                if (fabs(sum) > 1.0E-6) {
                    throw CanteraError("VCS_SOLVE::vcs_basopt", "we have a prob");
                }
            }
        }
        plogf("   ---               largest error in Stoich coeff = %g at rxn = %d ", sumMax, iMax);
        plogf("%-10.10s", m_speciesName[m_indexRxnToSpecies[iMax]]);
        plogf(" element = %d ", jMax);
        plogf("%-5.5s", m_elementName[jMax]);
        plogf("\n");
        plogf("   ");
        for (size_t i=0; i<77; i++) {
            plogf("-");
        }
        plogf("\n");
    }

    // EVALUATE DELTA N VALUES
    //
    // Evaluate the change in gas and liquid total moles due to reaction
    // vectors, DNG and DNL.

    //  Zero out the change of Phase Moles array
    m_deltaMolNumPhase.zero();
    m_phaseParticipation.zero();

    // Loop over each reaction, creating the change in Phase Moles array,
    // m_deltaMolNumPhase(iphase,irxn), and the phase participation array,
    // PhaseParticipation[irxn][iphase]
    for (size_t irxn = 0; irxn < m_numRxnTot; ++irxn) {
        double* scrxn_ptr = m_stoichCoeffRxnMatrix.ptrColumn(irxn);
        size_t kspec = m_indexRxnToSpecies[irxn];
        size_t iph = m_phaseID[kspec];
        m_deltaMolNumPhase(iph,irxn) = 1.0;
        m_phaseParticipation(iph,irxn)++;
        for (size_t j = 0; j < ncTrial; ++j) {
            iph = m_phaseID[j];
            if (fabs(scrxn_ptr[j]) <= 1.0e-6) {
                scrxn_ptr[j] = 0.0;
            } else {
                m_deltaMolNumPhase(iph,irxn) += scrxn_ptr[j];
                m_phaseParticipation(iph,irxn)++;
            }
        }
    }

L_CLEANUP:
    ;
    double tsecond = tickTock.secondsWC();
    m_VCount->Time_basopt += tsecond;
    m_VCount->Basis_Opts++;
    return VCS_SUCCESS;
}

size_t VCS_SOLVE::vcs_basisOptMax(const double* const molNum, const size_t j,
                                  const size_t n)
{
    // The factors of 1.01 and 1.001 are placed in this routine for a purpose.
    // The purpose is to ensure that roundoff errors don't influence major
    // decisions. This means that the optimized and non-optimized versions of
    // the code remain close to each other.
    //
    // (we try to avoid the logic:   a = b
    //                               if (a > b) { do this }
    //                               else       { do something else }
    // because roundoff error makes a difference in the inequality evaluation)
    //
    // Mole numbers are frequently equal to each other in equilibrium problems
    // due to constraints. Swaps are only done if there are a 1% difference in
    // the mole numbers. Of course this logic isn't foolproof.
    size_t largest = j;
    double big = molNum[j] * m_spSize[j] * 1.01;
    if (m_spSize[j] <= 0.0) {
        throw CanteraError("VCS_SOLVE::vcs_basisOptMax",
                           "spSize is nonpositive");
    }
    for (size_t i = j + 1; i < n; ++i) {
        if (m_spSize[i] <= 0.0) {
            throw CanteraError("VCS_SOLVE::vcs_basisOptMax",
                               "spSize is nonpositive");
        }
        bool doSwap = false;
        if (m_SSPhase[j]) {
            doSwap = (molNum[i] * m_spSize[i]) > big;
            if (!m_SSPhase[i] && doSwap) {
                doSwap = molNum[i] > (molNum[largest] * 1.001);
            }
        } else {
            if (m_SSPhase[i]) {
                doSwap = (molNum[i] * m_spSize[i]) > big;
                if (!doSwap) {
                    doSwap = molNum[i] > (molNum[largest] * 1.001);
                }
            } else {
                doSwap = (molNum[i] * m_spSize[i]) > big;
            }
        }
        if (doSwap) {
            largest = i;
            big = molNum[i] * m_spSize[i] * 1.01;
        }
    }
    return largest;
}

int VCS_SOLVE::vcs_species_type(const size_t kspec) const
{
    // Treat special cases first
    if (m_speciesUnknownType[kspec] == VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
        return VCS_SPECIES_INTERFACIALVOLTAGE;
    }

    size_t iph = m_phaseID[kspec];
    int irxn = int(kspec) - int(m_numComponents);
    int phaseExist = m_VolPhaseList[iph]->exists();

    // Treat zeroed out species first
    if (m_molNumSpecies_old[kspec] <= 0.0) {
        if (m_tPhaseMoles_old[iph] <= 0.0 && !m_SSPhase[kspec]) {
            return VCS_SPECIES_ZEROEDMS;
        }

        // see if the species has an element which is so low that species will
        // always be zero
        for (size_t j = 0; j < m_nelem; ++j) {
            if (m_elType[j] == VCS_ELEM_TYPE_ABSPOS) {
                double atomComp = m_formulaMatrix(kspec,j);
                if (atomComp > 0.0) {
                    double maxPermissible = m_elemAbundancesGoal[j] / atomComp;
                    if (maxPermissible < VCS_DELETE_MINORSPECIES_CUTOFF) {
                        if (m_debug_print_lvl >= 2) {
                            plogf("   ---   %s can not be nonzero because"
                                  " needed element %s is zero\n",
                                  m_speciesName[kspec], m_elementName[j]);
                        }
                        if (m_SSPhase[kspec]) {
                            return VCS_SPECIES_ZEROEDSS;
                        } else {
                            return VCS_SPECIES_STOICHZERO;
                        }
                    }
                }
            }
        }

        // The Gibbs free energy for this species is such that it will pop back
        // into existence. An exception to this is if a needed regular element
        // is also zeroed out. Then, don't pop the phase or the species back
        // into existence.
        if (irxn >= 0) {
            for (size_t j = 0; j < m_numComponents; ++j) {
                double stoicC = m_stoichCoeffRxnMatrix(j,irxn);
                if (stoicC != 0.0) {
                    double negChangeComp = - stoicC;
                    if (negChangeComp > 0.0) {
                        if (m_molNumSpecies_old[j] < 1.0E-60) {
                            if (m_debug_print_lvl >= 2) {
                                plogf("   ---   %s is prevented from popping into existence because"
                                      " a needed component to be consumed, %s, has a zero mole number\n",
                                      m_speciesName[kspec], m_speciesName[j]);
                            }
                            if (m_SSPhase[kspec]) {
                                return VCS_SPECIES_ZEROEDSS;
                            } else {
                                return VCS_SPECIES_STOICHZERO;
                            }
                        }
                    } else if (negChangeComp < 0.0) {
                        if (m_VolPhaseList[m_phaseID[j]]->exists() <= 0) {
                            if (m_debug_print_lvl >= 2) {
                                plogf("   ---   %s is prevented from popping into existence because"
                                      " a needed component %s is in a zeroed-phase that would be "
                                      "popped into existence at the same time\n",
                                      m_speciesName[kspec], m_speciesName[j]);
                            }
                            if (m_SSPhase[kspec]) {
                                return VCS_SPECIES_ZEROEDSS;
                            } else {
                                return VCS_SPECIES_STOICHZERO;
                            }
                        }
                    }
                }
            }
        }

        if (irxn >= 0 && m_deltaGRxn_old[irxn] >= 0.0) {
            // We are here when the species is or should remain zeroed out
            if (m_SSPhase[kspec]) {
                return VCS_SPECIES_ZEROEDSS;
            } else {
                if (phaseExist >= VCS_PHASE_EXIST_YES) {
                    return VCS_SPECIES_ACTIVEBUTZERO;
                } else if (phaseExist == VCS_PHASE_EXIST_ZEROEDPHASE) {
                    return VCS_SPECIES_ZEROEDPHASE;
                } else {
                    return VCS_SPECIES_ZEROEDMS;
                }
            }
        }

        // If the current phase already exists,
        if (m_tPhaseMoles_old[iph] > 0.0) {
            if (m_SSPhase[kspec]) {
                return VCS_SPECIES_MAJOR;
            } else {
                return VCS_SPECIES_ACTIVEBUTZERO;
            }
        }

        // The Gibbs free energy for this species is such that it will pop back
        // into existence. Set it to a major species in anticipation. Note, if
        // we had an estimate for the emerging mole fraction of the species in
        // the phase, we could do better here.
        if (m_tPhaseMoles_old[iph] <= 0.0) {
            if (m_SSPhase[kspec]) {
                return VCS_SPECIES_MAJOR;
            } else {
                return VCS_SPECIES_ZEROEDMS;
            }
        }
    }

    // Treat species with non-zero mole numbers next

    // Always treat species in single species phases as majors if the phase
    // exists.
    if (m_SSPhase[kspec]) {
        return VCS_SPECIES_MAJOR;
    }

    // Check to see whether the current species is a major component of its
    // phase. If it is, it is a major component. This is consistent with the
    // above rule about single species phases. A major component i.e., a species
    // with a high mole fraction) in any phase is always treated as a major
    // species
    if (m_molNumSpecies_old[kspec] > (m_tPhaseMoles_old[iph] * 0.001)) {
        return VCS_SPECIES_MAJOR;
    }

    // Main check in the loop: Check to see if there is a component with a mole
    // number that is within a factor of 100 of the current species. If there is
    // and that component is not part of a single species phase and shares a
    // non-zero stoichiometric coefficient, then the current species is a major
    // species.
    if (irxn < 0) {
        return VCS_SPECIES_MAJOR;
    } else {
        double szAdj = m_scSize[irxn] * std::sqrt((double)m_numRxnTot);
        for (size_t k = 0; k < m_numComponents; ++k) {
            if (!m_SSPhase[k] && m_stoichCoeffRxnMatrix(k,irxn) != 0.0 && m_molNumSpecies_old[kspec] * szAdj >= m_molNumSpecies_old[k] * 0.01) {
                return VCS_SPECIES_MAJOR;
            }
        }
    }
    return VCS_SPECIES_MINOR;
}

void VCS_SOLVE::vcs_dfe(const int stateCalc,
                        const int ll, const size_t lbot, const size_t ltop)
{
    double* tPhMoles_ptr=0;
    double* actCoeff_ptr=0;
    double* feSpecies=0;
    double* molNum=0;
    if (stateCalc == VCS_STATECALC_OLD) {
        feSpecies = &m_feSpecies_old[0];
        tPhMoles_ptr = &m_tPhaseMoles_old[0];
        actCoeff_ptr = &m_actCoeffSpecies_old[0];
        molNum = &m_molNumSpecies_old[0];
    } else if (stateCalc == VCS_STATECALC_NEW) {
        feSpecies = &m_feSpecies_new[0];
        tPhMoles_ptr = &m_tPhaseMoles_new[0];
        actCoeff_ptr = &m_actCoeffSpecies_new[0];
        molNum = &m_molNumSpecies_new[0];
    } else {
        throw CanteraError("VCS_SOLVE::vcs_dfe",
                "Subroutine vcs_dfe called with bad stateCalc value: {}", stateCalc);
    }

    if (m_debug_print_lvl >= 2) {
        if (ll == 0) {
            if (lbot != 0) {
                plogf("   --- Subroutine vcs_dfe called for one species: ");
                plogf("%-12.12s", m_speciesName[lbot]);
            } else {
                plogf("   --- Subroutine vcs_dfe called for all species");
            }
        } else if (ll > 0) {
            plogf("   --- Subroutine vcs_dfe called for components and minors");
        } else {
            plogf("   --- Subroutine vcs_dfe called for components and majors");
        }
        if (stateCalc == VCS_STATECALC_NEW) {
            plogf(" using tentative solution");
        }
        writelogendl();
    }

    double* tlogMoles = &m_TmpPhase[0];

    // Might as well recalculate the phase mole vector and compare to the stored
    // one. They should be correct.
    double* tPhInertMoles = &TPhInertMoles[0];
    for (size_t iph = 0; iph < m_numPhases; iph++) {
        tlogMoles[iph] = tPhInertMoles[iph];

    }
    for (size_t kspec = 0; kspec < m_nsp; kspec++) {
        if (m_speciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
            size_t iph = m_phaseID[kspec];
            tlogMoles[iph] += molNum[kspec];
        }
    }
    for (size_t iph = 0; iph < m_numPhases; iph++) {
        AssertThrowMsg(vcs_doubleEqual(tlogMoles[iph], tPhMoles_ptr[iph]),
            "VCS_SOLVE::vcs_dfe",
            "phase Moles may be off, iph = {}, {} {}",
            iph, tlogMoles[iph], tPhMoles_ptr[iph]);
    }
    m_TmpPhase.assign(m_TmpPhase.size(), 0.0);
    for (size_t iph = 0; iph < m_numPhases; iph++) {
        if (tPhMoles_ptr[iph] > 0.0) {
            tlogMoles[iph] = log(tPhMoles_ptr[iph]);
        }
    }

    size_t l1, l2;
    if (ll != 0) {
        l1 = lbot;
        l2 = m_numComponents;
    } else {
        l1 = lbot;
        l2 = ltop;
    }

    // Calculate activity coefficients for all phases that are not current. Here
    // we also trigger an update check for each VolPhase to see if its mole
    // numbers are current with vcs
    for (size_t iphase = 0; iphase < m_numPhases; iphase++) {
        vcs_VolPhase* Vphase = m_VolPhaseList[iphase].get();
        Vphase->updateFromVCS_MoleNumbers(stateCalc);
        if (!Vphase->m_singleSpecies) {
            Vphase->sendToVCS_ActCoeff(stateCalc, &actCoeff_ptr[0]);
        }
        m_phasePhi[iphase] = Vphase->electricPotential();
    }

    // ALL SPECIES, OR COMPONENTS
    //
    // Do all of the species when LL = 0. Then we are done for the routine When
    // LL ne 0., just do the initial components. We will then finish up below
    // with loops over either the major noncomponent species or the minor
    // noncomponent species.
    for (size_t kspec = l1; kspec < l2; ++kspec) {
        size_t iphase = m_phaseID[kspec];
        if (m_speciesUnknownType[kspec] == VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
            AssertThrowMsg(molNum[kspec] == m_phasePhi[iphase], "VCS_SOLVE::vcs_dfe",
                "We have an inconsistency!");
            AssertThrowMsg(m_chargeSpecies[kspec] == -1.0, "VCS_SOLVE::vcs_dfe",
                "We have an unexpected situation!");
            feSpecies[kspec] = m_SSfeSpecies[kspec]
                               + m_chargeSpecies[kspec] * m_Faraday_dim * m_phasePhi[iphase];
        } else {
            if (m_SSPhase[kspec]) {
                feSpecies[kspec] = m_SSfeSpecies[kspec]
                                   + m_chargeSpecies[kspec] * m_Faraday_dim * m_phasePhi[iphase];
            } else if ((m_speciesStatus[kspec] == VCS_SPECIES_ZEROEDMS) ||
                       (m_speciesStatus[kspec] == VCS_SPECIES_ZEROEDPHASE)) {
                feSpecies[kspec] = m_SSfeSpecies[kspec] - m_lnMnaughtSpecies[kspec]
                                   + m_chargeSpecies[kspec] * m_Faraday_dim * m_phasePhi[iphase];
            } else {
                if (molNum[kspec] <= VCS_DELETE_MINORSPECIES_CUTOFF) {
                    size_t iph = m_phaseID[kspec];
                    if (tPhMoles_ptr[iph] > 0.0) {
                        feSpecies[kspec] = m_SSfeSpecies[kspec]
                                           + log(actCoeff_ptr[kspec] * VCS_DELETE_MINORSPECIES_CUTOFF)
                                           - tlogMoles[m_phaseID[kspec]] - m_lnMnaughtSpecies[kspec]
                                           + m_chargeSpecies[kspec] * m_Faraday_dim * m_phasePhi[iphase];
                    } else {
                        feSpecies[kspec] = m_SSfeSpecies[kspec] - m_lnMnaughtSpecies[kspec]
                                           + m_chargeSpecies[kspec] * m_Faraday_dim * m_phasePhi[iphase];
                    }
                } else {
                    feSpecies[kspec] = m_SSfeSpecies[kspec]
                                       + log(actCoeff_ptr[kspec] * molNum[kspec])
                                       - tlogMoles[m_phaseID[kspec]] - m_lnMnaughtSpecies[kspec]
                                       + m_chargeSpecies[kspec] * m_Faraday_dim * m_phasePhi[iphase];
                }
            }
        }
    }

    if (ll < 0) {
        // MAJORS ONLY
        for (size_t irxn = 0; irxn < m_numRxnRdc; ++irxn) {
            size_t kspec = m_indexRxnToSpecies[irxn];
            if (m_speciesStatus[kspec] != VCS_SPECIES_MINOR) {
                size_t iphase = m_phaseID[kspec];
                if (m_speciesUnknownType[kspec] == VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
                    AssertThrowMsg(molNum[kspec] == m_phasePhi[iphase], "VCS_SOLVE::vcs_dfe",
                                   "We have an inconsistency!");
                    AssertThrowMsg(m_chargeSpecies[kspec] == -1.0, "VCS_SOLVE::vcs_dfe",
                                   "We have an unexpected situation!");
                    feSpecies[kspec] = m_SSfeSpecies[kspec]
                                       + m_chargeSpecies[kspec] * m_Faraday_dim * m_phasePhi[iphase];
                } else {
                    if (m_SSPhase[kspec]) {
                        feSpecies[kspec] = m_SSfeSpecies[kspec]
                                           + m_chargeSpecies[kspec] * m_Faraday_dim * m_phasePhi[iphase];
                    } else if ((m_speciesStatus[kspec] == VCS_SPECIES_ZEROEDMS) ||
                               (m_speciesStatus[kspec] == VCS_SPECIES_ZEROEDPHASE)) {
                        feSpecies[kspec] = m_SSfeSpecies[kspec] - m_lnMnaughtSpecies[kspec]
                                           + m_chargeSpecies[kspec] * m_Faraday_dim * m_phasePhi[iphase];
                    } else {
                        if (molNum[kspec] <= VCS_DELETE_MINORSPECIES_CUTOFF) {
                            size_t iph = m_phaseID[kspec];
                            if (tPhMoles_ptr[iph] > 0.0) {
                                feSpecies[kspec] = m_SSfeSpecies[kspec]
                                                   + log(actCoeff_ptr[kspec] * VCS_DELETE_MINORSPECIES_CUTOFF)
                                                   - tlogMoles[m_phaseID[kspec]] - m_lnMnaughtSpecies[kspec]
                                                   + m_chargeSpecies[kspec] * m_Faraday_dim * m_phasePhi[iphase];
                            } else {
                                feSpecies[kspec] = m_SSfeSpecies[kspec] - m_lnMnaughtSpecies[kspec]
                                                   + m_chargeSpecies[kspec] * m_Faraday_dim * m_phasePhi[iphase];
                            }
                        } else {
                            feSpecies[kspec] = m_SSfeSpecies[kspec]
                                               + log(actCoeff_ptr[kspec] * molNum[kspec])
                                               - tlogMoles[m_phaseID[kspec]] - m_lnMnaughtSpecies[kspec]
                                               + m_chargeSpecies[kspec] * m_Faraday_dim * m_phasePhi[iphase];
                        }
                    }
                }
            }
        }
    } else if (ll > 0) {
        // MINORS ONLY
        for (size_t irxn = 0; irxn < m_numRxnRdc; ++irxn) {
            size_t kspec = m_indexRxnToSpecies[irxn];
            if (m_speciesStatus[kspec] == VCS_SPECIES_MINOR) {
                size_t iphase = m_phaseID[kspec];
                if (m_speciesUnknownType[kspec] == VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
                    AssertThrowMsg(molNum[kspec] == m_phasePhi[iphase], "VCS_SOLVE::vcs_dfe",
                        "We have an inconsistency!");
                    AssertThrowMsg(m_chargeSpecies[kspec] == -1.0, "VCS_SOLVE::vcs_dfe",
                        "We have an unexpected situation!");
                    feSpecies[kspec] = m_SSfeSpecies[kspec]
                                       + m_chargeSpecies[kspec] * m_Faraday_dim * m_phasePhi[iphase];
                } else {
                    if (m_SSPhase[kspec]) {
                        feSpecies[kspec] = m_SSfeSpecies[kspec]
                                           + m_chargeSpecies[kspec] * m_Faraday_dim * m_phasePhi[iphase];
                    } else if ((m_speciesStatus[kspec] == VCS_SPECIES_ZEROEDMS) ||
                               (m_speciesStatus[kspec] == VCS_SPECIES_ZEROEDPHASE)) {
                        feSpecies[kspec] = m_SSfeSpecies[kspec] - m_lnMnaughtSpecies[kspec]
                                           + m_chargeSpecies[kspec] * m_Faraday_dim * m_phasePhi[iphase];
                    } else {
                        if (molNum[kspec] <= VCS_DELETE_MINORSPECIES_CUTOFF) {
                            size_t iph = m_phaseID[kspec];
                            if (tPhMoles_ptr[iph] > 0.0) {
                                feSpecies[kspec] = m_SSfeSpecies[kspec]
                                                   + log(actCoeff_ptr[kspec] * VCS_DELETE_MINORSPECIES_CUTOFF)
                                                   - tlogMoles[m_phaseID[kspec]] - m_lnMnaughtSpecies[kspec]
                                                   + m_chargeSpecies[kspec] * m_Faraday_dim * m_phasePhi[iphase];
                            } else {
                                feSpecies[kspec] = m_SSfeSpecies[kspec] - m_lnMnaughtSpecies[kspec]
                                                   + m_chargeSpecies[kspec] * m_Faraday_dim * m_phasePhi[iphase];
                            }
                        } else {
                            feSpecies[kspec] = m_SSfeSpecies[kspec]
                                               + log(actCoeff_ptr[kspec] * molNum[kspec])
                                               - tlogMoles[m_phaseID[kspec]] - m_lnMnaughtSpecies[kspec]
                                               + m_chargeSpecies[kspec] * m_Faraday_dim * m_phasePhi[iphase];
                        }
                    }
                }
            }
        }
    }
}

double VCS_SOLVE::l2normdg(double dgLocal[]) const
{
    if (m_numRxnRdc <= 0) {
        return 0.0;
    }
    double tmp = 0;
    for (size_t irxn = 0; irxn < m_numRxnRdc; ++irxn) {
        size_t kspec = irxn + m_numComponents;
        if (m_speciesStatus[kspec] == VCS_SPECIES_MAJOR || m_speciesStatus[kspec] == VCS_SPECIES_MINOR ||
                dgLocal[irxn] < 0.0) {
            if (m_speciesStatus[kspec] != VCS_SPECIES_ZEROEDMS) {
                tmp += dgLocal[irxn] * dgLocal[irxn];
            }
        }
    }
    return std::sqrt(tmp / m_numRxnRdc);
}

double VCS_SOLVE::vcs_tmoles()
{
    for (size_t i = 0; i < m_numPhases; i++) {
        m_tPhaseMoles_old[i] = TPhInertMoles[i];
    }
    for (size_t i = 0; i < m_nsp; i++) {
        if (m_speciesUnknownType[i] == VCS_SPECIES_TYPE_MOLNUM) {
            m_tPhaseMoles_old[m_phaseID[i]] += m_molNumSpecies_old[i];
        }
    }
    double sum = 0.0;
    for (size_t i = 0; i < m_numPhases; i++) {
        sum += m_tPhaseMoles_old[i];
        vcs_VolPhase* Vphase = m_VolPhaseList[i].get();
        if (m_tPhaseMoles_old[i] == 0.0) {
            Vphase->setTotalMoles(0.0);
        } else {
            Vphase->setTotalMoles(m_tPhaseMoles_old[i]);
        }
    }
    m_totalMolNum = sum;
    return m_totalMolNum;
}

void VCS_SOLVE::check_tmoles() const
{
    double sum = 0.0;
    for (size_t i = 0; i < m_numPhases; i++) {
        double m_tPhaseMoles_old_a = TPhInertMoles[i];

        for (size_t k = 0; k < m_nsp; k++) {
            if (m_speciesUnknownType[k] == VCS_SPECIES_TYPE_MOLNUM && m_phaseID[k] == i) {
                m_tPhaseMoles_old_a += m_molNumSpecies_old[k];
            }
        }
        sum += m_tPhaseMoles_old_a;

        double denom = m_tPhaseMoles_old[i]+ m_tPhaseMoles_old_a + 1.0E-19;
        if (!vcs_doubleEqual(m_tPhaseMoles_old[i]/denom, m_tPhaseMoles_old_a/denom)) {
            plogf("check_tmoles: we have found a problem with phase  %d: %20.15g, %20.15g\n",
                  i, m_tPhaseMoles_old[i], m_tPhaseMoles_old_a);
        }
    }
}

void VCS_SOLVE::vcs_updateVP(const int vcsState)
{
    for (size_t i = 0; i < m_numPhases; i++) {
        vcs_VolPhase* Vphase = m_VolPhaseList[i].get();
        if (vcsState == VCS_STATECALC_OLD) {
            Vphase->setMolesFromVCSCheck(VCS_STATECALC_OLD,
                                         &m_molNumSpecies_old[0],
                                         &m_tPhaseMoles_old[0]);
        } else if (vcsState == VCS_STATECALC_NEW) {
            Vphase->setMolesFromVCSCheck(VCS_STATECALC_NEW,
                                         &m_molNumSpecies_new[0],
                                         &m_tPhaseMoles_new[0]);
        } else {
            throw CanteraError("VCS_SOLVE::vcs_updateVP",
                               "wrong stateCalc value: {}", vcsState);
        }
    }
}

bool VCS_SOLVE::vcs_evaluate_speciesType()
{
    m_numRxnMinorZeroed = 0;
    if (m_debug_print_lvl >= 2) {
        plogf("  --- Species Status decision is reevaluated: All species are minor except for:\n");
    } else if (m_debug_print_lvl >= 5) {
        plogf("  --- Species Status decision is reevaluated\n");
    }
    for (size_t kspec = 0; kspec < m_nsp; ++kspec) {
        m_speciesStatus[kspec] = vcs_species_type(kspec);
        if (m_debug_print_lvl >= 5) {
            plogf("  --- %-16s: ", m_speciesName[kspec]);
            if (kspec < m_numComponents) {
                plogf("(COMP) ");
            } else {
                plogf("       ");
            }
            plogf(" %10.3g ", m_molNumSpecies_old[kspec]);
            const char* sString = vcs_speciesType_string(m_speciesStatus[kspec], 100);
            plogf("%s\n", sString);
        } else if (m_debug_print_lvl >= 2) {
            if (m_speciesStatus[kspec] != VCS_SPECIES_MINOR) {
                switch (m_speciesStatus[kspec]) {
                case VCS_SPECIES_COMPONENT:
                    break;
                case VCS_SPECIES_MAJOR:
                    plogf("  ---      Major Species          : %-s\n", m_speciesName[kspec]);
                    break;
                case VCS_SPECIES_ZEROEDPHASE:
                    plogf("  ---      Purposely Zeroed-Phase Species (not in problem): %-s\n",
                          m_speciesName[kspec]);
                    break;
                case VCS_SPECIES_ZEROEDMS:
                    plogf("  ---      Zeroed-MS Phase Species: %-s\n", m_speciesName[kspec]);
                    break;
                case VCS_SPECIES_ZEROEDSS:
                    plogf("  ---      Zeroed-SS Phase Species: %-s\n", m_speciesName[kspec]);
                    break;
                case VCS_SPECIES_DELETED:
                    plogf("  ---      Deleted-Small Species  : %-s\n", m_speciesName[kspec]);
                    break;
                case VCS_SPECIES_ACTIVEBUTZERO:
                    plogf("  ---      Zeroed Species in an active MS phase (tmp): %-s\n",
                          m_speciesName[kspec]);
                    break;
                case VCS_SPECIES_STOICHZERO:
                    plogf("  ---     Zeroed Species in an active MS phase (Stoich Constraint): %-s\n",
                          m_speciesName[kspec]);
                    break;
                case VCS_SPECIES_INTERFACIALVOLTAGE:
                    plogf("  ---      InterfaceVoltage Species: %-s\n", m_speciesName[kspec]);
                    break;
                default:
                    throw CanteraError("VCS_SOLVE::vcs_evaluate_speciesType",
                                       "Unknown type: {}", m_speciesStatus[kspec]);
                }
            }
        }
        if (kspec >= m_numComponents && m_speciesStatus[kspec] != VCS_SPECIES_MAJOR) {
            ++m_numRxnMinorZeroed;
        }
    }
    debuglog("  ---\n", m_debug_print_lvl >= 2);
    return (m_numRxnMinorZeroed >= m_numRxnRdc);
}

void VCS_SOLVE::vcs_deltag(const int L, const bool doDeleted,
                           const int vcsState, const bool alterZeroedPhases)
{
    size_t irxnl = m_numRxnRdc;
    if (doDeleted) {
        irxnl = m_numRxnTot;
    }

    double* deltaGRxn;
    double* feSpecies;
    double* molNumSpecies;
    double* actCoeffSpecies;
    if (vcsState == VCS_STATECALC_NEW) {
        deltaGRxn = &m_deltaGRxn_new[0];
        feSpecies = &m_feSpecies_new[0];
        molNumSpecies = &m_molNumSpecies_new[0];
        actCoeffSpecies = &m_actCoeffSpecies_new[0];
    } else if (vcsState == VCS_STATECALC_OLD) {
        deltaGRxn = &m_deltaGRxn_old[0];
        feSpecies = &m_feSpecies_old[0];
        molNumSpecies = &m_molNumSpecies_old[0];
        actCoeffSpecies = &m_actCoeffSpecies_old[0];
    } else {
        throw CanteraError("VCS_SOLVE::vcs_deltag", "bad vcsState");
    }

    if (m_debug_print_lvl >= 2) {
        plogf("   --- Subroutine vcs_deltag called for ");
        if (L < 0) {
            plogf("major noncomponents\n");
        } else if (L == 0) {
            plogf("all noncomponents\n");
        } else {
            plogf("minor noncomponents\n");
        }
    }

    if (L < 0) {
        // MAJORS and ZEROED SPECIES ONLY
        for (size_t irxn = 0; irxn < m_numRxnRdc; ++irxn) {
            size_t kspec = irxn + m_numComponents;
            if (m_speciesStatus[kspec] != VCS_SPECIES_MINOR) {
                int icase = 0;
                deltaGRxn[irxn] = feSpecies[m_indexRxnToSpecies[irxn]];
                double* dtmp_ptr = m_stoichCoeffRxnMatrix.ptrColumn(irxn);
                for (kspec = 0; kspec < m_numComponents; ++kspec) {
                    deltaGRxn[irxn] += dtmp_ptr[kspec] * feSpecies[kspec];
                    if (molNumSpecies[kspec] < VCS_DELETE_MINORSPECIES_CUTOFF && dtmp_ptr[kspec] < 0.0) {
                        icase = 1;
                    }
                }
                if (icase) {
                    deltaGRxn[irxn] = std::max(0.0, deltaGRxn[irxn]);
                }
            }
        }
    } else if (L == 0) {
        // ALL REACTIONS
        for (size_t irxn = 0; irxn < irxnl; ++irxn) {
            int icase = 0;
            deltaGRxn[irxn] = feSpecies[m_indexRxnToSpecies[irxn]];
            double* dtmp_ptr = m_stoichCoeffRxnMatrix.ptrColumn(irxn);
            for (size_t kspec = 0; kspec < m_numComponents; ++kspec) {
                deltaGRxn[irxn] += dtmp_ptr[kspec] * feSpecies[kspec];
                if (molNumSpecies[kspec] < VCS_DELETE_MINORSPECIES_CUTOFF &&
                        dtmp_ptr[kspec] < 0.0) {
                    icase = 1;
                }
            }
            if (icase) {
                deltaGRxn[irxn] = std::max(0.0, deltaGRxn[irxn]);
            }
        }
    } else {
        // MINORS AND ZEROED SPECIES
        for (size_t irxn = 0; irxn < m_numRxnRdc; ++irxn) {
            size_t kspec = irxn + m_numComponents;
            if (m_speciesStatus[kspec] <= VCS_SPECIES_MINOR) {
                int icase = 0;
                deltaGRxn[irxn] = feSpecies[m_indexRxnToSpecies[irxn]];
                double* dtmp_ptr = m_stoichCoeffRxnMatrix.ptrColumn(irxn);
                for (kspec = 0; kspec < m_numComponents; ++kspec) {
                    deltaGRxn[irxn] += dtmp_ptr[kspec] * feSpecies[kspec];
                    if (m_molNumSpecies_old[kspec] < VCS_DELETE_MINORSPECIES_CUTOFF &&
                            dtmp_ptr[kspec] < 0.0) {
                        icase = 1;
                    }
                }
                if (icase) {
                    deltaGRxn[irxn] = std::max(0.0, deltaGRxn[irxn]);
                }
            }
        }
    }

    /* **** MULTISPECIES PHASES WITH ZERO MOLES ******** */
    //
    // Massage the free energies for species with zero mole fractions in
    // multispecies phases.  This section implements the Equation 3.8-5 in Smith
    // and Missen, p.59. A multispecies phase will exist iff
    //
    //     1 < sum_i(exp(-dg_i)/AC_i)
    //
    // If DG is negative then that species wants to be reintroduced into the
    // calculation. For small dg_i, the expression below becomes:
    //
    //      1 - sum_i(exp(-dg_i)/AC_i) ~ sum_i((dg_i-1)/AC_i)  + 1
    //
    // So, what we are doing here is equalizing all DG's in a multispecies phase
    // whose total mole number has already been zeroed out. It must have to do
    // with the case where a complete multispecies phase is currently zeroed
    // out. In that case, when one species in that phase has a negative DG, then
    // the phase should kick in. This code section will cause that to happen,
    // because a negative DG will dominate the calculation of SDEL. Then, DG(I)
    // for all species in that phase will be forced to be equal and negative.
    // Thus, all species in that phase will come into being at the same time.
    //
    // HKM -> The ratio of mole fractions at the reinstatement time should be
    // equal to the normalized weighting of exp(-dg_i) / AC_i. This should be
    // implemented.
    //
    // HKM -> There is circular logic here. ActCoeff depends on the mole
    // fractions of a phase that does not exist. In actuality the proto-mole
    // fractions should be selected from the solution of a nonlinear problem
    // with NsPhase unknowns
    //
    //     X_i = exp(-dg[irxn]) / ActCoeff_i / denom
    //
    // where
    //
    //     denom = sum_i[  exp(-dg[irxn]) / ActCoeff_i  ]
    //
    // This can probably be solved by successive iteration. This should be
    // implemented.
    if (alterZeroedPhases && false) {
        for (size_t iph = 0; iph < m_numPhases; iph++) {
            bool lneed = false;
            vcs_VolPhase* Vphase = m_VolPhaseList[iph].get();
            if (! Vphase->m_singleSpecies) {
                double sum = 0.0;
                for (size_t k = 0; k < Vphase->nSpecies(); k++) {
                    size_t kspec = Vphase->spGlobalIndexVCS(k);
                    if (m_speciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
                        sum += molNumSpecies[kspec];
                    }
                    if (sum > 0.0) {
                        break;
                    }
                }
                if (sum == 0.0) {
                    lneed = true;
                }
            }

            if (lneed) {
                double poly = 0.0;
                for (size_t k = 0; k < Vphase->nSpecies(); k++) {
                    size_t kspec = Vphase->spGlobalIndexVCS(k);
                    // We may need to look at deltaGRxn for components!
                    if (kspec >= m_numComponents) {
                        size_t irxn = kspec - m_numComponents;
                        deltaGRxn[irxn] = clip(deltaGRxn[irxn], -50.0, 50.0);
                        poly += exp(-deltaGRxn[irxn])/actCoeffSpecies[kspec];
                    }
                }

                // Calculate deltaGRxn[] for each species in a zeroed
                // multispecies phase. All of the m_deltaGRxn_new[]'s will be
                // equal. If deltaGRxn[] is negative, then the phase will come
                // back into existence.
                for (size_t k = 0; k < Vphase->nSpecies(); k++) {
                    size_t kspec = Vphase->spGlobalIndexVCS(k);
                    if (kspec >= m_numComponents) {
                        size_t irxn = kspec - m_numComponents;
                        deltaGRxn[irxn] = 1.0 - poly;
                    }
                }
            }
        }
    }
}

void VCS_SOLVE::vcs_printDeltaG(const int stateCalc)
{
    double* deltaGRxn = &m_deltaGRxn_old[0];
    double* feSpecies = &m_feSpecies_old[0];
    double* molNumSpecies = &m_molNumSpecies_old[0];
    const double* tPhMoles_ptr = &m_tPhaseMoles_old[0];
    const double* actCoeff_ptr = &m_actCoeffSpecies_old[0];
    if (stateCalc == VCS_STATECALC_NEW) {
        deltaGRxn = &m_deltaGRxn_new[0];
        feSpecies = &m_feSpecies_new[0];
        molNumSpecies = &m_molNumSpecies_new[0];
        actCoeff_ptr = &m_actCoeffSpecies_new[0];
        tPhMoles_ptr = &m_tPhaseMoles_new[0];
    }
    double RT = m_temperature * GasConstant;
    bool zeroedPhase = false;
    if (m_debug_print_lvl >= 2) {
        plogf("   --- DELTA_G TABLE  Components:");
        for (size_t j = 0; j < m_numComponents; j++) {
            plogf("     %3d  ", j);
        }
        plogf("\n   ---          Components Moles:");
        for (size_t j = 0; j < m_numComponents; j++) {
            plogf("%10.3g", m_molNumSpecies_old[j]);
        }
        plogf("\n   ---   NonComponent|   Moles  |       ");
        for (size_t j = 0; j < m_numComponents; j++) {
            plogf("%-10.10s", m_speciesName[j]);
        }
        plogf("\n");
        for (size_t i = 0; i < m_numRxnTot; i++) {
            plogf("   --- %3d ", m_indexRxnToSpecies[i]);
            plogf("%-10.10s", m_speciesName[m_indexRxnToSpecies[i]]);
            if (m_speciesUnknownType[m_indexRxnToSpecies[i]] == VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
                plogf("|   NA     |");
            } else {
                plogf("|%10.3g|", m_molNumSpecies_old[m_indexRxnToSpecies[i]]);
            }
            for (size_t j = 0; j < m_numComponents; j++) {
                plogf("     %6.2f", m_stoichCoeffRxnMatrix(j,i));
            }
            plogf("\n");
        }
        plogf("   ");
        for (int i=0; i<77; i++) {
            plogf("-");
        }
        plogf("\n");
    }

    writelog("   --- DeltaG Table (J/kmol) Name       PhID   MoleNum      MolFR  "
           "  ElectrChemStar ElectrChem    DeltaGStar   DeltaG(Pred)  Stability\n");
    writelog("   ");
    writeline('-', 132);

    for (size_t kspec = 0; kspec < m_nsp; kspec++) {
        size_t irxn = npos;
        if (kspec >= m_numComponents) {
            irxn = kspec - m_numComponents;
        }
        double mfValue = 1.0;
        size_t iphase = m_phaseID[kspec];
        const vcs_VolPhase* Vphase = m_VolPhaseList[iphase].get();
        if ((m_speciesStatus[kspec] == VCS_SPECIES_ZEROEDMS) ||
                (m_speciesStatus[kspec] == VCS_SPECIES_ZEROEDPHASE) ||
                (m_speciesStatus[kspec] == VCS_SPECIES_ZEROEDSS)) {
            zeroedPhase = true;
        } else {
            zeroedPhase = false;
        }
        if (tPhMoles_ptr[iphase] > 0.0) {
            if (molNumSpecies[kspec] <= VCS_DELETE_MINORSPECIES_CUTOFF) {
                mfValue = VCS_DELETE_MINORSPECIES_CUTOFF / tPhMoles_ptr[iphase];
            } else {
                mfValue = molNumSpecies[kspec] / tPhMoles_ptr[iphase];
            }
        } else {
            size_t klocal = m_speciesLocalPhaseIndex[kspec];
            mfValue = Vphase->moleFraction(klocal);
        }
        if (zeroedPhase) {
            writelog("   --- ** zp *** ");
        } else {
            writelog("   ---           ");
        }
        double feFull = feSpecies[kspec];
        if ((m_speciesStatus[kspec] == VCS_SPECIES_ZEROEDMS) ||
                (m_speciesStatus[kspec] == VCS_SPECIES_ZEROEDPHASE)) {
            feFull += log(actCoeff_ptr[kspec]) + log(mfValue);
        }
        writelogf("%-24.24s", m_speciesName[kspec]);
        writelogf(" %3d", iphase);
        if (m_speciesUnknownType[kspec] == VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
            writelog("    NA       ");
        } else {
            writelogf(" % -12.4e", molNumSpecies[kspec]);
        }
        writelogf(" % -12.4e", mfValue);
        writelogf(" % -12.4e", feSpecies[kspec] * RT);
        writelogf(" % -12.4e", feFull * RT);
        if (irxn != npos) {
            writelogf(" % -12.4e", deltaGRxn[irxn] * RT);
            writelogf(" % -12.4e", (deltaGRxn[irxn] + feFull - feSpecies[kspec]) * RT);

            if (deltaGRxn[irxn] < 0.0) {
                if (molNumSpecies[kspec] > 0.0) {
                    writelog("   growing");
                } else {
                    writelog("    stable");
                }
            } else if (deltaGRxn[irxn] > 0.0) {
                if (molNumSpecies[kspec] > 0.0) {
                    writelog(" shrinking");
                } else {
                    writelog("  unstable");
                }
            } else {
                writelog(" balanced");
            }
        }
        writelog(" \n");
    }
    writelog("   ");
    writeline('-', 132);
}

void VCS_SOLVE::vcs_switch_pos(const bool ifunc, const size_t k1, const size_t k2)
{
    if (k1 == k2) {
        return;
    }
    if (k1 >= m_nsp || k2 >= m_nsp) {
        plogf("vcs_switch_pos: ifunc = 0: inappropriate args: %d %d\n",
              k1, k2);
    }

    // Handle the index pointer in the phase structures first
    vcs_VolPhase* pv1 = m_VolPhaseList[m_phaseID[k1]].get();
    vcs_VolPhase* pv2 = m_VolPhaseList[m_phaseID[k2]].get();
    size_t kp1 = m_speciesLocalPhaseIndex[k1];
    size_t kp2 = m_speciesLocalPhaseIndex[k2];
    AssertThrowMsg(pv1->spGlobalIndexVCS(kp1) == k1, "VCS_SOLVE::vcs_switch_pos",
        "Indexing error");
    AssertThrowMsg(pv2->spGlobalIndexVCS(kp2) == k2, "VCS_SOLVE::vcs_switch_pos",
        "Indexing error");
    pv1->setSpGlobalIndexVCS(kp1, k2);
    pv2->setSpGlobalIndexVCS(kp2, k1);
    std::swap(m_speciesName[k1], m_speciesName[k2]);
    std::swap(m_molNumSpecies_old[k1], m_molNumSpecies_old[k2]);
    std::swap(m_speciesUnknownType[k1], m_speciesUnknownType[k2]);
    std::swap(m_molNumSpecies_new[k1], m_molNumSpecies_new[k2]);
    std::swap(m_SSfeSpecies[k1], m_SSfeSpecies[k2]);
    std::swap(m_spSize[k1], m_spSize[k2]);
    std::swap(m_deltaMolNumSpecies[k1], m_deltaMolNumSpecies[k2]);
    std::swap(m_feSpecies_old[k1], m_feSpecies_old[k2]);
    std::swap(m_feSpecies_new[k1], m_feSpecies_new[k2]);
    std::swap(m_SSPhase[k1], m_SSPhase[k2]);
    std::swap(m_phaseID[k1], m_phaseID[k2]);
    std::swap(m_speciesMapIndex[k1], m_speciesMapIndex[k2]);
    std::swap(m_speciesLocalPhaseIndex[k1], m_speciesLocalPhaseIndex[k2]);
    std::swap(m_actConventionSpecies[k1], m_actConventionSpecies[k2]);
    std::swap(m_lnMnaughtSpecies[k1], m_lnMnaughtSpecies[k2]);
    std::swap(m_actCoeffSpecies_new[k1], m_actCoeffSpecies_new[k2]);
    std::swap(m_actCoeffSpecies_old[k1], m_actCoeffSpecies_old[k2]);
    std::swap(m_wtSpecies[k1], m_wtSpecies[k2]);
    std::swap(m_chargeSpecies[k1], m_chargeSpecies[k2]);
    std::swap(m_speciesThermoList[k1], m_speciesThermoList[k2]);
    std::swap(m_PMVolumeSpecies[k1], m_PMVolumeSpecies[k2]);

    for (size_t j = 0; j < m_nelem; ++j) {
        std::swap(m_formulaMatrix(k1,j), m_formulaMatrix(k2,j));
    }
    if (m_useActCoeffJac && k1 != k2) {
        for (size_t i = 0; i < m_nsp; i++) {
            std::swap(m_np_dLnActCoeffdMolNum(k1,i), m_np_dLnActCoeffdMolNum(k2,i));
        }
        for (size_t i = 0; i < m_nsp; i++) {
            std::swap(m_np_dLnActCoeffdMolNum(i,k1), m_np_dLnActCoeffdMolNum(i,k2));
        }
    }
    std::swap(m_speciesStatus[k1], m_speciesStatus[k2]);

    // Handle the index pointer in the phase structures
    if (ifunc) {
        // Find the Rxn indices corresponding to the two species
        size_t i1 = k1 - m_numComponents;
        size_t i2 = k2 - m_numComponents;
        if (i1 > m_numRxnTot || i2 >= m_numRxnTot) {
            plogf("switch_pos: ifunc = 1: inappropriate noncomp values: %d %d\n",
                  i1 , i2);
        }
        for (size_t j = 0; j < m_numComponents; ++j) {
            std::swap(m_stoichCoeffRxnMatrix(j,i1), m_stoichCoeffRxnMatrix(j,i2));
        }
        std::swap(m_scSize[i1], m_scSize[i2]);
        for (size_t iph = 0; iph < m_numPhases; iph++) {
            std::swap(m_deltaMolNumPhase(iph,i1), m_deltaMolNumPhase(iph,i2));
            std::swap(m_phaseParticipation(iph,i1),
                      m_phaseParticipation(iph,i2));
        }
        std::swap(m_deltaGRxn_new[i1], m_deltaGRxn_new[i2]);
        std::swap(m_deltaGRxn_old[i1], m_deltaGRxn_old[i2]);
        std::swap(m_deltaGRxn_tmp[i1], m_deltaGRxn_tmp[i2]);
    }
}

void VCS_SOLVE::vcs_setFlagsVolPhases(const bool upToDate, const int stateCalc)
{
    if (!upToDate) {
        for (size_t iph = 0; iph < m_numPhases; iph++) {
            m_VolPhaseList[iph]->setMolesOutOfDate(stateCalc);
        }
    } else {
        for (size_t iph = 0; iph < m_numPhases; iph++) {
            m_VolPhaseList[iph]->setMolesCurrent(stateCalc);
        }
    }
}

void VCS_SOLVE::vcs_setFlagsVolPhase(const size_t iph, const bool upToDate,
                                     const int stateCalc)
{
    if (!upToDate) {
        m_VolPhaseList[iph]->setMolesOutOfDate(stateCalc);
    } else {
        m_VolPhaseList[iph]->setMolesCurrent(stateCalc);
    }
}

void VCS_SOLVE::vcs_updateMolNumVolPhases(const int stateCalc)
{
    for (size_t iph = 0; iph < m_numPhases; iph++) {
        m_VolPhaseList[iph]->updateFromVCS_MoleNumbers(stateCalc);
    }
}

}
