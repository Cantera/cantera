/**
 * @file vcs_rxnadj.cpp
 *  Routines for carrying out various adjustments to the reaction steps
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/equil/vcs_solve.h"
#include "cantera/equil/vcs_VolPhase.h"
#include "cantera/base/ctexceptions.h"

#include <cstdio>

namespace Cantera
{

size_t VCS_SOLVE::vcs_RxnStepSizes(int& forceComponentCalc, size_t& kSpecial)
{
    size_t iphDel = npos;
    size_t k = 0;
    std::string ANOTE;
    if (m_debug_print_lvl >= 2) {
        plogf("   ");
        for (int j = 0; j < 82; j++) {
            plogf("-");
        }
        plogf("\n");
        plogf("   --- Subroutine vcs_RxnStepSizes called - Details:\n");
        plogf("   ");
        for (int j = 0; j < 82; j++) {
            plogf("-");
        }
        plogf("\n");
        plogf("   --- Species        KMoles     Rxn_Adjustment    DeltaG"
              "   | Comment\n");
    }

    // We update the matrix dlnActCoeffdmolNumber[][] at the top of the loop,
    // when necessary
    if (m_useActCoeffJac) {
        vcs_CalcLnActCoeffJac(&m_molNumSpecies_old[0]);
    }

    // LOOP OVER THE FORMATION REACTIONS
    for (size_t irxn = 0; irxn < m_numRxnRdc; ++irxn) {
        ANOTE = "Normal Calc";

        size_t kspec = m_indexRxnToSpecies[irxn];
        if (m_speciesStatus[kspec] == VCS_SPECIES_ZEROEDPHASE) {
            m_deltaMolNumSpecies[kspec] = 0.0;
            ANOTE = "ZeroedPhase: Phase is artificially zeroed";
        } else if (m_speciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
            if (m_molNumSpecies_old[kspec] == 0.0 && (!m_SSPhase[kspec])) {
                // MULTISPECIES PHASE WITH total moles equal to zero
                //
                // If dg[irxn] is negative, then the multispecies phase should
                // come alive again. Add a small positive step size to make it
                // come alive.
                if (m_deltaGRxn_new[irxn] < -1.0e-4) {
                    // First decide if this species is part of a multiphase that
                    // is nontrivial in size.
                    size_t iph = m_phaseID[kspec];
                    double tphmoles = m_tPhaseMoles_old[iph];
                    double trphmoles = tphmoles / m_totalMolNum;
                    vcs_VolPhase* Vphase = m_VolPhaseList[iph].get();
                    if (Vphase->exists() && (trphmoles > VCS_DELETE_PHASE_CUTOFF)) {
                        m_deltaMolNumSpecies[kspec] = m_totalMolNum * VCS_SMALL_MULTIPHASE_SPECIES;
                        if (m_speciesStatus[kspec] == VCS_SPECIES_STOICHZERO) {
                            m_deltaMolNumSpecies[kspec] = 0.0;
                            ANOTE = fmt::sprintf("MultSpec (%s): Species not born due to STOICH/PHASEPOP even though DG = %11.3E",
                                vcs_speciesType_string(m_speciesStatus[kspec], 15), m_deltaGRxn_new[irxn]);
                        } else {
                            m_deltaMolNumSpecies[kspec] = m_totalMolNum * VCS_SMALL_MULTIPHASE_SPECIES * 10.0;
                            ANOTE = fmt::sprintf("MultSpec (%s): small species born again DG = %11.3E",
                                vcs_speciesType_string(m_speciesStatus[kspec], 15), m_deltaGRxn_new[irxn]);
                        }
                    } else {
                        ANOTE = fmt::sprintf("MultSpec (%s):still dead, no phase pop, even though DG = %11.3E",
                            vcs_speciesType_string(m_speciesStatus[kspec], 15), m_deltaGRxn_new[irxn]);
                        m_deltaMolNumSpecies[kspec] = 0.0;
                        if (Vphase->exists() > 0 && trphmoles > 0.0) {
                            m_deltaMolNumSpecies[kspec] = m_totalMolNum * VCS_SMALL_MULTIPHASE_SPECIES * 10.;
                            ANOTE = fmt::sprintf("MultSpec (%s): birthed species because it was zero in a small existing phase with DG = %11.3E",
                                vcs_speciesType_string(m_speciesStatus[kspec], 15), m_deltaGRxn_new[irxn]);
                        }
                    }
                } else {
                    ANOTE = fmt::sprintf("MultSpec (%s): still dead DG = %11.3E",
                        vcs_speciesType_string(m_speciesStatus[kspec], 15), m_deltaGRxn_new[irxn]);
                    m_deltaMolNumSpecies[kspec] = 0.0;
                }
            } else {
                // REGULAR PROCESSING
                //
                // First take care of cases where we want to bail out. Don't
                // bother if superconvergence has already been achieved in this
                // mode.
                if (fabs(m_deltaGRxn_new[irxn]) <= m_tolmaj2) {
                    ANOTE = fmt::sprintf("Skipped: superconverged DG = %11.3E", m_deltaGRxn_new[irxn]);
                    if (m_debug_print_lvl >= 2) {
                        plogf("   --- %-12.12s", m_speciesName[kspec]);
                        plogf("  %12.4E %12.4E %12.4E | %s\n",
                              m_molNumSpecies_old[kspec], m_deltaMolNumSpecies[kspec],
                              m_deltaGRxn_new[irxn], ANOTE);
                    }
                    continue;
                }

                // Don't calculate for minor or nonexistent species if their
                // values are to be decreasing anyway.
                if ((m_speciesStatus[kspec] != VCS_SPECIES_MAJOR) && (m_deltaGRxn_new[irxn] >= 0.0)) {
                    ANOTE = fmt::sprintf("Skipped: IC = %3d and DG >0: %11.3E",
                        m_speciesStatus[kspec], m_deltaGRxn_new[irxn]);
                    if (m_debug_print_lvl >= 2) {
                        plogf("   --- %-12.12s", m_speciesName[kspec]);
                        plogf("  %12.4E %12.4E %12.4E | %s\n",
                              m_molNumSpecies_old[kspec], m_deltaMolNumSpecies[kspec],
                              m_deltaGRxn_new[irxn], ANOTE);
                    }
                    continue;
                }

                // Start of the regular processing
                double s;
                if (m_SSPhase[kspec]) {
                    s = 0.0;
                } else {
                    s = 1.0 / m_molNumSpecies_old[kspec];
                }
                for (size_t j = 0; j < m_numComponents; ++j) {
                    if (!m_SSPhase[j] && m_molNumSpecies_old[j] > 0.0) {
                        s += pow(m_stoichCoeffRxnMatrix(j,irxn), 2) / m_molNumSpecies_old[j];
                    }
                }
                for (size_t j = 0; j < m_numPhases; j++) {
                    vcs_VolPhase* Vphase = m_VolPhaseList[j].get();
                    if (!Vphase->m_singleSpecies && m_tPhaseMoles_old[j] > 0.0) {
                        s -= pow(m_deltaMolNumPhase(j,irxn), 2) / m_tPhaseMoles_old[j];
                    }
                }
                if (s != 0.0) {
                    // Take into account of the derivatives of the activity
                    // coefficients with respect to the mole numbers, even in
                    // our diagonal approximation.
                    if (m_useActCoeffJac) {
                        double s_old = s;
                        s = vcs_Hessian_diag_adj(irxn, s_old);
                        ANOTE = fmt::sprintf("Normal calc: diag adjusted from %g "
                            "to %g due to act coeff", s_old, s);
                    }

                    m_deltaMolNumSpecies[kspec] = -m_deltaGRxn_new[irxn] / s;
                    // New section to do damping of the m_deltaMolNumSpecies[]
                    for (size_t j = 0; j < m_numComponents; ++j) {
                        double stoicC = m_stoichCoeffRxnMatrix(j,irxn);
                        if (stoicC != 0.0) {
                            double negChangeComp = -stoicC * m_deltaMolNumSpecies[kspec];
                            if (negChangeComp > m_molNumSpecies_old[j]) {
                                if (m_molNumSpecies_old[j] > 0.0) {
                                    ANOTE = fmt::sprintf("Delta damped from %g "
                                        "to %g due to component %d (%10s) going neg", m_deltaMolNumSpecies[kspec],
                                        -m_molNumSpecies_old[j] / stoicC, j, m_speciesName[j]);
                                    m_deltaMolNumSpecies[kspec] = -m_molNumSpecies_old[j] / stoicC;
                                } else {
                                    ANOTE = fmt::sprintf("Delta damped from %g "
                                        "to %g due to component %d (%10s) zero", m_deltaMolNumSpecies[kspec],
                                        -m_molNumSpecies_old[j] / stoicC, j, m_speciesName[j]);
                                    m_deltaMolNumSpecies[kspec] = 0.0;
                                }
                            }
                        }
                    }
                    // Implement a damping term that limits m_deltaMolNumSpecies
                    // to the size of the mole number
                    if (-m_deltaMolNumSpecies[kspec] > m_molNumSpecies_old[kspec]) {
                        ANOTE = fmt::sprintf("Delta damped from %g "
                            "to %g due to %s going negative", m_deltaMolNumSpecies[kspec], -m_molNumSpecies_old[kspec],
                            m_speciesName[kspec]);
                        m_deltaMolNumSpecies[kspec] = -m_molNumSpecies_old[kspec];
                    }
                } else {
                    // REACTION IS ENTIRELY AMONGST SINGLE SPECIES PHASES.
                    // DELETE ONE OF THE PHASES AND RECOMPUTE BASIS.
                    //
                    // Either the species L will disappear or one of the
                    // component single species phases will disappear. The sign
                    // of DG(I) will indicate which way the reaction will go.
                    // Then, we need to follow the reaction to see which species
                    // will zero out first. The species to be zeroed out will be
                    // "k".
                    double dss;
                    if (m_deltaGRxn_new[irxn] > 0.0) {
                        dss = m_molNumSpecies_old[kspec];
                        k = kspec;
                        for (size_t j = 0; j < m_numComponents; ++j) {
                            if (m_stoichCoeffRxnMatrix(j,irxn) > 0.0) {
                                double xx = m_molNumSpecies_old[j] / m_stoichCoeffRxnMatrix(j,irxn);
                                if (xx < dss) {
                                    dss = xx;
                                    k = j;
                                }
                            }
                        }
                        dss = -dss;
                    } else {
                        dss = 1.0e10;
                        for (size_t j = 0; j < m_numComponents; ++j) {
                            if (m_stoichCoeffRxnMatrix(j,irxn) < 0.0) {
                                double xx = -m_molNumSpecies_old[j] / m_stoichCoeffRxnMatrix(j,irxn);
                                if (xx < dss) {
                                    dss = xx;
                                    k = j;
                                }
                            }
                        }
                    }

                    // Here we adjust the mole fractions according to DSS and
                    // the stoichiometric array to take into account that we are
                    // eliminating the kth species. DSS contains the amount of
                    // moles of the kth species that needs to be added back into
                    // the component species.
                    if (dss != 0.0) {
                        if ((k == kspec) && (m_SSPhase[kspec] != 1)) {
                            // Found out that we can be in this spot, when
                            // components of multispecies phases are zeroed,
                            // leaving noncomponent species of the same phase
                            // having all of the mole numbers of that phases. it
                            // seems that we can suggest a zero of the species
                            // and the code will recover.
                            ANOTE = fmt::sprintf("Delta damped from %g to %g due to delete %s", m_deltaMolNumSpecies[kspec],
                                -m_molNumSpecies_old[kspec], m_speciesName[kspec]);
                            m_deltaMolNumSpecies[kspec] = -m_molNumSpecies_old[kspec];
                            if (m_debug_print_lvl >= 2) {
                                plogf("   --- %-12.12s", m_speciesName[kspec]);
                                plogf("  %12.4E %12.4E %12.4E | %s\n",
                                      m_molNumSpecies_old[kspec], m_deltaMolNumSpecies[kspec],
                                      m_deltaGRxn_new[irxn], ANOTE);
                            }
                            continue;
                        }

                        // Delete the single species phase
                        for (size_t j = 0; j < m_nsp; j++) {
                            m_deltaMolNumSpecies[j] = 0.0;
                        }
                        m_deltaMolNumSpecies[kspec] = dss;
                        for (size_t j = 0; j < m_numComponents; ++j) {
                            m_deltaMolNumSpecies[j] = dss * m_stoichCoeffRxnMatrix(j,irxn);
                        }

                        iphDel = m_phaseID[k];
                        kSpecial = k;

                        if (k != kspec) {
                            ANOTE = fmt::sprintf("Delete component SS phase %d named %s - SS phases only",
                                iphDel, m_speciesName[k]);
                        } else {
                            ANOTE = fmt::sprintf("Delete this SS phase %d - SS components only", iphDel);
                        }
                        if (m_debug_print_lvl >= 2) {
                            plogf("   --- %-12.12s", m_speciesName[kspec]);
                            plogf("  %12.4E %12.4E %12.4E | %s\n",
                                  m_molNumSpecies_old[kspec], m_deltaMolNumSpecies[kspec],
                                  m_deltaGRxn_new[irxn], ANOTE);
                            plogf("   --- vcs_RxnStepSizes Special section to set up to delete %s\n",
                                  m_speciesName[k]);
                        }
                        if (k != kspec) {
                            forceComponentCalc = 1;
                            debuglog("   ---   Force a component recalculation\n\n", m_debug_print_lvl >= 2);
                        }
                        if (m_debug_print_lvl >= 2) {
                            plogf("   ");
                            writeline('-', 82);
                        }
                        return iphDel;
                    }
                }
            } // End of regular processing
            if (m_debug_print_lvl >= 2) {
                plogf("   --- %-12.12s", m_speciesName[kspec]);
                plogf("  %12.4E %12.4E %12.4E | %s\n",
                      m_molNumSpecies_old[kspec], m_deltaMolNumSpecies[kspec],
                      m_deltaGRxn_new[irxn], ANOTE);
            }
        } // End of loop over m_speciesUnknownType
    } // End of loop over non-component stoichiometric formation reactions
    if (m_debug_print_lvl >= 2) {
        plogf("   ");
        writeline('-', 82);
    }
    return iphDel;
}

double VCS_SOLVE::vcs_Hessian_diag_adj(size_t irxn, double hessianDiag_Ideal)
{
    double diag = hessianDiag_Ideal;
    double hessActCoef = vcs_Hessian_actCoeff_diag(irxn);
    if (hessianDiag_Ideal <= 0.0) {
        throw CanteraError("VCS_SOLVE::vcs_Hessian_diag_adj",
                           "We shouldn't be here");
    }
    if (hessActCoef >= 0.0) {
        diag += hessActCoef;
    } else if (fabs(hessActCoef) < 0.6666 * hessianDiag_Ideal) {
        diag += hessActCoef;
    } else {
        diag -= 0.6666 * hessianDiag_Ideal;
    }
    return diag;
}

double VCS_SOLVE::vcs_Hessian_actCoeff_diag(size_t irxn)
{
    size_t kspec = m_indexRxnToSpecies[irxn];
    size_t kph = m_phaseID[kspec];
    double np_kspec = std::max(m_tPhaseMoles_old[kph], 1e-13);
    double* sc_irxn = m_stoichCoeffRxnMatrix.ptrColumn(irxn);

    // First the diagonal term of the Jacobian
    double s = m_np_dLnActCoeffdMolNum(kspec,kspec) / np_kspec;

    // Next, the other terms. Note this only a loop over the components So, it's
    // not too expensive to calculate.
    for (size_t j = 0; j < m_numComponents; j++) {
        if (!m_SSPhase[j]) {
            for (size_t k = 0; k < m_numComponents; ++k) {
                if (m_phaseID[k] == m_phaseID[j]) {
                    double np = m_tPhaseMoles_old[m_phaseID[k]];
                    if (np > 0.0) {
                        s += sc_irxn[k] * sc_irxn[j] * m_np_dLnActCoeffdMolNum(j,k) / np;
                    }
                }
            }
            if (kph == m_phaseID[j]) {
                s += sc_irxn[j] * (m_np_dLnActCoeffdMolNum(j,kspec) + m_np_dLnActCoeffdMolNum(kspec,j)) / np_kspec;
            }
        }
    }
    return s;
}

void VCS_SOLVE::vcs_CalcLnActCoeffJac(const double* const moleSpeciesVCS)
{
    // Loop over all of the phases in the problem
    for (size_t iphase = 0; iphase < m_numPhases; iphase++) {
        vcs_VolPhase* Vphase = m_VolPhaseList[iphase].get();

        // We don't need to call single species phases;
        if (!Vphase->m_singleSpecies && !Vphase->isIdealSoln()) {
            // update the mole numbers
            Vphase->setMolesFromVCS(VCS_STATECALC_OLD, moleSpeciesVCS);

            // Download the resulting calculation into the full vector. This
            // scatter calculation is carried out in the vcs_VolPhase object.
            Vphase->sendToVCS_LnActCoeffJac(m_np_dLnActCoeffdMolNum);
        }
    }
}

}
