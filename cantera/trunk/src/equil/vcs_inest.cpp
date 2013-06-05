/**
 * @file vcs_inest.cpp
 *   Implementation methods for obtaining a good initial guess
 */
/*
 * Copyright (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#include "cantera/equil/vcs_solve.h"
#include "cantera/equil/vcs_internal.h"
#include "cantera/equil/vcs_VolPhase.h"

#include "cantera/base/clockWC.h"

namespace VCSnonideal
{

static char pprefix[20] = "   --- vcs_inest: ";

void VCS_SOLVE::vcs_inest(double* const aw, double* const sa, double* const sm,
                          double* const ss, double test)
{
    size_t lt, ikl, kspec, iph, irxn;
    double s;
    double s1 = 0.0;
    double xl, par;
    bool finished, conv;
    size_t nspecies = m_numSpeciesTot;
    size_t nrxn = m_numRxnTot;

    // double *molNum   = VCS_DATA_PTR(m_molNumSpecies_old);
    double TMolesMultiphase;
    double* xtphMax = VCS_DATA_PTR(m_TmpPhase);
    double* xtphMin = VCS_DATA_PTR(m_TmpPhase2);

    ikl = 0;
    lt = 0;


    /*
     *       CALL ROUTINE TO SOLVE MAX(CC*molNum) SUCH THAT AX*molNum = BB
     *           AND molNum(I) .GE. 0.0
     *
     *   Note, both of these programs do this.
     */
    vcs_setMolesLinProg();

#ifdef DEBUG_MODE
    if (m_debug_print_lvl >= 2) {
        plogf("%s Mole Numbers returned from linear programming (vcs_inest initial guess):\n",
              pprefix);
        plogf("%s     SPECIES          MOLE_NUMBER      -SS_ChemPotential\n", pprefix);
        for (kspec = 0; kspec < nspecies; ++kspec) {
            plogf("%s     ", pprefix);
            plogf("%-12.12s", m_speciesName[kspec].c_str());
            plogf(" %15.5g  %12.3g\n", m_molNumSpecies_old[kspec], -m_SSfeSpecies[kspec]);
        }
        plogf("%s Element Abundance Agreement returned from linear "
              "programming (vcs_inest initial guess):",
              pprefix);
        plogendl();
        plogf("%s     Element           Goal         Actual\n", pprefix);
        int jj = 0;
        for (size_t j = 0; j < m_numElemConstraints; j++) {
            if (m_elementActive[j]) {
                double tmp = 0.0;
                for (kspec = 0; kspec < nspecies; ++kspec) {
                    tmp +=  m_formulaMatrix[j][kspec] * m_molNumSpecies_old[kspec];
                }
                plogf("%s     ", pprefix);
                plogf("   %-9.9s", (m_elementName[j]).c_str());
                plogf(" %12.3g %12.3g\n", m_elemAbundancesGoal[j], tmp);
                jj++;
            }
        }
        plogendl();
    }
#endif


    /*
     *     Make sure all species have positive definite mole numbers
     *     Set voltages to zero for now, until we figure out what to do
     */
    vcs_dzero(VCS_DATA_PTR(m_deltaMolNumSpecies), nspecies);
    for (kspec = 0; kspec < nspecies; ++kspec) {
        iph = m_phaseID[kspec];
        if (m_speciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
            if (m_molNumSpecies_old[kspec] <= 0.0) {
                /*
                 * HKM Should eventually include logic here for non SS phases
                 */
                if (!m_SSPhase[kspec]) {
                    m_molNumSpecies_old[kspec] = 1.0e-30;
                }
            }
        } else {
            m_molNumSpecies_old[kspec] = 0.0;
        }
    }

    /*
     *      Now find the optimized basis that spans the stoichiometric
     *      coefficient matrix
     */
    (void) vcs_basopt(false, aw, sa, sm, ss, test, &conv);

    /* ***************************************************************** */
    /* **** CALCULATE TOTAL MOLES,                    ****************** */
    /* **** CHEMICAL POTENTIALS OF BASIS              ****************** */
    /* ***************************************************************** */
    /*
     * Calculate TMoles and m_tPhaseMoles_old[]
     */
    vcs_tmoles();
    /*
     * m_tPhaseMoles_new[] will consist of just the component moles
     */
    for (iph = 0; iph < m_numPhases; iph++) {
        m_tPhaseMoles_new[iph] = TPhInertMoles[iph] + 1.0E-20;
    }
    for (kspec = 0; kspec < m_numComponents; ++kspec) {
        if (m_speciesUnknownType[kspec] == VCS_SPECIES_TYPE_MOLNUM) {
            m_tPhaseMoles_new[m_phaseID[kspec]] += m_molNumSpecies_old[kspec];
        }
    }
    TMolesMultiphase = 0.0;
    for (iph = 0; iph < m_numPhases; iph++) {
        if (! m_VolPhaseList[iph]->m_singleSpecies) {
            TMolesMultiphase += m_tPhaseMoles_new[iph];
        }
    }
    vcs_dcopy(VCS_DATA_PTR(m_molNumSpecies_new), VCS_DATA_PTR(m_molNumSpecies_old),  nspecies);
    for (kspec = 0; kspec < m_numComponents; ++kspec) {
        if (m_speciesUnknownType[kspec] != VCS_SPECIES_TYPE_MOLNUM) {
            m_molNumSpecies_new[kspec] = 0.0;
        }
    }
    vcs_dcopy(VCS_DATA_PTR(m_feSpecies_new), VCS_DATA_PTR(m_SSfeSpecies),
              nspecies);

    for (kspec = 0; kspec < m_numComponents; ++kspec) {
        if (m_speciesUnknownType[kspec] == VCS_SPECIES_TYPE_MOLNUM) {
            if (! m_SSPhase[kspec]) {
                iph = m_phaseID[kspec];
                m_feSpecies_new[kspec] += log(m_molNumSpecies_new[kspec] / m_tPhaseMoles_old[iph]);
            }
        } else {
            m_molNumSpecies_new[kspec] = 0.0;
        }
    }
    vcs_deltag(0, true, VCS_STATECALC_NEW);
#ifdef DEBUG_MODE
    if (m_debug_print_lvl >= 2) {
        for (kspec = 0; kspec < nspecies; ++kspec) {
            plogf("%s", pprefix);
            plogf("%-12.12s", m_speciesName[kspec].c_str());
            if (kspec < m_numComponents)
                plogf("fe* = %15.5g ff = %15.5g\n", m_feSpecies_new[kspec],
                      m_SSfeSpecies[kspec]);
            else
                plogf("fe* = %15.5g ff = %15.5g dg* = %15.5g\n",
                      m_feSpecies_new[kspec], m_SSfeSpecies[kspec], m_deltaGRxn_new[kspec-m_numComponents]);
        }
    }
#endif
    /* ********************************************************** */
    /* **** ESTIMATE REACTION ADJUSTMENTS *********************** */
    /* ********************************************************** */
    vcs_dzero(VCS_DATA_PTR(m_deltaPhaseMoles), m_numPhases);
    for (iph = 0; iph < m_numPhases; iph++) {
        xtphMax[iph] = log(m_tPhaseMoles_new[iph] * 1.0E32);
        xtphMin[iph] = log(m_tPhaseMoles_new[iph] * 1.0E-32);
    }
    for (irxn = 0; irxn < nrxn; ++irxn) {
        kspec = m_indexRxnToSpecies[irxn];
        /*
         * For single species phases, we will not estimate the
         * mole numbers. If the phase exists, it stays. If it
         * doesn't exist in the estimate, it doesn't come into
         * existence here.
         */
        if (! m_SSPhase[kspec]) {
            iph = m_phaseID[kspec];
            if (m_deltaGRxn_new[irxn] > xtphMax[iph]) {
                m_deltaGRxn_new[irxn] = 0.8 * xtphMax[iph];
            }
            if (m_deltaGRxn_new[irxn] < xtphMin[iph]) {
                m_deltaGRxn_new[irxn] = 0.8 * xtphMin[iph];
            }
            /*
             *   HKM -> The TMolesMultiphase is a change of mine.
             *          It more evenly distributes the initial moles amongst
             *          multiple multispecies phases according to the
             *          relative values of the standard state free energies.
             *          There is no change for problems with one multispecies
             *          phase.
             *            It cut diamond4.vin iterations down from 62 to 14.
             */
            m_deltaMolNumSpecies[kspec] = 0.5 * (m_tPhaseMoles_new[iph] + TMolesMultiphase)
                                          * exp(-m_deltaGRxn_new[irxn]);

            for (size_t k = 0; k < m_numComponents; ++k) {
                m_deltaMolNumSpecies[k] += m_stoichCoeffRxnMatrix[irxn][k] * m_deltaMolNumSpecies[kspec];
            }

            for (iph = 0; iph < m_numPhases; iph++) {
                m_deltaPhaseMoles[iph] += m_deltaMolNumPhase[irxn][iph] * m_deltaMolNumSpecies[kspec];
            }
        }
    }
#ifdef DEBUG_MODE
    if (m_debug_print_lvl >= 2) {
        for (kspec = 0; kspec < nspecies; ++kspec) {
            if (m_speciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
                plogf("%sdirection (", pprefix);
                plogf("%-12.12s", m_speciesName[kspec].c_str());
                plogf(") = %g", m_deltaMolNumSpecies[kspec]);
                if (m_SSPhase[kspec]) {
                    if (m_molNumSpecies_old[kspec] > 0.0) {
                        plogf(" (ssPhase exists at w = %g moles)", m_molNumSpecies_old[kspec]);
                    } else {
                        plogf(" (ssPhase doesn't exist -> stability not checked)");
                    }
                }
                plogendl();
            }
        }
    }
#endif
    /* *********************************************************** */
    /* **** KEEP COMPONENT SPECIES POSITIVE ********************** */
    /* *********************************************************** */
    par = 0.5;
    for (kspec = 0; kspec < m_numComponents; ++kspec) {
        if (m_speciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
            if (par < -m_deltaMolNumSpecies[kspec] / m_molNumSpecies_new[kspec]) {
                par = -m_deltaMolNumSpecies[kspec] / m_molNumSpecies_new[kspec];
            }
        }
    }
    par = 1. / par;
    if (par <= 1.0 && par > 0.0) {
        par *= 0.8;
    } else {
        par = 1.0;
    }
    /* ******************************************** */
    /* **** CALCULATE NEW MOLE NUMBERS ************ */
    /* ******************************************** */
    finished = false;
    do {
        for (kspec = 0; kspec < m_numComponents; ++kspec) {
            if (m_speciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
                m_molNumSpecies_old[kspec] = m_molNumSpecies_new[kspec] + par * m_deltaMolNumSpecies[kspec];
            } else {
                m_deltaMolNumSpecies[kspec] = 0.0;
            }
        }
        for (kspec = m_numComponents; kspec < nspecies; ++kspec) {
            if (m_speciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
                if (m_deltaMolNumSpecies[kspec] != 0.0) {
                    m_molNumSpecies_old[kspec] = m_deltaMolNumSpecies[kspec] * par;
                }
            }
        }
        /*
         * We have a new w[] estimate, go get the
         * TMoles and m_tPhaseMoles_old[] values
         */
        vcs_tmoles();
        if (lt > 0) {
            goto finished;
        }
        /* ******************************************* */
        /* **** CONVERGENCE FORCING SECTION ********** */
        /* ******************************************* */
        vcs_setFlagsVolPhases(false, VCS_STATECALC_OLD);
        vcs_dfe(VCS_STATECALC_OLD, 0, 0, nspecies);
        for (kspec = 0, s = 0.0; kspec < nspecies; ++kspec) {
            s += m_deltaMolNumSpecies[kspec] * m_feSpecies_old[kspec];
        }
        if (s == 0.0) {
            finished = true;
            continue;
        }
        if (s < 0.0) {
            if (ikl == 0) {
                finished = true;
                continue;
            }
        }
        /* ***************************************** */
        /* *** TRY HALF STEP SIZE ****************** */
        /* ***************************************** */
        if (ikl == 0) {
            s1 = s;
            par *= 0.5;
            ikl = 1;
            continue;
        }
        /* **************************************************** */
        /* **** FIT PARABOLA THROUGH HALF AND FULL STEPS ****** */
        /* **************************************************** */
        xl = (1.0 - s / (s1 - s)) * 0.5;
        if (xl < 0.0) {
            /* *************************************************** */
            /* *** POOR DIRECTION, REDUCE STEP SIZE TO 0.2 ******* */
            /* *************************************************** */
            par *= 0.2;
        } else {
            if (xl > 1.0) {
                /* *************************************************** */
                /* **** TOO BIG A STEP, TAKE ORIGINAL FULL STEP ****** */
                /* *************************************************** */
                par *= 2.0;
            } else {
                /* *************************************************** */
                /* **** ACCEPT RESULTS OF FORCER ********************* */
                /* *************************************************** */
                par = par * 2.0 * xl;
            }
        }
        lt = 1;
    } while (!finished);
finished:
    ;
#ifdef DEBUG_MODE
    if (m_debug_print_lvl >= 2) {
        plogf("%s     Final Mole Numbers produced by inest:\n",
              pprefix);
        plogf("%s     SPECIES      MOLE_NUMBER\n", pprefix);
        for (kspec = 0; kspec < nspecies; ++kspec) {
            plogf("%s     ", pprefix);
            plogf("%-12.12s", m_speciesName[kspec].c_str());
            plogf(" %g", m_molNumSpecies_old[kspec]);
            plogendl();
        }
    }
#endif
}

int VCS_SOLVE::vcs_inest_TP()
{
    int retn = 0;
    double    test;
    Cantera::clockWC tickTock;
    test = -1.0E20;

    if (m_doEstimateEquil > 0) {
        /*
         *  Calculate the elemental abundances
         */
        vcs_elab();
        if (vcs_elabcheck(0)) {
#ifdef DEBUG_MODE
            if (m_debug_print_lvl >= 2) {
                plogf("%s Initial guess passed element abundances on input\n", pprefix);
                plogf("%s m_doEstimateEquil = 1 so will use the input mole "
                      "numbers as estimates", pprefix);
                plogendl();
            }
#endif
            return retn;
#ifdef DEBUG_MODE
        } else {
            if (m_debug_print_lvl >= 2) {
                plogf("%s Initial guess failed element abundances on input\n", pprefix);
                plogf("%s m_doEstimateEquil = 1 so will discard input "
                      "mole numbers and find our own estimate", pprefix);
                plogendl();
            }
#endif
        }
    }

    /*
     *  Malloc temporary space for usage in this routine and in
     *  subroutines
     *        sm[ne*ne]
     *        ss[ne]
     *        sa[ne]
     *        aw[m]
     */
    std::vector<double> sm(m_numElemConstraints*m_numElemConstraints, 0.0);
    std::vector<double> ss(m_numElemConstraints, 0.0);
    std::vector<double> sa(m_numElemConstraints, 0.0);
    std::vector<double> aw(m_numSpeciesTot+ m_numElemConstraints, 0.0);
    /*
     *  Go get the estimate of the solution
     */
#ifdef DEBUG_MODE
    if (m_debug_print_lvl >= 2) {
        plogf("%sGo find an initial estimate for the equilibrium problem",
              pprefix);
        plogendl();
    }
#endif
    vcs_inest(VCS_DATA_PTR(aw), VCS_DATA_PTR(sa), VCS_DATA_PTR(sm),
              VCS_DATA_PTR(ss), test);
    /*
     *  Calculate the elemental abundances
     */
    vcs_elab();

    /*
     *      If we still fail to achieve the correct elemental abundances,
     *      try to fix the problem again by calling the main elemental abundances
     *      fixer routine, used in the main program. This
     *      attempts to tweak the mole numbers of the component species to
     *      satisfy the element abundance constraints.
     *
     *       Note: We won't do this unless we have to since it involves inverting
     *             a matrix.
     */
    bool rangeCheck = vcs_elabcheck(1);
    if (!vcs_elabcheck(0)) {
        if (DEBUG_MODE_ENABLED && m_debug_print_lvl >= 2) {
            plogf("%sInitial guess failed element abundances\n", pprefix);
            plogf("%sCall vcs_elcorr to attempt fix", pprefix);
            plogendl();
        }
        vcs_elcorr(VCS_DATA_PTR(sm), VCS_DATA_PTR(aw));
        rangeCheck  = vcs_elabcheck(1);
        if (!vcs_elabcheck(0)) {
            plogf("%sInitial guess still fails element abundance equations\n",
                  pprefix);
            plogf("%s - Inability to ever satisfy element abundance "
                  "constraints is probable", pprefix);
            plogendl();
            retn = -1;
        } else {
            if (DEBUG_MODE_ENABLED && m_debug_print_lvl >= 2) {
                if (rangeCheck) {
                    plogf("%sInitial guess now satisfies element abundances", pprefix);
                    plogendl();
                } else {
                    plogf("%sElement Abundances RANGE ERROR\n", pprefix);
                    plogf("%s - Initial guess satisfies NC=%d element abundances, "
                          "BUT not NE=%d element abundances", pprefix,
                          m_numComponents, m_numElemConstraints);
                    plogendl();
                }
            }
        }
    } else {
        if (DEBUG_MODE_ENABLED && m_debug_print_lvl >= 2) {
            if (rangeCheck) {
                plogf("%sInitial guess satisfies element abundances", pprefix);
                plogendl();
            } else {
                plogf("%sElement Abundances RANGE ERROR\n", pprefix);
                plogf("%s - Initial guess satisfies NC=%d element abundances, "
                      "BUT not NE=%d element abundances", pprefix,
                      m_numComponents, m_numElemConstraints);
                plogendl();
            }
        }
    }

#ifdef DEBUG_MODE
    if (m_debug_print_lvl >= 2) {
        plogf("%sTotal Dimensionless Gibbs Free Energy = %15.7E", pprefix,
              vcs_Total_Gibbs(VCS_DATA_PTR(m_molNumSpecies_old), VCS_DATA_PTR(m_feSpecies_new),
                              VCS_DATA_PTR(m_tPhaseMoles_old)));
        plogendl();
    }
#endif

    /*
     * Record time
     */
    double tsecond = tickTock.secondsWC();
    m_VCount->T_Time_inest += tsecond;
    (m_VCount->T_Calls_Inest)++;
    return retn;
}

}
