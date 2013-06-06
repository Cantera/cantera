//! @file vcs_report.cpp
/*
 * Copyright (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#include "cantera/equil/vcs_solve.h"
#include "cantera/equil/vcs_internal.h"
#include "cantera/equil/vcs_VolPhase.h"

namespace VCSnonideal
{
static void print_space(int num)
{
    for (int j = 0; j < num; j++) {
        plogf(" ");
    }
}

static void print_line(const std::string& schar, size_t num)
{
    for (size_t j = 0; j < num; j++) {
        plogf("%s", schar.c_str());
    }
    plogf("\n");
}

int VCS_SOLVE::vcs_report(int iconv)
{
    bool printActualMoles = true, inertYes = false;
    size_t i, j, l, k, kspec;
    size_t nspecies = m_numSpeciesTot;
    double  g;

    char originalUnitsState = m_unitsState;


    std::vector<size_t> sortindex(nspecies,0);
    std::vector<double> xy(nspecies,0.0);

    /* ************************************************************** */
    /* **** SORT DEPENDENT SPECIES IN DECREASING ORDER OF MOLES ***** */
    /* ************************************************************** */

    for (i = 0; i < nspecies; ++i) {
        sortindex[i] = i;
        xy[i] = m_molNumSpecies_old[i];
    }
    /*
     *       Sort the XY vector, the mole fraction vector,
     *       and the sort index vector, sortindex, according to
     *       the magnitude of the mole fraction vector.
     */
    for (l = m_numComponents; l < m_numSpeciesRdc; ++l) {
        k = vcs_optMax(VCS_DATA_PTR(xy), 0, l, m_numSpeciesRdc);
        if (k != l) {
            std::swap(xy[k], xy[l]);
            std::swap(sortindex[k], sortindex[l]);
        }
    }

    /*
     *  Decide whether we have to nondimensionalize the equations.
     *     -> For the printouts from this routine, we will use nondimensional
     *        representations. This may be expanded in the future.
     */
    if (m_unitsState == VCS_DIMENSIONAL_G) {
        vcs_nondim_TP();
    }
    double molScale = 1.0;
    if (printActualMoles) {
        molScale = m_totalMoleScale;
    }

    vcs_setFlagsVolPhases(false, VCS_STATECALC_OLD);
    vcs_dfe(VCS_STATECALC_OLD, 0, 0, m_numSpeciesTot);
    /* ******************************************************** */
    /* *** PRINT OUT RESULTS ********************************** */
    /* ******************************************************** */

    plogf("\n\n\n\n");
    print_line("-", 80);
    print_line("-", 80);
    plogf("\t\t VCS_TP REPORT\n");
    print_line("-", 80);
    print_line("-", 80);
    if (iconv < 0) {
        plogf(" ERROR: CONVERGENCE CRITERION NOT SATISFIED.\n");
    } else if (iconv == 1) {
        plogf(" RANGE SPACE ERROR: Equilibrium Found but not all Element Abundances are Satisfied\n");
    }
    /*
     *   Calculate some quantities that may need updating
     */
    vcs_tmoles();
    m_totalVol = vcs_VolTotal(m_temperature, m_pressurePA,
                              VCS_DATA_PTR(m_molNumSpecies_old), VCS_DATA_PTR(m_PMVolumeSpecies));

    plogf("\t\tTemperature  = %15.2g Kelvin\n", m_temperature);
    plogf("\t\tPressure     = %15.5g Pa \n", m_pressurePA);
    plogf("\t\ttotal Volume = %15.5g m**3\n", m_totalVol * molScale);
    if (!printActualMoles) {
        plogf("\t\tMole Scale = %15.5g kmol (all mole numbers and volumes are scaled by this value)\n",
              molScale);
    }

    /*
     * -------- TABLE OF SPECIES IN DECREASING MOLE NUMBERS --------------
     */
    plogf("\n\n");
    print_line("-", 80);
    plogf(" Species                 Equilibrium kmoles   ");
    plogf("Mole Fraction    ChemPot/RT    SpecUnkType\n");
    print_line("-", 80);
    for (i = 0; i < m_numComponents; ++i) {
        plogf(" %-12.12s", m_speciesName[i].c_str());
        print_space(13);
        plogf("%14.7E     %14.7E    %12.4E", m_molNumSpecies_old[i] * molScale,
              m_molNumSpecies_new[i] * molScale, m_feSpecies_old[i]);
        plogf("   %3d", m_speciesUnknownType[i]);
        plogf("\n");
    }
    for (i = m_numComponents; i < m_numSpeciesRdc; ++i) {
        l = sortindex[i];
        plogf(" %-12.12s", m_speciesName[l].c_str());
        print_space(13);

        if (m_speciesUnknownType[l] == VCS_SPECIES_TYPE_MOLNUM) {
            plogf("%14.7E     %14.7E    %12.4E", m_molNumSpecies_old[l] * molScale,
                  m_molNumSpecies_new[l] * molScale, m_feSpecies_old[l]);
            plogf("  KMolNum ");
        } else if (m_speciesUnknownType[l] == VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
            plogf("        NA         %14.7E    %12.4E", 1.0, m_feSpecies_old[l]);
            plogf("   Voltage = %14.7E", m_molNumSpecies_old[l] * molScale);
        } else {
            plogf("we have a problem\n");
            exit(EXIT_FAILURE);
        }
        plogf("\n");
    }
    for (i = 0; i < m_numPhases; i++) {
        if (TPhInertMoles[i] > 0.0) {
            inertYes = true;
            if (i == 0) {
                plogf(" Inert Gas Species        ");
            } else {
                plogf(" Inert Species in phase %16s ",
                      (m_VolPhaseList[i])->PhaseName.c_str());
            }
            plogf("%14.7E     %14.7E    %12.4E\n", TPhInertMoles[i] * molScale,
                  TPhInertMoles[i] /  m_tPhaseMoles_old[i], 0.0);
        }
    }
    if (m_numSpeciesRdc != nspecies) {
        plogf("\n SPECIES WITH LESS THAN 1.0E-32 KMOLES:\n\n");
        for (kspec = m_numSpeciesRdc; kspec < nspecies; ++kspec) {
            plogf(" %-12.12s", m_speciesName[kspec].c_str());
            // Note m_deltaGRxn_new[] stores in kspec slot not irxn slot, after solve
            plogf("             %14.7E     %14.7E    %12.4E",
                  m_molNumSpecies_old[kspec]*molScale,
                  m_molNumSpecies_new[kspec]*molScale, m_deltaGRxn_new[kspec]);
            if (m_speciesUnknownType[i] == VCS_SPECIES_TYPE_MOLNUM) {
                plogf("  KMol_Num");
            } else if (m_speciesUnknownType[i] == VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
                plogf("   Voltage");
            } else {
                plogf("   Unknown");
            }

            plogf("\n");
        }
    }
    print_line("-", 80);
    plogf("\n");

    /*
     * ---------- TABLE OF SPECIES FORMATION REACTIONS ------------------
     */
    plogf("\n");
    print_line("-", m_numComponents*10 + 45);
    plogf("               |ComponentID|");
    for (j = 0; j < m_numComponents; j++) {
        plogf("        %3d", j);
    }
    plogf(" |           |\n");
    plogf("               | Components|");
    for (j = 0; j < m_numComponents; j++) {
        plogf(" %10.10s", m_speciesName[j].c_str());
    }
    plogf(" |           |\n");
    plogf(" NonComponent  |   Moles   |");
    for (j = 0; j < m_numComponents; j++) {
        plogf(" %10.3g", m_molNumSpecies_old[j] * molScale);
    }
    plogf(" | DG/RT Rxn |\n");
    print_line("-", m_numComponents*10 + 45);
    for (size_t irxn = 0; irxn < m_numRxnTot; irxn++) {
        size_t kspec = m_indexRxnToSpecies[irxn];
        plogf(" %3d ", kspec);
        plogf("%-10.10s", m_speciesName[kspec].c_str());
        plogf("|%10.3g |", m_molNumSpecies_old[kspec]*molScale);
        for (j = 0; j < m_numComponents; j++) {
            plogf("     %6.2f", m_stoichCoeffRxnMatrix[irxn][j]);
        }
        plogf(" |%10.3g |", m_deltaGRxn_new[irxn]);
        plogf("\n");
    }
    print_line("-", m_numComponents*10 + 45);
    plogf("\n");

    /*
     * ------------------ TABLE OF PHASE INFORMATION ---------------------
     */
    std::vector<double> gaPhase(m_numElemConstraints, 0.0);
    std::vector<double> gaTPhase(m_numElemConstraints, 0.0);
    double totalMoles = 0.0;
    double gibbsPhase = 0.0;
    double gibbsTotal = 0.0;
    plogf("\n\n");
    plogf("\n");
    print_line("-", m_numElemConstraints*10 + 58);
    plogf("                  | ElementID |");
    for (j = 0; j < m_numElemConstraints; j++) {
        plogf("        %3d", j);
    }
    plogf(" |                     |\n");
    plogf("                  | Element   |");
    for (j = 0; j < m_numElemConstraints; j++) {
        plogf(" %10.10s", (m_elementName[j]).c_str());
    }
    plogf(" |                     |\n");
    plogf("    PhaseName     |KMolTarget |");
    for (j = 0; j < m_numElemConstraints; j++) {
        plogf(" %10.3g", m_elemAbundancesGoal[j]);
    }
    plogf(" |     Gibbs Total     |\n");
    print_line("-", m_numElemConstraints*10 + 58);
    for (size_t iphase = 0; iphase < m_numPhases; iphase++) {
        plogf(" %3d ", iphase);
        vcs_VolPhase* VPhase = m_VolPhaseList[iphase];
        plogf("%-12.12s |",VPhase->PhaseName.c_str());
        plogf("%10.3e |", m_tPhaseMoles_old[iphase]*molScale);
        totalMoles +=  m_tPhaseMoles_old[iphase];
        if (m_tPhaseMoles_old[iphase] != VPhase->totalMoles()) {
            if (! vcs_doubleEqual(m_tPhaseMoles_old[iphase], VPhase->totalMoles())) {
                plogf("We have a problem\n");
                exit(EXIT_FAILURE);
            }
        }
        vcs_elabPhase(iphase, VCS_DATA_PTR(gaPhase));
        for (j = 0; j < m_numElemConstraints; j++) {
            plogf(" %10.3g", gaPhase[j]);
            gaTPhase[j] += gaPhase[j];
        }
        gibbsPhase = vcs_GibbsPhase(iphase, VCS_DATA_PTR(m_molNumSpecies_old),
                                    VCS_DATA_PTR(m_feSpecies_old));
        gibbsTotal += gibbsPhase;
        plogf(" | %18.11E |\n", gibbsPhase);
    }
    print_line("-", m_numElemConstraints*10 + 58);
    plogf("    TOTAL         |%10.3e |", totalMoles);
    for (j = 0; j < m_numElemConstraints; j++) {
        plogf(" %10.3g", gaTPhase[j]);
    }
    plogf(" | %18.11E |\n", gibbsTotal);

    print_line("-", m_numElemConstraints*10 + 58);
    plogf("\n");

    /*
     * ----------- GLOBAL SATISFACTION INFORMATION -----------------------
     */

    /*
     *    Calculate the total dimensionless Gibbs Free Energy
     *     -> Inert species are handled as if they had a standard free
     *        energy of zero
     */

    g = vcs_Total_Gibbs(VCS_DATA_PTR(m_molNumSpecies_old), VCS_DATA_PTR(m_feSpecies_old),
                        VCS_DATA_PTR(m_tPhaseMoles_old));
    plogf("\n\tTotal Dimensionless Gibbs Free Energy = G/RT = %15.7E\n", g);
    if (inertYes) {
        plogf("\t\t(Inert species have standard free energy of zero)\n");
    }

    plogf("\nElemental Abundances (kmol): ");
    plogf("         Actual                    Target         Type      ElActive\n");
    for (i = 0; i < m_numElemConstraints; ++i) {
        print_space(26);
        plogf("%-2.2s", (m_elementName[i]).c_str());
        plogf("%20.12E  %20.12E", m_elemAbundances[i]*molScale, m_elemAbundancesGoal[i]*molScale);
        plogf("   %3d     %3d\n", m_elType[i], m_elementActive[i]);
    }
    plogf("\n");

    /*
     * ------------------ TABLE OF SPECIES CHEM POTS ---------------------
     */
    plogf("\n");
    print_line("-", 93);
    plogf("Chemical Potentials of the Species: (dimensionless)\n");

    double rt = vcs_nondimMult_TP(m_VCS_UnitsFormat, m_temperature);
    plogf("\t\t(RT = %g ", rt);
    vcs_printChemPotUnits(m_VCS_UnitsFormat);
    plogf(")\n");
    plogf("    Name        TKMoles     StandStateChemPot   "
          "   ln(AC)       ln(X_i)      |   F z_i phi   |    ChemPot    | (-lnMnaught)");
    plogf("|  (MolNum ChemPot)|");
    plogf("\n");
    print_line("-", 147);
    for (i = 0; i < nspecies; ++i) {
        l = sortindex[i];
        size_t pid = m_phaseID[l];
        plogf(" %-12.12s", m_speciesName[l].c_str());
        plogf(" %14.7E ", m_molNumSpecies_old[l]*molScale);
        plogf("%14.7E  ", m_SSfeSpecies[l]);
        plogf("%14.7E  ", log(m_actCoeffSpecies_old[l]));
        double tpmoles = m_tPhaseMoles_old[pid];
        double phi = m_phasePhi[pid];
        double eContrib = phi * m_chargeSpecies[l] * m_Faraday_dim;
        double lx = 0.0;
        if (m_speciesUnknownType[l] == VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
            lx = 0.0;
        } else {
            if (tpmoles > 0.0 && m_molNumSpecies_old[l] > 0.0) {
                double tmp = std::max(VCS_DELETE_MINORSPECIES_CUTOFF, m_molNumSpecies_old[l]);
                lx = log(tmp) - log(tpmoles);
            } else {
                lx = m_feSpecies_old[l] - m_SSfeSpecies[l]
                     - log(m_actCoeffSpecies_old[l]) + m_lnMnaughtSpecies[l];
            }
        }
        plogf("%14.7E  |", lx);
        plogf("%14.7E | ", eContrib);
        double tmp = m_SSfeSpecies[l] + log(m_actCoeffSpecies_old[l])
                     + lx - m_lnMnaughtSpecies[l] + eContrib;
        if (fabs(m_feSpecies_old[l] - tmp) > 1.0E-7) {
            plogf("\n\t\twe have a problem - doesn't add up\n");
            exit(EXIT_FAILURE);
        }
        plogf(" %12.4E |", m_feSpecies_old[l]);
        if (m_lnMnaughtSpecies[l] != 0.0) {
            plogf("(%11.5E)", - m_lnMnaughtSpecies[l]);
        } else {
            plogf("             ");
        }

        plogf("|  %20.9E |", m_feSpecies_old[l] * m_molNumSpecies_old[l] * molScale);
        plogf("\n");
    }
    for (i = 0; i < 125; i++) {
        plogf(" ");
    }
    plogf(" %20.9E\n", g);
    print_line("-", 147);

    /*
     * ------------- TABLE OF SOLUTION COUNTERS --------------------------
     */
    plogf("\n");
    plogf("\nCounters:         Iterations          Time (seconds)\n");
    if (m_timing_print_lvl > 0) {
        plogf("    vcs_basopt:   %5d             %11.5E\n",
              m_VCount->Basis_Opts, m_VCount->Time_basopt);
        plogf("    vcs_TP:       %5d             %11.5E\n",
              m_VCount->Its, m_VCount->Time_vcs_TP);
    } else {
        plogf("    vcs_basopt:   %5d             %11s\n",
              m_VCount->Basis_Opts,"    NA     ");
        plogf("    vcs_TP:       %5d             %11s\n",
              m_VCount->Its,"    NA     ");
    }
    print_line("-", 80);
    print_line("-", 80);

    /*
     *   Set the Units state of the system back to where it was when we
     *   entered the program.
     */
    if (originalUnitsState != m_unitsState) {
        if (originalUnitsState == VCS_DIMENSIONAL_G) {
            vcs_redim_TP();
        } else {
            vcs_nondim_TP();
        }
    }
    /*
     *   Return a successful completion flag
     */
    return VCS_SUCCESS;
}

void VCS_SOLVE::vcs_TCounters_report(int timing_print_lvl)
{
    plogf("\nTCounters:   Num_Calls   Total_Its       Total_Time (seconds)\n");
    if (timing_print_lvl > 0) {
        plogf("    vcs_basopt:   %5d      %5d         %11.5E\n",
              m_VCount->T_Basis_Opts, m_VCount->T_Basis_Opts,
              m_VCount->T_Time_basopt);
        plogf("    vcs_TP:       %5d      %5d         %11.5E\n",
              m_VCount->T_Calls_vcs_TP, m_VCount->T_Its,
              m_VCount->T_Time_vcs_TP);
        plogf("    vcs_inest:    %5d                    %11.5E\n",
              m_VCount->T_Calls_Inest,  m_VCount->T_Time_inest);
        plogf("    vcs_TotalTime:                         %11.5E\n",
              m_VCount->T_Time_vcs);
    } else {
        plogf("    vcs_basopt:   %5d      %5d         %11s\n",
              m_VCount->T_Basis_Opts, m_VCount->T_Basis_Opts,"    NA     ");
        plogf("    vcs_TP:       %5d      %5d         %11s\n",
              m_VCount->T_Calls_vcs_TP, m_VCount->T_Its,"    NA     ");
        plogf("    vcs_inest:    %5d                    %11s\n",
              m_VCount->T_Calls_Inest, "    NA     ");
        plogf("    vcs_TotalTime:                         %11s\n",
              "    NA     ");
    }
}

}
