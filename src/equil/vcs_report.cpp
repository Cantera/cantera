//! @file vcs_report.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/equil/vcs_solve.h"
#include "cantera/equil/vcs_VolPhase.h"
#include "cantera/base/ctexceptions.h"

namespace Cantera
{
int VCS_SOLVE::vcs_report(int iconv)
{
    bool inertYes = false;

    // SORT DEPENDENT SPECIES IN DECREASING ORDER OF MOLES
    std::vector<std::pair<double, size_t>> x_order;
    for (size_t i = 0; i < m_nsp; i++) {
        x_order.push_back({-m_molNumSpecies_old[i], i});
    }
    std::sort(x_order.begin() + m_numComponents,
              x_order.begin() + m_numSpeciesRdc);

    vcs_setFlagsVolPhases(false, VCS_STATECALC_OLD);
    vcs_dfe(VCS_STATECALC_OLD, 0, 0, m_nsp);

    // PRINT OUT RESULTS
    plogf("\n\n\n\n");
    writeline('-', 80);
    writeline('-', 80);
    plogf("\t\t VCS_TP REPORT\n");
    writeline('-', 80);
    writeline('-', 80);
    if (iconv < 0) {
        plogf(" ERROR: CONVERGENCE CRITERION NOT SATISFIED.\n");
    } else if (iconv == 1) {
        plogf(" RANGE SPACE ERROR: Equilibrium Found but not all Element Abundances are Satisfied\n");
    }

    // Calculate some quantities that may need updating
    vcs_tmoles();
    m_totalVol = vcs_VolTotal(m_temperature, m_pressurePA,
                              &m_molNumSpecies_old[0], &m_PMVolumeSpecies[0]);

    plogf("\t\tTemperature  = %15.2g Kelvin\n", m_temperature);
    plogf("\t\tPressure     = %15.5g Pa \n", m_pressurePA);
    plogf("\t\ttotal Volume = %15.5g m**3\n", m_totalVol);

    // TABLE OF SPECIES IN DECREASING MOLE NUMBERS
    plogf("\n\n");
    writeline('-', 80);
    plogf(" Species                 Equilibrium kmoles   ");
    plogf("Mole Fraction    ChemPot/RT    SpecUnkType\n");
    writeline('-', 80);
    for (size_t i = 0; i < m_numComponents; ++i) {
        plogf(" %-12.12s", m_speciesName[i]);
        writeline(' ', 13, false);
        plogf("%14.7E     %14.7E    %12.4E", m_molNumSpecies_old[i],
              m_molNumSpecies_new[i], m_feSpecies_old[i]);
        plogf("   %3d", m_speciesUnknownType[i]);
        plogf("\n");
    }
    for (size_t i = m_numComponents; i < m_numSpeciesRdc; ++i) {
        size_t j = x_order[i].second;
        plogf(" %-12.12s", m_speciesName[j]);
        writeline(' ', 13, false);

        if (m_speciesUnknownType[j] == VCS_SPECIES_TYPE_MOLNUM) {
            plogf("%14.7E     %14.7E    %12.4E", m_molNumSpecies_old[j],
                  m_molNumSpecies_new[j], m_feSpecies_old[j]);
            plogf("  KMolNum ");
        } else if (m_speciesUnknownType[j] == VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
            plogf("        NA         %14.7E    %12.4E", 1.0, m_feSpecies_old[j]);
            plogf("   Voltage = %14.7E", m_molNumSpecies_old[j]);
        } else {
            throw CanteraError("VCS_SOLVE::vcs_report", "we have a problem");
        }
        plogf("\n");
    }
    for (size_t i = 0; i < m_numPhases; i++) {
        if (TPhInertMoles[i] > 0.0) {
            inertYes = true;
            if (i == 0) {
                plogf(" Inert Gas Species        ");
            } else {
                plogf(" Inert Species in phase %16s ",
                      m_VolPhaseList[i]->PhaseName);
            }
            plogf("%14.7E     %14.7E    %12.4E\n", TPhInertMoles[i],
                  TPhInertMoles[i] / m_tPhaseMoles_old[i], 0.0);
        }
    }
    if (m_numSpeciesRdc != m_nsp) {
        plogf("\n SPECIES WITH LESS THAN 1.0E-32 KMOLES:\n\n");
        for (size_t kspec = m_numSpeciesRdc; kspec < m_nsp; ++kspec) {
            plogf(" %-12.12s", m_speciesName[kspec]);
            // Note m_deltaGRxn_new[] stores in kspec slot not irxn slot, after solve
            plogf("             %14.7E     %14.7E    %12.4E",
                  m_molNumSpecies_old[kspec],
                  m_molNumSpecies_new[kspec], m_deltaGRxn_new[kspec]);
            if (m_speciesUnknownType[kspec] == VCS_SPECIES_TYPE_MOLNUM) {
                plogf("  KMol_Num");
            } else if (m_speciesUnknownType[kspec] == VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
                plogf("   Voltage");
            } else {
                plogf("   Unknown");
            }
            plogf("\n");
        }
    }
    writeline('-', 80);
    plogf("\n");

    // TABLE OF SPECIES FORMATION REACTIONS
    writeline('-', m_numComponents*10 + 45, true, true);
    plogf("               |ComponentID|");
    for (size_t j = 0; j < m_numComponents; j++) {
        plogf("        %3d", j);
    }
    plogf(" |           |\n");
    plogf("               | Components|");
    for (size_t j = 0; j < m_numComponents; j++) {
        plogf(" %10.10s", m_speciesName[j]);
    }
    plogf(" |           |\n");
    plogf(" NonComponent  |   Moles   |");
    for (size_t j = 0; j < m_numComponents; j++) {
        plogf(" %10.3g", m_molNumSpecies_old[j]);
    }
    plogf(" | DG/RT Rxn |\n");
    writeline('-', m_numComponents*10 + 45);
    for (size_t irxn = 0; irxn < m_numRxnTot; irxn++) {
        size_t kspec = m_indexRxnToSpecies[irxn];
        plogf(" %3d ", kspec);
        plogf("%-10.10s", m_speciesName[kspec]);
        plogf("|%10.3g |", m_molNumSpecies_old[kspec]);
        for (size_t j = 0; j < m_numComponents; j++) {
            plogf("     %6.2f", m_stoichCoeffRxnMatrix(j,irxn));
        }
        plogf(" |%10.3g |", m_deltaGRxn_new[irxn]);
        plogf("\n");
    }
    writeline('-', m_numComponents*10 + 45);
    plogf("\n");

    // TABLE OF PHASE INFORMATION
    vector_fp gaPhase(m_nelem, 0.0);
    vector_fp gaTPhase(m_nelem, 0.0);
    double totalMoles = 0.0;
    double gibbsPhase = 0.0;
    double gibbsTotal = 0.0;
    plogf("\n\n");
    plogf("\n");
    writeline('-', m_nelem*10 + 58);
    plogf("                  | ElementID |");
    for (size_t j = 0; j < m_nelem; j++) {
        plogf("        %3d", j);
    }
    plogf(" |                     |\n");
    plogf("                  | Element   |");
    for (size_t j = 0; j < m_nelem; j++) {
        plogf(" %10.10s", m_elementName[j]);
    }
    plogf(" |                     |\n");
    plogf("    PhaseName     |KMolTarget |");
    for (size_t j = 0; j < m_nelem; j++) {
        plogf(" %10.3g", m_elemAbundancesGoal[j]);
    }
    plogf(" |     Gibbs Total     |\n");
    writeline('-', m_nelem*10 + 58);
    for (size_t iphase = 0; iphase < m_numPhases; iphase++) {
        plogf(" %3d ", iphase);
        vcs_VolPhase* VPhase = m_VolPhaseList[iphase].get();
        plogf("%-12.12s |",VPhase->PhaseName);
        plogf("%10.3e |", m_tPhaseMoles_old[iphase]);
        totalMoles += m_tPhaseMoles_old[iphase];
        if (m_tPhaseMoles_old[iphase] != VPhase->totalMoles() &&
            !vcs_doubleEqual(m_tPhaseMoles_old[iphase], VPhase->totalMoles())) {
            throw CanteraError("VCS_SOLVE::vcs_report", "we have a problem");
        }
        vcs_elabPhase(iphase, &gaPhase[0]);
        for (size_t j = 0; j < m_nelem; j++) {
            plogf(" %10.3g", gaPhase[j]);
            gaTPhase[j] += gaPhase[j];
        }
        gibbsPhase = vcs_GibbsPhase(iphase, &m_molNumSpecies_old[0],
                                    &m_feSpecies_old[0]);
        gibbsTotal += gibbsPhase;
        plogf(" | %18.11E |\n", gibbsPhase);
    }
    writeline('-', m_nelem*10 + 58);
    plogf("    TOTAL         |%10.3e |", totalMoles);
    for (size_t j = 0; j < m_nelem; j++) {
        plogf(" %10.3g", gaTPhase[j]);
    }
    plogf(" | %18.11E |\n", gibbsTotal);

    writeline('-', m_nelem*10 + 58);
    plogf("\n");

    // GLOBAL SATISFACTION INFORMATION

    // Calculate the total dimensionless Gibbs Free Energy. Inert species are
    // handled as if they had a standard free energy of zero
    double g = vcs_Total_Gibbs(&m_molNumSpecies_old[0], &m_feSpecies_old[0],
                               &m_tPhaseMoles_old[0]);
    plogf("\n\tTotal Dimensionless Gibbs Free Energy = G/RT = %15.7E\n", g);
    if (inertYes) {
        plogf("\t\t(Inert species have standard free energy of zero)\n");
    }

    plogf("\nElemental Abundances (kmol): ");
    plogf("         Actual                    Target         Type      ElActive\n");
    for (size_t i = 0; i < m_nelem; ++i) {
        writeline(' ', 26, false);
        plogf("%-2.2s", m_elementName[i]);
        plogf("%20.12E  %20.12E", m_elemAbundances[i], m_elemAbundancesGoal[i]);
        plogf("   %3d     %3d\n", m_elType[i], m_elementActive[i]);
    }
    plogf("\n");

    // TABLE OF SPECIES CHEM POTS
    writeline('-', 93, true, true);
    plogf("Chemical Potentials of the Species: (dimensionless)\n");

    plogf("\t\t(RT = %g J/kmol)\n", GasConstant * m_temperature);
    plogf("    Name        TKMoles     StandStateChemPot   "
          "   ln(AC)       ln(X_i)      |   F z_i phi   |    ChemPot    | (-lnMnaught)");
    plogf("|  (MolNum ChemPot)|");
    writeline('-', 147, true, true);
    for (size_t i = 0; i < m_nsp; ++i) {
        size_t j = x_order[i].second;
        size_t pid = m_phaseID[j];
        plogf(" %-12.12s", m_speciesName[j]);
        plogf(" %14.7E ", m_molNumSpecies_old[j]);
        plogf("%14.7E  ", m_SSfeSpecies[j]);
        plogf("%14.7E  ", log(m_actCoeffSpecies_old[j]));
        double tpmoles = m_tPhaseMoles_old[pid];
        double phi = m_phasePhi[pid];
        double eContrib = phi * m_chargeSpecies[j] * m_Faraday_dim;
        double lx = 0.0;
        if (m_speciesUnknownType[j] == VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
            lx = 0.0;
        } else {
            if (tpmoles > 0.0 && m_molNumSpecies_old[j] > 0.0) {
                double tmp = std::max(VCS_DELETE_MINORSPECIES_CUTOFF, m_molNumSpecies_old[j]);
                lx = log(tmp) - log(tpmoles);
            } else {
                lx = m_feSpecies_old[j] - m_SSfeSpecies[j]
                     - log(m_actCoeffSpecies_old[j]) + m_lnMnaughtSpecies[j];
            }
        }
        plogf("%14.7E  |", lx);
        plogf("%14.7E | ", eContrib);
        double tmp = m_SSfeSpecies[j] + log(m_actCoeffSpecies_old[j])
                     + lx - m_lnMnaughtSpecies[j] + eContrib;
        if (fabs(m_feSpecies_old[j] - tmp) > 1.0E-7) {
            throw CanteraError("VCS_SOLVE::vcs_report",
                               "we have a problem - doesn't add up");
        }
        plogf(" %12.4E |", m_feSpecies_old[j]);
        if (m_lnMnaughtSpecies[j] != 0.0) {
            plogf("(%11.5E)", - m_lnMnaughtSpecies[j]);
        } else {
            plogf("             ");
        }

        plogf("|  %20.9E |", m_feSpecies_old[j] * m_molNumSpecies_old[j]);
        plogf("\n");
    }
    for (size_t i = 0; i < 125; i++) {
        plogf(" ");
    }
    plogf(" %20.9E\n", g);
    writeline('-', 147);

    // TABLE OF SOLUTION COUNTERS
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
    writeline('-', 80);
    writeline('-', 80);

    // Return a successful completion flag
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
              m_VCount->T_Calls_Inest, m_VCount->T_Time_inest);
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
