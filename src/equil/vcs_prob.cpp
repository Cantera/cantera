/**
 * @file vcs_prob.cpp
 *  Implementation for the Interface class for the vcs thermo
 *  equilibrium solver package,
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/equil/vcs_VolPhase.h"
#include "cantera/equil/vcs_species_thermo.h"
#include "cantera/equil/vcs_internal.h"
#include "cantera/equil/vcs_defs.h"
#include "cantera/equil/vcs_solve.h"
#include "cantera/thermo/MolalityVPSSTP.h"
#include "cantera/equil/MultiPhase.h"

#include <cstdio>

using namespace std;

namespace Cantera
{

void VCS_SOLVE::prob_report(int print_lvl)
{
    m_printLvl = print_lvl;

    // Printout the species information: PhaseID's and mole nums
    if (m_printLvl > 0) {
        writeline('=', 80, true, true);
        writeline('=', 20, false);
        plogf(" VCS_PROB: PROBLEM STATEMENT ");
        writeline('=', 31);
        writeline('=', 80);
        plogf("\n");
        plogf("\tSolve a constant T, P problem:\n");
        plogf("\t\tT    = %g K\n", m_temperature);
        double pres_atm = m_pressurePA / 1.01325E5;

        plogf("\t\tPres = %g atm\n", pres_atm);
        plogf("\n");
        plogf("             Phase IDs of species\n");
        plogf("            species     phaseID        phaseName   ");
        plogf(" Initial_Estimated_Moles   Species_Type\n");
        for (size_t i = 0; i < m_nsp; i++) {
            vcs_VolPhase* Vphase = m_VolPhaseList[m_phaseID[i]].get();
            plogf("%16s      %5d   %16s", m_mix->speciesName(i), m_phaseID[i],
                  Vphase->PhaseName);
            if (m_doEstimateEquil >= 0) {
                plogf("             %-10.5g", m_molNumSpecies_old[i]);
            } else {
                plogf("                N/A");
            }
            if (m_speciesUnknownType[i] == VCS_SPECIES_TYPE_MOLNUM) {
                plogf("                 Mol_Num");
            } else if (m_speciesUnknownType[i] == VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
                plogf("                 Voltage");
            } else {
                plogf("                        ");
            }
            plogf("\n");
        }

        // Printout of the Phase structure information
        writeline('-', 80, true, true);
        plogf("             Information about phases\n");
        plogf("  PhaseName    PhaseNum SingSpec  GasPhase   "
              " EqnState    NumSpec");
        plogf("  TMolesInert      TKmoles\n");

        for (size_t iphase = 0; iphase < m_numPhases; iphase++) {
            vcs_VolPhase* Vphase = m_VolPhaseList[iphase].get();
            plogf("%16s %5d %5d %8d ", Vphase->PhaseName,
                  Vphase->VP_ID_, Vphase->m_singleSpecies, Vphase->m_gasPhase);
            plogf("%16s %8d %16e ", Vphase->eos_name(),
                  Vphase->nSpecies(), Vphase->totalMolesInert());
            if (m_doEstimateEquil >= 0) {
                plogf("%16e\n", Vphase->totalMoles());
            } else {
                plogf("   N/A\n");
            }
        }

        plogf("\nElemental Abundances:    ");
        plogf("         Target_kmol    ElemType ElActive\n");
        for (size_t i = 0; i < m_nelem; ++i) {
            writeline(' ', 26, false);
            plogf("%-2.2s", m_elementName[i]);
            plogf("%20.12E  ", m_elemAbundancesGoal[i]);
            plogf("%3d       %3d\n", m_elType[i], m_elementActive[i]);
        }

        plogf("\nChemical Potentials:  (J/kmol)\n");
        plogf("             Species       (phase)    "
              "    SS0ChemPot       StarChemPot\n");
        for (size_t iphase = 0; iphase < m_numPhases; iphase++) {
            vcs_VolPhase* Vphase = m_VolPhaseList[iphase].get();
            Vphase->setState_TP(m_temperature, m_pressurePA);
            for (size_t kindex = 0; kindex < Vphase->nSpecies(); kindex++) {
                size_t kglob = Vphase->spGlobalIndexVCS(kindex);
                plogf("%16s ", m_mix->speciesName(kglob));
                if (kindex == 0) {
                    plogf("%16s", Vphase->PhaseName);
                } else {
                    plogf("                ");
                }

                plogf("%16g   %16g\n", Vphase->G0_calc_one(kindex),
                      Vphase->GStar_calc_one(kindex));
            }
        }
        writeline('=', 80, true, true);
        writeline('=', 20, false);
        plogf(" VCS_PROB: END OF PROBLEM STATEMENT ");
        writeline('=', 24);
        writeline('=', 80);
        plogf("\n");
    }
}

void VCS_SOLVE::addPhaseElements(vcs_VolPhase* volPhase)
{
    size_t neVP = volPhase->nElemConstraints();
    // Loop through the elements in the vol phase object
    for (size_t eVP = 0; eVP < neVP; eVP++) {
        size_t foundPos = npos;
        std::string enVP = volPhase->elementName(eVP);

        // Search for matches with the existing elements. If found, then fill in
        // the entry in the global mapping array.
        for (size_t e = 0; e < m_nelem; e++) {
            std::string en = m_elementName[e];
            if (!strcmp(enVP.c_str(), en.c_str())) {
                volPhase->setElemGlobalIndex(eVP, e);
                foundPos = e;
            }
        }
        if (foundPos == npos) {
            int elType = volPhase->elementType(eVP);
            int elactive = volPhase->elementActive(eVP);
            size_t e = addElement(enVP.c_str(), elType, elactive);
            volPhase->setElemGlobalIndex(eVP, e);
        }
    }
}

size_t VCS_SOLVE::addElement(const char* elNameNew, int elType, int elactive)
{
    if (!elNameNew) {
        throw CanteraError("VCS_SOLVE::addElement",
                           "error: element must have a name");
    }
    m_nelem++;
    m_numComponents++;

    m_formulaMatrix.resize(m_nsp, m_nelem, 0.0);
    m_stoichCoeffRxnMatrix.resize(m_nelem, m_nsp, 0.0);
    m_elType.push_back(elType);
    m_elementActive.push_back(elactive);
    m_elemAbundances.push_back(0.0);
    m_elemAbundancesGoal.push_back(0.0);
    m_elementMapIndex.push_back(0);
    m_elementName.push_back(elNameNew);
    return m_nelem - 1;
}

size_t VCS_SOLVE::addOnePhaseSpecies(vcs_VolPhase* volPhase, size_t k, size_t kT)
{
    if (kT > m_nsp) {
        // Need to expand the number of species here
        throw CanteraError("VCS_SOLVE::addOnePhaseSpecies", "Shouldn't be here");
    }
    const Array2D& fm = volPhase->getFormulaMatrix();
    for (size_t eVP = 0; eVP < volPhase->nElemConstraints(); eVP++) {
        size_t e = volPhase->elemGlobalIndex(eVP);
        AssertThrowMsg(e != npos, "VCS_PROB::addOnePhaseSpecies",
                       "element not found");
        m_formulaMatrix(kT,e) = fm(k,eVP);
    }

    // Tell the phase object about the current position of the species within
    // the global species vector
    volPhase->setSpGlobalIndexVCS(k, kT);
    return kT;
}

}
