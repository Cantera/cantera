/**
 * @file vcs_prob.cpp
 *  Implementation for the Interface class for the vcs thermo
 *  equilibrium solver package,
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

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

void VCS_SOLVE::set_gai()
{
    gai.assign(gai.size(), 0.0);
    for (size_t j = 0; j < ne; j++) {
        for (size_t kspec = 0; kspec < nspecies; kspec++) {
            if (SpeciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
                gai[j] += FormulaMatrix(kspec,j) * w[kspec];
            }
        }
    }
}

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
        if (prob_type == 0) {
            plogf("\tSolve a constant T, P problem:\n");
            plogf("\t\tT    = %g K\n", T);
            double pres_atm = PresPA / 1.01325E5;

            plogf("\t\tPres = %g atm\n", pres_atm);
        } else {
            throw CanteraError("VCS_PROB::prob_report", "Unknown problem type");
        }
        plogf("\n");
        plogf("             Phase IDs of species\n");
        plogf("            species     phaseID        phaseName   ");
        plogf(" Initial_Estimated_Moles   Species_Type\n");
        for (size_t i = 0; i < nspecies; i++) {
            vcs_VolPhase* Vphase = VPhaseList[PhaseID[i]];
            plogf("%16s      %5d   %16s", m_mix->speciesName(i), PhaseID[i],
                  Vphase->PhaseName);
            if (iest >= 0) {
                plogf("             %-10.5g", w[i]);
            } else {
                plogf("                N/A");
            }
            if (SpeciesUnknownType[i] == VCS_SPECIES_TYPE_MOLNUM) {
                plogf("                 Mol_Num");
            } else if (SpeciesUnknownType[i] == VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
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

        for (size_t iphase = 0; iphase < NPhase; iphase++) {
            vcs_VolPhase* Vphase = VPhaseList[iphase];
            plogf("%16s %5d %5d %8d ", Vphase->PhaseName,
                  Vphase->VP_ID_, Vphase->m_singleSpecies, Vphase->m_gasPhase);
            plogf("%16s %8d %16e ", Vphase->eos_name(),
                  Vphase->nSpecies(), Vphase->totalMolesInert());
            if (iest >= 0) {
                plogf("%16e\n", Vphase->totalMoles());
            } else {
                plogf("   N/A\n");
            }
        }

        plogf("\nElemental Abundances:    ");
        plogf("         Target_kmol    ElemType ElActive\n");
        for (size_t i = 0; i < ne; ++i) {
            writeline(' ', 26, false);
            plogf("%-2.2s", m_elementName[i]);
            plogf("%20.12E  ", gai[i]);
            plogf("%3d       %3d\n", m_elType[i], ElActive[i]);
        }

        plogf("\nChemical Potentials:  (J/kmol)\n");
        plogf("             Species       (phase)    "
              "    SS0ChemPot       StarChemPot\n");
        for (size_t iphase = 0; iphase < NPhase; iphase++) {
            vcs_VolPhase* Vphase = VPhaseList[iphase];
            Vphase->setState_TP(T, PresPA);
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
        for (size_t e = 0; e < ne; e++) {
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
    ne++;
    m_numElemConstraints++;
    m_numComponents++;

    gai.push_back(0.0);
    FormulaMatrix.resize(NSPECIES0, ne, 0.0);
    m_formulaMatrix.resize(NSPECIES0, ne);
    m_stoichCoeffRxnMatrix.resize(ne, NSPECIES0, 0.0);
    m_elType.push_back(elType);
    ElActive.push_back(elactive);
    m_elementActive.push_back(elactive);
    m_elemAbundances.push_back(0.0);
    m_elemAbundancesGoal.push_back(0.0);
    m_elementMapIndex.push_back(0);
    m_elementName.push_back(elNameNew);
    NE0 = ne;
    return ne - 1;
}

size_t VCS_SOLVE::addOnePhaseSpecies(vcs_VolPhase* volPhase, size_t k, size_t kT)
{
    if (kT > nspecies) {
        // Need to expand the number of species here
        throw CanteraError("VCS_PROB::addOnePhaseSpecies", "Shouldn't be here");
    }
    const Array2D& fm = volPhase->getFormulaMatrix();
    for (size_t eVP = 0; eVP < volPhase->nElemConstraints(); eVP++) {
        size_t e = volPhase->elemGlobalIndex(eVP);
        AssertThrowMsg(e != npos, "VCS_PROB::addOnePhaseSpecies",
                       "element not found");
        FormulaMatrix(kT,e) = fm(k,eVP);
    }

    // Tell the phase object about the current position of the species within
    // the global species vector
    volPhase->setSpGlobalIndexVCS(k, kT);
    return kT;
}

}
