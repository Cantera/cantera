/**
 *  @file vcs_prep.cpp
 *    This file contains some prepatory functions.
 */

/*
 * Copyright (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#include "cantera/equil/vcs_solve.h"
#include "cantera/equil/vcs_prob.h"
#include "cantera/equil/vcs_VolPhase.h"

namespace Cantera
{
void VCS_SOLVE::vcs_SSPhase()
{
    std::vector<int> numPhSpecies(m_numPhases, 0);

    for (size_t kspec = 0; kspec < m_numSpeciesTot; ++kspec) {
        numPhSpecies[m_phaseID[kspec]]++;
    }
    /*
     *           Handle the special case of a single species in a phase that
     *           has been earmarked as a multispecies phase.
     *           Treat that species as a single-species phase
     */
    for (size_t iph = 0; iph < m_numPhases; iph++) {
        vcs_VolPhase* Vphase = m_VolPhaseList[iph];
        Vphase->m_singleSpecies = false;
        if (TPhInertMoles[iph] > 0.0) {
            Vphase->setExistence(2);
        }
        if (numPhSpecies[iph] <= 1) {
            if (TPhInertMoles[iph] == 0.0) {
                Vphase->m_singleSpecies = true;
            }
        }
    }

    /*
     *  Fill in some useful arrays here that have to do with the
     *  static information concerning the phase ID of species.
     *       SSPhase = Boolean indicating whether a species is in a
     *                 single species phase or not.
     */
    for (size_t kspec = 0; kspec < m_numSpeciesTot; kspec++) {
        size_t iph = m_phaseID[kspec];
        vcs_VolPhase* Vphase = m_VolPhaseList[iph];
        if (Vphase->m_singleSpecies) {
            m_SSPhase[kspec] = true;
        } else {
            m_SSPhase[kspec] = false;
        }
    }
}

int VCS_SOLVE::vcs_prep_oneTime(int printLvl)
{
    int retn = VCS_SUCCESS;
    m_debug_print_lvl = printLvl;

    /*
     *  Calculate the Single Species status of phases
     *  Also calculate the number of species per phase
     */
    vcs_SSPhase();

    /*
     *      Set an initial estimate for the number of noncomponent species
     *      equal to nspecies - nelements. This may be changed below
     */
    if (m_numElemConstraints > m_numSpeciesTot) {
        m_numRxnTot = 0;
    } else {
        m_numRxnTot = m_numSpeciesTot - m_numElemConstraints;
    }
    m_numRxnRdc = m_numRxnTot;
    m_numSpeciesRdc = m_numSpeciesTot;
    for (size_t i = 0; i < m_numRxnRdc; ++i) {
        m_indexRxnToSpecies[i] = m_numElemConstraints + i;
    }

    for (size_t kspec = 0; kspec < m_numSpeciesTot; ++kspec) {
        size_t pID = m_phaseID[kspec];
        size_t spPhIndex = m_speciesLocalPhaseIndex[kspec];
        vcs_VolPhase* vPhase =  m_VolPhaseList[pID];
        vcs_SpeciesProperties* spProp = vPhase->speciesProperty(spPhIndex);
        double sz = 0.0;
        size_t eSize = spProp->FormulaMatrixCol.size();
        for (size_t e = 0; e < eSize; e++) {
            sz += fabs(spProp->FormulaMatrixCol[e]);
        }
        if (sz > 0.0) {
            m_spSize[kspec] = sz;
        } else {
            m_spSize[kspec] = 1.0;
        }
    }

    /* ***************************************************** */
    /* **** DETERMINE THE NUMBER OF COMPONENTS ************* */
    /* ***************************************************** */

    /*
     *       Obtain a valid estimate of the mole fraction. This will
     *       be used as an initial ordering vector for prioritizing
     *       which species are defined as components.
     *
     *       If a mole number estimate was supplied from the
     *       input file, use that mole number estimate.
     *
     *       If a solution estimate wasn't supplied from the input file,
     *       supply an initial estimate for the mole fractions
     *       based on the relative reverse ordering of the
     *       chemical potentials.
     *
     *       For voltage unknowns, set these to zero for the moment.
     */
    double test = -1.0e-10;
    bool modifiedSoln = false;
    if (m_doEstimateEquil < 0) {
        double sum  = 0.0;
        for (size_t kspec = 0; kspec < m_numSpeciesTot; ++kspec) {
            if (m_speciesUnknownType[kspec] == VCS_SPECIES_TYPE_MOLNUM) {
                sum += fabs(m_molNumSpecies_old[kspec]);
            }
        }
        if (fabs(sum) < 1.0E-6) {
            modifiedSoln = true;
            double pres = (m_pressurePA <= 0.0) ? 1.01325E5 : m_pressurePA;
            retn = vcs_evalSS_TP(0, 0, m_temperature, pres);
            for (size_t kspec = 0; kspec < m_numSpeciesTot; ++kspec) {
                if (m_speciesUnknownType[kspec] == VCS_SPECIES_TYPE_MOLNUM) {
                    m_molNumSpecies_old[kspec] = - m_SSfeSpecies[kspec];
                } else {
                    m_molNumSpecies_old[kspec] = 0.0;
                }
            }
        }
        test = -1.0e20;
    }

    /*
     *      NC = number of components is in the vcs.h common block
     *     This call to BASOPT doesn't calculate the stoichiometric
     *     reaction matrix.
     */
    std::vector<double> awSpace(m_numSpeciesTot + (m_numElemConstraints + 2)*(m_numElemConstraints), 0.0);
    double* aw = VCS_DATA_PTR(awSpace);
    if (aw == NULL) {
        plogf("vcs_prep_oneTime: failed to get memory: global bailout\n");
        return VCS_NOMEMORY;
    }
    double* sa = aw + m_numSpeciesTot;
    double* sm = sa + m_numElemConstraints;
    double* ss = sm + (m_numElemConstraints)*(m_numElemConstraints);
    bool conv;
    retn = vcs_basopt(true, aw, sa, sm, ss, test, &conv);
    if (retn != VCS_SUCCESS) {
        plogf("vcs_prep_oneTime:");
        plogf(" Determination of number of components failed: %d\n",
              retn);
        plogf("          Global Bailout!\n");
        return retn;
    }

    if (m_numSpeciesTot >= m_numComponents) {
        m_numRxnTot = m_numRxnRdc = m_numSpeciesTot - m_numComponents;
        for (size_t i = 0; i < m_numRxnRdc; ++i) {
            m_indexRxnToSpecies[i] = m_numComponents + i;
        }
    } else {
        m_numRxnTot = m_numRxnRdc = 0;
    }

    /*
     *   The elements might need to be rearranged.
     */
    awSpace.resize(m_numElemConstraints + (m_numElemConstraints + 2)*(m_numElemConstraints), 0.0);
    aw = VCS_DATA_PTR(awSpace);
    sa = aw + m_numElemConstraints;
    sm = sa + m_numElemConstraints;
    ss = sm + (m_numElemConstraints)*(m_numElemConstraints);
    retn = vcs_elem_rearrange(aw, sa, sm, ss);
    if (retn != VCS_SUCCESS) {
        plogf("vcs_prep_oneTime:");
        plogf(" Determination of element reordering failed: %d\n",
              retn);
        plogf("          Global Bailout!\n");
        return retn;
    }

    // If we mucked up the solution unknowns because they were all
    // zero to start with, set them back to zero here
    if (modifiedSoln) {
        for (size_t kspec = 0; kspec < m_numSpeciesTot; ++kspec) {
            m_molNumSpecies_old[kspec] = 0.0;
        }
    }
    return VCS_SUCCESS;
}

int VCS_SOLVE::vcs_prep()
{
    /*
     *        Initialize various arrays in the data to zero
     */
    m_feSpecies_old.assign(m_feSpecies_old.size(), 0.0);
    m_feSpecies_new.assign(m_feSpecies_new.size(), 0.0);
    m_molNumSpecies_new.assign(m_molNumSpecies_new.size(), 0.0);
    m_deltaMolNumPhase.zero();
    m_phaseParticipation.zero();
    m_deltaPhaseMoles.assign(m_deltaPhaseMoles.size(), 0.0);
    m_tPhaseMoles_new.assign(m_tPhaseMoles_new.size(), 0.0);
    /*
     *   Calculate the total number of moles in all phases.
     */
    vcs_tmoles();
    return VCS_SUCCESS;
}

bool VCS_SOLVE::vcs_wellPosed(VCS_PROB* vprob)
{
    double sum = 0.0;
    for (size_t e = 0; e < vprob->ne; e++) {
        sum = sum + vprob->gai[e];
    }
    if (sum < 1.0E-20) {
        plogf("vcs_wellPosed: Element abundance is close to zero\n");
        return false;
    }
    return true;
}

}
