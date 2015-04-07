//! @file vcs_solve_phaseStability.cpp

/*
 * Copyright (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#include "cantera/equil/vcs_solve.h"
#include "cantera/equil/vcs_internal.h"
#include "cantera/equil/vcs_VolPhase.h"
#include "cantera/equil/vcs_species_thermo.h"
#include "cantera/equil/vcs_prob.h"

#include "cantera/base/clockWC.h"

using namespace std;

namespace VCSnonideal
{

int VCS_SOLVE::vcs_PS(VCS_PROB* vprob, int iphase, int printLvl, double& feStable)
{

    /*
     *  ifunc determines the problem type
     */
    int ifunc = 0;
    int iStab = 0;

    /*
     *         This function is called to create the private data
     *         using the public data.
     */
    size_t nspecies0 = vprob->nspecies + 10;
    size_t nelements0 = vprob->ne;
    size_t nphase0 = vprob->NPhase;

    vcs_initSizes(nspecies0, nelements0, nphase0);


    if (ifunc < 0 || ifunc > 2) {
        plogf("vcs: Unrecognized value of ifunc, %d: bailing!\n",
              ifunc);
        return VCS_PUB_BAD;
    }

    /*
     *         This function is called to copy the public data
     *         and the current problem specification
     *         into the current object's data structure.
     */
    int retn = vcs_prob_specifyFully(vprob);
    if (retn != 0) {
        plogf("vcs_pub_to_priv returned a bad status, %d: bailing!\n",
              retn);
        return retn;
    }
    /*
     *        Prep the problem data
     *           - adjust the identity of any phases
     *           - determine the number of components in the problem
     */
    retn = vcs_prep_oneTime(printLvl);
    if (retn != 0) {
        plogf("vcs_prep_oneTime returned a bad status, %d: bailing!\n",
              retn);
        return retn;
    }


    /*
     *         This function is called to copy the current problem
     *         into the current object's data structure.
     */
    retn = vcs_prob_specify(vprob);
    if (retn != 0) {
        plogf("vcs_prob_specify returned a bad status, %d: bailing!\n",
              retn);
        return retn;
    }


    /*
     *        Prep the problem data for this particular instantiation of
     *        the problem
     */
    retn = vcs_prep();
    if (retn != VCS_SUCCESS) {
        plogf("vcs_prep returned a bad status, %d: bailing!\n", retn);
        return retn;
    }
    /*
     *        Check to see if the current problem is well posed.
     */
    if (!vcs_wellPosed(vprob)) {
        plogf("vcs has determined the problem is not well posed: Bailing\n");
        return VCS_PUB_BAD;
    }

    /*
     *        Store the temperature and pressure in the private global variables
     */
    m_temperature = vprob->T;
    m_pressurePA = vprob->PresPA;
    /*
     *        Evaluate the standard state free energies
     *        at the current temperatures and pressures.
     */
    vcs_evalSS_TP(printLvl, printLvl, m_temperature, m_pressurePA);

    /*
     *        Prepare the problem data:
     *          ->nondimensionalize the free energies using
     *            the divisor, R * T
     */
    vcs_nondim_TP();
    /*
     * Prep the fe field
     */
    vcs_fePrep_TP();

    /*
     *        Solve the problem at a fixed Temperature and Pressure
     * (all information concerning Temperature and Pressure has already
     *  been derived. The free energies are now in dimensionless form.)
     */
    iStab  = vcs_solve_phaseStability(iphase, ifunc, feStable, printLvl);


    /*
     *        Redimensionalize the free energies using
     *        the reverse of vcs_nondim to add back units.
     */
    vcs_redim_TP();

    /*
    vcs_VolPhase *Vphase = m_VolPhaseList[iphase];

    std::vector<double> mfPop = Vphase->moleFractions();
    int nsp = Vphase->nSpecies();

    vcs_VolPhase *VPphase = vprob->VPhaseList[iphase];
    int kstart = Vphase->spGlobalIndexVCS(0);
    for (int k = 0; k < nsp; k++) {
      vprob->mf[kstart + k] = mfPop[k];
    }
    VPphase->setMoleFractionsState(Vphase->totalMoles(),
                   VCS_DATA_PTR(Vphase->moleFractions()),
                   VCS_STATECALC_TMP);
    */
    vcs_prob_update(vprob);
    /*
     *        Return the convergence success flag.
     */
    return iStab;
}

int VCS_SOLVE::vcs_solve_phaseStability(const int iph, const int ifunc,
                                        double& funcVal,
                                        int printLvl)
{
    double test = -1.0E-10;
    bool usedZeroedSpecies;
    // std::vector<size_t> phasePopPhaseIDs(0);
    int iStab = 0;

    std::vector<double> sm(m_numElemConstraints*m_numElemConstraints, 0.0);
    std::vector<double> ss(m_numElemConstraints, 0.0);
    std::vector<double> sa(m_numElemConstraints, 0.0);

    std::vector<double> aw(m_numSpeciesTot, 0.0);
    std::vector<double> wx(m_numElemConstraints, 0.0);


    vcs_basopt(false, VCS_DATA_PTR(aw), VCS_DATA_PTR(sa), VCS_DATA_PTR(sm), VCS_DATA_PTR(ss),
               test, &usedZeroedSpecies);
    vcs_evaluate_speciesType();

    vcs_dfe(VCS_STATECALC_OLD, 0, 0, m_numSpeciesRdc);
    if (printLvl > 3) {
        vcs_printSpeciesChemPot(VCS_STATECALC_OLD);
    }
    vcs_deltag(0, true, VCS_STATECALC_OLD);

    if (printLvl > 3) {
        vcs_printDeltaG(VCS_STATECALC_OLD);
    }
    vcs_dcopy(VCS_DATA_PTR(m_deltaGRxn_Deficient), VCS_DATA_PTR(m_deltaGRxn_old), m_numRxnRdc);
    // phasePopPhaseIDs.clear();
    // vcs_popPhaseID(phasePopPhaseIDs);
    funcVal = vcs_phaseStabilityTest(iph);
    if (funcVal > 0.0) {
        iStab = 1;
    } else {
        iStab = 0;
    }

    return iStab;
}

}
