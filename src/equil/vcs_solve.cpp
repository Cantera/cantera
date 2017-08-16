/*!
 * @file vcs_solve.cpp Implementation file for the internal class that holds
 *     the problem.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/equil/vcs_solve.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/stringUtils.h"
#include "cantera/equil/vcs_VolPhase.h"
#include "cantera/equil/vcs_species_thermo.h"
#include "cantera/base/clockWC.h"
#include "cantera/equil/MultiPhase.h"

using namespace std;

namespace Cantera
{

int vcs_timing_print_lvl = 1;

VCS_SOLVE::VCS_SOLVE(size_t nspecies_, size_t nelements_, size_t nphase_) :
    prob_type(VCS_PROBTYPE_TP),
    nspecies(nspecies_),
    ne(0),
    NE0(0),
    NPhase(nphase_),
    T(298.15),
    PresPA(1.0),
    Vol(0.0),
    // The default is to not expect an initial estimate  of the species
    // concentrations
    iest(-1),
    tolmaj(1.0E-8),
    tolmin(1.0E-6),
    m_Iterations(0),
    m_NumBasisOptimizations(0),
    m_printLvl(0),
    vcs_debug_print_lvl(0),
    NSPECIES0(nspecies_),
    NPHASE0(nphase_),
    m_numSpeciesTot(nspecies_),
    m_numElemConstraints(0),
    m_numComponents(0),
    m_numRxnTot(0),
    m_numSpeciesRdc(0),
    m_numRxnRdc(0),
    m_numRxnMinorZeroed(0),
    m_numPhases(0),
    m_doEstimateEquil(0),
    m_totalMolNum(0.0),
    m_temperature(0.0),
    m_pressurePA(0.0),
    m_tolmaj(0.0),
    m_tolmin(0.0),
    m_tolmaj2(0.0),
    m_tolmin2(0.0),
    m_unitsState(VCS_DIMENSIONAL_G),
    m_totalMoleScale(1.0),
    m_useActCoeffJac(0),
    m_totalVol(0.0),
    m_Faraday_dim(ElectronCharge * Avogadro),
    m_VCount(0),
    m_debug_print_lvl(0),
    m_timing_print_lvl(1)
{
    m_gibbsSpecies.resize(nspecies, 0.0);
    w.resize(nspecies, 0.0);
    mf.resize(nspecies, 0.0);
    SpeciesUnknownType.resize(nspecies, VCS_SPECIES_TYPE_MOLNUM);
    VolPM.resize(nspecies, 0.0);
    PhaseID.resize(nspecies, npos);
    WtSpecies.resize(nspecies, 0.0);
    Charge.resize(nspecies, 0.0);
    SpeciesThermo.resize(nspecies,0);
    for (size_t kspec = 0; kspec < nspecies; kspec++) {
        SpeciesThermo[kspec] = new VCS_SPECIES_THERMO(0, 0);
    }
    VPhaseList.resize(NPhase, 0);
    for (size_t iphase = 0; iphase < NPhase; iphase++) {
        VPhaseList[iphase] = new vcs_VolPhase();
    }

    string ser = "VCS_SOLVE: ERROR:\n\t";
    if (nspecies_ <= 0) {
        plogf("%s Number of species is nonpositive\n", ser);
        throw CanteraError("VCS_SOLVE()", ser +
                           " Number of species is nonpositive\n");
    }
    if (nphase_ <= 0) {
        plogf("%s Number of phases is nonpositive\n", ser);
        throw CanteraError("VCS_SOLVE()", ser +
                           " Number of species is nonpositive\n");
    }

    /*
     * We will initialize sc[] to note the fact that it needs to be
     * filled with meaningful information.
     */
    m_scSize.resize(NSPECIES0, 0.0);
    m_spSize.resize(NSPECIES0, 1.0);
    m_SSfeSpecies.resize(NSPECIES0, 0.0);
    m_feSpecies_new.resize(NSPECIES0, 0.0);
    m_molNumSpecies_old.resize(NSPECIES0, 0.0);
    m_speciesUnknownType.resize(NSPECIES0, VCS_SPECIES_TYPE_MOLNUM);
    m_deltaMolNumPhase.resize(NPHASE0, NSPECIES0, 0.0);
    m_phaseParticipation.resize(NPHASE0, NSPECIES0, 0);
    m_phasePhi.resize(NPHASE0, 0.0);
    m_molNumSpecies_new.resize(NSPECIES0, 0.0);
    m_deltaGRxn_new.resize(NSPECIES0, 0.0);
    m_deltaGRxn_old.resize(NSPECIES0, 0.0);
    m_deltaGRxn_Deficient.resize(NSPECIES0, 0.0);
    m_deltaGRxn_tmp.resize(NSPECIES0, 0.0);
    m_deltaMolNumSpecies.resize(NSPECIES0, 0.0);
    m_feSpecies_old.resize(NSPECIES0, 0.0);
    m_tPhaseMoles_old.resize(NPHASE0, 0.0);
    m_tPhaseMoles_new.resize(NPHASE0, 0.0);
    m_deltaPhaseMoles.resize(NPHASE0, 0.0);
    m_TmpPhase.resize(NPHASE0, 0.0);
    m_TmpPhase2.resize(NPHASE0, 0.0);
    TPhInertMoles.resize(NPHASE0, 0.0);

    // ind[] is an index variable that keep track of solution vector rotations.
    m_speciesMapIndex.resize(NSPECIES0, 0);
    m_speciesLocalPhaseIndex.resize(NSPECIES0, 0);

    // ir[] is an index vector that keeps track of the irxn to species mapping.
    // We can't fill it in until we know the number of c components in the
    // problem
    m_indexRxnToSpecies.resize(NSPECIES0, 0);

    // Initialize all species to be major species
    m_speciesStatus.resize(NSPECIES0, 1);

    m_SSPhase.resize(2*NSPECIES0, 0);
    m_phaseID.resize(NSPECIES0, 0);
    m_speciesName.resize(NSPECIES0);

    // space for activity coefficients for all species. Set it equal to one.
    m_actConventionSpecies.resize(NSPECIES0, 0);
    m_phaseActConvention.resize(NPHASE0, 0);
    m_lnMnaughtSpecies.resize(NSPECIES0, 0.0);
    m_actCoeffSpecies_new.resize(NSPECIES0, 1.0);
    m_actCoeffSpecies_old.resize(NSPECIES0, 1.0);
    m_wtSpecies.resize(NSPECIES0, 0.0);
    m_chargeSpecies.resize(NSPECIES0, 0.0);
    m_speciesThermoList.resize(NSPECIES0, (VCS_SPECIES_THERMO*)0);

    // Phase Info
    m_VolPhaseList.resize(NPHASE0, 0);
    for (size_t iph = 0; iph < NPHASE0; iph++) {
        m_VolPhaseList[iph] = new vcs_VolPhase(this);
    }

    // For Future expansion
    m_useActCoeffJac = true;
    if (m_useActCoeffJac) {
        m_np_dLnActCoeffdMolNum.resize(NSPECIES0, NSPECIES0, 0.0);
    }

    m_PMVolumeSpecies.resize(NSPECIES0, 0.0);

    // counters kept within vcs
    m_VCount = new VCS_COUNTERS();
    vcs_counters_init(1);

    if (vcs_timing_print_lvl == 0) {
        m_timing_print_lvl = 0;
    }

    return;
}

VCS_SOLVE::~VCS_SOLVE()
{
    vcs_delete_memory();
}

void VCS_SOLVE::vcs_delete_memory()
{
    for (size_t j = 0; j < m_VolPhaseList.size(); j++) {
        delete m_VolPhaseList[j];
        m_VolPhaseList[j] = 0;
    }

    for (size_t j = 0; j < m_speciesThermoList.size(); j++) {
        delete m_speciesThermoList[j];
        m_speciesThermoList[j] = 0;
    }

    for (size_t i = 0; i < SpeciesThermo.size(); i++) {
        delete SpeciesThermo[i];
        SpeciesThermo[i] = 0;
    }
    for (size_t iph = 0; iph < VPhaseList.size(); iph++) {
        delete VPhaseList[iph];
        VPhaseList[iph] = 0;
    }

    delete m_VCount;
    m_VCount = 0;

    NSPECIES0 = 0;
    NPHASE0 = 0;
    m_numElemConstraints = 0;
    m_numComponents = 0;

}

int VCS_SOLVE::vcs(int ifunc, int ipr, int ip1, int maxit)
{
    int retn = 0, iconv = 0;
    clockWC tickTock;

    int iprintTime = std::max(ipr, ip1);
    //iprintTime = std::min(iprintTime, m_timing_print_lvl);
    // writelog("iprintTime = {}; m_timing_print_lvl = {}; ipr = {}; ip1 = {}\n",
    //     iprintTime, m_timing_print_lvl, ipr, ip1);

    if (ifunc > 2) {
        plogf("vcs: Unrecognized value of ifunc, %d: bailing!\n",
              ifunc);
        return VCS_PUB_BAD;
    }

    if (ifunc == 0) {
        // This function is called to copy the public data and the current
        // problem specification into the current object's data structure.
        retn = vcs_prob_specifyFully();
        if (retn != 0) {
            plogf("vcs_pub_to_priv returned a bad status, %d: bailing!\n",
                  retn);
            return retn;
        }
        // Prep the problem data
        //    - adjust the identity of any phases
        //    - determine the number of components in the problem
        retn = vcs_prep_oneTime(ip1);
        if (retn != 0) {
            plogf("vcs_prep_oneTime returned a bad status, %d: bailing!\n",
                  retn);
            return retn;
        }
    }
    if (ifunc == 1) {
        // This function is called to copy the current problem into the current
        // object's data structure.
        retn = vcs_prob_specify();
        if (retn != 0) {
            plogf("vcs_prob_specify returned a bad status, %d: bailing!\n",
                  retn);
            return retn;
        }
    }
    if (ifunc != 2) {
        // Prep the problem data for this particular instantiation of
        // the problem
        retn = vcs_prep();
        if (retn != VCS_SUCCESS) {
            plogf("vcs_prep returned a bad status, %d: bailing!\n", retn);
            return retn;
        }

        // Check to see if the current problem is well posed.
        if (!vcs_wellPosed()) {
            plogf("vcs has determined the problem is not well posed: Bailing\n");
            return VCS_PUB_BAD;
        }

        // Once we have defined the global internal data structure defining the
        // problem, then we go ahead and solve the problem.
        //
        // (right now, all we do is solve fixed T, P problems. Methods for other
        // problem types will go in at this level. For example, solving for
        // fixed T, V problems will involve a 2x2 Newton's method, using loops
        // over vcs_TP() to calculate the residual and Jacobian)
        iconv = vcs_TP(ipr, ip1, maxit, T, PresPA);

        // If requested to print anything out, go ahead and do so;
        if (ipr > 0) {
            vcs_report(iconv);
        }

        vcs_prob_update();
    }

    // Report on the time if requested to do so
    double te = tickTock.secondsWC();
    m_VCount->T_Time_vcs += te;
    if (iprintTime > 0) {
        vcs_TCounters_report(m_timing_print_lvl);
    }

    // Now, destroy the private data, if requested to do so
    // FILL IN
    if (iconv < 0) {
        plogf("ERROR: FAILURE its = %d!\n", m_VCount->Its);
    } else if (iconv == 1) {
        plogf("WARNING: RANGE SPACE ERROR encountered\n");
    }
    return iconv;
}

int VCS_SOLVE::vcs_prob_specifyFully()
{
    const char* ser =
        "vcs_pub_to_priv ERROR :ill defined interface -> bailout:\n\t";

    // First Check to see whether we have room for the current problem size
    if (NSPECIES0 < nspecies) {
        plogf("%sPrivate Data is dimensioned too small\n", ser);
        return VCS_PUB_BAD;
    }
    size_t nph = NPhase;
    if (NPHASE0 < nph) {
        plogf("%sPrivate Data is dimensioned too small\n", ser);
        return VCS_PUB_BAD;
    }

    // OK, We have room. Now, transfer the integer numbers
    m_numSpeciesTot = nspecies;
    m_numSpeciesRdc = m_numSpeciesTot;

    // m_numRxnTot = number of noncomponents, also equal to the number of
    // reactions. Note, it's possible that the number of elements is greater
    // than the number of species. In that case set the number of reactions to
    // zero.
    if (ne > nspecies) {
        m_numRxnTot = 0;
    } else {
        m_numRxnTot = nspecies - ne;
    }
    m_numRxnRdc = m_numRxnTot;

    // number of minor species rxn -> all species rxn are major at the start.
    m_numRxnMinorZeroed = 0;

    // NPhase = number of phases
    m_numPhases = nph;

    m_debug_print_lvl = vcs_debug_print_lvl;

    // FormulaMatrix[] -> Copy the formula matrix over
    for (size_t i = 0; i < nspecies; i++) {
        bool nonzero = false;
        for (size_t j = 0; j < ne; j++) {
            if (FormulaMatrix(i,j) != 0.0) {
                nonzero = true;
            }
            m_formulaMatrix(i,j) = FormulaMatrix(i,j);
        }
        if (!nonzero) {
            plogf("vcs_prob_specifyFully:: species %d %s has a zero formula matrix!\n", i,
                  m_speciesName[i]);
            return VCS_PUB_BAD;
        }
    }

    // Copy over the species molecular weights
    m_wtSpecies = WtSpecies;

    // Copy over the charges
    m_chargeSpecies = Charge;

    // Copy the VCS_SPECIES_THERMO structures
    for (size_t kspec = 0; kspec < nspecies; kspec++) {
        delete m_speciesThermoList[kspec];
        VCS_SPECIES_THERMO* spf = SpeciesThermo[kspec];
        m_speciesThermoList[kspec] = spf->duplMyselfAsVCS_SPECIES_THERMO();
        if (m_speciesThermoList[kspec] == NULL) {
            plogf(" duplMyselfAsVCS_SPECIES_THERMO returned an error!\n");
            return VCS_PUB_BAD;
        }
    }

    // Copy the species unknown type
    m_speciesUnknownType = SpeciesUnknownType;

    // iest => Do we have an initial estimate of the species mole numbers ?
    m_doEstimateEquil = iest;

    // w[] -> Copy the equilibrium mole number estimate if it exists.
    if (w.size() != 0) {
        m_molNumSpecies_old = w;
    } else {
        m_doEstimateEquil = -1;
        m_molNumSpecies_old.assign(m_molNumSpecies_old.size(), 0.0);
    }

    // Formulate the Goal Element Abundance Vector
    if (gai.size() != 0) {
        for (size_t i = 0; i < ne; i++) {
            m_elemAbundancesGoal[i] = gai[i];
            if (m_elType[i] == VCS_ELEM_TYPE_LATTICERATIO && m_elemAbundancesGoal[i] < 1.0E-10) {
                m_elemAbundancesGoal[i] = 0.0;
            }
        }
    } else {
        if (m_doEstimateEquil == 0) {
            double sum = 0;
            for (size_t j = 0; j < ne; j++) {
                m_elemAbundancesGoal[j] = 0.0;
                for (size_t kspec = 0; kspec < nspecies; kspec++) {
                    if (m_speciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
                        sum += m_molNumSpecies_old[kspec];
                        m_elemAbundancesGoal[j] += m_formulaMatrix(kspec,j) * m_molNumSpecies_old[kspec];
                    }
                }
                if (m_elType[j] == VCS_ELEM_TYPE_LATTICERATIO && m_elemAbundancesGoal[j] < 1.0E-10 * sum) {
                    m_elemAbundancesGoal[j] = 0.0;
                }
            }
        } else {
            plogf("%sElement Abundances, m_elemAbundancesGoal[], not specified\n", ser);
            return VCS_PUB_BAD;
        }
    }

    // zero out values that will be filled in later
    //
    // TPhMoles[]      -> Untouched here. These will be filled in vcs_prep.c
    // TPhMoles1[]
    // DelTPhMoles[]
    //
    // T, Pres, copy over here
    if (T > 0.0) {
        m_temperature = T;
    } else {
        m_temperature = 293.15;
    }
    if (PresPA > 0.0) {
        m_pressurePA = PresPA;
    } else {
        m_pressurePA = OneAtm;
    }

    // TPhInertMoles[] -> must be copied over here
    for (size_t iph = 0; iph < nph; iph++) {
        vcs_VolPhase* Vphase = VPhaseList[iph];
        TPhInertMoles[iph] = Vphase->totalMolesInert();
    }

    // tolerance requirements -> copy them over here and later
    m_tolmaj = tolmaj;
    m_tolmin = tolmin;
    m_tolmaj2 = 0.01 * m_tolmaj;
    m_tolmin2 = 0.01 * m_tolmin;

    // m_speciesIndexVector[] is an index variable that keep track of solution
    // vector rotations.
    for (size_t i = 0; i < nspecies; i++) {
        m_speciesMapIndex[i] = i;
    }

    // IndEl[] is an index variable that keep track of element vector rotations.
    for (size_t i = 0; i < ne; i++) {
        m_elementMapIndex[i] = i;
    }

    // Define all species to be major species, initially.
    for (size_t i = 0; i < nspecies; i++) {
        m_speciesStatus[i] = VCS_SPECIES_MAJOR;
    }

    // PhaseID: Fill in the species to phase mapping. Check for bad values at
    // the same time.
    if (PhaseID.size() != 0) {
        std::vector<size_t> numPhSp(nph, 0);
        for (size_t kspec = 0; kspec < nspecies; kspec++) {
            size_t iph = PhaseID[kspec];
            if (iph >= nph) {
                plogf("%sSpecies to Phase Mapping, PhaseID, has a bad value\n",
                      ser);
                plogf("\tPhaseID[%d] = %d\n", kspec, iph);
                plogf("\tAllowed values: 0 to %d\n", nph - 1);
                return VCS_PUB_BAD;
            }
            m_phaseID[kspec] = PhaseID[kspec];
            m_speciesLocalPhaseIndex[kspec] = numPhSp[iph];
            numPhSp[iph]++;
        }
        for (size_t iph = 0; iph < nph; iph++) {
            vcs_VolPhase* Vphase = VPhaseList[iph];
            if (numPhSp[iph] != Vphase->nSpecies()) {
                plogf("%sNumber of species in phase %d, %s, doesn't match (%d != %d) [vphase = %d]\n",
                      ser, iph, Vphase->PhaseName, numPhSp[iph], Vphase->nSpecies(), (size_t) Vphase);
                return VCS_PUB_BAD;
            }
        }
    } else {
        if (m_numPhases == 1) {
            for (size_t kspec = 0; kspec < nspecies; kspec++) {
                m_phaseID[kspec] = 0;
                m_speciesLocalPhaseIndex[kspec] = kspec;
            }
        } else {
            plogf("%sSpecies to Phase Mapping, PhaseID, is not defined\n", ser);
            return VCS_PUB_BAD;
        }
    }

    for (size_t i = 0; i < ne; i++) {
        if (m_elType[i] == VCS_ELEM_TYPE_CHARGENEUTRALITY) {
            if (m_elemAbundancesGoal[i] != 0.0) {
                if (fabs(m_elemAbundancesGoal[i]) > 1.0E-9) {
                    throw CanteraError("VCS_SOLVE::vcs_prob_specifyFully",
                            "Charge neutrality condition {} is signicantly "
                            "nonzero, {}. Giving up",
                            m_elementName[i], m_elemAbundancesGoal[i]);
                } else {
                    if (m_debug_print_lvl >= 2) {
                        plogf("Charge neutrality condition %s not zero, %g. Setting it zero\n",
                              m_elementName[i], m_elemAbundancesGoal[i]);
                    }
                    m_elemAbundancesGoal[i] = 0.0;
                }

            }
        }
    }

    // Copy over the species names
    for (size_t i = 0; i < nspecies; i++) {
        m_speciesName[i] = m_mix->speciesName(i);
    }

    // Copy over all of the phase information. Use the object's assignment
    // operator
    for (size_t iph = 0; iph < nph; iph++) {
        *m_VolPhaseList[iph] = *VPhaseList[iph];

        // Fix up the species thermo pointer in the vcs_SpeciesThermo object. It
        // should point to the species thermo pointer in the private data space.
        vcs_VolPhase* Vphase = m_VolPhaseList[iph];
        for (size_t k = 0; k < Vphase->nSpecies(); k++) {
            vcs_SpeciesProperties* sProp = Vphase->speciesProperty(k);
            size_t kT = Vphase->spGlobalIndexVCS(k);
            sProp->SpeciesThermo = m_speciesThermoList[kT];
        }
    }

    // Specify the Activity Convention information
    for (size_t iph = 0; iph < nph; iph++) {
        vcs_VolPhase* Vphase = m_VolPhaseList[iph];
        m_phaseActConvention[iph] = Vphase->p_activityConvention;
        if (Vphase->p_activityConvention != 0) {
            // We assume here that species 0 is the solvent. The solvent isn't
            // on a unity activity basis The activity for the solvent assumes
            // that the it goes to one as the species mole fraction goes to one;
            // i.e., it's really on a molarity framework. So
            // SpecLnMnaught[iSolvent] = 0.0, and the loop below starts at 1,
            // not 0.
            size_t iSolvent = Vphase->spGlobalIndexVCS(0);
            double mnaught = m_wtSpecies[iSolvent] / 1000.;
            for (size_t k = 1; k < Vphase->nSpecies(); k++) {
                size_t kspec = Vphase->spGlobalIndexVCS(k);
                m_actConventionSpecies[kspec] = Vphase->p_activityConvention;
                m_lnMnaughtSpecies[kspec] = log(mnaught);
            }
        }
    }

    // Copy the title info
    if (Title.size() == 0) {
        m_title = "Unspecified Problem Title";
    } else {
        m_title = Title;
    }

    // Copy the volume info
    m_totalVol = Vol;
    if (m_PMVolumeSpecies.size() != 0) {
        m_PMVolumeSpecies = VolPM;
    }

    // Return the success flag
    return VCS_SUCCESS;
}

int VCS_SOLVE::vcs_prob_specify()
{
    string yo("vcs_prob_specify ERROR: ");
    int retn = VCS_SUCCESS;

    m_temperature = T;
    m_pressurePA = PresPA;
    m_doEstimateEquil = iest;
    m_totalVol = Vol;
    m_tolmaj = tolmaj;
    m_tolmin = tolmin;
    m_tolmaj2 = 0.01 * m_tolmaj;
    m_tolmin2 = 0.01 * m_tolmin;

    for (size_t kspec = 0; kspec < m_numSpeciesTot; ++kspec) {
        size_t k = m_speciesMapIndex[kspec];
        m_molNumSpecies_old[kspec] = w[k];
        m_molNumSpecies_new[kspec] = mf[k];
        m_feSpecies_old[kspec] = m_gibbsSpecies[k];
    }

    // Transfer the element abundance goals to the solve object
    for (size_t i = 0; i < m_numElemConstraints; i++) {
        size_t j = m_elementMapIndex[i];
        m_elemAbundancesGoal[i] = gai[j];
    }

    // Try to do the best job at guessing at the title
    if (Title.size() == 0) {
        if (m_title.size() == 0) {
            m_title = "Unspecified Problem Title";
        }
    } else {
        m_title = Title;
    }

    // Copy over the phase information. For each entry in the phase structure,
    // determine if that entry can change from its initial value Either copy
    // over the new value or create an error condition.
    bool status_change = false;
    for (size_t iph = 0; iph < m_numPhases; iph++) {
        vcs_VolPhase* vPhase = m_VolPhaseList[iph];
        vcs_VolPhase* pub_phase_ptr = VPhaseList[iph];

        if (vPhase->VP_ID_ != pub_phase_ptr->VP_ID_) {
            plogf("%sPhase numbers have changed:%d %d\n",
                  yo, vPhase->VP_ID_, pub_phase_ptr->VP_ID_);
            retn = VCS_PUB_BAD;
        }

        if (vPhase->m_singleSpecies != pub_phase_ptr->m_singleSpecies) {
            plogf("%sSingleSpecies value have changed:%d %d\n",
                  yo, vPhase->m_singleSpecies, pub_phase_ptr->m_singleSpecies);
            retn = VCS_PUB_BAD;
        }

        if (vPhase->m_gasPhase != pub_phase_ptr->m_gasPhase) {
            plogf("%sGasPhase value have changed:%d %d\n",
                  yo, vPhase->m_gasPhase, pub_phase_ptr->m_gasPhase);
            retn = VCS_PUB_BAD;
        }

        vPhase->m_eqnState = pub_phase_ptr->m_eqnState;

        if (vPhase->nSpecies() != pub_phase_ptr->nSpecies()) {
            plogf("%sNVolSpecies value have changed:%d %d\n",
                  yo, vPhase->nSpecies(), pub_phase_ptr->nSpecies());
            retn = VCS_PUB_BAD;
        }

        if (vPhase->PhaseName != pub_phase_ptr->PhaseName) {
            plogf("%sPhaseName value have changed:%s %s\n",
                  yo, vPhase->PhaseName, pub_phase_ptr->PhaseName);
            retn = VCS_PUB_BAD;
        }

        if (vPhase->totalMolesInert() != pub_phase_ptr->totalMolesInert()) {
            status_change = true;
        }

        // Copy over the number of inert moles if it has changed.
        TPhInertMoles[iph] = pub_phase_ptr->totalMolesInert();
        vPhase->setTotalMolesInert(pub_phase_ptr->totalMolesInert());
        if (TPhInertMoles[iph] > 0.0) {
            vPhase->setExistence(2);
            vPhase->m_singleSpecies = false;
        }

        // Copy over the interfacial potential
        double phi = pub_phase_ptr->electricPotential();
        vPhase->setElectricPotential(phi);
    }

    if (status_change) {
        vcs_SSPhase();
    }

    // Calculate the total number of moles in all phases.
    vcs_tmoles();
    return retn;
}

int VCS_SOLVE::vcs_prob_update()
{
    size_t k1 = 0;
    vcs_tmoles();
    m_totalVol = vcs_VolTotal(m_temperature, m_pressurePA,
                              &m_molNumSpecies_old[0], &m_PMVolumeSpecies[0]);

    for (size_t i = 0; i < m_numSpeciesTot; ++i) {
        // Find the index of I in the index vector, m_speciesIndexVector[]. Call
        // it K1 and continue.
        for (size_t j = 0; j < m_numSpeciesTot; ++j) {
            k1 = j;
            if (m_speciesMapIndex[j] == i) {
                break;
            }
        }

        // Switch the species data back from K1 into I
        if (SpeciesUnknownType[i] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
            w[i] = m_molNumSpecies_old[k1];
        } else {
            w[i] = 0.0;
        }
        m_gibbsSpecies[i] = m_feSpecies_old[k1];
        VolPM[i] = m_PMVolumeSpecies[k1];
    }

    T = m_temperature;
    PresPA = m_pressurePA;
    Vol = m_totalVol;
    size_t kT = 0;
    for (size_t iph = 0; iph < NPhase; iph++) {
        vcs_VolPhase* pubPhase = VPhaseList[iph];
        vcs_VolPhase* vPhase = m_VolPhaseList[iph];
        pubPhase->setTotalMolesInert(vPhase->totalMolesInert());
        pubPhase->setTotalMoles(vPhase->totalMoles());
        pubPhase->setElectricPotential(vPhase->electricPotential());
        double sumMoles = pubPhase->totalMolesInert();
        pubPhase->setMoleFractionsState(vPhase->totalMoles(),
                                        &vPhase->moleFractions()[0],
                                        VCS_STATECALC_TMP);
        const vector_fp & mfVector = pubPhase->moleFractions();
        for (size_t k = 0; k < pubPhase->nSpecies(); k++) {
            kT = pubPhase->spGlobalIndexVCS(k);
            mf[kT] = mfVector[k];
            if (pubPhase->phiVarIndex() == k) {
                k1 = vPhase->spGlobalIndexVCS(k);
                double tmp = m_molNumSpecies_old[k1];
                if (! vcs_doubleEqual(pubPhase->electricPotential() , tmp)) {
                    throw CanteraError("VCS_SOLVE::vcs_prob_update",
                            "We have an inconsistency in voltage, {} {}",
                            pubPhase->electricPotential(), tmp);
                }
            }

            if (! vcs_doubleEqual(mf[kT], vPhase->molefraction(k))) {
                throw CanteraError("VCS_SOLVE::vcs_prob_update",
                        "We have an inconsistency in mole fraction, {} {}",
                        mf[kT], vPhase->molefraction(k));
            }
            if (pubPhase->speciesUnknownType(k) != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
                sumMoles += w[kT];
            }
        }
        if (! vcs_doubleEqual(sumMoles, vPhase->totalMoles())) {
            throw CanteraError("VCS_SOLVE::vcs_prob_update",
                            "We have an inconsistency in total moles, {} {}",
                            sumMoles, pubPhase->totalMoles());
        }
    }

    m_Iterations = m_VCount->Its;
    m_NumBasisOptimizations = m_VCount->Basis_Opts;
    return VCS_SUCCESS;
}

void VCS_SOLVE::vcs_counters_init(int ifunc)
{
    m_VCount->Its = 0;
    m_VCount->Basis_Opts = 0;
    m_VCount->Time_vcs_TP = 0.0;
    m_VCount->Time_basopt = 0.0;
    if (ifunc) {
        m_VCount->T_Its = 0;
        m_VCount->T_Basis_Opts = 0;
        m_VCount->T_Calls_Inest = 0;
        m_VCount->T_Calls_vcs_TP = 0;
        m_VCount->T_Time_vcs_TP = 0.0;
        m_VCount->T_Time_basopt = 0.0;
        m_VCount->T_Time_inest = 0.0;
        m_VCount->T_Time_vcs = 0.0;
    }
}

double VCS_SOLVE::vcs_VolTotal(const double tkelvin, const double pres,
                               const double w[], double volPM[])
{
    double VolTot = 0.0;
    for (size_t iphase = 0; iphase < m_numPhases; iphase++) {
        vcs_VolPhase* Vphase = m_VolPhaseList[iphase];
        Vphase->setState_TP(tkelvin, pres);
        Vphase->setMolesFromVCS(VCS_STATECALC_OLD, w);
        double Volp = Vphase->sendToVCS_VolPM(volPM);
        VolTot += Volp;
    }
    return VolTot;
}

void VCS_SOLVE::disableTiming() {
    vcs_timing_print_lvl = 0;
}

}
