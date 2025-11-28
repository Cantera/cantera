/*!
 * @file vcs_solve.cpp Implementation file for the internal class that holds
 *     the problem.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/equil/vcs_solve.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/stringUtils.h"
#include "cantera/equil/vcs_VolPhase.h"
#include "cantera/equil/MultiPhase.h"
#include "cantera/thermo/speciesThermoTypes.h"
#include "cantera/thermo/ThermoPhase.h"

using namespace std;

namespace Cantera
{

namespace {

void printProgress(const vector<string>& spName, const vector<double>& soln,
                   const vector<double>& ff)
{
    double sum = 0.0;
    plogf(" --- Summary of current progress:\n");
    plogf(" ---                   Name           Moles  -       SSGibbs \n");
    plogf(" -------------------------------------------------------------------------------------\n");
    for (size_t k = 0; k < soln.size(); k++) {
        plogf(" ---      %20s %12.4g  - %12.4g\n", spName[k], soln[k], ff[k]);
        sum += soln[k] * ff[k];
    }
    plogf(" ---  Total sum to be minimized = %g\n", sum);
}

} // anonymous namespace

VCS_SOLVE::VCS_SOLVE(MultiPhase* mphase, int printLvl) :
    m_printLvl(printLvl),
    m_mix(mphase),
    m_nsp(mphase->nSpecies()),
    m_numSpeciesRdc(mphase->nSpecies()),
    m_numPhases(mphase->nPhases()),
    m_temperature(mphase->temperature()),
    m_pressurePA(mphase->pressure()),
    m_Faraday_dim(Faraday / (m_temperature * GasConstant))
{
    string ser = "VCS_SOLVE: ERROR:\n\t";
    if (m_nsp <= 0) {
        plogf("%s Number of species is nonpositive\n", ser);
        throw CanteraError("VCS_SOLVE::VCS_SOLVE", ser +
                           " Number of species is nonpositive\n");
    }
    if (m_numPhases <= 0) {
        plogf("%s Number of phases is nonpositive\n", ser);
        throw CanteraError("VCS_SOLVE::VCS_SOLVE", ser +
                           " Number of phases is nonpositive\n");
    }

    /*
     * We will initialize sc[] to note the fact that it needs to be
     * filled with meaningful information.
     */
    m_scSize.resize(m_nsp, 0.0);
    m_spSize.resize(m_nsp, 1.0);
    m_SSfeSpecies.resize(m_nsp, 0.0);
    m_feSpecies_new.resize(m_nsp, 0.0);
    m_molNumSpecies_old.resize(m_nsp, 0.0);
    m_speciesUnknownType.resize(m_nsp, VCS_SPECIES_TYPE_MOLNUM);
    m_deltaMolNumPhase.resize(m_numPhases, m_nsp, 0.0);
    m_phaseParticipation.resize(m_numPhases, m_nsp, 0);
    m_phasePhi.resize(m_numPhases, 0.0);
    m_molNumSpecies_new.resize(m_nsp, 0.0);
    m_deltaGRxn_new.resize(m_nsp, 0.0);
    m_deltaGRxn_old.resize(m_nsp, 0.0);
    m_deltaGRxn_Deficient.resize(m_nsp, 0.0);
    m_deltaGRxn_tmp.resize(m_nsp, 0.0);
    m_deltaMolNumSpecies.resize(m_nsp, 0.0);
    m_feSpecies_old.resize(m_nsp, 0.0);
    m_tPhaseMoles_old.resize(m_numPhases, 0.0);
    m_tPhaseMoles_new.resize(m_numPhases, 0.0);
    m_deltaPhaseMoles.resize(m_numPhases, 0.0);
    m_TmpPhase.resize(m_numPhases, 0.0);
    m_TmpPhase2.resize(m_numPhases, 0.0);

    // ind[] is an index variable that keep track of solution vector rotations.
    m_speciesMapIndex.resize(m_nsp, 0);
    m_speciesLocalPhaseIndex.resize(m_nsp, 0);

    // ir[] is an index vector that keeps track of the irxn to species mapping.
    // We can't fill it in until we know the number of c components in the
    // problem
    m_indexRxnToSpecies.resize(m_nsp, 0);

    // Initialize all species to be major species
    m_speciesStatus.resize(m_nsp, VCS_SPECIES_MAJOR);

    m_SSPhase.resize(2*m_nsp, 0);
    m_phaseID.resize(m_nsp, 0);
    m_speciesName.resize(m_nsp);

    // space for activity coefficients for all species. Set it equal to one.
    m_actConventionSpecies.resize(m_nsp, 0);
    m_phaseActConvention.resize(m_numPhases, 0);
    m_lnMnaughtSpecies.resize(m_nsp, 0.0);
    m_actCoeffSpecies_new.resize(m_nsp, 1.0);
    m_actCoeffSpecies_old.resize(m_nsp, 1.0);
    m_wtSpecies.resize(m_nsp, 0.0);
    m_chargeSpecies.resize(m_nsp, 0.0);

    // Phase Info
    m_VolPhaseList.resize(m_numPhases);
    for (size_t iph = 0; iph < m_numPhases; iph++) {
        m_VolPhaseList[iph] = make_unique<vcs_VolPhase>(this);
    }

    // For Future expansion
    m_useActCoeffJac = true;
    if (m_useActCoeffJac) {
        m_np_dLnActCoeffdMolNum.resize(m_nsp, m_nsp, 0.0);
    }

    m_PMVolumeSpecies.resize(m_nsp, 0.0);

    // counters kept within vcs
    m_VCount = new VCS_COUNTERS();
    vcs_counters_init(1);

    // Loop over the phases, transferring pertinent information
    int kT = 0;
    for (size_t iphase = 0; iphase < m_numPhases; iphase++) {
        // Get the ThermoPhase object - assume volume phase
        ThermoPhase* tPhase = &mphase->phase(iphase);
        size_t nelem = tPhase->nElements();

        // Query Cantera for the equation of state type of the current phase.
        string eos = tPhase->type();

        // Find out the number of species in the phase
        size_t nSpPhase = tPhase->nSpecies();
        // Find out the name of the phase
        string phaseName = tPhase->name();

        // Call the basic vcs_VolPhase creation routine.
        // Properties set here:
        //    ->PhaseNum = phase number in the thermo problem
        //    ->NumSpecies = number of species in the phase
        //    ->PhaseName  = Name of the phase
        vcs_VolPhase* VolPhase = m_VolPhaseList[iphase].get();
        VolPhase->resize(iphase, nSpPhase, nelem, phaseName.c_str());

        // Tell the vcs_VolPhase pointer about cantera
        VolPhase->setPtrThermoPhase(tPhase);
        VolPhase->setTotalMoles(0.0);

        // Set the electric potential of the volume phase from the
        // ThermoPhase object's value.
        VolPhase->setElectricPotential(tPhase->electricPotential());

        // Query the ThermoPhase object to find out what convention
        // it uses for the specification of activity and Standard State.
        VolPhase->p_activityConvention = tPhase->activityConvention();

        // Assign the value of eqn of state. Handle conflicts here.
        if (eos == "ideal-gas") {
            VolPhase->m_eqnState = VCS_EOS_IDEAL_GAS;
        } else if (eos == "fixed-stoichiometry") {
            VolPhase->m_eqnState = VCS_EOS_STOICH_SUB;
        } else if (eos == "ideal-condensed") {
            VolPhase->m_eqnState = VCS_EOS_IDEAL_SOLN;
        } else if (tPhase->nDim() != 3) {
            throw CanteraError("VCS_SOLVE::VCS_SOLVE",
                               "Surface/edge phase not handled yet.");
        } else {
            if (m_printLvl > 1) {
                writelog("Unknown Cantera EOS to VCSnonideal: '{}'\n", eos);
            }
            VolPhase->m_eqnState = VCS_EOS_UNK_CANTERA;
        }

        // Transfer all of the element information from the ThermoPhase object
        // to the vcs_VolPhase object. Also decide whether we need a new charge
        // neutrality element in the phase to enforce a charge neutrality
        // constraint. We also decide whether this is a single species phase
        // with the voltage being the independent variable setting the chemical
        // potential of the electrons.
        VolPhase->transferElementsFM(tPhase);

        // Combine the element information in the vcs_VolPhase
        // object into the vprob object.
        addPhaseElements(VolPhase);
        VolPhase->setState_TP(m_temperature, m_pressurePA);

        // Loop through each species in the current phase
        for (size_t k = 0; k < nSpPhase; k++) {
            // Obtain the molecular weight of the species from the
            // ThermoPhase object
            m_wtSpecies[kT] = tPhase->molecularWeight(k);

            // Obtain the charges of the species from the ThermoPhase object
            m_chargeSpecies[kT] = tPhase->charge(k);

            // Set the phaseid of the species
            m_phaseID[kT] = iphase;

            // Transfer the type of unknown
            m_speciesUnknownType[kT] = VolPhase->speciesUnknownType(k);
            if (m_speciesUnknownType[kT] == VCS_SPECIES_TYPE_MOLNUM) {
                // Set the initial number of kmoles of the species
                m_molNumSpecies_old[kT] = mphase->speciesMoles(kT);
             } else if (m_speciesUnknownType[kT] == VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
                m_molNumSpecies_old[kT] = tPhase->electricPotential();
             } else {
               throw CanteraError("VCS_SOLVE::VCS_SOLVE",
                                  "Unknown species type: {}", m_speciesUnknownType[kT]);
             }

            // Transfer the species information from the
            // volPhase structure to the VPROB structure
            // This includes:
            //      FormulaMatrix[][]
            //      VolPhase->IndSpecies[]
            addOnePhaseSpecies(VolPhase, k, kT);

            kT++;
        }

        VolPhase->setMolesFromVCS(VCS_STATECALC_OLD, &m_molNumSpecies_old[0]);
    }

    // Work arrays used by vcs_basopt
    m_sm.assign(m_nelem * m_nelem, 0.0);
    m_ss.assign(m_nelem, 0.0);
    m_sa.assign(m_nelem, 0.0);
    m_aw.assign(m_nsp, 0.0);
    m_wx.assign(m_nelem, 0.0);

    // Transfer initial element abundances based on the species mole numbers
    for (size_t j = 0; j < m_nelem; j++) {
        for (size_t kspec = 0; kspec < m_nsp; kspec++) {
            if (m_speciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
                m_elemAbundancesGoal[j] += m_formulaMatrix(kspec,j) * m_molNumSpecies_old[kspec];
            }
        }
    }

    // Printout the species information: PhaseID's and mole nums
    if (m_printLvl > 1) {
        writeline('=', 80, true, true);
        writeline('=', 16, false);
        plogf(" Cantera_to_vprob: START OF PROBLEM STATEMENT ");
        writeline('=', 20);
        writeline('=', 80);
        plogf("             Phase IDs of species\n");
        plogf("            species     phaseID        phaseName   ");
        plogf(" Initial_Estimated_kMols\n");
        for (size_t i = 0; i < m_nsp; i++) {
            size_t iphase = m_phaseID[i];

            vcs_VolPhase* VolPhase = m_VolPhaseList[iphase].get();
            plogf("%16s      %5d   %16s", mphase->speciesName(i).c_str(), iphase,
                  VolPhase->PhaseName.c_str());
            if (m_speciesUnknownType[i] == VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
                plogf("     Volts = %-10.5g\n", m_molNumSpecies_old[i]);
            } else {
                plogf("             %-10.5g\n", m_molNumSpecies_old[i]);
            }
        }

        // Printout of the Phase structure information
        writeline('-', 80, true, true);
        plogf("             Information about phases\n");
        plogf("  PhaseName      PhaseNum SingSpec   EqnState        NumSpec     "
              "Tmoles(kmol)\n");

        for (size_t iphase = 0; iphase < m_numPhases; iphase++) {
            vcs_VolPhase* VolPhase = m_VolPhaseList[iphase].get();
            plogf("%16s %8d %8d %16s %8d ", VolPhase->PhaseName.c_str(),
                  VolPhase->VP_ID_, VolPhase->m_singleSpecies,
                  VolPhase->eos_name(), VolPhase->nSpecies());
            plogf("%16e\n", VolPhase->totalMoles());
        }

        writeline('=', 80, true, true);
        writeline('=', 16, false);
        plogf(" Cantera_to_vprob: END OF PROBLEM STATEMENT ");
        writeline('=', 20);
        writeline('=', 80);
        plogf("\n");
    }

    // m_speciesIndexVector[] is an index variable that keep track of solution
    // vector rotations.
    for (size_t i = 0; i < m_nsp; i++) {
        m_speciesMapIndex[i] = i;
    }

    // IndEl[] is an index variable that keep track of element vector rotations.
    for (size_t i = 0; i < m_nelem; i++) {
        m_elementMapIndex[i] = i;
    }

    // Fill in the species to phase mapping. Check for bad values at the same
    // time.
    vector<size_t> numPhSp(m_numPhases, 0);
    for (size_t kspec = 0; kspec < m_nsp; kspec++) {
        size_t iph = m_phaseID[kspec];
        if (iph >= m_numPhases) {
            throw CanteraError("VCS_SOLVE::VCS_SOLVE",
                "Species to Phase Mapping, PhaseID, has a bad value\n"
                "\tm_phaseID[{}] = {}\n"
                "Allowed values: 0 to {}", kspec, iph, m_numPhases - 1);
        }
        m_phaseID[kspec] = m_phaseID[kspec];
        m_speciesLocalPhaseIndex[kspec] = numPhSp[iph];
        numPhSp[iph]++;
    }
    for (size_t iph = 0; iph < m_numPhases; iph++) {
        vcs_VolPhase* Vphase = m_VolPhaseList[iph].get();
        if (numPhSp[iph] != Vphase->nSpecies()) {
            throw CanteraError("VCS_SOLVE::VCS_SOLVE",
                "Number of species in phase {}, {}, doesn't match ({} != {}) [vphase = {}]",
                ser, iph, Vphase->PhaseName, numPhSp[iph], Vphase->nSpecies(), (size_t) Vphase);
        }
    }

    for (size_t i = 0; i < m_nelem; i++) {
        if (m_elType[i] == VCS_ELEM_TYPE_CHARGENEUTRALITY) {
            if (m_elemAbundancesGoal[i] != 0.0) {
                if (fabs(m_elemAbundancesGoal[i]) > 1.0E-9) {
                    throw CanteraError("VCS_SOLVE::VCS_SOLVE",
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
    for (size_t i = 0; i < m_nsp; i++) {
        m_speciesName[i] = m_mix->speciesName(i);
    }

    // Specify the Activity Convention information
    for (size_t iph = 0; iph < m_numPhases; iph++) {
        vcs_VolPhase* Vphase = m_VolPhaseList[iph].get();
        m_phaseActConvention[iph] = Vphase->p_activityConvention;
        if (Vphase->p_activityConvention != 0) {
            // We assume here that species 0 is the solvent. The solvent isn't
            // on a unity activity basis The activity for the solvent assumes
            // that the it goes to one as the species mole fraction goes to one;
            // that is, it's really on a molarity framework. So
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

}

VCS_SOLVE::~VCS_SOLVE()
{
    delete m_VCount;
}

bool VCS_SOLVE::vcs_popPhasePossible(const size_t iphasePop) const
{
    vcs_VolPhase* Vphase = m_VolPhaseList[iphasePop].get();
    AssertThrowMsg(!Vphase->exists(), "VCS_SOLVE::vcs_popPhasePossible",
                   "called for a phase that exists!");

    // Loop through all of the species in the phase. We say the phase can be
    // popped, if there is one species in the phase that can be popped. This
    // does not mean that the phase will be popped or that it leads to a lower
    // Gibbs free energy.
    for (size_t k = 0; k < Vphase->nSpecies(); k++) {
        size_t kspec = Vphase->spGlobalIndexVCS(k);
        AssertThrowMsg(m_molNumSpecies_old[kspec] <= 0.0,
                       "VCS_SOLVE::vcs_popPhasePossible",
                       "we shouldn't be here {}: {} > 0.0", kspec,
                       m_molNumSpecies_old[kspec]);
        size_t irxn = kspec - m_numComponents;
        if (kspec >= m_numComponents) {
            bool iPopPossible = true;

            // Note one case is if the component is a member of the popping
            // phase. This component will be zeroed and the logic here will
            // negate the current species from causing a positive if this
            // component is consumed.
            for (size_t j = 0; j < m_numComponents; ++j) {
                if (m_elType[j] == VCS_ELEM_TYPE_ABSPOS) {
                    double stoicC = m_stoichCoeffRxnMatrix(j,irxn);
                    if (stoicC != 0.0) {
                        double negChangeComp = - stoicC;
                        if (negChangeComp > 0.0) {
                            // If there is no component to give, then the
                            // species can't be created
                            if (m_molNumSpecies_old[j] <= VCS_DELETE_ELEMENTABS_CUTOFF*0.5) {
                                iPopPossible = false;
                            }
                        }
                    }
                }
            }
            // We are here when the species can be popped because all its needed
            // components have positive mole numbers
            if (iPopPossible) {
                return true;
            }
        } else {
            // We are here when the species, k, in the phase is a component. Its
            // mole number is zero. We loop through the regular reaction looking
            // for a reaction that can pop the component.
            for (size_t jrxn = 0; jrxn < m_numRxnRdc; jrxn++) {
                bool foundJrxn = false;
                // First, if the component is a product of the reaction
                if (m_stoichCoeffRxnMatrix(kspec,jrxn) > 0.0) {
                    foundJrxn = true;
                    // We can do the reaction if all other reactant components
                    // have positive mole fractions
                    for (size_t kcomp = 0; kcomp < m_numComponents; kcomp++) {
                        if (m_stoichCoeffRxnMatrix(kcomp,jrxn) < 0.0 && m_molNumSpecies_old[kcomp] <= VCS_DELETE_ELEMENTABS_CUTOFF*0.5) {
                            foundJrxn = false;
                        }
                    }
                    if (foundJrxn) {
                        return true;
                    }
                } else if (m_stoichCoeffRxnMatrix(kspec,jrxn) < 0.0) {
                    // Second we are here if the component is a reactant in the
                    // reaction, and the reaction goes backwards.
                    foundJrxn = true;
                    size_t jspec = jrxn + m_numComponents;
                    if (m_molNumSpecies_old[jspec] <= VCS_DELETE_ELEMENTABS_CUTOFF*0.5) {
                        foundJrxn = false;
                        continue;
                    }
                    // We can do the backwards reaction if all of the product
                    // components species are positive
                    for (size_t kcomp = 0; kcomp < m_numComponents; kcomp++) {
                        if (m_stoichCoeffRxnMatrix(kcomp,jrxn) > 0.0 && m_molNumSpecies_old[kcomp] <= VCS_DELETE_ELEMENTABS_CUTOFF*0.5) {
                            foundJrxn = false;
                        }
                    }
                    if (foundJrxn) {
                        return true;
                    }
                }
            }
        }
    }
    return false;
}

size_t VCS_SOLVE::vcs_popPhaseID(vector<size_t> & phasePopPhaseIDs)
{
    size_t iphasePop = npos;
    double FephaseMax = -1.0E30;
    double Fephase = -1.0E30;

    string note;
    if (m_debug_print_lvl >= 2) {
        plogf("   --- vcs_popPhaseID() called\n");
        plogf("   ---   Phase                 Status       F_e        MoleNum\n");
        plogf("   --------------------------------------------------------------------------\n");
    }
    for (size_t iph = 0; iph < m_numPhases; iph++) {
        vcs_VolPhase* Vphase = m_VolPhaseList[iph].get();
        int existence = Vphase->exists();
        note = "";
        if (existence > 0) {
            if (m_debug_print_lvl >= 2) {
                plogf("  ---    %18s %5d           NA       %11.3e\n",
                      Vphase->PhaseName, existence, m_tPhaseMoles_old[iph]);
            }
        } else {
            if (Vphase->m_singleSpecies) {
                // Single Phase Stability Resolution
                size_t kspec = Vphase->spGlobalIndexVCS(0);
                size_t irxn = kspec - m_numComponents;
                if (irxn > m_deltaGRxn_old.size()) {
                    throw CanteraError("VCS_SOLVE::vcs_popPhaseID",
                        "Index out of bounds due to logic error.");
                }
                double deltaGRxn = m_deltaGRxn_old[irxn];
                Fephase = exp(-deltaGRxn) - 1.0;
                if (Fephase > 0.0) {
                    note = "(ready to be birthed)";
                    if (Fephase > FephaseMax) {
                        iphasePop = iph;
                        FephaseMax = Fephase;
                        note = "(chosen to be birthed)";
                    }
                }
                if (Fephase < 0.0) {
                    note = "(not stable)";
                    AssertThrowMsg(m_tPhaseMoles_old[iph] <= 0.0,
                        "VCS_SOLVE::vcs_popPhaseID", "shouldn't be here");
                }

                if (m_debug_print_lvl >= 2) {
                    plogf("  ---    %18s %5d %10.3g %10.3g  %s\n",
                          Vphase->PhaseName, existence, Fephase,
                          m_tPhaseMoles_old[iph], note.c_str());
                }
            } else if (vcs_popPhasePossible(iph)) {
                // MultiSpecies Phase Stability Resolution
                Fephase = vcs_phaseStabilityTest(iph);
                if (Fephase > 0.0) {
                    if (Fephase > FephaseMax) {
                        iphasePop = iph;
                        FephaseMax = Fephase;
                    }
                } else {
                    FephaseMax = std::max(FephaseMax, Fephase);
                }
                if (m_debug_print_lvl >= 2) {
                    plogf("  ---    %18s %5d  %11.3g %11.3g\n",
                            Vphase->PhaseName, existence, Fephase,
                            m_tPhaseMoles_old[iph]);
                }
            } else if (m_debug_print_lvl >= 2) {
                plogf("  ---    %18s %5d   blocked  %11.3g\n",
                        Vphase->PhaseName,
                        existence, m_tPhaseMoles_old[iph]);
            }
        }
    }
    phasePopPhaseIDs.resize(0);
    if (iphasePop != npos) {
        phasePopPhaseIDs.push_back(iphasePop);
    }

    // Insert logic here to figure out if phase pops are linked together. Only
    // do one linked pop at a time.
    if (m_debug_print_lvl >= 2) {
        plogf("   ---------------------------------------------------------------------\n");
    }
    return iphasePop;
}

int VCS_SOLVE::vcs_popPhaseRxnStepSizes(const size_t iphasePop)
{
    vcs_VolPhase* Vphase = m_VolPhaseList[iphasePop].get();
    // Identify the first species in the phase
    size_t kspec = Vphase->spGlobalIndexVCS(0);
    // Identify the formation reaction for that species
    size_t irxn = kspec - m_numComponents;
    vector<size_t> creationGlobalRxnNumbers;

    // Calculate the initial moles of the phase being born.
    //   Here we set it to 10x of the value which would cause the phase to be
    //   zeroed out within the algorithm.  We may later adjust the value.
    double tPhaseMoles = 10. * m_totalMolNum * VCS_DELETE_PHASE_CUTOFF;

    AssertThrowMsg(!Vphase->exists(), "VCS_SOLVE::vcs_popPhaseRxnStepSizes",
                   "called for a phase that exists!");
    if (m_debug_print_lvl >= 2) {
        plogf("  ---  vcs_popPhaseRxnStepSizes() called to pop phase %s %d into existence\n",
              Vphase->PhaseName, iphasePop);
    }
    // Section for a single-species phase
    if (Vphase->m_singleSpecies) {
        double s = 0.0;
        for (size_t j = 0; j < m_numComponents; ++j) {
            if (!m_SSPhase[j] && m_molNumSpecies_old[j] > 0.0) {
                s += pow(m_stoichCoeffRxnMatrix(j,irxn), 2) / m_molNumSpecies_old[j];
            }
        }
        for (size_t j = 0; j < m_numPhases; j++) {
            Vphase = m_VolPhaseList[j].get();
            if (! Vphase->m_singleSpecies && m_tPhaseMoles_old[j] > 0.0) {
                s -= pow(m_deltaMolNumPhase(j,irxn), 2) / m_tPhaseMoles_old[j];
            }
        }
        if (s != 0.0) {
            s = vcs_Hessian_diag_adj(irxn, s);
            m_deltaMolNumSpecies[kspec] = -m_deltaGRxn_new[irxn] / s;
        } else {
            // Ok, s is equal to zero. We can not apply a sophisticated theory
            // to birth the phase. Just pick a small delta and go with it.
            m_deltaMolNumSpecies[kspec] = tPhaseMoles;
        }

        // section to do damping of the m_deltaMolNumSpecies[]
        for (size_t j = 0; j < m_numComponents; ++j) {
            double stoicC = m_stoichCoeffRxnMatrix(j,irxn);
            if (stoicC != 0.0 && m_elType[j] == VCS_ELEM_TYPE_ABSPOS) {
                double negChangeComp = - stoicC * m_deltaMolNumSpecies[kspec];
                if (negChangeComp > m_molNumSpecies_old[j]) {
                    if (m_molNumSpecies_old[j] > 0.0) {
                        m_deltaMolNumSpecies[kspec] = - 0.5 * m_molNumSpecies_old[j] / stoicC;
                    } else {
                        m_deltaMolNumSpecies[kspec] = 0.0;
                    }
                }
            }
        }
        // Implement a damping term that limits m_deltaMolNumSpecies to the size
        // of the mole number
        if (-m_deltaMolNumSpecies[kspec] > m_molNumSpecies_old[kspec]) {
            m_deltaMolNumSpecies[kspec] = -m_molNumSpecies_old[kspec];
        }
    } else {
        vector<double> fracDelta(Vphase->nSpecies());
        vector<double> X_est(Vphase->nSpecies());
        fracDelta = Vphase->creationMoleNumbers(creationGlobalRxnNumbers);

        double sumFrac = 0.0;
        for (size_t k = 0; k < Vphase->nSpecies(); k++) {
            sumFrac += fracDelta[k];
        }
        for (size_t k = 0; k < Vphase->nSpecies(); k++) {
            X_est[k] = fracDelta[k] / sumFrac;
        }

        double deltaMolNumPhase = tPhaseMoles;
        double damp = 1.0;
        m_deltaGRxn_tmp = m_molNumSpecies_old;
        double* molNumSpecies_tmp = m_deltaGRxn_tmp.data();

        for (size_t k = 0; k < Vphase->nSpecies(); k++) {
            kspec = Vphase->spGlobalIndexVCS(k);
            double delmol = deltaMolNumPhase * X_est[k];
            if (kspec >= m_numComponents) {
                irxn = kspec - m_numComponents;
                if (irxn > m_stoichCoeffRxnMatrix.nColumns()) {
                    throw CanteraError("VCS_SOLVE::vcs_popPhaseRxnStepSizes",
                        "Index out of bounds due to logic error.");
                }
                for (size_t j = 0; j < m_numComponents; ++j) {
                    double stoicC = m_stoichCoeffRxnMatrix(j,irxn);
                    if (stoicC != 0.0 && m_elType[j] == VCS_ELEM_TYPE_ABSPOS) {
                        molNumSpecies_tmp[j] += stoicC * delmol;
                    }
                }
            }
        }

        double ratioComp = 0.0;
        for (size_t j = 0; j < m_numComponents; ++j) {
            double deltaJ = m_molNumSpecies_old[j] - molNumSpecies_tmp[j];
            if (molNumSpecies_tmp[j] < 0.0) {
                ratioComp = 1.0;
                if (deltaJ > 0.0) {
                    double delta0 = m_molNumSpecies_old[j];
                    damp = std::min(damp, delta0 / deltaJ * 0.9);
                }
            } else {
                if (m_elType[j] == VCS_ELEM_TYPE_ABSPOS) {
                    size_t jph = m_phaseID[j];
                    if ((jph != iphasePop) && (!m_SSPhase[j])) {
                        double fdeltaJ = fabs(deltaJ);
                        if (m_molNumSpecies_old[j] > 0.0) {
                            ratioComp = std::max(ratioComp, fdeltaJ/ m_molNumSpecies_old[j]);
                        }
                    }
                }
            }
        }

        // We may have greatly underestimated the deltaMoles for the phase pop
        // Here we create a damp > 1 to account for this possibility. We adjust
        // upwards to make sure that a component in an existing multispecies
        // phase is modified by a factor of 1/1000.
        if (ratioComp > 1.0E-30 && ratioComp < 0.001) {
            damp = 0.001 / ratioComp;
        }
        if (damp <= 1.0E-6) {
            return 3;
        }

        for (size_t k = 0; k < Vphase->nSpecies(); k++) {
            kspec = Vphase->spGlobalIndexVCS(k);
            if (kspec < m_numComponents) {
                m_speciesStatus[kspec] = VCS_SPECIES_COMPONENT;
            } else {
                m_deltaMolNumSpecies[kspec] = deltaMolNumPhase * X_est[k] * damp;
                m_speciesStatus[kspec] = (X_est[k] > 1.0E-3) ? VCS_SPECIES_MAJOR
                                                             : VCS_SPECIES_MINOR;
            }
        }
    }
    return 0;
}

size_t VCS_SOLVE::vcs_RxnStepSizes(int& forceComponentCalc, size_t& kSpecial)
{
    size_t iphDel = npos;
    size_t k = 0;
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
        plogf("   --- Species        KMoles     Rxn_Adjustment    DeltaG\n");
    }

    // We update the matrix dlnActCoeffdmolNumber[][] at the top of the loop,
    // when necessary
    if (m_useActCoeffJac) {
        vcs_CalcLnActCoeffJac(&m_molNumSpecies_old[0]);
    }

    // LOOP OVER THE FORMATION REACTIONS
    for (size_t irxn = 0; irxn < m_numRxnRdc; ++irxn) {

        size_t kspec = m_indexRxnToSpecies[irxn];
        if (m_speciesStatus[kspec] == VCS_SPECIES_ZEROEDPHASE) {
            m_deltaMolNumSpecies[kspec] = 0.0;
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
                        } else {
                            m_deltaMolNumSpecies[kspec] = m_totalMolNum * VCS_SMALL_MULTIPHASE_SPECIES * 10.0;
                        }
                    } else {
                        m_deltaMolNumSpecies[kspec] = 0.0;
                        if (Vphase->exists() > 0 && trphmoles > 0.0) {
                            m_deltaMolNumSpecies[kspec] = m_totalMolNum * VCS_SMALL_MULTIPHASE_SPECIES * 10.;
                        }
                    }
                } else {
                    m_deltaMolNumSpecies[kspec] = 0.0;
                }
            } else {
                // REGULAR PROCESSING
                //
                // First take care of cases where we want to bail out. Don't
                // bother if superconvergence has already been achieved in this
                // mode.
                if (fabs(m_deltaGRxn_new[irxn]) <= m_tolmaj2) {
                    if (m_debug_print_lvl >= 2) {
                        plogf("   --- %-12.12s", m_speciesName[kspec]);
                        plogf("  %12.4E %12.4E %12.4E\n",
                              m_molNumSpecies_old[kspec], m_deltaMolNumSpecies[kspec],
                              m_deltaGRxn_new[irxn]);
                    }
                    continue;
                }

                // Don't calculate for minor or nonexistent species if their
                // values are to be decreasing anyway.
                if ((m_speciesStatus[kspec] != VCS_SPECIES_MAJOR) && (m_deltaGRxn_new[irxn] >= 0.0)) {
                    if (m_debug_print_lvl >= 2) {
                        plogf("   --- %-12.12s", m_speciesName[kspec]);
                        plogf("  %12.4E %12.4E %12.4E\n",
                              m_molNumSpecies_old[kspec], m_deltaMolNumSpecies[kspec],
                              m_deltaGRxn_new[irxn]);
                    }
                    continue;
                }

                // Start of the regular processing
                double s = (m_SSPhase[kspec]) ? 0.0 : 1.0 / m_molNumSpecies_old[kspec];
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
                    }

                    m_deltaMolNumSpecies[kspec] = -m_deltaGRxn_new[irxn] / s;
                    // New section to do damping of the m_deltaMolNumSpecies[]
                    for (size_t j = 0; j < m_numComponents; ++j) {
                        double stoicC = m_stoichCoeffRxnMatrix(j,irxn);
                        if (stoicC != 0.0) {
                            double negChangeComp = -stoicC * m_deltaMolNumSpecies[kspec];
                            if (negChangeComp > m_molNumSpecies_old[j]) {
                                if (m_molNumSpecies_old[j] > 0.0) {
                                    m_deltaMolNumSpecies[kspec] = -m_molNumSpecies_old[j] / stoicC;
                                } else {
                                    m_deltaMolNumSpecies[kspec] = 0.0;
                                }
                            }
                        }
                    }
                    // Implement a damping term that limits m_deltaMolNumSpecies
                    // to the size of the mole number
                    if (-m_deltaMolNumSpecies[kspec] > m_molNumSpecies_old[kspec]) {
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
                            m_deltaMolNumSpecies[kspec] = -m_molNumSpecies_old[kspec];
                            if (m_debug_print_lvl >= 2) {
                                plogf("   --- %-12.12s", m_speciesName[kspec]);
                                plogf("  %12.4E %12.4E %12.4E\n",
                                      m_molNumSpecies_old[kspec], m_deltaMolNumSpecies[kspec],
                                      m_deltaGRxn_new[irxn]);
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
                        } else {
                        }
                        if (m_debug_print_lvl >= 2) {
                            plogf("   --- %-12.12s", m_speciesName[kspec]);
                            plogf("  %12.4E %12.4E %12.4E\n",
                                  m_molNumSpecies_old[kspec], m_deltaMolNumSpecies[kspec],
                                  m_deltaGRxn_new[irxn]);
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
                plogf("  %12.4E %12.4E %12.4E\n",
                      m_molNumSpecies_old[kspec], m_deltaMolNumSpecies[kspec],
                      m_deltaGRxn_new[irxn]);
            }
        } // End of loop over m_speciesUnknownType
    } // End of loop over non-component stoichiometric formation reactions
    if (m_debug_print_lvl >= 2) {
        plogf("   ");
        writeline('-', 82);
    }
    return iphDel;
}

double VCS_SOLVE::vcs_phaseStabilityTest(const size_t iph)
{
    // We will use the _new state calc here
    vcs_VolPhase* Vphase = m_VolPhaseList[iph].get();
    const size_t nsp = Vphase->nSpecies();
    int minNumberIterations = 3;
    if (nsp <= 1) {
        minNumberIterations = 1;
    }

    // We will do a full Newton calculation later, but for now, ...
    bool doSuccessiveSubstitution = true;
    double funcPhaseStability;
    vector<double> X_est(nsp, 0.0);
    vector<double> delFrac(nsp, 0.0);
    vector<double> E_phi(nsp, 0.0);
    vector<double> fracDelta_old(nsp, 0.0);
    vector<double> fracDelta_raw(nsp, 0.0);
    vector<size_t> creationGlobalRxnNumbers(nsp, npos);
    m_deltaGRxn_Deficient = m_deltaGRxn_old;
    vector<double> feSpecies_Deficient = m_feSpecies_old;

    // get the activity coefficients
    Vphase->sendToVCS_ActCoeff(VCS_STATECALC_OLD, &m_actCoeffSpecies_new[0]);

    // Get the stored estimate for the composition of the phase if
    // it gets created
    vector<double> fracDelta_new = Vphase->creationMoleNumbers(creationGlobalRxnNumbers);

    vector<size_t> componentList;
    for (size_t k = 0; k < nsp; k++) {
        size_t kspec = Vphase->spGlobalIndexVCS(k);
        if (kspec < m_numComponents) {
            componentList.push_back(k);
        }
    }

    double normUpdate = 0.1 * vcs_l2norm(fracDelta_new);
    double damp = 1.0E-2;

    if (doSuccessiveSubstitution) {
        int KP = 0;
        if (m_debug_print_lvl >= 2) {
            plogf("   --- vcs_phaseStabilityTest() called\n");
            plogf("   ---  Its   X_old[%2d]  FracDel_old[%2d]  deltaF[%2d] FracDel_new[%2d]"
                  "  normUpdate     damp     FuncPhaseStability\n", KP, KP, KP, KP);
            plogf("   --------------------------------------------------------------"
                  "--------------------------------------------------------\n");
        } else if (m_debug_print_lvl == 1) {
            plogf("   --- vcs_phaseStabilityTest() called for phase %d\n", iph);
        }

        for (size_t k = 0; k < nsp; k++) {
            if (fracDelta_new[k] < 1.0E-13) {
               fracDelta_new[k] =1.0E-13;
            }
        }
        bool converged = false;
        double dirProd = 0.0;
        for (int its = 0; its < 200 && (!converged); its++) {
            double dampOld = damp;
            double normUpdateOld = normUpdate;
            fracDelta_old = fracDelta_new;
            double dirProdOld = dirProd;

            // Given a set of fracDelta's, we calculate the fracDelta's
            // for the component species, if any
            for (size_t i = 0; i < componentList.size(); i++) {
                size_t kc = componentList[i];
                size_t kc_spec = Vphase->spGlobalIndexVCS(kc);
                fracDelta_old[kc] = 0.0;
                for (size_t k = 0; k < nsp; k++) {
                    size_t kspec = Vphase->spGlobalIndexVCS(k);
                    if (kspec >= m_numComponents) {
                        size_t irxn = kspec - m_numComponents;
                        fracDelta_old[kc] += m_stoichCoeffRxnMatrix(kc_spec,irxn) * fracDelta_old[k];
                    }
                }
            }

            // Now, calculate the predicted mole fractions, X_est[k]
            double sumFrac = 0.0;
            for (size_t k = 0; k < nsp; k++) {
                sumFrac += fracDelta_old[k];
            }
            // Necessary because this can be identically zero. -> we need to fix
            // this algorithm!
            if (sumFrac <= 0.0) {
                sumFrac = 1.0;
            }
            double sum_Xcomp = 0.0;
            for (size_t k = 0; k < nsp; k++) {
                X_est[k] = fracDelta_old[k] / sumFrac;
                if (Vphase->spGlobalIndexVCS(k) < m_numComponents) {
                    sum_Xcomp += X_est[k];
                }
            }

            // Feed the newly formed estimate of the mole fractions back into the
            // ThermoPhase object
            Vphase->setMoleFractionsState(0.0, &X_est[0], VCS_STATECALC_PHASESTABILITY);

            // get the activity coefficients
            Vphase->sendToVCS_ActCoeff(VCS_STATECALC_OLD, &m_actCoeffSpecies_new[0]);

            // First calculate altered chemical potentials for component species
            // belonging to this phase.
            for (size_t i = 0; i < componentList.size(); i++) {
                size_t kc = componentList[i];
                size_t kc_spec = Vphase->spGlobalIndexVCS(kc);
                double moleFrac = std::max(X_est[kc], VCS_DELETE_MINORSPECIES_CUTOFF);
                feSpecies_Deficient[kc_spec] = m_feSpecies_old[kc_spec]
                    + log(m_actCoeffSpecies_new[kc_spec] * moleFrac);
            }

            for (size_t i = 0; i < componentList.size(); i++) {
                size_t kc_spec = Vphase->spGlobalIndexVCS(componentList[i]);
                for (size_t k = 0; k < Vphase->nSpecies(); k++) {
                    size_t kspec = Vphase->spGlobalIndexVCS(k);
                    if (kspec >= m_numComponents) {
                        size_t irxn = kspec - m_numComponents;
                        if (i == 0) {
                            m_deltaGRxn_Deficient[irxn] = m_deltaGRxn_old[irxn];
                        }
                        if (m_stoichCoeffRxnMatrix(kc_spec,irxn) != 0.0) {
                            m_deltaGRxn_Deficient[irxn] +=
                                m_stoichCoeffRxnMatrix(kc_spec,irxn) * (feSpecies_Deficient[kc_spec]- m_feSpecies_old[kc_spec]);
                        }
                    }
                }
            }

            // Calculate the E_phi's
            double sum = 0.0;
            funcPhaseStability = sum_Xcomp - 1.0;
            for (size_t k = 0; k < nsp; k++) {
                size_t kspec = Vphase->spGlobalIndexVCS(k);
                if (kspec >= m_numComponents) {
                    size_t irxn = kspec - m_numComponents;
                    double deltaGRxn = clip(m_deltaGRxn_Deficient[irxn], -50.0, 50.0);
                    E_phi[k] = std::exp(-deltaGRxn) / m_actCoeffSpecies_new[kspec];
                    sum += E_phi[k];
                    funcPhaseStability += E_phi[k];
                } else {
                    E_phi[k] = 0.0;
                }
            }

            // Calculate the raw estimate of the new fracs
            for (size_t k = 0; k < nsp; k++) {
                size_t kspec = Vphase->spGlobalIndexVCS(k);
                double b = E_phi[k] / sum * (1.0 - sum_Xcomp);
                if (kspec >= m_numComponents) {
                    fracDelta_raw[k] = b;
                }
            }

            // Given a set of fracDelta's, we calculate the fracDelta's
            // for the component species, if any
            for (size_t i = 0; i < componentList.size(); i++) {
                size_t kc = componentList[i];
                size_t kc_spec = Vphase->spGlobalIndexVCS(kc);
                fracDelta_raw[kc] = 0.0;
                for (size_t k = 0; k < nsp; k++) {
                    size_t kspec = Vphase->spGlobalIndexVCS(k);
                    if (kspec >= m_numComponents) {
                        size_t irxn = kspec - m_numComponents;
                        fracDelta_raw[kc] += m_stoichCoeffRxnMatrix(kc_spec,irxn) * fracDelta_raw[k];
                    }
                }
            }

            // Now possibly dampen the estimate.
            for (size_t k = 0; k < nsp; k++) {
                delFrac[k] = fracDelta_raw[k] - fracDelta_old[k];
            }
            normUpdate = vcs_l2norm(delFrac);

            dirProd = 0.0;
            for (size_t k = 0; k < nsp; k++) {
                dirProd += fracDelta_old[k] * delFrac[k];
            }
            bool crossedSign = false;
            if (dirProd * dirProdOld < 0.0) {
                crossedSign = true;
            }

            damp = 0.5;
            if (dampOld < 0.25) {
                damp = 2.0 * dampOld;
            }
            if (crossedSign) {
                if (normUpdate *1.5 > normUpdateOld) {
                    damp = 0.5 * dampOld;
                } else if (normUpdate *2.0 > normUpdateOld) {
                    damp = 0.8 * dampOld;
                }
            } else {
                if (normUpdate > normUpdateOld * 2.0) {
                    damp = 0.6 * dampOld;
                } else if (normUpdate > normUpdateOld * 1.2) {
                    damp = 0.9 * dampOld;
                }
            }

            for (size_t k = 0; k < nsp; k++) {
                if (fabs(damp * delFrac[k]) > 0.3*fabs(fracDelta_old[k])) {
                    damp = std::max(0.3*fabs(fracDelta_old[k]) / fabs(delFrac[k]),
                                    1.0E-8/fabs(delFrac[k]));
                }
                if (delFrac[k] < 0.0 && 2.0 * damp * (-delFrac[k]) > fracDelta_old[k]) {
                    damp = fracDelta_old[k] / (2.0 * -delFrac[k]);
                }
                if (delFrac[k] > 0.0 && 2.0 * damp * delFrac[k] > fracDelta_old[k]) {
                    damp = fracDelta_old[k] / (2.0 * delFrac[k]);
                }
            }
            damp = std::max(damp, 0.000001);
            for (size_t k = 0; k < nsp; k++) {
                fracDelta_new[k] = fracDelta_old[k] + damp * delFrac[k];
            }

            if (m_debug_print_lvl >= 2) {
                plogf("  --- %3d %12g %12g %12g %12g %12g %12g %12g\n", its, X_est[KP], fracDelta_old[KP],
                      delFrac[KP], fracDelta_new[KP], normUpdate, damp, funcPhaseStability);
            }

            if (normUpdate < 1.0E-5 * damp) {
                converged = true;
                if (its < minNumberIterations) {
                    converged = false;
                }
            }
        }

        if (converged) {
            // Save the final optimized stated back into the VolPhase object for later use
            Vphase->setMoleFractionsState(0.0, &X_est[0], VCS_STATECALC_PHASESTABILITY);

            // Save fracDelta for later use to initialize the problem better
            // @todo  creationGlobalRxnNumbers needs to be calculated here and stored.
            Vphase->setCreationMoleNumbers(&fracDelta_new[0], creationGlobalRxnNumbers);
        }
    } else {
        throw CanteraError("VCS_SOLVE::vcs_phaseStabilityTest", "not done yet");
    }
    if (m_debug_print_lvl >= 2) {
        plogf("  ------------------------------------------------------------"
              "-------------------------------------------------------------\n");
    } else if (m_debug_print_lvl == 1) {
        if (funcPhaseStability > 0.0) {
            plogf("  --- phase %d with func = %g is to be born\n", iph, funcPhaseStability);
        } else {
            plogf("  --- phase %d with func = %g stays dead\n", iph, funcPhaseStability);
        }
    }
    return funcPhaseStability;
}

int VCS_SOLVE::vcs_evalSS_TP(int ipr, int ip1, double Temp, double pres)
{
    for (size_t iph = 0; iph < m_numPhases; iph++) {
        vcs_VolPhase* vph = m_VolPhaseList[iph].get();
        vph->setState_TP(m_temperature, m_pressurePA);
        vph->sendToVCS_GStar(&m_SSfeSpecies[0]);
    }
    for (size_t k = 0; k < m_nsp; k++) {
        m_SSfeSpecies[k] /= GasConstant * m_temperature;
    }

    return VCS_SUCCESS;
}

void VCS_SOLVE::vcs_prep(int printLvl)
{
    m_debug_print_lvl = printLvl;

    // Calculate the Single Species status of phases
    // Also calculate the number of species per phase
    vcs_SSPhase();

    // Set an initial estimate for the number of noncomponent species equal to
    // nspecies - nelements. This may be changed below
    m_numRxnTot = (m_nelem > m_nsp) ? 0 : m_nsp - m_nelem;
    m_numRxnRdc = m_numRxnTot;
    m_numSpeciesRdc = m_nsp;
    for (size_t i = 0; i < m_numRxnRdc; ++i) {
        m_indexRxnToSpecies[i] = m_nelem + i;
    }

    for (size_t kspec = 0; kspec < m_nsp; ++kspec) {
        double sz = 0.0;
        for (size_t e = 0; e < m_nelem; e++) {
            sz += fabs(m_formulaMatrix(kspec, e));
        }
        m_spSize[kspec] = (sz > 0.0) ? sz : 1.0;
    }

    // DETERMINE THE NUMBER OF COMPONENTS
    //
    // Obtain a valid estimate of the mole fraction. This will be used as an
    // initial ordering vector for prioritizing which species are defined as
    // components.
    //
    // If a mole number estimate was supplied from the input file, use that mole
    // number estimate.
    //
    // If a solution estimate wasn't supplied from the input file, supply an
    // initial estimate for the mole fractions based on the relative reverse
    // ordering of the chemical potentials.
    //
    // For voltage unknowns, set these to zero for the moment.
    double test = -1.0e-10;
    bool modifiedSoln = false;
    if (m_doEstimateEquil < 0) {
        double sum = 0.0;
        for (size_t kspec = 0; kspec < m_nsp; ++kspec) {
            if (m_speciesUnknownType[kspec] == VCS_SPECIES_TYPE_MOLNUM) {
                sum += fabs(m_molNumSpecies_old[kspec]);
            }
        }
        if (fabs(sum) < 1.0E-6) {
            modifiedSoln = true;
            double pres = (m_pressurePA <= 0.0) ? 1.01325E5 : m_pressurePA;
            vcs_evalSS_TP(0, 0, m_temperature, pres);
            for (size_t kspec = 0; kspec < m_nsp; ++kspec) {
                if (m_speciesUnknownType[kspec] == VCS_SPECIES_TYPE_MOLNUM) {
                    m_molNumSpecies_old[kspec] = - m_SSfeSpecies[kspec];
                } else {
                    m_molNumSpecies_old[kspec] = 0.0;
                }
            }
        }
        test = -1.0e20;
    }

    // NC = number of components is in the vcs.h common block. This call to
    // BASOPT doesn't calculate the stoichiometric reaction matrix.
    vcs_basopt(true, test);

    if (m_nsp >= m_numComponents) {
        m_numRxnTot = m_numRxnRdc = m_nsp - m_numComponents;
        for (size_t i = 0; i < m_numRxnRdc; ++i) {
            m_indexRxnToSpecies[i] = m_numComponents + i;
        }
    } else {
        m_numRxnTot = m_numRxnRdc = 0;
    }

    // The elements might need to be rearranged.
    vcs_elem_rearrange();

    // If we mucked up the solution unknowns because they were all
    // zero to start with, set them back to zero here
    if (modifiedSoln) {
        for (size_t kspec = 0; kspec < m_nsp; ++kspec) {
            m_molNumSpecies_old[kspec] = 0.0;
        }
    }

    // Initialize various arrays in the data to zero
    m_feSpecies_old.assign(m_feSpecies_old.size(), 0.0);
    m_feSpecies_new.assign(m_feSpecies_new.size(), 0.0);
    m_molNumSpecies_new.assign(m_molNumSpecies_new.size(), 0.0);
    m_deltaMolNumPhase.zero();
    m_phaseParticipation.zero();
    m_deltaPhaseMoles.assign(m_deltaPhaseMoles.size(), 0.0);
    m_tPhaseMoles_new.assign(m_tPhaseMoles_new.size(), 0.0);

    // Calculate the total number of moles in all phases.
    vcs_tmoles();

    // Check to see if the current problem is well posed.
    double sum = 0.0;
    for (size_t e = 0; e < m_mix->nElements(); e++) {
        sum += m_mix->elementMoles(e);
    }
    if (sum < 1.0E-20) {
        throw CanteraError("VCS_SOLVE::vcs_prep", "The problem is not well posed"
                           " because the total number of element moles is zero.");
    }
}

void VCS_SOLVE::vcs_elem_rearrange()
{
    vector<double> awSpace(m_nsp + (m_nelem + 2)*(m_nelem), 0.0);
    vector<double> aw(m_nelem), sa(m_nelem), ss(m_nelem);
    vector<double> sm(m_nelem*m_nelem);

    size_t ncomponents = m_numComponents;
    if (m_debug_print_lvl >= 2) {
        plogf("   ");
        writeline('-', 77);
        plogf("   --- Subroutine elem_rearrange() called to ");
        plogf("check stoich. coefficient matrix\n");
        plogf("   ---    and to rearrange the element ordering once\n");
    }

    // Use a temporary work array for the element numbers
    // Also make sure the value of test is unique.
    bool lindep = true;
    double test = -1.0E10;
    while (lindep) {
        lindep = false;
        for (size_t i = 0; i < m_nelem; ++i) {
            test -= 1.0;
            aw[i] = m_elemAbundancesGoal[i];
            if (test == aw[i]) {
                lindep = true;
            }
        }
    }

    // Top of a loop of some sort based on the index JR. JR is the current
    // number independent elements found.
    size_t jr = 0;
    while (jr < ncomponents) {
        size_t k;

        // Top of another loop point based on finding a linearly independent
        // species
        while (true) {
            // Search the remaining part of the mole fraction vector, AW, for
            // the largest remaining species. Return its identity in K.
            k = m_nelem;
            for (size_t ielem = jr; ielem < m_nelem; ielem++) {
                if (m_elementActive[ielem] && aw[ielem] != test) {
                    k = ielem;
                    break;
                }
            }
            if (k == m_nelem) {
                throw CanteraError("VCS_SOLVE::vcs_elem_rearrange",
                        "Shouldn't be here. Algorithm misfired.");
            }

            // Assign a large negative number to the element that we have just
            // found, in order to take it out of further consideration.
            aw[k] = test;

            // CHECK LINEAR INDEPENDENCE OF CURRENT FORMULA MATRIX LINE WITH
            // PREVIOUS LINES OF THE FORMULA MATRIX
            //
            // Modified Gram-Schmidt Method, p. 202 Dalquist QR factorization of
            // a matrix without row pivoting.
            size_t jl = jr;

            // Fill in the row for the current element, k, under consideration
            // The row will contain the Formula matrix value for that element
            // from the current component.
            for (size_t j = 0; j < ncomponents; ++j) {
                sm[j + jr*ncomponents] = m_formulaMatrix(j,k);
            }
            if (jl > 0) {
                // Compute the coefficients of JA column of the the upper
                // triangular R matrix, SS(J) = R_J_JR (this is slightly
                // different than Dalquist) R_JA_JA = 1
                for (size_t j = 0; j < jl; ++j) {
                    ss[j] = 0.0;
                    for (size_t i = 0; i < ncomponents; ++i) {
                        ss[j] += sm[i + jr*ncomponents] * sm[i + j*ncomponents];
                    }
                    ss[j] /= sa[j];
                }

                // Now make the new column, (*,JR), orthogonal to the previous
                // columns
                for (size_t j = 0; j < jl; ++j) {
                    for (size_t i = 0; i < ncomponents; ++i) {
                        sm[i + jr*ncomponents] -= ss[j] * sm[i + j*ncomponents];
                    }
                }
            }

            // Find the new length of the new column in Q. It will be used in
            // the denominator in future row calcs.
            sa[jr] = 0.0;
            for (size_t ml = 0; ml < ncomponents; ++ml) {
                sa[jr] += pow(sm[ml + jr*ncomponents], 2);
            }
            // IF NORM OF NEW ROW  .LT. 1E-6 REJECT
            if (sa[jr] > 1.0e-6) {
                break;
            }
        }
        // REARRANGE THE DATA
        if (jr != k) {
            if (m_debug_print_lvl >= 2) {
                plogf("   ---   %-2.2s(%9.2g) replaces %-2.2s(%9.2g) as element %3d\n",
                    m_elementName[k], m_elemAbundancesGoal[k],
                    m_elementName[jr], m_elemAbundancesGoal[jr], jr);
            }
            vcs_switch_elem_pos(jr, k);
            std::swap(aw[jr], aw[k]);
        }

        // If we haven't found enough components, go back and find some more.
        jr++;
    }
}

void VCS_SOLVE::vcs_switch_elem_pos(size_t ipos, size_t jpos)
{
    if (ipos == jpos) {
        return;
    }
    AssertThrowMsg(ipos < m_nelem && jpos < m_nelem,
                   "vcs_switch_elem_pos",
                   "inappropriate args: {} {}", ipos, jpos);

    // Change the element Global Index list in each vcs_VolPhase object
    // to reflect the switch in the element positions.
    for (size_t iph = 0; iph < m_numPhases; iph++) {
        vcs_VolPhase* volPhase = m_VolPhaseList[iph].get();
        for (size_t e = 0; e < volPhase->nElemConstraints(); e++) {
            if (volPhase->elemGlobalIndex(e) == ipos) {
                volPhase->setElemGlobalIndex(e, jpos);
            }
            if (volPhase->elemGlobalIndex(e) == jpos) {
                volPhase->setElemGlobalIndex(e, ipos);
            }
        }
    }
    std::swap(m_elemAbundancesGoal[ipos], m_elemAbundancesGoal[jpos]);
    std::swap(m_elemAbundances[ipos], m_elemAbundances[jpos]);
    std::swap(m_elementMapIndex[ipos], m_elementMapIndex[jpos]);
    std::swap(m_elType[ipos], m_elType[jpos]);
    std::swap(m_elementActive[ipos], m_elementActive[jpos]);
    for (size_t j = 0; j < m_nsp; ++j) {
        std::swap(m_formulaMatrix(j,ipos), m_formulaMatrix(j,jpos));
    }
    std::swap(m_elementName[ipos], m_elementName[jpos]);
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
    for (auto& Vphase : m_VolPhaseList) {
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

void VCS_SOLVE::vcs_report(int iconv)
{
    // SORT DEPENDENT SPECIES IN DECREASING ORDER OF MOLES
    vector<pair<double, size_t>> x_order;
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
    double totalVolume = 0.0;
    for (size_t iphase = 0; iphase < m_numPhases; iphase++) {
        vcs_VolPhase* Vphase = m_VolPhaseList[iphase].get();
        Vphase->setState_TP(m_temperature, m_pressurePA);
        Vphase->setMolesFromVCS(VCS_STATECALC_OLD, m_molNumSpecies_old.data());
        double Volp = Vphase->sendToVCS_VolPM(m_PMVolumeSpecies.data());
        totalVolume += Volp;
    }

    plogf("\t\tTemperature  = %15.2g Kelvin\n", m_temperature);
    plogf("\t\tPressure     = %15.5g Pa \n", m_pressurePA);
    plogf("\t\ttotal Volume = %15.5g m**3\n", totalVolume);

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
    vector<double> gaTPhase(m_nelem, 0.0);
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
        // Compute the elemental abundances for each element in this phase
        for (size_t j = 0; j < m_nelem; ++j) {
            double abundance_j = 0.0;
            for (size_t i = 0; i < m_nsp; ++i) {
                if (m_speciesUnknownType[i] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE
                    && m_phaseID[i] == iphase)
                {
                    abundance_j += m_formulaMatrix(i,j) * m_molNumSpecies_old[i];
                }
            }
            plogf(" %10.3g", abundance_j);
            gaTPhase[j] += abundance_j;
        }
        gibbsPhase = 0.0;
        for (size_t kspec = 0; kspec < m_numSpeciesRdc; ++kspec) {
            if (m_phaseID[kspec] == iphase && m_speciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
                gibbsPhase += m_molNumSpecies_old[kspec] * m_feSpecies_old[kspec];
            }
        }
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

    // Calculate the total dimensionless Gibbs Free Energy.
    double g = vcs_Total_Gibbs(&m_molNumSpecies_old[0], &m_feSpecies_old[0],
                               &m_tPhaseMoles_old[0]);
    plogf("\n\tTotal Dimensionless Gibbs Free Energy = G/RT = %15.7E\n", g);
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
    plogf("\nCounters:         Iterations\n");
    plogf("    vcs_basopt:   %5d\n", m_VCount->Basis_Opts);
    plogf("    vcs_TP:       %5d\n", m_VCount->Its);
    writeline('-', 80);
    writeline('-', 80);
}

void VCS_SOLVE::vcs_elab()
{
    for (size_t j = 0; j < m_nelem; ++j) {
        m_elemAbundances[j] = 0.0;
        for (size_t i = 0; i < m_nsp; ++i) {
            if (m_speciesUnknownType[i] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
                m_elemAbundances[j] += m_formulaMatrix(i,j) * m_molNumSpecies_old[i];
            }
        }
    }
}

bool VCS_SOLVE::vcs_elabcheck(int ibound)
{
    size_t top = (ibound) ? m_nelem : m_numComponents;
    for (size_t i = 0; i < top; ++i) {
        // Require 12 digits of accuracy on non-zero constraints.
        if (m_elementActive[i] && fabs(m_elemAbundances[i] - m_elemAbundancesGoal[i]) > fabs(m_elemAbundancesGoal[i]) * 1.0e-12) {
            // This logic is for charge neutrality condition
            if (m_elType[i] == VCS_ELEM_TYPE_CHARGENEUTRALITY &&
                    m_elemAbundancesGoal[i] != 0.0) {
                throw CanteraError("VCS_SOLVE::vcs_elabcheck",
                                   "Problem with charge neutrality condition");
            }
            if (m_elemAbundancesGoal[i] == 0.0 || (m_elType[i] == VCS_ELEM_TYPE_ELECTRONCHARGE)) {
                double scale = VCS_DELETE_MINORSPECIES_CUTOFF;

                // Find out if the constraint is a multisign constraint. If it
                // is, then we have to worry about roundoff error in the
                // addition of terms. We are limited to 13 digits of finite
                // arithmetic accuracy.
                bool multisign = false;
                for (size_t kspec = 0; kspec < m_nsp; kspec++) {
                    double eval = m_formulaMatrix(kspec,i);
                    if (eval < 0.0) {
                        multisign = true;
                    }
                    if (eval != 0.0) {
                        scale = std::max(scale, fabs(eval * m_molNumSpecies_old[kspec]));
                    }
                }
                if (multisign) {
                    if (fabs(m_elemAbundances[i] - m_elemAbundancesGoal[i]) > 1e-11 * scale) {
                        return false;
                    }
                } else {
                    if (fabs(m_elemAbundances[i] - m_elemAbundancesGoal[i]) > VCS_DELETE_MINORSPECIES_CUTOFF) {
                        return false;
                    }
                }
            } else {
                return false;
            }
        }
    }
    return true;
}

int VCS_SOLVE::vcs_elcorr(double aa[], double x[])
{
    int retn = 0;

    vector<double> ga_save(m_elemAbundances);
    if (m_debug_print_lvl >= 2) {
        plogf("   --- vcsc_elcorr: Element abundances correction routine");
        if (m_nelem != m_numComponents) {
            plogf(" (m_numComponents != m_nelem)");
        }
        plogf("\n");
    }

    for (size_t i = 0; i < m_nelem; ++i) {
        x[i] = m_elemAbundances[i] - m_elemAbundancesGoal[i];
    }
    double l2before = 0.0;
    for (size_t i = 0; i < m_nelem; ++i) {
        l2before += x[i] * x[i];
    }
    l2before = sqrt(l2before/m_nelem);

    // Special section to take out single species, single component,
    // moles. These are species which have non-zero entries in the
    // formula matrix, and no other species have zero values either.
    bool changed = false;
    for (size_t i = 0; i < m_nelem; ++i) {
        int numNonZero = 0;
        bool multisign = false;
        for (size_t kspec = 0; kspec < m_nsp; kspec++) {
            if (m_speciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
                double eval = m_formulaMatrix(kspec,i);
                if (eval < 0.0) {
                    multisign = true;
                }
                if (eval != 0.0) {
                    numNonZero++;
                }
            }
        }
        if (!multisign) {
            if (numNonZero < 2) {
                for (size_t kspec = 0; kspec < m_nsp; kspec++) {
                    if (m_speciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
                        double eval = m_formulaMatrix(kspec,i);
                        if (eval > 0.0) {
                            m_molNumSpecies_old[kspec] = m_elemAbundancesGoal[i] / eval;
                            changed = true;
                        }
                    }
                }
            } else {
                int numCompNonZero = 0;
                size_t compID = npos;
                for (size_t kspec = 0; kspec < m_numComponents; kspec++) {
                    if (m_speciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
                        double eval = m_formulaMatrix(kspec,i);
                        if (eval > 0.0) {
                            compID = kspec;
                            numCompNonZero++;
                        }
                    }
                }
                if (numCompNonZero == 1) {
                    double diff = m_elemAbundancesGoal[i];
                    for (size_t kspec = m_numComponents; kspec < m_nsp; kspec++) {
                        if (m_speciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
                            double eval = m_formulaMatrix(kspec,i);
                            diff -= eval * m_molNumSpecies_old[kspec];
                        }
                        m_molNumSpecies_old[compID] = std::max(0.0,diff/m_formulaMatrix(compID,i));
                        changed = true;
                    }
                }
            }
        }
    }
    if (changed) {
        vcs_elab();
    }

    // Section to check for maximum bounds errors on all species due to
    // elements. This may only be tried on element types which are
    // VCS_ELEM_TYPE_ABSPOS. This is because no other species may have a
    // negative number of these.
    //
    // Note, also we can do this over ne, the number of elements, not just the
    // number of components.
    changed = false;
    for (size_t i = 0; i < m_nelem; ++i) {
        int elType = m_elType[i];
        if (elType == VCS_ELEM_TYPE_ABSPOS) {
            for (size_t kspec = 0; kspec < m_nsp; kspec++) {
                if (m_speciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
                    double atomComp = m_formulaMatrix(kspec,i);
                    if (atomComp > 0.0) {
                        double maxPermissible = m_elemAbundancesGoal[i] / atomComp;
                        if (m_molNumSpecies_old[kspec] > maxPermissible) {
                            if (m_debug_print_lvl >= 3) {
                                plogf("  ---  vcs_elcorr: Reduced species %s from %g to %g "
                                      "due to %s max bounds constraint\n",
                                      m_speciesName[kspec], m_molNumSpecies_old[kspec],
                                      maxPermissible, m_elementName[i]);
                            }
                            m_molNumSpecies_old[kspec] = maxPermissible;
                            changed = true;
                            if (m_molNumSpecies_old[kspec] < VCS_DELETE_MINORSPECIES_CUTOFF) {
                                m_molNumSpecies_old[kspec] = 0.0;
                                if (m_SSPhase[kspec]) {
                                    m_speciesStatus[kspec] = VCS_SPECIES_ZEROEDSS;
                                } else {
                                    m_speciesStatus[kspec] = VCS_SPECIES_ACTIVEBUTZERO;
                                }
                                if (m_debug_print_lvl >= 2) {
                                    plogf("  ---  vcs_elcorr: Zeroed species %s and changed "
                                          "status to %d due to max bounds constraint\n",
                                          m_speciesName[kspec], m_speciesStatus[kspec]);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Recalculate the element abundances if something has changed.
    if (changed) {
        vcs_elab();
    }

    // Ok, do the general case. Linear algebra problem is of length nc, not ne,
    // as there may be degenerate rows when nc .ne. ne.
    DenseMatrix A(m_numComponents, m_numComponents);
    for (size_t i = 0; i < m_numComponents; ++i) {
        x[i] = m_elemAbundances[i] - m_elemAbundancesGoal[i];
        if (fabs(x[i]) > 1.0E-13) {
            retn = 1;
        }
        for (size_t j = 0; j < m_numComponents; ++j) {
            A(j, i) = - m_formulaMatrix(i,j);
        }
    }

    solve(A, x, 1, m_nelem);

    // Now apply the new direction without creating negative species.
    double par = 0.5;
    for (size_t i = 0; i < m_numComponents; ++i) {
        if (m_molNumSpecies_old[i] > 0.0) {
            par = std::max(par, -x[i] / m_molNumSpecies_old[i]);
        }
    }
    par = std::min(par, 100.0);
    par = 1.0 / par;
    if (par < 1.0 && par > 0.0) {
        retn = 2;
        par *= 0.9999;
    } else {
        par = 1.0;
    }
    for (size_t i = 0; i < m_numComponents; ++i) {
        double tmp = m_molNumSpecies_old[i] + par * x[i];
        if (tmp > 0.0) {
            m_molNumSpecies_old[i] = tmp;
        } else if (m_SSPhase[i]) {
            m_molNumSpecies_old[i] = 0.0;
        } else {
            m_molNumSpecies_old[i] *= 0.0001;
        }
    }

    // We have changed the element abundances. Calculate them again
    vcs_elab();

    // We have changed the total moles in each phase. Calculate them again
    vcs_tmoles();

    // Try some ad hoc procedures for fixing the problem
    if (retn >= 2) {
        // First find a species whose adjustment is a win-win situation.
        for (size_t kspec = 0; kspec < m_nsp; kspec++) {
            if (m_speciesUnknownType[kspec] == VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
                continue;
            }
            double saveDir = 0.0;
            bool goodSpec = true;
            for (size_t i = 0; i < m_numComponents; ++i) {
                double dir = m_formulaMatrix(kspec,i) * (m_elemAbundancesGoal[i] - m_elemAbundances[i]);
                if (fabs(dir) > 1.0E-10) {
                    if (dir > 0.0) {
                        if (saveDir < 0.0) {
                            goodSpec = false;
                            break;
                        }
                    } else if (saveDir > 0.0) {
                        goodSpec = false;
                        break;
                    }
                    saveDir = dir;
                } else if (m_formulaMatrix(kspec,i) != 0.0) {
                    goodSpec = false;
                    break;
                }
            }
            if (goodSpec) {
                int its = 0;
                double xx = 0.0;
                for (size_t i = 0; i < m_numComponents; ++i) {
                    if (m_formulaMatrix(kspec,i) != 0.0) {
                        xx += (m_elemAbundancesGoal[i] - m_elemAbundances[i]) / m_formulaMatrix(kspec,i);
                        its++;
                    }
                }
                if (its > 0) {
                    xx /= its;
                }
                m_molNumSpecies_old[kspec] += xx;
                m_molNumSpecies_old[kspec] = std::max(m_molNumSpecies_old[kspec], 1.0E-10);

                // If we are dealing with a deleted species, then we need to
                // reinsert it into the active list.
                if (kspec >= m_numSpeciesRdc) {
                    vcs_reinsert_deleted(kspec);
                    m_molNumSpecies_old[m_numSpeciesRdc - 1] = xx;
                    vcs_elab();
                    goto L_CLEANUP;
                }
                vcs_elab();
            }
        }
    }
    if (vcs_elabcheck(0)) {
        retn = 1;
        goto L_CLEANUP;
    }

    for (size_t i = 0; i < m_nelem; ++i) {
        if (m_elType[i] == VCS_ELEM_TYPE_CHARGENEUTRALITY ||
                (m_elType[i] == VCS_ELEM_TYPE_ABSPOS && m_elemAbundancesGoal[i] == 0.0)) {
            for (size_t kspec = 0; kspec < m_numSpeciesRdc; kspec++) {
                if (m_elemAbundances[i] > 0.0 && m_formulaMatrix(kspec,i) < 0.0) {
                    m_molNumSpecies_old[kspec] -= m_elemAbundances[i] / m_formulaMatrix(kspec,i);
                    m_molNumSpecies_old[kspec] = std::max(m_molNumSpecies_old[kspec], 0.0);
                    vcs_elab();
                    break;
                }
                if (m_elemAbundances[i] < 0.0 && m_formulaMatrix(kspec,i) > 0.0) {
                    m_molNumSpecies_old[kspec] -= m_elemAbundances[i] / m_formulaMatrix(kspec,i);
                    m_molNumSpecies_old[kspec] = std::max(m_molNumSpecies_old[kspec], 0.0);
                    vcs_elab();
                    break;
                }
            }
        }
    }
    if (vcs_elabcheck(1)) {
        retn = 1;
        goto L_CLEANUP;
    }

    // For electron charges element types, we try positive deltas in the species
    // concentrations to match the desired electron charge exactly.
    for (size_t i = 0; i < m_nelem; ++i) {
        double dev = m_elemAbundancesGoal[i] - m_elemAbundances[i];
        if (m_elType[i] == VCS_ELEM_TYPE_ELECTRONCHARGE && (fabs(dev) > 1.0E-300)) {
            bool useZeroed = true;
            for (size_t kspec = 0; kspec < m_numSpeciesRdc; kspec++) {
                if (dev < 0.0) {
                    if (m_formulaMatrix(kspec,i) < 0.0 && m_molNumSpecies_old[kspec] > 0.0) {
                        useZeroed = false;
                    }
                } else {
                    if (m_formulaMatrix(kspec,i) > 0.0 && m_molNumSpecies_old[kspec] > 0.0) {
                        useZeroed = false;
                    }
                }
            }
            for (size_t kspec = 0; kspec < m_numSpeciesRdc; kspec++) {
                if (m_molNumSpecies_old[kspec] > 0.0 || useZeroed) {
                    if (dev < 0.0 && m_formulaMatrix(kspec,i) < 0.0) {
                        double delta = dev / m_formulaMatrix(kspec,i);
                        m_molNumSpecies_old[kspec] += delta;
                        m_molNumSpecies_old[kspec] = std::max(m_molNumSpecies_old[kspec], 0.0);
                        vcs_elab();
                        break;
                    }
                    if (dev > 0.0 && m_formulaMatrix(kspec,i) > 0.0) {
                        double delta = dev / m_formulaMatrix(kspec,i);
                        m_molNumSpecies_old[kspec] += delta;
                        m_molNumSpecies_old[kspec] = std::max(m_molNumSpecies_old[kspec], 0.0);
                        vcs_elab();
                        break;
                    }
                }
            }
        }
    }
    if (vcs_elabcheck(1)) {
        retn = 1;
        goto L_CLEANUP;
    }

L_CLEANUP:
    ;
    vcs_tmoles();
    double l2after = 0.0;
    for (size_t i = 0; i < m_nelem; ++i) {
        l2after += pow(m_elemAbundances[i] - m_elemAbundancesGoal[i], 2);
    }
    l2after = sqrt(l2after/m_nelem);
    if (m_debug_print_lvl >= 2) {
        plogf("   ---    Elem_Abund:  Correct             Initial  "
              "              Final\n");
        for (size_t i = 0; i < m_nelem; ++i) {
            plogf("   ---       ");
            plogf("%-2.2s", m_elementName[i]);
            plogf(" %20.12E %20.12E %20.12E\n", m_elemAbundancesGoal[i], ga_save[i], m_elemAbundances[i]);
        }
        plogf("   ---            Diff_Norm:         %20.12E %20.12E\n",
              l2before, l2after);
    }
    return retn;
}

int VCS_SOLVE::vcs_setMolesLinProg()
{
    double test = -1.0E-10;

    if (m_debug_print_lvl >= 2) {
        plogf("   --- call setInitialMoles\n");
    }

    int iter = 0;
    bool abundancesOK = true;
    vector<double> sm(m_nelem * m_nelem, 0.0);
    vector<double> ss(m_nelem, 0.0);
    vector<double> sa(m_nelem, 0.0);
    vector<double> wx(m_nelem, 0.0);
    vector<double> aw(m_nsp, 0.0);

    for (size_t ik = 0; ik < m_nsp; ik++) {
        if (m_speciesUnknownType[ik] != VCS_SPECIES_INTERFACIALVOLTAGE) {
            m_molNumSpecies_old[ik] = max(0.0, m_molNumSpecies_old[ik]);
        }
    }

    if (m_debug_print_lvl >= 2) {
        printProgress(m_speciesName, m_molNumSpecies_old, m_SSfeSpecies);
    }

    bool redo = true;
    while (redo) {
        if (!vcs_elabcheck(0)) {
            if (m_debug_print_lvl >= 2) {
                plogf(" --- seMolesLinProg  Mole numbers failing element abundances\n");
                plogf(" --- seMolesLinProg  Call vcs_elcorr to attempt fix\n");
            }
            int retn = vcs_elcorr(&sm[0], &wx[0]);
            abundancesOK = (retn < 2);
        } else {
            abundancesOK = true;
        }

        // Now find the optimized basis that spans the stoichiometric
        // coefficient matrix, based on the current composition,
        // m_molNumSpecies_old[] We also calculate sc[][], the reaction matrix.
        vcs_basopt(false, test);

        if (m_debug_print_lvl >= 2) {
            plogf("iteration %d\n", iter);
        }
        redo = false;
        iter++;
        if (iter > 15) {
            break;
        }

        // loop over all reactions
        for (size_t irxn = 0; irxn < m_numRxnTot; irxn++) {
            // dg_rt is the Delta_G / RT value for the reaction
            size_t ik = m_numComponents + irxn;
            double dg_rt = m_SSfeSpecies[ik];
            double dxi_min = 1.0e10;
            const double* sc_irxn = m_stoichCoeffRxnMatrix.ptrColumn(irxn);
            for (size_t jcomp = 0; jcomp < m_nelem; jcomp++) {
                dg_rt += m_SSfeSpecies[jcomp] * sc_irxn[jcomp];
            }
            // fwd or rev direction.
            //  idir > 0 implies increasing the current species
            //  idir < 0 implies decreasing the current species
            int idir = (dg_rt < 0.0 ? 1 : -1);
            if (idir < 0) {
                dxi_min = m_molNumSpecies_old[ik];
            }

            for (size_t jcomp = 0; jcomp < m_numComponents; jcomp++) {
                double nu = sc_irxn[jcomp];
                // set max change in progress variable by
                // non-negativity requirement
                if (nu*idir < 0) {
                    double delta_xi = fabs(m_molNumSpecies_old[jcomp]/nu);
                    // if a component has nearly zero moles, redo
                    // with a new set of components
                    if (!redo && delta_xi < 1.0e-10 && (m_molNumSpecies_old[ik] >= 1.0E-10)) {
                        if (m_debug_print_lvl >= 2) {
                            plogf("   --- Component too small: %s\n", m_speciesName[jcomp]);
                        }
                        redo = true;
                    }
                    dxi_min = std::min(dxi_min, delta_xi);
                }
            }

            // step the composition by dxi_min, check against zero, since
            // we are zeroing components and species on every step.
            // Redo the iteration, if a component went from positive to zero on this step.
            double dsLocal = idir*dxi_min;
            m_molNumSpecies_old[ik] += dsLocal;
            m_molNumSpecies_old[ik] = max(0.0, m_molNumSpecies_old[ik]);
            for (size_t jcomp = 0; jcomp < m_numComponents; jcomp++) {
                bool full = false;
                if (m_molNumSpecies_old[jcomp] > 1.0E-15) {
                    full = true;
                }
                m_molNumSpecies_old[jcomp] += sc_irxn[jcomp] * dsLocal;
                m_molNumSpecies_old[jcomp] = max(0.0, m_molNumSpecies_old[jcomp]);
                if (full && m_molNumSpecies_old[jcomp] < 1.0E-60) {
                    redo = true;
                }
            }
        }

        if (m_debug_print_lvl >= 2) {
            printProgress(m_speciesName, m_molNumSpecies_old, m_SSfeSpecies);
        }
    }

    if (m_debug_print_lvl == 1) {
        printProgress(m_speciesName, m_molNumSpecies_old, m_SSfeSpecies);
        plogf("   --- setInitialMoles end\n");
    }
    if (!abundancesOK) {
        return -1;
    } else if (iter > 15) {
        return 1;
    }
    return 0;
}

double VCS_SOLVE::vcs_Total_Gibbs(double* molesSp, double* chemPot,
                                  double* tPhMoles)
{
    double g = 0.0;

    for (size_t kspec = 0; kspec < m_numSpeciesRdc; ++kspec) {
        if (m_speciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
            g += molesSp[kspec] * chemPot[kspec];
        }
    }

    return g;
}

void VCS_SOLVE::vcs_prob_specifyFully()
{
    // Whether we have an estimate or not gets overwritten on
    // the call to the equilibrium solver.
    m_temperature = m_mix->temperature();
    m_pressurePA = m_mix->pressure();
    m_Faraday_dim = Faraday / (m_temperature * GasConstant);

    for (size_t iphase = 0; iphase < m_numPhases; iphase++) {
        vcs_VolPhase* volPhase = m_VolPhaseList[iphase].get();

        volPhase->setState_TP(m_temperature, m_pressurePA);

        // Loop through each species in the current phase
        size_t nSpPhase = m_mix->phase(iphase).nSpecies();
        if ((nSpPhase == 1) && (volPhase->phiVarIndex() == 0)) {
            volPhase->setExistence(VCS_PHASE_EXIST_ALWAYS);
        } else if (volPhase->totalMoles() > 0.0) {
            volPhase->setExistence(VCS_PHASE_EXIST_YES);
        } else {
            volPhase->setExistence(VCS_PHASE_EXIST_NO);
        }
    }

    // Printout the species information: PhaseID's and mole nums
    if (m_printLvl > 1) {
        writeline('=', 80, true, true);
        writeline('=', 20, false);
        plogf(" Cantera_to_vprob: START OF PROBLEM STATEMENT ");
        writeline('=', 20);
        writeline('=', 80);
        plogf("\n");
        plogf("             Phase IDs of species\n");
        plogf("            species     phaseID        phaseName   ");
        plogf(" Initial_Estimated_kMols\n");
        for (size_t i = 0; i < m_nsp; i++) {
            size_t iphase = m_phaseID[i];
            vcs_VolPhase* VolPhase = m_VolPhaseList[iphase].get();
            plogf("%16s      %5d   %16s", m_speciesName[i].c_str(), iphase,
                  VolPhase->PhaseName.c_str());
            if (m_speciesUnknownType[i] == VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
                plogf("     Volts = %-10.5g\n", m_molNumSpecies_old[i]);
            } else {
                plogf("             %-10.5g\n", m_molNumSpecies_old[i]);
            }
        }

        // Printout of the Phase structure information
        writeline('-', 80, true, true);
        plogf("             Information about phases\n");
        plogf("  PhaseName      PhaseNum SingSpec   EqnState        NumSpec     "
              "Tmoles(kmol)\n");

        for (size_t iphase = 0; iphase < m_numPhases; iphase++) {
            vcs_VolPhase* VolPhase = m_VolPhaseList[iphase].get();
            plogf("%16s %8d %8d %16s %8d ", VolPhase->PhaseName.c_str(),
                  VolPhase->VP_ID_, VolPhase->m_singleSpecies,
                  VolPhase->eos_name(), VolPhase->nSpecies());
            plogf("%16e\n", VolPhase->totalMoles());
        }

        writeline('=', 80, true, true);
        writeline('=', 20, false);
        plogf(" Cantera_to_vprob: END OF PROBLEM STATEMENT ");
        writeline('=', 20);
        writeline('=', 80);
        plogf("\n");
    }

    // m_numRxnTot = number of noncomponents, also equal to the number of
    // reactions. Note, it's possible that the number of elements is greater
    // than the number of species. In that case set the number of reactions to
    // zero.
    m_numRxnTot = (m_nelem > m_nsp) ? 0 : m_nsp - m_nelem;
    m_numRxnRdc = m_numRxnTot;
}

void VCS_SOLVE::vcs_inest()
{
    const char* pprefix = "   --- vcs_inest: ";
    if (m_doEstimateEquil > 0) {
        // Calculate the elemental abundances
        vcs_elab();
        if (vcs_elabcheck(0)) {
            if (m_debug_print_lvl >= 2) {
                plogf("%s Initial guess passed element abundances on input\n", pprefix);
                plogf("%s m_doEstimateEquil = 1 so will use the input mole "
                      "numbers as estimates\n", pprefix);
            }
            return;
        } else if (m_debug_print_lvl >= 2) {
            plogf("%s Initial guess failed element abundances on input\n", pprefix);
            plogf("%s m_doEstimateEquil = 1 so will discard input "
                  "mole numbers and find our own estimate\n", pprefix);
        }
    }

    // temporary space for usage in this routine and in subroutines
    vector<double> sm(m_nelem*m_nelem, 0.0);
    vector<double> ss(m_nelem, 0.0);
    vector<double> sa(m_nelem, 0.0);
    vector<double> aw(m_nsp + m_nelem, 0.0);

    // Go get the estimate of the solution
    if (m_debug_print_lvl >= 2) {
        plogf("%sGo find an initial estimate for the equilibrium problem\n",
              pprefix);
    }
    double test = -1.0E20;
    size_t nrxn = m_numRxnTot;

    // CALL ROUTINE TO SOLVE MAX(CC*molNum) SUCH THAT AX*molNum = BB AND
    // molNum(I) .GE. 0.0. Note, both of these programs do this.
    vcs_setMolesLinProg();

    if (m_debug_print_lvl >= 2) {
        plogf("%s Mole Numbers returned from linear programming (vcs_inest initial guess):\n",
              pprefix);
        plogf("%s     SPECIES          MOLE_NUMBER      -SS_ChemPotential\n", pprefix);
        for (size_t kspec = 0; kspec < m_nsp; ++kspec) {
            plogf("%s     ", pprefix);
            plogf("%-12.12s", m_speciesName[kspec]);
            plogf(" %15.5g  %12.3g\n", m_molNumSpecies_old[kspec], -m_SSfeSpecies[kspec]);
        }
        plogf("%s Element Abundance Agreement returned from linear "
              "programming (vcs_inest initial guess):\n", pprefix);
        plogf("%s     Element           Goal         Actual\n", pprefix);
        for (size_t j = 0; j < m_nelem; j++) {
            if (m_elementActive[j]) {
                double tmp = 0.0;
                for (size_t kspec = 0; kspec < m_nsp; ++kspec) {
                    tmp += m_formulaMatrix(kspec,j) * m_molNumSpecies_old[kspec];
                }
                plogf("%s     ", pprefix);
                plogf("   %-9.9s", m_elementName[j]);
                plogf(" %12.3g %12.3g\n", m_elemAbundancesGoal[j], tmp);
            }
        }
        writelogendl();
    }

    // Make sure all species have positive definite mole numbers Set voltages to
    // zero for now, until we figure out what to do
    m_deltaMolNumSpecies.assign(m_deltaMolNumSpecies.size(), 0.0);
    for (size_t kspec = 0; kspec < m_nsp; ++kspec) {
        if (m_speciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
            if (m_molNumSpecies_old[kspec] <= 0.0) {
                // HKM Should eventually include logic here for non SS phases
                if (!m_SSPhase[kspec]) {
                    m_molNumSpecies_old[kspec] = 1.0e-30;
                }
            }
        } else {
            m_molNumSpecies_old[kspec] = 0.0;
        }
    }

    // Now find the optimized basis that spans the stoichiometric coefficient matrix
    vcs_basopt(false, test);

    // CALCULATE TOTAL MOLES, CHEMICAL POTENTIALS OF BASIS

    // Calculate TMoles and m_tPhaseMoles_old[]
    vcs_tmoles();

    // m_tPhaseMoles_new[] will consist of just the component moles
    m_tPhaseMoles_new.assign(m_numPhases, 1.0e-20);
    for (size_t kspec = 0; kspec < m_numComponents; ++kspec) {
        if (m_speciesUnknownType[kspec] == VCS_SPECIES_TYPE_MOLNUM) {
            m_tPhaseMoles_new[m_phaseID[kspec]] += m_molNumSpecies_old[kspec];
        }
    }
    double TMolesMultiphase = 0.0;
    for (size_t iph = 0; iph < m_numPhases; iph++) {
        if (! m_VolPhaseList[iph]->m_singleSpecies) {
            TMolesMultiphase += m_tPhaseMoles_new[iph];
        }
    }
    m_molNumSpecies_new = m_molNumSpecies_old;
    for (size_t kspec = 0; kspec < m_numComponents; ++kspec) {
        if (m_speciesUnknownType[kspec] != VCS_SPECIES_TYPE_MOLNUM) {
            m_molNumSpecies_new[kspec] = 0.0;
        }
    }
    m_feSpecies_new = m_SSfeSpecies;

    for (size_t kspec = 0; kspec < m_numComponents; ++kspec) {
        if (m_speciesUnknownType[kspec] == VCS_SPECIES_TYPE_MOLNUM) {
            if (! m_SSPhase[kspec]) {
                size_t iph = m_phaseID[kspec];
                m_feSpecies_new[kspec] += log(m_molNumSpecies_new[kspec] / m_tPhaseMoles_old[iph]);
            }
        } else {
            m_molNumSpecies_new[kspec] = 0.0;
        }
    }
    vcs_deltag(0, true, VCS_STATECALC_NEW);
    if (m_debug_print_lvl >= 2) {
        for (size_t kspec = 0; kspec < m_nsp; ++kspec) {
            plogf("%s", pprefix);
            plogf("%-12.12s", m_speciesName[kspec]);
            if (kspec < m_numComponents) {
                plogf("fe* = %15.5g ff = %15.5g\n", m_feSpecies_new[kspec],
                      m_SSfeSpecies[kspec]);
            } else {
                plogf("fe* = %15.5g ff = %15.5g dg* = %15.5g\n",
                      m_feSpecies_new[kspec], m_SSfeSpecies[kspec], m_deltaGRxn_new[kspec-m_numComponents]);
            }
        }
    }

    // ESTIMATE REACTION ADJUSTMENTS
    vector<double>& xtphMax = m_TmpPhase;
    vector<double>& xtphMin = m_TmpPhase2;
    m_deltaPhaseMoles.assign(m_deltaPhaseMoles.size(), 0.0);
    for (size_t iph = 0; iph < m_numPhases; iph++) {
        xtphMax[iph] = log(m_tPhaseMoles_new[iph] * 1.0E32);
        xtphMin[iph] = log(m_tPhaseMoles_new[iph] * 1.0E-32);
    }
    for (size_t irxn = 0; irxn < nrxn; ++irxn) {
        size_t kspec = m_indexRxnToSpecies[irxn];

        // For single species phases, we will not estimate the mole numbers. If
        // the phase exists, it stays. If it doesn't exist in the estimate, it
        // doesn't come into existence here.
        if (! m_SSPhase[kspec]) {
            size_t iph = m_phaseID[kspec];
            if (m_deltaGRxn_new[irxn] > xtphMax[iph]) {
                m_deltaGRxn_new[irxn] = 0.8 * xtphMax[iph];
            }
            if (m_deltaGRxn_new[irxn] < xtphMin[iph]) {
                m_deltaGRxn_new[irxn] = 0.8 * xtphMin[iph];
            }

            // HKM -> The TMolesMultiphase is a change of mine. It more evenly
            // distributes the initial moles amongst multiple multispecies
            // phases according to the relative values of the standard state
            // free energies. There is no change for problems with one
            // multispecies phase. It cut diamond4.vin iterations down from 62
            // to 14.
            m_deltaMolNumSpecies[kspec] = 0.5 * (m_tPhaseMoles_new[iph] + TMolesMultiphase)
                                          * exp(-m_deltaGRxn_new[irxn]);

            for (size_t k = 0; k < m_numComponents; ++k) {
                m_deltaMolNumSpecies[k] += m_stoichCoeffRxnMatrix(k,irxn) * m_deltaMolNumSpecies[kspec];
            }

            for (iph = 0; iph < m_numPhases; iph++) {
                m_deltaPhaseMoles[iph] += m_deltaMolNumPhase(iph,irxn) * m_deltaMolNumSpecies[kspec];
            }
        }
    }
    if (m_debug_print_lvl >= 2) {
        for (size_t kspec = 0; kspec < m_nsp; ++kspec) {
            if (m_speciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
                plogf("%sdirection (", pprefix);
                plogf("%-12.12s", m_speciesName[kspec]);
                plogf(") = %g", m_deltaMolNumSpecies[kspec]);
                if (m_SSPhase[kspec]) {
                    if (m_molNumSpecies_old[kspec] > 0.0) {
                        plogf(" (ssPhase exists at w = %g moles)", m_molNumSpecies_old[kspec]);
                    } else {
                        plogf(" (ssPhase doesn't exist -> stability not checked)");
                    }
                }
                writelogendl();
            }
        }
    }

    // KEEP COMPONENT SPECIES POSITIVE
    double par = 0.5;
    for (size_t kspec = 0; kspec < m_numComponents; ++kspec) {
        if (m_speciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE &&
            par < -m_deltaMolNumSpecies[kspec] / m_molNumSpecies_new[kspec]) {
            par = -m_deltaMolNumSpecies[kspec] / m_molNumSpecies_new[kspec];
        }
    }
    par = 1. / par;
    if (par <= 1.0 && par > 0.0) {
        par *= 0.8;
    } else {
        par = 1.0;
    }

    // CALCULATE NEW MOLE NUMBERS
    size_t lt = 0;
    size_t ikl = 0;
    double s1 = 0.0;
    while (true) {
        for (size_t kspec = 0; kspec < m_numComponents; ++kspec) {
            if (m_speciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
                m_molNumSpecies_old[kspec] = m_molNumSpecies_new[kspec] + par * m_deltaMolNumSpecies[kspec];
            } else {
                m_deltaMolNumSpecies[kspec] = 0.0;
            }
        }
        for (size_t kspec = m_numComponents; kspec < m_nsp; ++kspec) {
            if (m_speciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE &&
                m_deltaMolNumSpecies[kspec] != 0.0) {
                m_molNumSpecies_old[kspec] = m_deltaMolNumSpecies[kspec] * par;
            }
        }

        // We have a new w[] estimate, go get the TMoles and m_tPhaseMoles_old[]
        // values
        vcs_tmoles();
        if (lt > 0) {
            break;
        }

        // CONVERGENCE FORCING SECTION
        vcs_setFlagsVolPhases(false, VCS_STATECALC_OLD);
        vcs_dfe(VCS_STATECALC_OLD, 0, 0, m_nsp);
        double s = 0.0;
        for (size_t kspec = 0; kspec < m_nsp; ++kspec) {
            s += m_deltaMolNumSpecies[kspec] * m_feSpecies_old[kspec];
        }
        if (s == 0.0) {
            break;
        }
        if (s < 0.0 && ikl == 0) {
            break;
        }

        // TRY HALF STEP SIZE
        if (ikl == 0) {
            s1 = s;
            par *= 0.5;
            ikl = 1;
            continue;
        }

        // FIT PARABOLA THROUGH HALF AND FULL STEPS
        double xl = (1.0 - s / (s1 - s)) * 0.5;
        if (xl < 0.0) {
            // POOR DIRECTION, REDUCE STEP SIZE TO 0.2
            par *= 0.2;
        } else {
            if (xl > 1.0) {
                // TOO BIG A STEP, TAKE ORIGINAL FULL STEP
                par *= 2.0;
            } else {
                // ACCEPT RESULTS OF FORCER
                par = par * 2.0 * xl;
            }
        }
        lt = 1;
    }

    if (m_debug_print_lvl >= 2) {
        plogf("%s     Final Mole Numbers produced by inest:\n",
              pprefix);
        plogf("%s     SPECIES      MOLE_NUMBER\n", pprefix);
        for (size_t kspec = 0; kspec < m_nsp; ++kspec) {
            plogf("%s     %-12.12s %g\n",
                pprefix, m_speciesName[kspec], m_molNumSpecies_old[kspec]);
        }
    }

    // Calculate the elemental abundances
    vcs_elab();

    // If we still fail to achieve the correct elemental abundances, try to fix
    // the problem again by calling the main elemental abundances fixer routine,
    // used in the main program. This attempts to tweak the mole numbers of the
    // component species to satisfy the element abundance constraints.
    //
    // Note: We won't do this unless we have to since it involves inverting a
    // matrix.
    bool rangeCheck = vcs_elabcheck(1);
    if (!vcs_elabcheck(0)) {
        if (m_debug_print_lvl >= 2) {
            plogf("%sInitial guess failed element abundances\n", pprefix);
            plogf("%sCall vcs_elcorr to attempt fix\n", pprefix);
        }
        vcs_elcorr(&sm[0], &aw[0]);
        rangeCheck = vcs_elabcheck(1);
        if (!vcs_elabcheck(0)) {
            throw CanteraError("VCS_SOLVE::vcs_inest",
                "Initial guess still fails element abundance equations\n"
                "Inability to ever satisfy element abundance constraints is probable");
        } else {
            if (m_debug_print_lvl >= 2) {
                if (rangeCheck) {
                    plogf("%sInitial guess now satisfies element abundances\n", pprefix);
                } else {
                    plogf("%sElement Abundances RANGE ERROR\n", pprefix);
                    plogf("%s - Initial guess satisfies NC=%d element abundances, "
                          "BUT not NE=%d element abundances\n", pprefix,
                          m_numComponents, m_nelem);
                }
            }
        }
    } else {
        if (m_debug_print_lvl >= 2) {
            if (rangeCheck) {
                plogf("%sInitial guess satisfies element abundances\n", pprefix);
            } else {
                plogf("%sElement Abundances RANGE ERROR\n", pprefix);
                plogf("%s - Initial guess satisfies NC=%d element abundances, "
                      "BUT not NE=%d element abundances\n", pprefix,
                      m_numComponents, m_nelem);
            }
        }
    }

    if (m_debug_print_lvl >= 2) {
        plogf("%sTotal Dimensionless Gibbs Free Energy = %15.7E\n", pprefix,
              vcs_Total_Gibbs(&m_molNumSpecies_old[0], &m_feSpecies_new[0],
                              &m_tPhaseMoles_old[0]));
    }

    m_VCount->T_Calls_Inest++;
}

void VCS_SOLVE::vcs_SSPhase()
{
    vector<int> numPhSpecies(m_numPhases, 0);
    for (size_t kspec = 0; kspec < m_nsp; ++kspec) {
        numPhSpecies[m_phaseID[kspec]]++;
    }

    // Handle the special case of a single species in a phase that has been
    // earmarked as a multispecies phase. Treat that species as a single-species
    // phase
    for (size_t iph = 0; iph < m_numPhases; iph++) {
        m_VolPhaseList[iph]->m_singleSpecies = (numPhSpecies[iph] <= 1);
    }

    // Fill in some useful arrays here that have to do with the static
    // information concerning the phase ID of species. SSPhase = Boolean
    // indicating whether a species is in a single species phase or not.
    for (size_t kspec = 0; kspec < m_nsp; kspec++) {
        m_SSPhase[kspec] = m_VolPhaseList[m_phaseID[kspec]]->m_singleSpecies;
    }
}

void VCS_SOLVE::vcs_counters_init(int ifunc)
{
    m_VCount->Its = 0;
    m_VCount->Basis_Opts = 0;
    if (ifunc) {
        m_VCount->T_Its = 0;
        m_VCount->T_Basis_Opts = 0;
        m_VCount->T_Calls_Inest = 0;
        m_VCount->T_Calls_vcs_TP = 0;
    }
}

void VCS_SOLVE::vcs_TCounters_report()
{
    plogf("\nTCounters:   Num_Calls   Total_Its\n");
    plogf("    vcs_basopt:   %5d      %5d\n",
          m_VCount->T_Basis_Opts, m_VCount->T_Basis_Opts);
    plogf("    vcs_TP:       %5d      %5d\n",
          m_VCount->T_Calls_vcs_TP, m_VCount->T_Its);
    plogf("    vcs_inest:    %5d\n", m_VCount->T_Calls_Inest);
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
        plogf("  PhaseName      PhaseNum SingSpec   EqnState        NumSpec     "
              "Tmoles(kmol)\n");

        for (size_t iphase = 0; iphase < m_numPhases; iphase++) {
            vcs_VolPhase* Vphase = m_VolPhaseList[iphase].get();
            plogf("%16s %8d %8d ", Vphase->PhaseName,
                  Vphase->VP_ID_, Vphase->m_singleSpecies);
            plogf("%16s %8d ", Vphase->eos_name(),
                  Vphase->nSpecies());
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
        string enVP = volPhase->elementName(eVP);

        // Search for matches with the existing elements. If found, then fill in
        // the entry in the global mapping array.
        for (size_t e = 0; e < m_nelem; e++) {
            string en = m_elementName[e];
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

}
