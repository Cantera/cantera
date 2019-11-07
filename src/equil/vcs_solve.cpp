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
#include "cantera/equil/vcs_species_thermo.h"
#include "cantera/base/clockWC.h"
#include "cantera/equil/MultiPhase.h"

using namespace std;

namespace Cantera
{

int vcs_timing_print_lvl = 1;

VCS_SOLVE::VCS_SOLVE(MultiPhase* mphase, int printLvl) :
    m_printLvl(printLvl),
    m_mix(mphase),
    m_nsp(mphase->nSpecies()),
    m_nelem(0),
    m_numComponents(0),
    m_numRxnTot(0),
    m_numSpeciesRdc(mphase->nSpecies()),
    m_numRxnRdc(0),
    m_numRxnMinorZeroed(0),
    m_numPhases(mphase->nPhases()),
    m_doEstimateEquil(-1),
    m_totalMolNum(0.0),
    m_temperature(mphase->temperature()),
    m_pressurePA(mphase->pressure()),
    m_tolmaj(1.0E-8),
    m_tolmin(1.0E-6),
    m_tolmaj2(1.0E-10),
    m_tolmin2(1.0E-8),
    m_useActCoeffJac(0),
    m_totalVol(mphase->volume()),
    m_Faraday_dim(Faraday / (m_temperature * GasConstant)),
    m_VCount(0),
    m_debug_print_lvl(0),
    m_timing_print_lvl(1)
{
    m_speciesThermoList.resize(m_nsp);
    for (size_t kspec = 0; kspec < m_nsp; kspec++) {
        m_speciesThermoList[kspec].reset(new VCS_SPECIES_THERMO());
    }

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
    TPhInertMoles.resize(m_numPhases, 0.0);

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
        m_VolPhaseList[iph].reset(new vcs_VolPhase(this));
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

    if (vcs_timing_print_lvl == 0) {
        m_timing_print_lvl = 0;
    }

    VCS_SPECIES_THERMO* ts_ptr = 0;

    // Loop over the phases, transferring pertinent information
    int kT = 0;
    for (size_t iphase = 0; iphase < m_numPhases; iphase++) {
        // Get the ThermoPhase object - assume volume phase
        ThermoPhase* tPhase = &mphase->phase(iphase);
        size_t nelem = tPhase->nElements();

        // Query Cantera for the equation of state type of the current phase.
        std::string eos = tPhase->type();
        bool gasPhase = (eos == "IdealGas");

        // Find out the number of species in the phase
        size_t nSpPhase = tPhase->nSpecies();
        // Find out the name of the phase
        string phaseName = tPhase->name();

        // Call the basic vcs_VolPhase creation routine.
        // Properties set here:
        //    ->PhaseNum = phase number in the thermo problem
        //    ->GasPhase = Boolean indicating whether it is a gas phase
        //    ->NumSpecies = number of species in the phase
        //    ->TMolesInert = Inerts in the phase = 0.0 for cantera
        //    ->PhaseName  = Name of the phase
        vcs_VolPhase* VolPhase = m_VolPhaseList[iphase].get();
        VolPhase->resize(iphase, nSpPhase, nelem, phaseName.c_str(), 0.0);
        VolPhase->m_gasPhase = gasPhase;

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
        if (eos == "IdealGas") {
            VolPhase->m_eqnState = VCS_EOS_IDEAL_GAS;
        } else if (eos == "ConstDensity") {
            VolPhase->m_eqnState = VCS_EOS_CONSTANT;
        } else if (eos == "StoichSubstance") {
            VolPhase->m_eqnState = VCS_EOS_STOICH_SUB;
        } else if (eos == "IdealSolidSoln") {
            VolPhase->m_eqnState = VCS_EOS_IDEAL_SOLN;
        } else if (eos == "Surf" || eos == "Edge") {
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

            // Get a pointer to the thermo object
            ts_ptr = m_speciesThermoList[kT].get();

            // Fill in the vcs_SpeciesProperty structure
            vcs_SpeciesProperties* sProp = VolPhase->speciesProperty(k);
            sProp->NumElements = m_nelem;
            sProp->SpName = mphase->speciesName(kT);
            sProp->SpeciesThermo = ts_ptr;
            sProp->WtSpecies = tPhase->molecularWeight(k);
            sProp->FormulaMatrixCol.resize(m_nelem, 0.0);
            for (size_t e = 0; e < m_nelem; e++) {
                sProp->FormulaMatrixCol[e] = m_formulaMatrix(kT,e);
            }
            sProp->Charge = tPhase->charge(k);
            sProp->SurfaceSpecies = false;
            sProp->VolPM = 0.0;

            // Transfer the thermo specification of the species
            //              vsolve->SpeciesThermo[]

            // Add lookback connectivity into the thermo object first
            ts_ptr->IndexPhase = iphase;
            ts_ptr->IndexSpeciesPhase = k;
            ts_ptr->OwningPhase = VolPhase;

            // get a reference to the Cantera species thermo.
            MultiSpeciesThermo& sp = tPhase->speciesThermo();

            int spType = sp.reportType(k);
            if (spType == SIMPLE) {
                double c[4];
                double minTemp, maxTemp, refPressure;
                sp.reportParams(k, spType, c, minTemp, maxTemp, refPressure);
                ts_ptr->SS0_Model = VCS_SS0_CONSTANT;
                ts_ptr->SS0_T0 = c[0];
                ts_ptr->SS0_H0 = c[1];
                ts_ptr->SS0_S0 = c[2];
                ts_ptr->SS0_Cp0 = c[3];
                if (gasPhase) {
                    ts_ptr->SSStar_Model = VCS_SSSTAR_IDEAL_GAS;
                    ts_ptr->SSStar_Vol_Model = VCS_SSVOL_IDEALGAS;
                } else {
                    ts_ptr->SSStar_Model = VCS_SSSTAR_CONSTANT;
                    ts_ptr->SSStar_Vol_Model = VCS_SSVOL_CONSTANT;
                }
            } else {
                if (m_printLvl > 2) {
                    plogf("vcs_Cantera_convert: Species Type %d not known \n",
                          spType);
                }
                ts_ptr->SS0_Model = VCS_SS0_NOTHANDLED;
                ts_ptr->SSStar_Model = VCS_SSSTAR_NOTHANDLED;
            }

            // Transfer the Volume Information -> NEEDS WORK
            if (gasPhase) {
                ts_ptr->SSStar_Vol_Model = VCS_SSVOL_IDEALGAS;
                ts_ptr->SSStar_Vol0 = 82.05 * 273.15 / 1.0;
            } else {
                vector_fp phaseTermCoeff(nSpPhase, 0.0);
                int nCoeff;
                tPhase->getParameters(nCoeff, &phaseTermCoeff[0]);
                ts_ptr->SSStar_Vol_Model = VCS_SSVOL_CONSTANT;
                ts_ptr->SSStar_Vol0 = phaseTermCoeff[k];
            }
            kT++;
        }

        VolPhase->setMolesFromVCS(VCS_STATECALC_OLD, &m_molNumSpecies_old[0]);

        // Now, calculate a sample naught Gibbs free energy calculation
        // at the specified temperature.
        for (size_t k = 0; k < nSpPhase; k++) {
            vcs_SpeciesProperties* sProp = VolPhase->speciesProperty(k);
            ts_ptr = sProp->SpeciesThermo;
            ts_ptr->SS0_feSave = VolPhase->G0_calc_one(k)/ GasConstant;
            ts_ptr->SS0_TSave = m_temperature;
        }
    }

    // Transfer initial element abundances based on the species mole numbers
    for (size_t j = 0; j < m_nelem; j++) {
        for (size_t kspec = 0; kspec < m_nsp; kspec++) {
            if (m_speciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
                m_elemAbundancesGoal[j] += m_formulaMatrix(kspec,j) * m_molNumSpecies_old[kspec];
            }
        }
        if (m_elType[j] == VCS_ELEM_TYPE_LATTICERATIO && m_elemAbundancesGoal[j] < 1.0E-10) {
            m_elemAbundancesGoal[j] = 0.0;
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
        plogf("  PhaseName    PhaseNum SingSpec GasPhase EqnState NumSpec");
        plogf("  TMolesInert       Tmoles(kmol)\n");

        for (size_t iphase = 0; iphase < m_numPhases; iphase++) {
            vcs_VolPhase* VolPhase = m_VolPhaseList[iphase].get();
            plogf("%16s %5d %5d %8d %16s %8d %16e ", VolPhase->PhaseName.c_str(),
                  VolPhase->VP_ID_, VolPhase->m_singleSpecies,
                  VolPhase->m_gasPhase, VolPhase->eos_name(),
                  VolPhase->nSpecies(), VolPhase->totalMolesInert());
            plogf("%16e\n", VolPhase->totalMoles());
        }

        writeline('=', 80, true, true);
        writeline('=', 16, false);
        plogf(" Cantera_to_vprob: END OF PROBLEM STATEMENT ");
        writeline('=', 20);
        writeline('=', 80);
        plogf("\n");
    }

    // TPhInertMoles[] -> must be copied over here
    for (size_t iph = 0; iph < m_numPhases; iph++) {
        vcs_VolPhase* Vphase = m_VolPhaseList[iph].get();
        TPhInertMoles[iph] = Vphase->totalMolesInert();
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
    std::vector<size_t> numPhSp(m_numPhases, 0);
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

}

VCS_SOLVE::~VCS_SOLVE()
{
    vcs_delete_memory();
}

void VCS_SOLVE::vcs_delete_memory()
{
    delete m_VCount;
    m_VCount = 0;

    m_nsp = 0;
    m_nelem = 0;
    m_numComponents = 0;

}

int VCS_SOLVE::vcs(int ipr, int ip1, int maxit)
{
    clockWC tickTock;

    // This function is called to copy the public data and the current
    // problem specification into the current object's data structure.
    vcs_prob_specifyFully();

    prob_report(m_printLvl);

    // Prep the problem data
    //    - adjust the identity of any phases
    //    - determine the number of components in the problem
    int retn = vcs_prep(ip1);
    if (retn != 0) {
        plogf("vcs_prep_oneTime returned a bad status, %d: bailing!\n",
              retn);
        return retn;
    }

    // Once we have defined the global internal data structure defining the
    // problem, then we go ahead and solve the problem.
    //
    // (right now, all we do is solve fixed T, P problems. Methods for other
    // problem types will go in at this level. For example, solving for
    // fixed T, V problems will involve a 2x2 Newton's method, using loops
    // over vcs_TP() to calculate the residual and Jacobian)
    int iconv = vcs_TP(ipr, ip1, maxit, m_temperature, m_pressurePA);

    // If requested to print anything out, go ahead and do so;
    if (ipr > 0) {
        vcs_report(iconv);
    }

    vcs_prob_update();

    // Report on the time if requested to do so
    double te = tickTock.secondsWC();
    m_VCount->T_Time_vcs += te;
    if (ipr > 0 || ip1 > 0) {
        vcs_TCounters_report(m_timing_print_lvl);
    }

    // FILL IN
    if (iconv < 0) {
        plogf("ERROR: FAILURE its = %d!\n", m_VCount->Its);
    } else if (iconv == 1) {
        plogf("WARNING: RANGE SPACE ERROR encountered\n");
    }
    return iconv;
}

void VCS_SOLVE::vcs_prob_specifyFully()
{
    // Whether we have an estimate or not gets overwritten on
    // the call to the equilibrium solver.
    m_temperature = m_mix->temperature();
    m_pressurePA = m_mix->pressure();
    m_Faraday_dim = Faraday / (m_temperature * GasConstant);
    m_totalVol = m_mix->volume();

    vector<size_t> invSpecies(m_nsp);
    for (size_t k = 0; k < m_nsp; k++) {
        invSpecies[m_speciesMapIndex[k]] = k;
    }

    for (size_t iphase = 0; iphase < m_numPhases; iphase++) {
        ThermoPhase* tPhase = &m_mix->phase(iphase);
        vcs_VolPhase* volPhase = m_VolPhaseList[iphase].get();

        volPhase->setState_TP(m_temperature, m_pressurePA);

        // Loop through each species in the current phase
        size_t nSpPhase = tPhase->nSpecies();
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
        plogf("  PhaseName    PhaseNum SingSpec GasPhase EqnState NumSpec");
        plogf("  TMolesInert       Tmoles(kmol)\n");

        for (size_t iphase = 0; iphase < m_numPhases; iphase++) {
            vcs_VolPhase* VolPhase = m_VolPhaseList[iphase].get();
            plogf("%16s %5d %5d %8d %16s %8d %16e ", VolPhase->PhaseName.c_str(),
                  VolPhase->VP_ID_, VolPhase->m_singleSpecies,
                  VolPhase->m_gasPhase, VolPhase->eos_name(),
                  VolPhase->nSpecies(), VolPhase->totalMolesInert());
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
    if (m_nelem > m_nsp) {
        m_numRxnTot = 0;
    } else {
        m_numRxnTot = m_nsp - m_nelem;
    }
    m_numRxnRdc = m_numRxnTot;
}

void VCS_SOLVE::vcs_prob_update()
{
    // Transfer the information back to the MultiPhase object. Note we don't
    // just call setMoles, because some multispecies solution phases may be
    // zeroed out, and that would cause a problem for that routine. Also, the
    // mole fractions of such zeroed out phases actually contain information
    // about likely reemergent states.
    m_mix->uploadMoleFractionsFromPhases();
    for (size_t ip = 0; ip < m_numPhases; ip++) {
        m_mix->setPhaseMoles(ip, m_VolPhaseList[ip]->totalMoles());
    }
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
        vcs_VolPhase* Vphase = m_VolPhaseList[iphase].get();
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
