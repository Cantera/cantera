/**
 * @file vcs_phaseStability.cpp
 *  Implementation class for functions associated with determining the stability of a phase
 *   (see Class \link VCSnonideal::VCS_SOLVE VCS_SOLVE\endlink and \ref equilfunctions ).
 */
#include "cantera/equil/vcs_solve.h"
#include "cantera/equil/vcs_internal.h"
#include "cantera/equil/vcs_species_thermo.h"
#include "cantera/equil/vcs_VolPhase.h"

#include <cstdio>
#include <algorithm>

using namespace std;

namespace VCSnonideal
{

bool VCS_SOLVE::vcs_popPhasePossible(const size_t iphasePop) const
{
    vcs_VolPhase* Vphase = m_VolPhaseList[iphasePop];

#ifdef DEBUG_MODE
    int existence = Vphase->exists();
    if (existence > 0) {
        printf("ERROR vcs_popPhasePossible called for a phase that exists!");
        std::exit(-1);
    }
#endif

    /*
     * Loop through all of the species in the phase. We say the phase
     * can be popped, if there is one species in the phase that can be
     * popped.
     */
    for (size_t k = 0; k < Vphase->nSpecies(); k++) {
        size_t kspec = Vphase->spGlobalIndexVCS(k);
#ifdef DEBUG_MODE
        if (m_molNumSpecies_old[kspec] > 0.0) {
            printf("ERROR vcs_popPhasePossible we shouldn't be here %lu %g > 0.0",
                   kspec, m_molNumSpecies_old[kspec]);
            exit(-1);
        }
#endif
        size_t irxn = kspec - m_numComponents;
        if (kspec >= m_numComponents) {
            bool iPopPossible = true;
            for (size_t j = 0; j < m_numComponents; ++j) {
                if (m_elType[j] == VCS_ELEM_TYPE_ABSPOS) {
                    double stoicC = m_stoichCoeffRxnMatrix[irxn][j];
                    if (stoicC != 0.0) {
                        double negChangeComp = - stoicC * 1.0;
                        if (negChangeComp > 0.0) {
                            // TODO: We may have to come up with a tolerance here
                            if (m_molNumSpecies_old[j] <= VCS_DELETE_ELEMENTABS_CUTOFF*0.5) {
                                iPopPossible = false;
                            }
                        }
                    }
                }
            }
            if (iPopPossible) {
                return true;
            }
        } else {
            /*
             * We are here when the species in the phase is a component. Its mole number is zero.
             * We loop through the regular reaction looking for a reaction that can pop the
             * component.
             */
            //printf("WE are here at new logic - CHECK\n");
            for (size_t jrxn = 0; jrxn < m_numRxnRdc; jrxn++) {
                bool foundJrxn = false;
                // First, if the component is a product of the reaction
                if (m_stoichCoeffRxnMatrix[jrxn][kspec] > 0.0) {
                    foundJrxn = true;
                    for (size_t kcomp = 0; kcomp < m_numComponents; kcomp++) {
                        if (m_stoichCoeffRxnMatrix[jrxn][kcomp] < 0.0) {
                            if (m_molNumSpecies_old[kcomp] <= VCS_DELETE_ELEMENTABS_CUTOFF*0.5) {
                                foundJrxn = false;
                            }
                        }
                    }
                    if (foundJrxn) {
                        //printf("We have found a component phase pop! CHECK1 \n");
                        return true;
                    }
                }
                // Second we are here if the component is a reactant in the reaction, and the reaction goes backwards.
                else if (m_stoichCoeffRxnMatrix[jrxn][kspec] < 0.0) {
                    foundJrxn = true;
                    size_t jspec = jrxn + m_numComponents;
                    if (m_molNumSpecies_old[jspec] <= VCS_DELETE_ELEMENTABS_CUTOFF*0.5) {
                        foundJrxn = false;
                        continue;
                    }
                    for (size_t kcomp = 0; kcomp < m_numComponents; kcomp++) {
                        if (m_stoichCoeffRxnMatrix[jrxn][kcomp] > 0.0) {
                            if (m_molNumSpecies_old[kcomp] <= VCS_DELETE_ELEMENTABS_CUTOFF*0.5) {
                                foundJrxn = false;
                            }
                        }
                    }
                    if (foundJrxn) {
                        //printf("We have found a component phase pop! CHECK2 \n");
                        return true;
                    }
                }
            }
        }
    }
    return false;
}

int  VCS_SOLVE::vcs_phasePopDeterminePossibleList()
{
    int nfound = 0;
    vcs_VolPhase* Vphase = 0;
    double stoicC;
    double molComp;
    std::vector<int> linkedPhases;
    phasePopProblemLists_.clear();

    /*
     *  This is a vector over each component.
     *  For zeroed components it lists the phases, which are currently zeroed,
     *     which have a species with a positive stoichiometric value wrt the component.
     *     Therefore, we could pop the component species and pop that phase at the same time
     *     if we considered no other factors than keeping the component mole number positve.
     *
     *     It does not count species with positive stoichiometric values if that species
     *     already has a positive mole number. The phase is already popped.
     */
    std::vector< std::vector<size_t> > zeroedComponentLinkedPhasePops(m_numComponents);
    /*
     *  The logic below calculates zeroedComponentLinkedPhasePops
     */
    for (size_t j = 0; j < m_numComponents; j++) {
        if (m_elType[j] == VCS_ELEM_TYPE_ABSPOS) {
            molComp = m_molNumSpecies_old[j];
            if (molComp <= 0.0) {
                std::vector<size_t> &jList = zeroedComponentLinkedPhasePops[j];
                size_t iph = m_phaseID[j];
                jList.push_back(iph);
                for (size_t irxn = 0; irxn < m_numRxnTot; irxn++) {
                    size_t kspec = irxn +  m_numComponents;
                    iph = m_phaseID[kspec];
                    Vphase = m_VolPhaseList[iph];
                    int existence = Vphase->exists();
                    if (existence < 0) {
                        stoicC = m_stoichCoeffRxnMatrix[irxn][j];
                        if (stoicC > 0.0) {
                            if (std::find(jList.begin(), jList.end(), iph) != jList.end()) {
                                jList.push_back(iph);
                            }
                        }
                    }
                }
            }
        }
    }
    /*
     *   This is a vector over each zeroed phase
     *   For zeroed phases, it lists the components, which are currently zeroed,
     *     which have a species with a negative stoichiometric value wrt one or more species in the phase.
     *     Cut out components which have a pos stoichiometric value with another species in the phase.
     */
    std::vector< std::vector<size_t> > zeroedPhaseLinkedZeroComponents(m_numPhases);
    /*
     *   The logic below calculates  zeroedPhaseLinkedZeroComponents
     */
    for (size_t iph = 0; iph < m_numPhases; iph++) {
        std::vector<size_t> &iphList = zeroedPhaseLinkedZeroComponents[iph];
        iphList.clear();
        Vphase = m_VolPhaseList[iph];
        int existence = Vphase->exists();
        if (existence < 0) {

            linkedPhases.clear();
            size_t nsp = Vphase->nSpecies();
            for (size_t k = 0; k < nsp; k++) {
                size_t kspec = Vphase->spGlobalIndexVCS(k);
                size_t irxn = kspec - m_numComponents;

                for (size_t j = 0; j < m_numComponents; j++) {
                    if (m_elType[j] == VCS_ELEM_TYPE_ABSPOS) {
                        molComp = m_molNumSpecies_old[j];
                        if (molComp <= 0.0) {
                            stoicC = m_stoichCoeffRxnMatrix[irxn][j];
                            if (stoicC < 0.0) {
                                bool foundPos = false;
                                for (size_t kk = 0; kk < nsp; kk++) {
                                    size_t kkspec  = Vphase->spGlobalIndexVCS(kk);
                                    if (kkspec >= m_numComponents) {
                                        size_t iirxn = kkspec - m_numComponents;
                                        if (m_stoichCoeffRxnMatrix[iirxn][j] > 0.0) {
                                            foundPos = true;
                                        }
                                    }
                                }
                                if (!foundPos) {
                                    if (std::find(iphList.begin(), iphList.end(), j) != iphList.end()) {
                                        iphList.push_back(j);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    /*
     *  Now fill in the   phasePopProblemLists_  list.
     *
     */
    for (size_t iph = 0; iph < m_numPhases; iph++) {
        Vphase = m_VolPhaseList[iph];
        int existence = Vphase->exists();
        if (existence < 0) {
            std::vector<size_t> &iphList = zeroedPhaseLinkedZeroComponents[iph];
            std::vector<size_t> popProblem(0);
            popProblem.push_back(iph);
            for (size_t i = 0; i < iphList.size(); i++) {
                size_t j = iphList[i];
                std::vector<size_t> &jList = zeroedComponentLinkedPhasePops[j];
                for (size_t jjl = 0; jjl < jList.size(); jjl++) {
                    size_t jph = jList[jjl];
                    if (std::find(popProblem.begin(), popProblem.end(), jph) != popProblem.end()) {
                        popProblem.push_back(jph);
                    }
                }
            }
            phasePopProblemLists_.push_back(popProblem);
        }
    }

    return nfound;
}

size_t VCS_SOLVE::vcs_popPhaseID(std::vector<size_t> & phasePopPhaseIDs)
{
    size_t iphasePop = npos;
    size_t irxn, kspec;
    doublereal FephaseMax = -1.0E30;
    doublereal Fephase = -1.0E30;
    vcs_VolPhase* Vphase = 0;


#ifdef DEBUG_MODE
    char anote[128];
    if (m_debug_print_lvl >= 2) {
        plogf("   --- vcs_popPhaseID() called\n");
        plogf("   ---   Phase                 Status       F_e        MoleNum\n");
        plogf("   --------------------------------------------------------------------------\n");
    }
#endif
    for (size_t iph = 0; iph < m_numPhases; iph++) {
        Vphase = m_VolPhaseList[iph];
        int existence = Vphase->exists();
#ifdef DEBUG_MODE
        strcpy(anote, "");
#endif
        if (existence > 0) {

#ifdef DEBUG_MODE
            if (m_debug_print_lvl >= 2) {
                plogf("  ---    %18s %5d           NA       %11.3e\n",
                      Vphase->PhaseName.c_str(),
                      existence,
                      m_tPhaseMoles_old[iph]);
            }
#endif
        } else {
            if (Vphase->m_singleSpecies) {
                /***********************************************************************
                 *
                 *  Single Phase Stability Resolution
                 *
                 ***********************************************************************/
                kspec = Vphase->spGlobalIndexVCS(0);
                irxn = kspec - m_numComponents;
                doublereal deltaGRxn = m_deltaGRxn_old[irxn];
                Fephase = exp(-deltaGRxn) - 1.0;
                if (Fephase > 0.0) {
#ifdef DEBUG_MODE
                    strcpy(anote," (ready to be birthed)");
#endif
                    if (Fephase > FephaseMax) {
                        iphasePop = iph;
                        FephaseMax = Fephase;
#ifdef DEBUG_MODE
                        strcpy(anote," (chosen to be birthed)");
#endif
                    }
                }
#ifdef DEBUG_MODE
                if (Fephase < 0.0) {
                    strcpy(anote," (not stable)");
                    if (m_tPhaseMoles_old[iph] > 0.0) {
                        printf("shouldn't be here\n");
                        exit(-1);
                    }
                }
#endif

#ifdef DEBUG_MODE
                if (m_debug_print_lvl >= 2) {
                    plogf("  ---    %18s %5d %10.3g %10.3g %s\n",
                          Vphase->PhaseName.c_str(),
                          existence, Fephase,
                          m_tPhaseMoles_old[iph], anote);
                }
#endif

            } else {
                /***********************************************************************
                 *
                 * MultiSpecies Phase Stability Resolution
                 *
                 ***********************************************************************/
                if (vcs_popPhasePossible(iph)) {
                    Fephase = vcs_phaseStabilityTest(iph);
                    if (Fephase > 0.0) {
                        if (Fephase > FephaseMax) {
                            iphasePop = iph;
                            FephaseMax = Fephase;
                        }
                    } else {
                        if (Fephase > FephaseMax) {
                            FephaseMax = Fephase;
                        }
                    }
#ifdef DEBUG_MODE
                    if (m_debug_print_lvl >= 2) {
                        plogf("  ---    %18s %5d  %11.3g %11.3g\n",
                              Vphase->PhaseName.c_str(),
                              existence, Fephase,
                              m_tPhaseMoles_old[iph]);
                    }
#endif
                } else {
#ifdef DEBUG_MODE
                    if (m_debug_print_lvl >= 2) {
                        plogf("  ---    %18s %5d   blocked  %11.3g\n",
                              Vphase->PhaseName.c_str(),
                              existence, m_tPhaseMoles_old[iph]);
                    }
#endif
                }
            }
        }
    }
    phasePopPhaseIDs.resize(0);
    if (iphasePop != npos) {
        phasePopPhaseIDs.push_back(iphasePop);
    }

    /*
     *   Insert logic here to figure out if phase pops are linked together. Only do one linked
     *   pop at a time.
     */

#ifdef DEBUG_MODE
    if (m_debug_print_lvl >= 2) {
        plogf("   ---------------------------------------------------------------------\n");
    }
#endif
    return iphasePop;
}

int VCS_SOLVE::vcs_popPhaseRxnStepSizes(const size_t iphasePop)
{
    vcs_VolPhase* Vphase = m_VolPhaseList[iphasePop];
    // Identify the first species in the phase
    size_t kspec = Vphase->spGlobalIndexVCS(0);
    // Identify the formation reaction for that species
    size_t irxn = kspec - m_numComponents;
    std::vector<size_t> creationGlobalRxnNumbers;

    doublereal s;
    // Calculate the initial moles of the phase being born.
    //   Here we set it to 10x of the value which would cause the phase to be
    //   zeroed out within the algorithm.  We may later adjust the value.
    doublereal tPhaseMoles = 10. * m_totalMolNum * VCS_DELETE_PHASE_CUTOFF;


#ifdef DEBUG_MODE
    int existence = Vphase->exists();
    if (existence > 0) {
        printf("ERROR vcs_popPhaseRxnStepSizes called for a phase that exists!");
        exit(-1);
    }
    char anote[256];
    if (m_debug_print_lvl >= 2) {
        plogf("  ---  vcs_popPhaseRxnStepSizes() called to pop phase %s %d into existence\n",
              Vphase->PhaseName.c_str(), iphasePop);
    }
#endif
    // Section for a single-species phase
    if (Vphase->m_singleSpecies) {
        s = 0.0;
        double* dnPhase_irxn = m_deltaMolNumPhase[irxn];
        for (size_t j = 0; j < m_numComponents; ++j) {
            if (!m_SSPhase[j]) {
                if (m_molNumSpecies_old[j] > 0.0) {
                    s += SQUARE(m_stoichCoeffRxnMatrix[irxn][j]) / m_molNumSpecies_old[j];
                }
            }
        }
        for (size_t j = 0; j < m_numPhases; j++) {
            Vphase = m_VolPhaseList[j];
            if (! Vphase->m_singleSpecies) {
                if (m_tPhaseMoles_old[j] > 0.0) {
                    s -= SQUARE(dnPhase_irxn[j]) / m_tPhaseMoles_old[j];
                }
            }
        }
        if (s != 0.0) {
            double s_old = s;
            s = vcs_Hessian_diag_adj(irxn, s_old);
#ifdef DEBUG_MODE
            if (s_old != s) {
                sprintf(anote, "Normal calc: diag adjusted from %g "
                        "to %g due to act coeff",  s_old, s);
            }
#endif
            m_deltaMolNumSpecies[kspec] = -m_deltaGRxn_new[irxn] / s;
        } else {
            // Ok, s is equal to zero. We can not apply a sophisticated theory
            // to birth the phase. Just pick a small delta and go with it.
            m_deltaMolNumSpecies[kspec] = tPhaseMoles;
        }

        /*
         * section to do damping of the m_deltaMolNumSpecies[]
         */
        for (size_t j = 0; j < m_numComponents; ++j) {
            double stoicC = m_stoichCoeffRxnMatrix[irxn][j];
            if (stoicC != 0.0) {
                if (m_elType[j] == VCS_ELEM_TYPE_ABSPOS) {
                    double negChangeComp = - stoicC * m_deltaMolNumSpecies[kspec];
                    if (negChangeComp > m_molNumSpecies_old[j]) {
                        if (m_molNumSpecies_old[j] > 0.0) {
#ifdef DEBUG_MODE
                            sprintf(anote, "Delta damped from %g "
                                    "to %g due to component %lu (%10s) going neg", m_deltaMolNumSpecies[kspec],
                                    -m_molNumSpecies_old[j]/stoicC, j,  m_speciesName[j].c_str());
#endif
                            m_deltaMolNumSpecies[kspec] = - 0.5 * m_molNumSpecies_old[j] / stoicC;
                        } else {
#ifdef DEBUG_MODE
                            sprintf(anote, "Delta damped from %g "
                                    "to %g due to component %lu (%10s) zero", m_deltaMolNumSpecies[kspec],
                                    -m_molNumSpecies_old[j]/stoicC, j,  m_speciesName[j].c_str());
#endif
                            m_deltaMolNumSpecies[kspec] = 0.0;
                        }
                    }
                }
            }
        }
        // Implement a damping term that limits m_deltaMolNumSpecies to the size of the mole number
        if (-m_deltaMolNumSpecies[kspec] > m_molNumSpecies_old[kspec]) {
#ifdef DEBUG_MODE
            sprintf(anote, "Delta damped from %g "
                    "to %g due to %s going negative", m_deltaMolNumSpecies[kspec],
                    -m_molNumSpecies_old[kspec],  m_speciesName[kspec].c_str());
#endif
            m_deltaMolNumSpecies[kspec] = -m_molNumSpecies_old[kspec];
        }


    } else {
        vector<doublereal> fracDelta(Vphase->nSpecies());
        vector<doublereal> X_est(Vphase->nSpecies());
        fracDelta = Vphase->creationMoleNumbers(creationGlobalRxnNumbers);

        double sumFrac = 0.0;
        for (size_t k = 0; k < Vphase->nSpecies(); k++) {
            sumFrac += fracDelta[k];
        }
        for (size_t k = 0; k < Vphase->nSpecies(); k++) {
            X_est[k] = fracDelta[k] / sumFrac;
        }

        doublereal deltaMolNumPhase = tPhaseMoles;
        doublereal damp = 1.0;
        m_deltaGRxn_tmp = m_molNumSpecies_old;
        double* molNumSpecies_tmp = DATA_PTR(m_deltaGRxn_tmp);


        for (size_t k = 0; k < Vphase->nSpecies(); k++) {
            kspec = Vphase->spGlobalIndexVCS(k);
            double delmol =  deltaMolNumPhase * X_est[k];
            if (kspec >= m_numComponents) {
                irxn = kspec - m_numComponents;
                for (size_t j = 0; j < m_numComponents; ++j) {
                    double stoicC = m_stoichCoeffRxnMatrix[irxn][j];
                    if (stoicC != 0.0) {
                        if (m_elType[j] == VCS_ELEM_TYPE_ABSPOS) {
                            molNumSpecies_tmp[j] +=  stoicC * delmol;
                        }
                    }
                }
            }
        }

        doublereal ratioComp = 0.0;
        for (size_t j = 0; j < m_numComponents; ++j) {
            double deltaJ = m_molNumSpecies_old[j] - molNumSpecies_tmp[j];
            if (molNumSpecies_tmp[j] < 0.0) {
                ratioComp = 1.0;
                if (deltaJ > 0.0) {
                    double delta0 = m_molNumSpecies_old[j];
                    double dampj = delta0 / deltaJ * 0.9;
                    if (dampj < damp) {
                        damp = dampj;
                    }
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
        // Here we create a damp > 1 to account for this possibility.
        // We adjust upwards to make sure that a component in an existing multispecies
        // phase is modified by a factor of 1/1000.
        if (ratioComp > 1.0E-30) {
            if (ratioComp < 0.001) {
                damp = 0.001 / ratioComp;
            }
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
                if (X_est[k] > 1.0E-3) {
                    m_speciesStatus[kspec] = VCS_SPECIES_MAJOR;
                } else {
                    m_speciesStatus[kspec] = VCS_SPECIES_MINOR;
                }
            }
        }

    }

    return 0;
}

double VCS_SOLVE::vcs_phaseStabilityTest(const size_t iph)
{
    /*
     * We will use the _new state calc here
     */
    size_t kspec, irxn, k, i, kc, kc_spec;
    vcs_VolPhase* Vphase = m_VolPhaseList[iph];
    const size_t nsp = Vphase->nSpecies();
    doublereal deltaGRxn;
    int minNumberIterations = 3;
    if (nsp <= 1) {
	minNumberIterations = 1;
    }

    // We will do a full newton calculation later, but for now, ...
    bool doSuccessiveSubstitution = true;
    double funcPhaseStability;
    vector<doublereal> X_est(nsp, 0.0);
    vector<doublereal> delFrac(nsp, 0.0);
    vector<doublereal> E_phi(nsp, 0.0);
    vector<doublereal> fracDelta_new(nsp, 0.0);
    vector<doublereal> fracDelta_old(nsp, 0.0);
    vector<doublereal> fracDelta_raw(nsp, 0.0);
    vector<size_t> creationGlobalRxnNumbers(nsp, npos);
    vcs_dcopy(VCS_DATA_PTR(m_deltaGRxn_Deficient), VCS_DATA_PTR(m_deltaGRxn_old), m_numRxnRdc);

    vector<doublereal> m_feSpecies_Deficient(m_numComponents, 0.0);
    doublereal damp = 1.0;
    doublereal dampOld = 1.0;
    doublereal normUpdate = 1.0;
    doublereal normUpdateOld = 1.0;
    doublereal sum = 0.0;
    doublereal dirProd = 0.0;
    doublereal dirProdOld = 0.0;

    // get the activity coefficients
    Vphase->sendToVCS_ActCoeff(VCS_STATECALC_OLD, VCS_DATA_PTR(m_actCoeffSpecies_new));

    // Get the stored estimate for the composition of the phase if
    // it gets created
    fracDelta_new = Vphase->creationMoleNumbers(creationGlobalRxnNumbers);


    std::vector<size_t> componentList;

    for (k = 0; k < nsp; k++) {
        kspec = Vphase->spGlobalIndexVCS(k);
        if (kspec < m_numComponents) {
            componentList.push_back(k);
        }
    }

    for (k = 0; k < m_numComponents; k++) {
        m_feSpecies_Deficient[k]  = m_feSpecies_old[k];
    }
    normUpdate = 0.1 * vcs_l2norm(fracDelta_new);
    damp = 1.0E-2;

    if (doSuccessiveSubstitution) {

#ifdef DEBUG_MODE
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
#endif

        for (k = 0; k < nsp; k++) {
            if (fracDelta_new[k] < 1.0E-13) {
                fracDelta_new[k] = 1.0E-13;
            }
        }
        bool converged = false;
        for (int its = 0; its < 200  && (!converged); its++) {

            dampOld = damp;
            normUpdateOld = normUpdate;
            fracDelta_old = fracDelta_new;
            dirProdOld = dirProd;

            // Given a set of fracDelta's, we calculate the fracDelta's
            // for the component species, if any
            for (i = 0; i < componentList.size(); i++) {
                kc = componentList[i];
                kc_spec = Vphase->spGlobalIndexVCS(kc);
                fracDelta_old[kc] = 0.0;
                for (k = 0; k <  nsp; k++) {
                    kspec = Vphase->spGlobalIndexVCS(k);
                    if (kspec >= m_numComponents) {
                        irxn = kspec - m_numComponents;
                        fracDelta_old[kc] += m_stoichCoeffRxnMatrix[irxn][kc_spec] *  fracDelta_old[k];
                    }
                }
            }

            // Now, calculate the predicted mole fractions, X_est[k]
            double sumFrac = 0.0;
            for (k = 0; k < nsp; k++) {
                sumFrac += fracDelta_old[k];
            }
            // Necessary because this can be identically zero. -> we need to fix this algorithm!
            if (sumFrac <= 0.0) {
                sumFrac = 1.0;
            }
            double sum_Xcomp = 0.0;
            for (k = 0; k < nsp; k++) {
                X_est[k] = fracDelta_old[k] / sumFrac;
                kc_spec = Vphase->spGlobalIndexVCS(k);
                if (kc_spec < m_numComponents) {
                    sum_Xcomp += X_est[k];
                }
            }


            /*
             * Feed the newly formed estimate of the mole fractions back into the
             * ThermoPhase object
             */
            Vphase->setMoleFractionsState(0.0, VCS_DATA_PTR(X_est), VCS_STATECALC_PHASESTABILITY);

            /*
             *   get the activity coefficients
             */
            Vphase->sendToVCS_ActCoeff(VCS_STATECALC_OLD, VCS_DATA_PTR(m_actCoeffSpecies_new));

            /*
             * First calculate altered chemical potentials for component species
             * belonging to this phase.
             */
            for (i = 0; i < componentList.size(); i++) {
                kc = componentList[i];
                kc_spec = Vphase->spGlobalIndexVCS(kc);
                if (X_est[kc] > VCS_DELETE_MINORSPECIES_CUTOFF) {
                    m_feSpecies_Deficient[kc_spec] = m_feSpecies_old[kc_spec]
                                                     + log(m_actCoeffSpecies_new[kc_spec] * X_est[kc]);
                } else {
                    m_feSpecies_Deficient[kc_spec] = m_feSpecies_old[kc_spec]
                                                     + log(m_actCoeffSpecies_new[kc_spec] * VCS_DELETE_MINORSPECIES_CUTOFF);
                }
            }

            for (i = 0; i < componentList.size(); i++) {
                kc = componentList[i];
                kc_spec = Vphase->spGlobalIndexVCS(kc);

                for (k = 0; k <  Vphase->nSpecies(); k++) {
                    kspec = Vphase->spGlobalIndexVCS(k);
                    if (kspec >= m_numComponents) {
                        irxn = kspec - m_numComponents;
                        if (i == 0) {
                            m_deltaGRxn_Deficient[irxn] = m_deltaGRxn_old[irxn];
                        }
                        double* dtmp_ptr = m_stoichCoeffRxnMatrix[irxn];
                        if (dtmp_ptr[kc_spec] != 0.0) {
                            m_deltaGRxn_Deficient[irxn] +=
                                dtmp_ptr[kc_spec] * (m_feSpecies_Deficient[kc_spec]- m_feSpecies_old[kc_spec]);
                        }
                    }

                }
            }

            /*
             *  Calculate the E_phi's
             */
            sum = 0.0;
            funcPhaseStability = sum_Xcomp - 1.0;
            for (k = 0; k < nsp; k++) {
                kspec = Vphase->spGlobalIndexVCS(k);
                if (kspec >= m_numComponents) {
                    irxn = kspec - m_numComponents;
                    deltaGRxn = m_deltaGRxn_Deficient[irxn];
                    if (deltaGRxn >  50.0) {
                        deltaGRxn =  50.0;
                    }
                    if (deltaGRxn < -50.0) {
                        deltaGRxn = -50.0;
                    }
                    E_phi[k] = std::exp(-deltaGRxn) / m_actCoeffSpecies_new[kspec];
                    sum +=  E_phi[k];
                    funcPhaseStability += E_phi[k];
                } else {
                    E_phi[k] = 0.0;
                }
            }

            /*
             * Calculate the raw estimate of the new fracs
             */
            for (k = 0; k < nsp; k++) {
                kspec = Vphase->spGlobalIndexVCS(k);
                double b =  E_phi[k] / sum * (1.0 - sum_Xcomp);
                if (kspec >= m_numComponents) {
                    irxn = kspec - m_numComponents;
                    fracDelta_raw[k] = b;
                }
            }


            // Given a set of fracDelta's, we calculate the fracDelta's
            // for the component species, if any
            for (i = 0; i < componentList.size(); i++) {
                kc = componentList[i];
                kc_spec = Vphase->spGlobalIndexVCS(kc);
                fracDelta_raw[kc] = 0.0;
                for (k = 0; k < nsp; k++) {
                    kspec = Vphase->spGlobalIndexVCS(k);
                    if (kspec >= m_numComponents) {
                        irxn = kspec - m_numComponents;
                        fracDelta_raw[kc] += m_stoichCoeffRxnMatrix[irxn][kc_spec] * fracDelta_raw[k];
                    }
                }
            }



            /*
             * Now possibly dampen the estimate.
             */
            doublereal sumADel = 0.0;
            for (k = 0; k < nsp; k++) {
                delFrac[k] = fracDelta_raw[k] - fracDelta_old[k];
                sumADel += fabs(delFrac[k]);
            }
            normUpdate = vcs_l2norm(delFrac);

            dirProd = 0.0;
            for (k = 0; k < nsp; k++) {
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

            for (k = 0; k < nsp; k++) {
                if (fabs(damp * delFrac[k]) > 0.3*fabs(fracDelta_old[k])) {
                    damp = std::max(0.3*fabs(fracDelta_old[k]) / fabs(delFrac[k]),
                                    1.0E-8/fabs(delFrac[k]));
                }
                if (delFrac[k] < 0.0) {
                    if (2.0 * damp * (-delFrac[k]) > fracDelta_old[k]) {
                        damp = fracDelta_old[k] / (2.0 * (-delFrac[k]));
                    }
                }
                if (delFrac[k] > 0.0) {
                    if (2.0 * damp * delFrac[k] > fracDelta_old[k]) {
                        damp = fracDelta_old[k] / (2.0 * delFrac[k]);
                    }
                }
            }
            if (damp < 0.000001) {
                damp = 0.000001;
            }

            for (k = 0; k < nsp; k++) {
                fracDelta_new[k] = fracDelta_old[k] + damp * (delFrac[k]);
            }

#ifdef DEBUG_MODE
            if (m_debug_print_lvl >= 2) {
                plogf("  --- %3d %12g %12g %12g %12g %12g %12g %12g\n", its, X_est[KP], fracDelta_old[KP],
                      delFrac[KP], fracDelta_new[KP], normUpdate, damp, funcPhaseStability);
            }
#endif

            if (normUpdate < 1.0E-5 * damp) {
                converged = true;
		if (its < minNumberIterations) {
		    converged = false;
		}
            }

        }

        if (converged) {
	    /*
	     *  Save the final optimized stated back into the VolPhase object for later use
	     */
            Vphase->setMoleFractionsState(0.0, VCS_DATA_PTR(X_est), VCS_STATECALC_PHASESTABILITY);
	    /*
	     * Save fracDelta for later use to initialize the problem better
	     *  @TODO  creationGlobalRxnNumbers needs to be calculated here and storred.
	     */
            Vphase->setCreationMoleNumbers(VCS_DATA_PTR(fracDelta_new), creationGlobalRxnNumbers);
        }


    } else {
        printf("not done yet\n");
        exit(-1);
    }
#ifdef DEBUG_MODE
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
#endif
    return funcPhaseStability;
}

}
