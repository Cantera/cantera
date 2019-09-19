/**
 * @file vcs_phaseStability.cpp
 *  Implementation class for functions associated with determining the stability of a phase
 *   (see Class \link Cantera::VCS_SOLVE VCS_SOLVE\endlink and \ref equilfunctions ).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/equil/vcs_solve.h"
#include "cantera/equil/vcs_VolPhase.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ctexceptions.h"

using namespace std;

namespace Cantera
{

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

size_t VCS_SOLVE::vcs_popPhaseID(std::vector<size_t> & phasePopPhaseIDs)
{
    size_t iphasePop = npos;
    doublereal FephaseMax = -1.0E30;
    doublereal Fephase = -1.0E30;

    char anote[128];
    if (m_debug_print_lvl >= 2) {
        plogf("   --- vcs_popPhaseID() called\n");
        plogf("   ---   Phase                 Status       F_e        MoleNum\n");
        plogf("   --------------------------------------------------------------------------\n");
    }
    for (size_t iph = 0; iph < m_numPhases; iph++) {
        vcs_VolPhase* Vphase = m_VolPhaseList[iph].get();
        int existence = Vphase->exists();
        strcpy(anote, "");
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
                doublereal deltaGRxn = m_deltaGRxn_old[irxn];
                Fephase = exp(-deltaGRxn) - 1.0;
                if (Fephase > 0.0) {
                    strcpy(anote," (ready to be birthed)");
                    if (Fephase > FephaseMax) {
                        iphasePop = iph;
                        FephaseMax = Fephase;
                        strcpy(anote," (chosen to be birthed)");
                    }
                }
                if (Fephase < 0.0) {
                    strcpy(anote," (not stable)");
                    AssertThrowMsg(m_tPhaseMoles_old[iph] <= 0.0,
                        "VCS_SOLVE::vcs_popPhaseID", "shouldn't be here");
                }

                if (m_debug_print_lvl >= 2) {
                    plogf("  ---    %18s %5d %10.3g %10.3g %s\n",
                          Vphase->PhaseName, existence, Fephase,
                          m_tPhaseMoles_old[iph], anote);
                }
            } else {
                // MultiSpecies Phase Stability Resolution
                if (vcs_popPhasePossible(iph)) {
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
                } else {
                    if (m_debug_print_lvl >= 2) {
                        plogf("  ---    %18s %5d   blocked  %11.3g\n",
                              Vphase->PhaseName,
                              existence, m_tPhaseMoles_old[iph]);
                    }
                }
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
    std::vector<size_t> creationGlobalRxnNumbers;

    // Calculate the initial moles of the phase being born.
    //   Here we set it to 10x of the value which would cause the phase to be
    //   zeroed out within the algorithm.  We may later adjust the value.
    doublereal tPhaseMoles = 10. * m_totalMolNum * VCS_DELETE_PHASE_CUTOFF;

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
            double s_old = s;
            s = vcs_Hessian_diag_adj(irxn, s_old);
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
        vector_fp fracDelta(Vphase->nSpecies());
        vector_fp X_est(Vphase->nSpecies());
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

        doublereal ratioComp = 0.0;
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
    vector_fp X_est(nsp, 0.0);
    vector_fp delFrac(nsp, 0.0);
    vector_fp E_phi(nsp, 0.0);
    vector_fp fracDelta_old(nsp, 0.0);
    vector_fp fracDelta_raw(nsp, 0.0);
    vector<size_t> creationGlobalRxnNumbers(nsp, npos);
    m_deltaGRxn_Deficient = m_deltaGRxn_old;
    vector_fp feSpecies_Deficient = m_feSpecies_old;

    // get the activity coefficients
    Vphase->sendToVCS_ActCoeff(VCS_STATECALC_OLD, &m_actCoeffSpecies_new[0]);

    // Get the stored estimate for the composition of the phase if
    // it gets created
    vector_fp fracDelta_new = Vphase->creationMoleNumbers(creationGlobalRxnNumbers);

    std::vector<size_t> componentList;
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
                if (X_est[kc] > VCS_DELETE_MINORSPECIES_CUTOFF) {
                    feSpecies_Deficient[kc_spec] = m_feSpecies_old[kc_spec]
                                                     + log(m_actCoeffSpecies_new[kc_spec] * X_est[kc]);
                } else {
                    feSpecies_Deficient[kc_spec] = m_feSpecies_old[kc_spec]
                                                     + log(m_actCoeffSpecies_new[kc_spec] * VCS_DELETE_MINORSPECIES_CUTOFF);
                }
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
            doublereal sumADel = 0.0;
            for (size_t k = 0; k < nsp; k++) {
                delFrac[k] = fracDelta_raw[k] - fracDelta_old[k];
                sumADel += fabs(delFrac[k]);
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
            // @TODO  creationGlobalRxnNumbers needs to be calculated here and stored.
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

}
