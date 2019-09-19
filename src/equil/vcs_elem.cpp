/**
 * @file vcs_elem.cpp
 *  This file contains the algorithm for checking the satisfaction of the
 *  element abundances constraints and the algorithm for fixing violations
 *  of the element abundances constraints.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/equil/vcs_solve.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/numerics/DenseMatrix.h"

namespace Cantera
{
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
    size_t top = m_numComponents;
    if (ibound) {
        top = m_nelem;
    }

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
                // For normal element balances, we require absolute compliance
                // even for ridiculously small numbers.
                if (m_elType[i] == VCS_ELEM_TYPE_ABSPOS) {
                    return false;
                } else {
                    return false;
                }
            }
        }
    }
    return true;
}

void VCS_SOLVE::vcs_elabPhase(size_t iphase, double* const elemAbundPhase)
{
    for (size_t j = 0; j < m_nelem; ++j) {
        elemAbundPhase[j] = 0.0;
        for (size_t i = 0; i < m_nsp; ++i) {
            if (m_speciesUnknownType[i] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE && m_phaseID[i] == iphase) {
                elemAbundPhase[j] += m_formulaMatrix(i,j) * m_molNumSpecies_old[i];
            }
        }
    }
}

int VCS_SOLVE::vcs_elcorr(double aa[], double x[])
{
    int retn = 0;

    vector_fp ga_save(m_elemAbundances);
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
        for (size_t i = 0; i < m_numComponents; ++i) {
            double tmp = m_molNumSpecies_old[i] + par * x[i];
            if (tmp > 0.0) {
                m_molNumSpecies_old[i] = tmp;
            } else {
                if (m_SSPhase[i]) {
                    m_molNumSpecies_old[i] = 0.0;
                }  else {
                    m_molNumSpecies_old[i] = m_molNumSpecies_old[i] * 0.0001;
                }
            }
        }
    } else {
        for (size_t i = 0; i < m_numComponents; ++i) {
            double tmp = m_molNumSpecies_old[i] + x[i];
            if (tmp > 0.0) {
                m_molNumSpecies_old[i] = tmp;
            } else {
                if (m_SSPhase[i]) {
                    m_molNumSpecies_old[i] = 0.0;
                }  else {
                    m_molNumSpecies_old[i] = m_molNumSpecies_old[i] * 0.0001;
                }
            }
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
                    } else {
                        if (saveDir > 0.0) {
                            goodSpec = false;
                            break;
                        }
                    }
                    saveDir = dir;
                } else {
                    if (m_formulaMatrix(kspec,i) != 0.) {
                        goodSpec = false;
                        break;
                    }
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

}
