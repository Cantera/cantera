/**
 * @file vcs_elem_rearrange.cpp
 *   Contains implementations for rearranging the element columns, and
 *   it contains the algorithm for choosing the rearrangement.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/equil/vcs_solve.h"
#include "cantera/equil/vcs_VolPhase.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ctexceptions.h"

namespace Cantera
{

int VCS_SOLVE::vcs_elem_rearrange(double* const aw, double* const sa,
                                  double* const sm, double* const ss)
{
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
    return VCS_SUCCESS;
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
}
