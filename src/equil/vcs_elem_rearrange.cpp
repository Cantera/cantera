/**
 * @file vcs_elem_rearrange.cpp
 *   Contains implementations for rearranging the element columns, and
 *   it contains the algorithm for choosing the rearrangement.
 */
/*
 * Copyright (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#include "cantera/equil/vcs_solve.h"
#include "cantera/equil/vcs_internal.h"
#include "cantera/equil/vcs_VolPhase.h"

namespace VCSnonideal
{

int VCS_SOLVE::vcs_elem_rearrange(double* const aw, double* const sa,
                                  double* const sm, double* const ss)
{
    size_t  j, k, l, i, jl, ml, jr, ielem;
    bool lindep;
    size_t ncomponents = m_numComponents;
    double test = -1.0E10;
#ifdef DEBUG_MODE
    if (m_debug_print_lvl >= 2) {
        plogf("   ");
        for (i=0; i<77; i++) {
            plogf("-");
        }
        plogf("\n");
        plogf("   --- Subroutine elem_rearrange() called to ");
        plogf("check stoich. coefficient matrix\n");
        plogf("   ---    and to rearrange the element ordering once");
        plogendl();
    }
#endif

    /*
     *        Use a temporary work array for the element numbers
     *        Also make sure the value of test is unique.
     */
    lindep = false;
    do {
        lindep = false;
        for (i = 0; i < m_numElemConstraints; ++i) {
            test -= 1.0;
            aw[i] = m_elemAbundancesGoal[i];
            if (test == aw[i]) {
                lindep = true;
            }
        }
    } while (lindep);

    /*
     *        Top of a loop of some sort based on the index JR. JR is the
     *       current number independent elements found.
     */
    jr = npos;
    do {
        ++jr;
        /*
         *     Top of another loop point based on finding a linearly
         *     independent species
         */
        do {
            /*
             *    Search the remaining part of the mole fraction vector, AW,
             *    for the largest remaining species. Return its identity in K.
             */
            k = m_numElemConstraints;
            for (ielem = jr; ielem < m_numElemConstraints; ielem++) {
                if (m_elementActive[ielem]) {
                    if (aw[ielem] != test) {
                        k = ielem;
                        break;
                    }
                }
            }
            if (k == m_numElemConstraints) {
                plogf("vcs_elem_rearrange::Shouldn't be here. Algorithm misfired.");
                plogendl();
                exit(EXIT_FAILURE);
            }

            /*
             *  Assign a large negative number to the element that we have
             *  just found, in order to take it out of further consideration.
             */
            aw[k] = test;

            /* *********************************************************** */
            /* **** CHECK LINEAR INDEPENDENCE OF CURRENT FORMULA MATRIX    */
            /* **** LINE WITH PREVIOUS LINES OF THE FORMULA MATRIX  ****** */
            /* *********************************************************** */
            /*
             *          Modified Gram-Schmidt Method, p. 202 Dalquist
             *          QR factorization of a matrix without row pivoting.
             */
            jl = jr;
            /*
             *   Fill in the row for the current element, k, under consideration
             *   The row will contain the Formula matrix value for that element
             *   from the current component.
             */
            for (j = 0; j < ncomponents; ++j) {
                sm[j + jr*ncomponents] = m_formulaMatrix[k][j];
            }
            if (jl > 0) {
                /*
                 *         Compute the coefficients of JA column of the
                 *         the upper triangular R matrix, SS(J) = R_J_JR
                 *         (this is slightly different than Dalquist)
                 *         R_JA_JA = 1
                 */
                for (j = 0; j < jl; ++j) {
                    ss[j] = 0.0;
                    for (i = 0; i < ncomponents; ++i) {
                        ss[j] += sm[i + jr*ncomponents] * sm[i + j*ncomponents];
                    }
                    ss[j] /= sa[j];
                }
                /*
                 *     Now make the new column, (*,JR), orthogonal to the
                 *     previous columns
                 */
                for (j = 0; j < jl; ++j) {
                    for (l = 0; l < ncomponents; ++l) {
                        sm[l + jr*ncomponents] -= ss[j] * sm[l + j*ncomponents];
                    }
                }
            }

            /*
             *        Find the new length of the new column in Q.
             *        It will be used in the denominator in future row calcs.
             */
            sa[jr] = 0.0;
            for (ml = 0; ml < ncomponents; ++ml) {
                sa[jr] += SQUARE(sm[ml + jr*ncomponents]);
            }
            /* **************************************************** */
            /* **** IF NORM OF NEW ROW  .LT. 1E-6 REJECT ********** */
            /* **************************************************** */
            if (sa[jr] < 1.0e-6) {
                lindep = true;
            } else {
                lindep = false;
            }
        } while (lindep);
        /* ****************************************** */
        /* **** REARRANGE THE DATA ****************** */
        /* ****************************************** */
        if (jr != k) {
#ifdef DEBUG_MODE
            if (m_debug_print_lvl >= 2) {
                plogf("   ---   ");
                plogf("%-2.2s", (m_elementName[k]).c_str());
                plogf("(%9.2g) replaces ", m_elemAbundancesGoal[k]);
                plogf("%-2.2s", (m_elementName[jr]).c_str());
                plogf("(%9.2g) as element %3d", m_elemAbundancesGoal[jr], jr);
                plogendl();
            }
#endif
            vcs_switch_elem_pos(jr, k);
            std::swap(aw[jr], aw[k]);
        }

        /*
         *      If we haven't found enough components, go back
         *      and find some more. (nc -1 is used below, because
         *      jr is counted from 0, via the C convention.
         */
    } while (jr < (ncomponents-1));
    return VCS_SUCCESS;
}

void VCS_SOLVE::vcs_switch_elem_pos(size_t ipos, size_t jpos)
{
    if (ipos == jpos) {
        return;
    }
    size_t j;
    vcs_VolPhase* volPhase;
#ifdef DEBUG_MODE
    if (ipos > (m_numElemConstraints - 1) ||
            jpos > (m_numElemConstraints - 1)) {
        plogf("vcs_switch_elem_pos: ifunc = 0: inappropriate args: %d %d\n",
              ipos, jpos);
        plogendl();
        exit(EXIT_FAILURE);
    }
#endif
    /*
     * Change the element Global Index list in each vcs_VolPhase object
     * to reflect the switch in the element positions.
     */
    for (size_t iph = 0; iph < m_numPhases; iph++) {
        volPhase = m_VolPhaseList[iph];
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
    for (j = 0; j < m_numSpeciesTot; ++j) {
        std::swap(m_formulaMatrix[ipos][j], m_formulaMatrix[jpos][j]);
    }
    std::swap(m_elementName[ipos], m_elementName[jpos]);
}
}
