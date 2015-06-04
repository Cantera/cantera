/**
 * @file vcs_rearrange.cpp
 *    implementation file for rearranging species.
 */
/*
 * Copyright (2007) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#include "cantera/equil/vcs_solve.h"

namespace Cantera
{
int VCS_SOLVE::vcs_rearrange()
{
    size_t k1 = 0;

    /* - Loop over all of the species */
    for (size_t i = 0; i < m_numSpeciesTot; ++i) {
        /*
         *         Find the index of I in the index vector m_speciesIndexVector[].
         *         Call it k1 and continue.
         */
        for (size_t j = 0; j < m_numSpeciesTot; ++j) {
            size_t l = m_speciesMapIndex[j];
            k1 = j;
            if (l == i) {
                break;
            }
        }
        /*
         * - Switch the species data back from k1 into i
         *     -> because we loop over all species, reaction data
         *        are now permanently hosed.
         */
        vcs_switch_pos(false, i, k1);
    }
    return 0;
}

}
