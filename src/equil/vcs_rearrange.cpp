/**
 * @file vcs_rearrange.cpp
 *    implementation file for rearranging species.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/equil/vcs_solve.h"

namespace Cantera
{
int VCS_SOLVE::vcs_rearrange()
{
    warn_deprecated("VCS_SOLVE::vcs_rearrange",
                    "Unused. To be removed after Cantera 2.3.");
    size_t k1 = 0;

    // Loop over all of the species
    for (size_t i = 0; i < m_numSpeciesTot; ++i) {
        // Find the index of I in the index vector m_speciesIndexVector[]. Call
        // it k1 and continue.
        for (size_t j = 0; j < m_numSpeciesTot; ++j) {
            size_t n = m_speciesMapIndex[j];
            k1 = j;
            if (n == i) {
                break;
            }
        }

        //Switch the species data back from k1 into i. because we loop over all
        //species, reaction data are now permanently hosed.
        vcs_switch_pos(false, i, k1);
    }
    return 0;
}

}
