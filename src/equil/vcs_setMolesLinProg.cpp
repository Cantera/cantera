/*!
 * @file vcs_setMolesLinProg.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/equil/vcs_solve.h"

using namespace std;

namespace Cantera
{

static void printProgress(const vector<string> &spName,
                          const vector_fp &soln,
                          const vector_fp &ff)
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

int VCS_SOLVE::vcs_setMolesLinProg()
{
    double test = -1.0E-10;

    if (m_debug_print_lvl >= 2) {
        plogf("   --- call setInitialMoles\n");
    }

    double dxi_min = 1.0e10;
    int retn;
    int iter = 0;
    bool abundancesOK = true;
    bool usedZeroedSpecies;
    vector_fp sm(m_nelem * m_nelem, 0.0);
    vector_fp ss(m_nelem, 0.0);
    vector_fp sa(m_nelem, 0.0);
    vector_fp wx(m_nelem, 0.0);
    vector_fp aw(m_nsp, 0.0);

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
            retn = vcs_elcorr(&sm[0], &wx[0]);
            if (retn >= 2) {
                abundancesOK = false;
            } else {
                abundancesOK = true;
            }
        } else {
            abundancesOK = true;
        }

        // Now find the optimized basis that spans the stoichiometric
        // coefficient matrix, based on the current composition,
        // m_molNumSpecies_old[] We also calculate sc[][], the reaction matrix.
        retn = vcs_basopt(false, &aw[0], &sa[0], &sm[0], &ss[0],
                          test, &usedZeroedSpecies);
        if (retn != VCS_SUCCESS) {
            return retn;
        }

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
            dxi_min = 1.0e10;
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
    retn = 0;
    if (!abundancesOK) {
        retn = -1;
    } else if (iter > 15) {
        retn = 1;
    }
    return retn;
}

}
