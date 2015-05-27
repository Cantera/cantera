/*!
 * @file vcs_setMolesLinProg.cpp
 */
/*
 * Copyright (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#include "cantera/equil/vcs_solve.h"

using namespace std;

namespace Cantera
{

static void printProgress(const vector<string> &spName,
                          const vector<double> &soln,
                          const vector<double> &ff)
{
    double sum = 0.0;
    plogf(" --- Summary of current progress:\n");
    plogf(" ---                   Name           Moles  -       SSGibbs \n");
    plogf(" -------------------------------------------------------------------------------------\n");
    for (size_t k = 0; k < soln.size(); k++) {
        plogf(" ---      %20s %12.4g  - %12.4g\n", spName[k].c_str(), soln[k], ff[k]);
        sum += soln[k] * ff[k];
    }
    plogf(" ---  Total sum to be minimized = %g\n", sum);
}

int VCS_SOLVE::vcs_setMolesLinProg()
{
    size_t ik, irxn;
    double test = -1.0E-10;

    if (DEBUG_MODE_ENABLED && m_debug_print_lvl >= 2) {
        plogf("   --- call setInitialMoles\n");
    }

    // m_mu are standard state chemical potentials
    //  Boolean on the end specifies standard chem potentials
    // m_mix->getValidChemPotentials(not_mu, DATA_PTR(m_mu), true);
    // -> This is already done coming into the routine.
    double dg_rt;

    int idir;
    double nu;
    double delta_xi, dxi_min = 1.0e10;
    bool redo = true;
    int retn;
    int iter = 0;
    bool abundancesOK = true;
    bool usedZeroedSpecies;

    std::vector<double> sm(m_numElemConstraints*m_numElemConstraints, 0.0);
    std::vector<double> ss(m_numElemConstraints, 0.0);
    std::vector<double> sa(m_numElemConstraints, 0.0);
    std::vector<double> wx(m_numElemConstraints, 0.0);
    std::vector<double> aw(m_numSpeciesTot, 0.0);

    for (ik = 0; ik <  m_numSpeciesTot; ik++) {
        if (m_speciesUnknownType[ik] !=  VCS_SPECIES_INTERFACIALVOLTAGE) {
            m_molNumSpecies_old[ik] = max(0.0, m_molNumSpecies_old[ik]);
        }
    }

    if (DEBUG_MODE_ENABLED && m_debug_print_lvl >= 2) {
        printProgress(m_speciesName, m_molNumSpecies_old, m_SSfeSpecies);
    }

    while (redo) {

        if (!vcs_elabcheck(0)) {
            if (DEBUG_MODE_ENABLED && m_debug_print_lvl >= 2) {
                plogf(" --- seMolesLinProg  Mole numbers failing element abundances\n");
                plogf(" --- seMolesLinProg  Call vcs_elcorr to attempt fix\n");
            }
            retn = vcs_elcorr(VCS_DATA_PTR(sm), VCS_DATA_PTR(wx));
            if (retn >= 2) {
                abundancesOK = false;
            } else {
                abundancesOK = true;
            }
        } else {
            abundancesOK = true;
        }
        /*
         *  Now find the optimized basis that spans the stoichiometric
         *  coefficient matrix, based on the current composition, m_molNumSpecies_old[]
         *  We also calculate sc[][], the reaction matrix.
         */
        retn = vcs_basopt(false, VCS_DATA_PTR(aw), VCS_DATA_PTR(sa),
                          VCS_DATA_PTR(sm), VCS_DATA_PTR(ss),
                          test, &usedZeroedSpecies);
        if (retn != VCS_SUCCESS) {
            return retn;
        }

        if (DEBUG_MODE_ENABLED && m_debug_print_lvl >= 2) {
            plogf("iteration %d\n", iter);
        }
        redo = false;
        iter++;
        if (iter > 15) {
            break;
        }

        // loop over all reactions
        for (irxn = 0; irxn < m_numRxnTot; irxn++) {

            // dg_rt is the Delta_G / RT value for the reaction
            ik = m_numComponents + irxn;
            dg_rt = m_SSfeSpecies[ik];
            dxi_min = 1.0e10;
            const double* sc_irxn = m_stoichCoeffRxnMatrix.ptrColumn(irxn);
            for (size_t jcomp = 0; jcomp < m_numElemConstraints; jcomp++) {
                dg_rt += m_SSfeSpecies[jcomp] * sc_irxn[jcomp];
            }
            // fwd or rev direction.
            //  idir > 0 implies increasing the current species
            //  idir < 0  implies decreasing the current species
            idir = (dg_rt < 0.0 ? 1 : -1);
            if (idir < 0) {
                dxi_min = m_molNumSpecies_old[ik];
            }

            for (size_t jcomp = 0; jcomp < m_numComponents; jcomp++) {
                nu = sc_irxn[jcomp];

                // set max change in progress variable by
                // non-negativity requirement
                if (nu*idir < 0) {
                    delta_xi = fabs(m_molNumSpecies_old[jcomp]/nu);
                    // if a component has nearly zero moles, redo
                    // with a new set of components
                    if (!redo) {
                        if (delta_xi < 1.0e-10 && (m_molNumSpecies_old[ik] >= 1.0E-10)) {
                            if (DEBUG_MODE_ENABLED && m_debug_print_lvl >= 2) {
                                plogf("   --- Component too small: %s\n", m_speciesName[jcomp].c_str());
                            }
                            redo = true;
                        }
                    }
                    dxi_min = std::min(dxi_min, delta_xi);
                }
            }

            // step the composition by dxi_min, check against zero, since
            // we are zeroing components and species on every step.
            // Redo the iteration, if a component went from positive to zero on this step.
            double dsLocal = idir*dxi_min;
            m_molNumSpecies_old[ik] += dsLocal;
            m_molNumSpecies_old[ik] = max(0.0,  m_molNumSpecies_old[ik]);
            for (size_t jcomp = 0; jcomp < m_numComponents; jcomp++) {
                bool full = false;
                if (m_molNumSpecies_old[jcomp] > 1.0E-15) {
                    full = true;
                }
                m_molNumSpecies_old[jcomp] += sc_irxn[jcomp] * dsLocal;
                m_molNumSpecies_old[jcomp] = max(0.0, m_molNumSpecies_old[jcomp]);
                if (full) {
                    if (m_molNumSpecies_old[jcomp] < 1.0E-60) {
                        redo = true;
                    }
                }
            }
        }

        if (DEBUG_MODE_ENABLED && m_debug_print_lvl >= 2) {
            printProgress(m_speciesName, m_molNumSpecies_old, m_SSfeSpecies);
        }
    }

    if (DEBUG_MODE_ENABLED && m_debug_print_lvl == 1) {
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
