/*!
 * @file vcs_setMolesLinProg.cpp
 *
 */
/*
 * Copyright (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#include "cantera/equil/vcs_internal.h"
#include "cantera/equil/vcs_VolPhase.h"
#include "vcs_species_thermo.h"
#include "cantera/equil/vcs_solve.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>

#ifdef hpux
#define dbocls_ dbocls
#endif
#ifdef DEBUG_MODE
//extern int vcs_debug_print_lvl;
#endif
#ifndef  MAX
#define  MAX(x,y) (( (x) > (y) ) ? (x) : (y))
#endif


extern "C" void dbocls_(double* W, int* MDW, int* MCON, int* MROWS,
                        int* NCOLS,
                        double* BL, double* BU, int* IND, int* IOPT,
                        double* X, double* RNORMC, double* RNORM,
                        int* MODE, double* RW, int* IW);

using namespace std;

namespace VCSnonideal
{

#ifdef DEBUG_MODE
static void printProgress(const vector<string> &spName,
                          const vector<double> &soln,
                          const vector<double> &ff)
{
    int nsp = soln.size();
    double sum = 0.0;
    plogf(" --- Summary of current progress:\n");
    plogf(" ---                   Name           Moles  -       SSGibbs \n");
    plogf(" -------------------------------------------------------------------------------------\n");
    for (int k = 0; k < nsp; k++) {
        plogf(" ---      %20s %12.4g  - %12.4g\n", spName[k].c_str(), soln[k], ff[k]);
        sum += soln[k] * ff[k];
    }
    plogf(" ---  Total sum to be minimized = %g\n", sum);
}
#endif

#ifdef ALTLINPROG
//! Estimate the initial mole numbers.
/*!
 *   This is done by running
 *   each reaction as far forward or backward as possible, subject
 *   to the constraint that all mole numbers remain
 *   non-negative. Reactions for which \f$ \Delta \mu^0 \f$ are
 *   positive are run in reverse, and ones for which it is negative
 *   are run in the forward direction. The end result is equivalent
 *   to solving the linear programming problem of minimizing the
 *   linear Gibbs function subject to the element and
 *   non-negativity constraints.
 */
int VCS_SOLVE::vcs_setMolesLinProg()
{
    size_t ik, irxn;
    double test = -1.0E-10;

#ifdef DEBUG_MODE
    std::string pprefix(" --- seMolesLinProg ");
    if (m_debug_print_lvl >= 2) {
        plogf("   --- call setInitialMoles\n");
    }
#endif


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
            m_molNumSpecies_old[ik] = MAX(0.0, m_molNumSpecies_old[ik]);
        }
    }

#ifdef DEBUG_MODE
    if (m_debug_print_lvl >= 2) {
        printProgress(m_speciesName, m_molNumSpecies_old, m_SSfeSpecies);
    }
#endif

    while (redo) {

        if (!vcs_elabcheck(0)) {
#ifdef DEBUG_MODE
            if (m_debug_print_lvl >= 2) {
                plogf("%s Mole numbers failing element abundances\n", pprefix.c_str());
                plogf("%sCall vcs_elcorr to attempt fix\n",          pprefix.c_str());
            }
#endif
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

#ifdef DEBUG_MODE
        if (m_debug_print_lvl >= 2) {
            plogf("iteration %d\n", iter);
        }
#endif
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
            const double* sc_irxn = m_stoichCoeffRxnMatrix[irxn];
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
#ifdef DEBUG_MODE
                            if (m_debug_print_lvl >= 2) {
                                plogf("   --- Component too small: %s\n", m_speciesName[jcomp].c_str());
                            }
#endif
                            redo = true;
                        }
                    }
                    if (delta_xi < dxi_min) {
                        dxi_min = delta_xi;
                    }
                }
            }

            // step the composition by dxi_min, check against zero, since
            // we are zeroing components and species on every step.
            // Redo the iteration, if a component went from positive to zero on this step.
            double dsLocal = idir*dxi_min;
            m_molNumSpecies_old[ik] += dsLocal;
            m_molNumSpecies_old[ik] = MAX(0.0,  m_molNumSpecies_old[ik]);
            for (size_t jcomp = 0; jcomp < m_numComponents; jcomp++) {
                bool full = false;
                if (m_molNumSpecies_old[jcomp] > 1.0E-15) {
                    full = true;
                }
                m_molNumSpecies_old[jcomp] += sc_irxn[jcomp] * dsLocal;
                m_molNumSpecies_old[jcomp] = MAX(0.0, m_molNumSpecies_old[jcomp]);
                if (full) {
                    if (m_molNumSpecies_old[jcomp] < 1.0E-60) {
                        redo = true;
                    }
                }
            }
        }

        // set the moles of the phase objects to match
        //  updateMixMoles();
        // Update the phase objects with the contents of the m_molNumSpecies_old vector
        // vcs_updateVP(0);
#ifdef DEBUG_MODE
        if (m_debug_print_lvl >= 2) {
            printProgress(m_speciesName, m_molNumSpecies_old, m_SSfeSpecies);
        }
#endif
    }

#ifdef DEBUG_MODE
    if (m_debug_print_lvl == 1) {
        printProgress(m_speciesName, m_molNumSpecies_old, m_SSfeSpecies);
        plogf("   --- setInitialMoles end\n");
    }
#endif
    retn = 0;
    if (!abundancesOK) {
        retn = -1;
    } else if (iter > 15) {
        retn = 1;
    }
    return retn;
}

#else // ALTLINPROG

int linprogmax(double* XMOLES, double* CC, double* AX, double* BB,
               size_t NE, size_t M, size_t NE0)

/*-----------------------------------------------------------------------
* Find XMOLES(I), i = 1, M such that
*   Maximize CC dot W,  subject to the NE constraints:
*
*             [AX] [XMOLES] = [BB]
*             and XMOLES(i) > 0
*
* Input
* ---------
*   AX(NE, M) - matrix of constraints  AX(I,J) = ax(i + j*ne0)
*   BB(NE)    - contraint values
*   CC(M)     -  Vector of "Good Values" to maximize
*
* Output
* ---------
*   XMOLES(M) - optimal value of XMOLES()
*----------------------------------------------------------------------*/
{
    int MROWS, MCON, NCOLS, NX, NI, MDW, i, j, MODE;
    double sum, F[1], RNORMC, RNORM, *W, *BL, *BU, *RW, *X;
    int*    IND, *IW, *IOPT;

    MROWS = 1;
    MCON = (int) NE;
    NCOLS = (int) M;
    MDW = MCON + NCOLS;
    NX = 0;
    NI = 0;

    sum = 0.0;
    for (i = 0; i < NCOLS; i++) {
        sum += fabs(CC[i]);
    }
    F[0] = sum * 1000.;
    if (F[0] <= 0.0) {
        F[0] = 1000.;
    }

    BL = (double*) malloc(2*(NCOLS+MCON)            * sizeof(double));
    BU = BL + (NCOLS+MCON);
    IND = (int*) malloc((NCOLS+MCON)              * sizeof(int));
    RW = (double*) malloc((6*NCOLS + 5*MCON)        * sizeof(double));
    IW = (int*) malloc((2*NCOLS + 2*MCON)        * sizeof(int));
    IOPT = (int*) malloc((17 + NI)                 * sizeof(int));
    X = (double*) malloc((2*(NCOLS+MCON) + 2 + NX) * sizeof(double));
    W = (double*) malloc((MDW*(NCOLS+MCON+1))      * sizeof(double));
    if (W == NULL) {
        plogf("linproxmax ERROR: can not malloc memory of size %d bytes\n",
              (int)((MDW*(NCOLS+MCON+1)) * sizeof(double)));
        if (BL != NULL) {
            free((void*) BL);
        }
        if (IND != NULL) {
            free((void*) IND);
        }
        if (RW != NULL) {
            free((void*) RW);
        }
        if (IW != NULL) {
            free((void*) IW);
        }
        if (IOPT != NULL) {
            free((void*) IOPT);
        }
        if (W != NULL) {
            free((void*) W);
        }
        return -1;
    }
    for (j = 0; j < MCON; j++) {
        for (i = 0; i < NCOLS; i++) {
            W[j + i*MDW] = AX[j + i*NE0];
        }
    }
    for (i = 0; i < NCOLS; i++) {
        W[MCON + i*MDW] = CC[i];
    }
    W[MCON + (NCOLS)*MDW] = F[0];
    IOPT[0] = 99;

    for (j = 0; j < NCOLS; j++) {
        IND[j] = 1;
        BL[j]  = 0.0;
        BU[j]  = 1.0e200;
    }
    for (j = 0; j < MCON; j++) {
        IND[j + NCOLS] = 3;
        BL[j + NCOLS]  = BB[j];
        BU[j + NCOLS]  = BL[j + NCOLS];
    }


    dbocls_(W, &MDW, &MCON, &MROWS, &NCOLS, BL, BU, IND, IOPT,
            X, &RNORMC, &RNORM, &MODE, RW, IW);
    if (MODE != 0) {
        plogf("Return from DBOCLS was not normal, MODE = %d\n", MODE);
        plogf("       refer to subroutine DBOCLS for resolution\n");
        plogf("       RNORMC = %g\n", RNORMC);
    }

    for (j = 0; j < NCOLS; j++) {
        XMOLES[j] = X[j];
    }
#ifdef DEBUG_MODE
    //sum = 0.0;
    //for (j = 0; j < NCOLS; j++) {
    //   sum += XMOLES[j] * CC[j];
    //}
    //if (vcs_debug_print_lvl >= 2) {
    //  plogf(" -- linmaxc: Final Maximized Value = %g\n", sum);
    //}
#endif

    free((void*)W);
    free((void*)BL);
    free((void*)IND);
    free((void*)RW);
    free((void*)IW);
    free((void*)IOPT);
    free((void*)X);

    return 0;
}
#endif // ALTLINPROG
}
