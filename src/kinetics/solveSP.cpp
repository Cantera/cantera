/*
 * @file: solveSP.cpp Implicit surface site concentration solver
 */
/*
 * Copyright 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "cantera/kinetics/solveSP.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/kinetics/ImplicitSurfChem.h"

#include <cstdio>

using namespace std;
namespace Cantera
{

/***************************************************************************
 *       STATIC ROUTINES DEFINED IN THIS FILE
 ***************************************************************************/

static doublereal calc_damping(doublereal* x, doublereal* dx, size_t dim, int*);
static doublereal calcWeightedNorm(const doublereal [], const doublereal dx[], size_t);

/***************************************************************************
 *  solveSP Class Definitions
 ***************************************************************************/

solveSP::solveSP(ImplicitSurfChem* surfChemPtr, int bulkFunc) :
    m_SurfChemPtr(surfChemPtr),
    m_objects(surfChemPtr->getObjects()),
    m_neq(0),
    m_bulkFunc(bulkFunc),
    m_numSurfPhases(0),
    m_numTotSurfSpecies(0),
    m_numBulkPhasesSS(0),
    m_numTotBulkSpeciesSS(0),
    m_atol(1.0E-15),
    m_rtol(1.0E-4),
    m_maxstep(1000),
    m_maxTotSpecies(0),
    m_ioflag(0)
{
    m_numSurfPhases = 0;
    size_t numPossibleSurfPhases = m_objects.size();
    for (size_t n = 0; n < numPossibleSurfPhases; n++) {
        InterfaceKinetics* m_kin = m_objects[n];
        size_t surfPhaseIndex = m_kin->surfacePhaseIndex();
        if (surfPhaseIndex != npos) {
            m_numSurfPhases++;
            m_indexKinObjSurfPhase.push_back(n);
            m_kinObjPhaseIDSurfPhase.push_back(surfPhaseIndex);
        } else {
            throw CanteraError("solveSP",
                               "InterfaceKinetics object has no surface phase");
        }
        ThermoPhase* tp = &(m_kin->thermo(surfPhaseIndex));
        SurfPhase* sp = dynamic_cast<SurfPhase*>(tp);
        if (!sp) {
            throw CanteraError("solveSP",
                               "Inconsistent ThermoPhase object within "
                               "InterfaceKinetics object");
        }

        m_ptrsSurfPhase.push_back(sp);
        size_t nsp = sp->nSpecies();
        m_nSpeciesSurfPhase.push_back(nsp);
        m_numTotSurfSpecies += nsp;

    }
    /*
     * We rely on ordering to figure things out
     */
    m_numBulkPhasesSS = 0;

    if (bulkFunc == BULK_DEPOSITION) {
        m_neq = m_numTotSurfSpecies + m_numTotBulkSpeciesSS;
    } else {
        m_neq = m_numTotSurfSpecies;
    }

    m_maxTotSpecies = 0;
    for (size_t n = 0; n < m_numSurfPhases; n++) {
        size_t tsp =  m_objects[n]->nTotalSpecies();
        m_maxTotSpecies = std::max(m_maxTotSpecies, tsp);
    }
    m_maxTotSpecies = std::max(m_maxTotSpecies, m_neq);


    m_netProductionRatesSave.resize(m_maxTotSpecies, 0.0);
    m_numEqn1.resize(m_maxTotSpecies, 0.0);
    m_numEqn2.resize(m_maxTotSpecies, 0.0);
    m_XMolKinSpecies.resize(m_maxTotSpecies, 0.0);
    m_CSolnSave.resize(m_neq, 0.0);

    m_spSurfLarge.resize(m_numSurfPhases, 0);

    m_kinSpecIndex.resize(m_numTotSurfSpecies + m_numTotBulkSpeciesSS, 0);
    m_kinObjIndex.resize(m_numTotSurfSpecies + m_numTotBulkSpeciesSS, 0);
    m_eqnIndexStartSolnPhase.resize(m_numSurfPhases + m_numBulkPhasesSS, 0);

    size_t kindexSP = 0;
    size_t isp, k, nsp, kstart;
    for (isp = 0; isp < m_numSurfPhases; isp++) {
        size_t iKinObject = m_indexKinObjSurfPhase[isp];
        InterfaceKinetics* m_kin = m_objects[iKinObject];
        size_t surfPhaseIndex = m_kinObjPhaseIDSurfPhase[isp];
        kstart = m_kin->kineticsSpeciesIndex(0, surfPhaseIndex);
        nsp = m_nSpeciesSurfPhase[isp];
        m_eqnIndexStartSolnPhase[isp] = kindexSP;
        for (k = 0; k < nsp; k++, kindexSP++) {
            m_kinSpecIndex[kindexSP] = kstart + k;
            m_kinObjIndex[kindexSP] = isp;
        }
    }

    // Dimension solution vector
    size_t dim1 = std::max<size_t>(1, m_neq);
    m_CSolnSP.resize(dim1, 0.0);
    m_CSolnSPInit.resize(dim1, 0.0);
    m_CSolnSPOld.resize(dim1, 0.0);
    m_wtResid.resize(dim1, 0.0);
    m_wtSpecies.resize(dim1, 0.0);
    m_resid.resize(dim1, 0.0);
    m_Jac.resize(dim1, dim1, 0.0);
}

int solveSP::solveSurfProb(int ifunc, doublereal time_scale, doublereal TKelvin,
                           doublereal PGas, doublereal reltol, doublereal abstol)
{
    doublereal EXTRA_ACCURACY = 0.001;
    if (ifunc == SFLUX_JACOBIAN) {
        EXTRA_ACCURACY *= 0.001;
    }
    int info = 0;
    int label_t=-1; /* Species IDs for time control */
    int label_d = -1; /* Species IDs for damping control */
    int        label_t_old=-1;
    doublereal     label_factor = 1.0;
    int iter=0; // iteration number on numlinear solver
    int iter_max=1000; // maximum number of nonlinear iterations
    doublereal deltaT = 1.0E-10; // Delta time step
    doublereal damp=1.0, tmp;
    //  Weighted L2 norm of the residual.  Currently, this is only
    //  used for IO purposes. It doesn't control convergence.
    doublereal  resid_norm;
    doublereal inv_t = 0.0;
    doublereal t_real = 0.0, update_norm = 1.0E6;

    bool do_time = false, not_converged = true;
    m_ioflag = std::min(m_ioflag, 1);

    /*
     *       Set the initial value of the do_time parameter
     */
    if (ifunc == SFLUX_INITIALIZE || ifunc == SFLUX_TRANSIENT) {
        do_time = true;
    }

    /*
     *    Store the initial guess for the surface problem in the soln vector,
     *  CSoln, and in an separate vector CSolnInit.
     */
    size_t loc = 0;
    for (size_t n = 0; n < m_numSurfPhases; n++) {
        SurfPhase* sf_ptr =  m_ptrsSurfPhase[n];
        sf_ptr->getConcentrations(DATA_PTR(m_numEqn1));
        size_t nsp = m_nSpeciesSurfPhase[n];
        for (size_t k = 0; k <nsp; k++) {
            m_CSolnSP[loc] = m_numEqn1[k];
            loc++;
        }
    }

    std::copy(m_CSolnSP.begin(), m_CSolnSP.end(), m_CSolnSPInit.begin());

    // Calculate the largest species in each phase
    evalSurfLarge(DATA_PTR(m_CSolnSP));

    if (m_ioflag) {
        print_header(m_ioflag, ifunc, time_scale, true, reltol, abstol);
    }

    /*
     *  Quick return when there isn't a surface problem to solve
     */
    if (m_neq == 0) {
        not_converged = false;
        update_norm = 0.0;
    }

    /* ------------------------------------------------------------------
     *                         Start of Newton's method
     * ------------------------------------------------------------------
     */
    while (not_converged && iter < iter_max) {
        iter++;
        /*
         *    Store previous iteration's solution in the old solution vector
         */
        std::copy(m_CSolnSP.begin(), m_CSolnSP.end(), m_CSolnSPOld.begin());

        /*
         * Evaluate the largest surface species for each surface phase every
         * 5 iterations.
         */
        if (iter%5 == 4) {
            evalSurfLarge(DATA_PTR(m_CSolnSP));
        }

        /*
         *    Calculate the value of the time step
         *       - heuristics to stop large oscillations in deltaT
         */
        if (do_time) {
            /* don't hurry increase in time step at the same time as damping */
            if (damp < 1.0) {
                label_factor = 1.0;
            }
            tmp = calc_t(DATA_PTR(m_netProductionRatesSave),
                         DATA_PTR(m_XMolKinSpecies),
                         &label_t, &label_t_old,  &label_factor, m_ioflag);
            if (iter < 10) {
                inv_t = tmp;
            } else if (tmp > 2.0*inv_t) {
                inv_t =  2.0*inv_t;
            } else {
                inv_t = tmp;
            }

            /*
             *   Check end condition
             */

            if (ifunc == SFLUX_TRANSIENT) {
                tmp = t_real + 1.0/inv_t;
                if (tmp > time_scale) {
                    inv_t = 1.0/(time_scale - t_real);
                }
            }
        } else {
            /* make steady state calc a step of 1 million seconds to
               prevent singular Jacobians for some pathological cases */
            inv_t = 1.0e-6;
        }
        deltaT = 1.0/inv_t;

        /*
         * Call the routine to numerically evaluation the Jacobian
         * and residual for the current iteration.
         */
        resjac_eval(m_Jac, DATA_PTR(m_resid), DATA_PTR(m_CSolnSP),
                    DATA_PTR(m_CSolnSPOld), do_time, deltaT);

        /*
         * Calculate the weights. Make sure the calculation is carried
         * out on the first iteration.
         */
        if (iter%4 == 1) {
            calcWeights(DATA_PTR(m_wtSpecies), DATA_PTR(m_wtResid),
                        m_Jac, DATA_PTR(m_CSolnSP), abstol, reltol);
        }

        /*
         *    Find the weighted norm of the residual
         */
        resid_norm = calcWeightedNorm(DATA_PTR(m_wtResid),
                                      DATA_PTR(m_resid), m_neq);

        /*
         *  Solve Linear system.  The solution is in resid[]
         */
        info = m_Jac.factor();
        if (info==0) {
            m_Jac.solve(&m_resid[0]);
        }
        /*
         *    Force convergence if residual is small to avoid
         *    "nan" results from the linear solve.
         */
        else {
            if (m_ioflag) {
                printf("solveSurfSS: Zero pivot, assuming converged: %g (%d)\n",
                       resid_norm, info);
            }
            for (size_t jcol = 0; jcol < m_neq; jcol++) {
                m_resid[jcol] = 0.0;
            }

            /* print out some helpful info */
            if (m_ioflag > 1) {
                printf("-----\n");
                printf("solveSurfProb: iter %d t_real %g delta_t %g\n\n",
                       iter,t_real, 1.0/inv_t);
                printf("solveSurfProb: init guess, current concentration,"
                       "and prod rate:\n");
                for (size_t jcol = 0; jcol < m_neq; jcol++) {
                    printf("\t%s  %g %g %g\n", int2str(jcol).c_str(),
                           m_CSolnSPInit[jcol], m_CSolnSP[jcol],
                           m_netProductionRatesSave[m_kinSpecIndex[jcol]]);
                }
                printf("-----\n");
            }
            if (do_time) {
                t_real += time_scale;
            }
        }

        /*
         *    Calculate the Damping factor needed to keep all unknowns
         *    between 0 and 1, and not allow too large a change (factor of 2)
         *    in any unknown.
         */

        damp = calc_damping(DATA_PTR(m_CSolnSP), DATA_PTR(m_resid), m_neq, &label_d);

        /*
         *    Calculate the weighted norm of the update vector
         *       Here, resid is the delta of the solution, in concentration
         *       units.
         */
        update_norm = calcWeightedNorm(DATA_PTR(m_wtSpecies),
                                       DATA_PTR(m_resid), m_neq);
        /*
         *    Update the solution vector and real time
         *    Crop the concentrations to zero.
         */
        for (size_t irow = 0; irow < m_neq; irow++) {
            m_CSolnSP[irow] -= damp * m_resid[irow];
        }
        for (size_t irow = 0; irow < m_neq; irow++) {
            m_CSolnSP[irow] = std::max(0.0, m_CSolnSP[irow]);
        }
        updateState(DATA_PTR(m_CSolnSP));

        if (do_time) {
            t_real += damp/inv_t;
        }

        if (m_ioflag) {
            printIteration(m_ioflag, damp, label_d, label_t,  inv_t, t_real, iter,
                           update_norm, resid_norm, do_time);
        }

        if (ifunc == SFLUX_TRANSIENT) {
            not_converged = (t_real < time_scale);
        } else {
            if (do_time) {
                if (t_real > time_scale ||
                        (resid_norm < 1.0e-7 &&
                         update_norm*time_scale/t_real < EXTRA_ACCURACY))  {
                    do_time = false;
                }
            } else {
                not_converged = ((update_norm > EXTRA_ACCURACY) ||
                                 (resid_norm  > EXTRA_ACCURACY));
            }
        }
    }  /* End of Newton's Method while statement */

    /*
     *  End Newton's method.  If not converged, print error message and
     *  recalculate sdot's at equal site fractions.
     */
    if (not_converged) {
        if (m_ioflag) {
            printf("#$#$#$# Error in solveSP $#$#$#$ \n");
            printf("Newton iter on surface species did not converge, "
                   "update_norm = %e \n", update_norm);
            printf("Continuing anyway\n");
        }
    }

    /*
     *        Decide on what to return in the solution vector
     *        - right now, will always return the last solution
     *          no matter how bad
     */
    if (m_ioflag) {
        fun_eval(DATA_PTR(m_resid), DATA_PTR(m_CSolnSP), DATA_PTR(m_CSolnSPOld),
                 false, deltaT);
        resid_norm = calcWeightedNorm(DATA_PTR(m_wtResid),
                                      DATA_PTR(m_resid), m_neq);
        printIteration(m_ioflag, damp, label_d, label_t,  inv_t, t_real, iter,
                       update_norm, resid_norm, do_time, true);
    }

    /*
     *        Return with the appropriate flag
     */
    if (update_norm > 1.0) {
        return -1;
    }
    return 1;
}

void solveSP::updateState(const doublereal* CSolnSP)
{
    size_t loc = 0;
    for (size_t n = 0; n < m_numSurfPhases; n++) {
        m_ptrsSurfPhase[n]->setConcentrations(CSolnSP + loc);
        loc += m_nSpeciesSurfPhase[n];
    }
}

void solveSP::updateMFSolnSP(doublereal* XMolSolnSP)
{
    for (size_t isp = 0; isp < m_numSurfPhases; isp++) {
        size_t keqnStart = m_eqnIndexStartSolnPhase[isp];
        m_ptrsSurfPhase[isp]->getMoleFractions(XMolSolnSP + keqnStart);
    }
}

void solveSP::updateMFKinSpecies(doublereal* XMolKinSpecies, int isp)
{
    InterfaceKinetics* m_kin = m_objects[isp];
    size_t nph = m_kin->nPhases();
    for (size_t iph = 0; iph < nph; iph++) {
        size_t ksi = m_kin->kineticsSpeciesIndex(0, iph);
        ThermoPhase& thref = m_kin->thermo(iph);
        thref.getMoleFractions(XMolKinSpecies + ksi);
    }
}

void solveSP::evalSurfLarge(const doublereal* CSolnSP)
{
    size_t kindexSP = 0;
    for (size_t isp = 0; isp < m_numSurfPhases; isp++) {
        size_t nsp = m_nSpeciesSurfPhase[isp];
        doublereal Clarge = CSolnSP[kindexSP];
        m_spSurfLarge[isp] = 0;
        kindexSP++;
        for (size_t k = 1; k < nsp; k++, kindexSP++) {
            if (CSolnSP[kindexSP] > Clarge) {
                Clarge = CSolnSP[kindexSP];
                m_spSurfLarge[isp] = k;
            }
        }
    }
}

void solveSP::fun_eval(doublereal* resid, const doublereal* CSoln,
                       const doublereal* CSolnOld,  const bool do_time,
                       const doublereal deltaT)
{
    size_t isp, nsp, kstart, k, kindexSP, kins, kspecial;
    doublereal lenScale = 1.0E-9;
    doublereal sd = 0.0;
    doublereal grRate;
    if (m_numSurfPhases > 0) {
        /*
         * update the surface concentrations with the input surface
         * concentration vector
         */
        updateState(CSoln);
        /*
         * Get the net production rates of all of the species in the
         * surface kinetics mechanism
         *
         * HKM Should do it here for all kinetics objects so that
         *     bulk will eventually work.
         */

        if (do_time) {
            kindexSP = 0;
            for (isp = 0; isp < m_numSurfPhases; isp++) {
                nsp = m_nSpeciesSurfPhase[isp];
                InterfaceKinetics* kinPtr = m_objects[isp];
                size_t surfIndex = kinPtr->surfacePhaseIndex();
                kstart = kinPtr->kineticsSpeciesIndex(0, surfIndex);
                kins = kindexSP;
                kinPtr->getNetProductionRates(DATA_PTR(m_netProductionRatesSave));
                for (k = 0; k < nsp; k++, kindexSP++) {
                    resid[kindexSP] =
                        (CSoln[kindexSP] - CSolnOld[kindexSP]) / deltaT
                        - m_netProductionRatesSave[kstart + k];
                }

                kspecial = kins + m_spSurfLarge[isp];
                sd = m_ptrsSurfPhase[isp]->siteDensity();
                resid[kspecial] = sd;
                for (k = 0; k < nsp; k++) {
                    resid[kspecial] -= CSoln[kins + k];
                }
            }
        } else {
            kindexSP = 0;
            for (isp = 0; isp < m_numSurfPhases; isp++) {
                nsp = m_nSpeciesSurfPhase[isp];
                InterfaceKinetics* kinPtr = m_objects[isp];
                size_t surfIndex = kinPtr->surfacePhaseIndex();
                kstart = kinPtr->kineticsSpeciesIndex(0, surfIndex);
                kins = kindexSP;
                kinPtr->getNetProductionRates(DATA_PTR(m_netProductionRatesSave));
                for (k = 0; k < nsp; k++, kindexSP++) {
                    resid[kindexSP] = - m_netProductionRatesSave[kstart + k];
                }
                kspecial = kins + m_spSurfLarge[isp];
                sd = m_ptrsSurfPhase[isp]->siteDensity();
                resid[kspecial] = sd;
                for (k = 0; k < nsp; k++) {
                    resid[kspecial] -= CSoln[kins + k];
                }
            }
        }

        if (m_bulkFunc == BULK_DEPOSITION) {
            kindexSP = m_numTotSurfSpecies;
            for (isp = 0; isp < m_numBulkPhasesSS; isp++) {
                doublereal* XBlk = DATA_PTR(m_numEqn1);
                nsp = m_nSpeciesSurfPhase[isp];
                size_t surfPhaseIndex = m_indexKinObjSurfPhase[isp];
                InterfaceKinetics* m_kin = m_objects[isp];
                grRate = 0.0;
                kstart = m_kin->kineticsSpeciesIndex(0, surfPhaseIndex);
                for (k = 0; k < nsp; k++) {
                    if (m_netProductionRatesSave[kstart + k] > 0.0) {
                        grRate += m_netProductionRatesSave[kstart + k];
                    }
                }
                resid[kindexSP] = m_bulkPhasePtrs[isp]->molarDensity();
                for (k = 0; k < nsp; k++) {
                    resid[kindexSP] -= CSoln[kindexSP + k];
                }
                if (grRate > 0.0) {
                    for (k = 1; k < nsp; k++) {
                        if (m_netProductionRatesSave[kstart + k] > 0.0) {
                            resid[kindexSP + k] = XBlk[k] * grRate
                                                  - m_netProductionRatesSave[kstart + k];
                        } else {
                            resid[kindexSP + k] = XBlk[k] * grRate;
                        }
                    }
                } else {
                    grRate = 1.0E-6;
                    grRate += fabs(m_netProductionRatesSave[kstart + k]);
                    for (k = 1; k < nsp; k++) {
                        resid[kindexSP + k] = grRate * (XBlk[k] - 1.0/nsp);
                    }
                }
                if (do_time) {
                    for (k = 1; k < nsp; k++) {
                        resid[kindexSP + k] +=
                            lenScale / deltaT *
                            (CSoln[kindexSP + k]- CSolnOld[kindexSP + k]);
                    }
                }
                kindexSP += nsp;
            }
        }
    }
}

void solveSP::resjac_eval(SquareMatrix& jac,
                          doublereal resid[], doublereal CSoln[],
                          const doublereal CSolnOld[], const bool do_time,
                          const doublereal deltaT)
{
    size_t kColIndex = 0, nsp, jsp, i, kCol;
    doublereal dc, cSave, sd;
    /*
     * Calculate the residual
     */
    fun_eval(resid, CSoln, CSolnOld, do_time, deltaT);
    /*
     * Now we will look over the columns perturbing each unknown.
     */
    for (jsp = 0; jsp < m_numSurfPhases; jsp++) {
        nsp = m_nSpeciesSurfPhase[jsp];
        sd = m_ptrsSurfPhase[jsp]->siteDensity();
        for (kCol = 0; kCol < nsp; kCol++) {
            cSave = CSoln[kColIndex];
            dc = std::max(1.0E-10 * sd, fabs(cSave) * 1.0E-7);
            CSoln[kColIndex] += dc;
            fun_eval(DATA_PTR(m_numEqn2), CSoln, CSolnOld, do_time, deltaT);
            for (i = 0; i < m_neq; i++) {
                jac(i, kColIndex) = (m_numEqn2[i] - resid[i])/dc;
            }
            CSoln[kColIndex] = cSave;
            kColIndex++;
        }
    }

    if (m_bulkFunc == BULK_DEPOSITION) {
        for (jsp = 0; jsp < m_numBulkPhasesSS; jsp++) {
            nsp = m_numBulkSpecies[jsp];
            sd = m_bulkPhasePtrs[jsp]->molarDensity();
            for (kCol = 0; kCol < nsp; kCol++) {
                cSave = CSoln[kColIndex];
                dc = std::max(1.0E-10 * sd, fabs(cSave) * 1.0E-7);
                CSoln[kColIndex] += dc;
                fun_eval(DATA_PTR(m_numEqn2), CSoln, CSolnOld, do_time, deltaT);
                for (i = 0; i < m_neq; i++) {
                    jac(i, kColIndex) = (m_numEqn2[i] - resid[i])/dc;
                }
                CSoln[kColIndex] = cSave;
                kColIndex++;
            }
        }
    }
}

/*!
 * This function calculates a damping factor for the Newton iteration update
 * vector, dxneg, to insure that all site and bulk fractions, x, remain
 * bounded between zero and one.
 *
 *      dxneg[] = negative of the update vector.
 *
 * The constant "APPROACH" sets the fraction of the distance to the boundary
 * that the step can take.  If the full step would not force any fraction
 * outside of 0-1, then Newton's method is allowed to operate normally.
 */
static doublereal calc_damping(doublereal x[], doublereal dxneg[], size_t dim, int* label)
{
    const doublereal APPROACH = 0.80;
    doublereal    damp = 1.0, xnew, xtop, xbot;
    static doublereal damp_old = 1.0;

    *label = -1;

    for (size_t i = 0; i < dim; i++) {

        /*
         * Calculate the new suggested new value of x[i]
         */

        xnew = x[i] - damp * dxneg[i];

        /*
         *  Calculate the allowed maximum and minimum values of x[i]
         *   - Only going to allow x[i] to converge to zero by a
         *     single order of magnitude at a time
         */

        xtop = 1.0 - 0.1*fabs(1.0-x[i]);
        xbot = fabs(x[i]*0.1) - 1.0e-16;
        if (xnew > xtop)  {
            damp = - APPROACH * (1.0 - x[i]) / dxneg[i];
            *label = int(i);
        } else if (xnew < xbot) {
            damp = APPROACH * x[i] / dxneg[i];
            *label = int(i);
        } else  if (xnew > 3.0*std::max(x[i], 1.0E-10)) {
            damp = - 2.0 * std::max(x[i], 1.0E-10) / dxneg[i];
            *label = int(i);
        }
    }
    damp = std::max(damp, 1e-2);
    /*
     * Only allow the damping parameter to increase by a factor of three each
     * iteration. Heuristic to avoid oscillations in the value of damp
     */

    if (damp > damp_old*3) {
        damp = damp_old*3;
        *label = -1;
    }

    /*
     *      Save old value of the damping parameter for use
     *      in subsequent calls.
     */

    damp_old = damp;
    return damp;

} /* calc_damping */

/*
 *    This function calculates the norm  of an update, dx[],
 *    based on the weighted values of x.
 */
static doublereal calcWeightedNorm(const doublereal wtX[], const doublereal dx[], size_t dim)
{
    doublereal norm = 0.0;
    doublereal tmp;
    if (dim == 0) {
        return 0.0;
    }
    for (size_t i = 0; i < dim; i++) {
        tmp = dx[i] / wtX[i];
        norm += tmp * tmp;
    }
    return sqrt(norm/dim);
}

void solveSP::calcWeights(doublereal wtSpecies[], doublereal wtResid[],
                          const Array2D& Jac, const doublereal CSoln[],
                          const doublereal abstol, const doublereal reltol)
{
    size_t k, jcol, kindex, isp, nsp;
    doublereal sd;
    /*
     * First calculate the weighting factor for the concentrations of
     * the surface species and bulk species.
     */
    kindex = 0;
    for (isp = 0; isp < m_numSurfPhases; isp++) {
        nsp = m_nSpeciesSurfPhase[isp];
        sd = m_ptrsSurfPhase[isp]->siteDensity();
        for (k = 0; k < nsp; k++, kindex++) {
            wtSpecies[kindex] = abstol * sd + reltol * fabs(CSoln[kindex]);
        }
    }
    if (m_bulkFunc == BULK_DEPOSITION) {
        for (isp = 0; isp < m_numBulkPhasesSS; isp++) {
            nsp = m_numBulkSpecies[isp];
            sd = m_bulkPhasePtrs[isp]->molarDensity();
            for (k = 0; k < nsp; k++, kindex++) {
                wtSpecies[kindex] = abstol * sd + reltol * fabs(CSoln[kindex]);
            }
        }
    }
    /*
     * Now do the residual Weights. Since we have the Jacobian, we
     * will use it to generate a number based on the what a significant
     * change in a solution variable does to each residual.
     * This is a row sum scale operation.
     */
    for (k = 0; k < m_neq; k++) {
        wtResid[k] = 0.0;
        for (jcol = 0; jcol < m_neq; jcol++) {
            wtResid[k] += fabs(Jac(k,jcol) * wtSpecies[jcol]);
        }
    }
}

doublereal solveSP::calc_t(doublereal netProdRateSolnSP[],
                          doublereal XMolSolnSP[],
                          int* label, int* label_old,
                          doublereal* label_factor, int ioflag)
{
    size_t k, isp, nsp, kstart;
    doublereal   inv_timeScale = 1.0E-10;
    doublereal sden, tmp;
    size_t kindexSP = 0;
    *label = 0;
    updateMFSolnSP(XMolSolnSP);
    for (isp = 0; isp < m_numSurfPhases; isp++) {
        nsp = m_nSpeciesSurfPhase[isp];

        // Get the interface kinetics associated with this surface
        InterfaceKinetics* m_kin = m_objects[isp];

        // Calculate the start of the species index for surfaces within
        // the InterfaceKinetics object
        size_t surfIndex = m_kin->surfacePhaseIndex();
        kstart = m_kin->kineticsSpeciesIndex(0, surfIndex);
        ThermoPhase& THref = m_kin->thermo(surfIndex);

        m_kin->getNetProductionRates(DATA_PTR(m_numEqn1));

        sden = THref.molarDensity();
        for (k = 0; k < nsp; k++, kindexSP++) {
            size_t kspindex = kstart + k;
            netProdRateSolnSP[kindexSP] = m_numEqn1[kspindex];
            if (XMolSolnSP[kindexSP] <= 1.0E-10) {
                tmp = 1.0E-10;
            } else {
                tmp = XMolSolnSP[kindexSP];
            }
            tmp *= sden;
            tmp = fabs(netProdRateSolnSP[kindexSP]/ tmp);
            if (netProdRateSolnSP[kindexSP]> 0.0) {
                tmp /= 100.;
            }
            if (tmp > inv_timeScale) {
                inv_timeScale = tmp;
                *label = int(kindexSP);
            }
        }
    }

    /*
     * Increase time step exponentially as same species repeatedly
     * controls time step
     */
    if (*label == *label_old) {
        *label_factor *= 1.5;
    } else {
        *label_old = *label;
        *label_factor = 1.0;
    }
    return inv_timeScale / *label_factor;
} /* calc_t */

void solveSP::print_header(int ioflag, int ifunc, doublereal time_scale,
                           int damping, doublereal reltol, doublereal abstol)
{
    if (ioflag) {
        printf("\n================================ SOLVESP CALL SETUP "
               "========================================\n");
        if (ifunc == SFLUX_INITIALIZE) {
            printf("\n  SOLVESP Called with Initialization turned on\n");
            printf("     Time scale input = %9.3e\n", time_scale);
        } else if (ifunc == SFLUX_RESIDUAL) {
            printf("\n   SOLVESP Called to calculate steady state residual\n");
            printf("           from a good initial guess\n");
        } else if (ifunc == SFLUX_JACOBIAN)  {
            printf("\n   SOLVESP Called to calculate steady state Jacobian\n");
            printf("           from a good initial guess\n");
        } else if (ifunc == SFLUX_TRANSIENT) {
            printf("\n   SOLVESP Called to integrate surface in time\n");
            printf("           for a total of %9.3e sec\n", time_scale);
        } else {
            throw CanteraError("solveSP::print_header",
                               "Unknown ifunc flag = " + int2str(ifunc));
        }

        if (m_bulkFunc == BULK_DEPOSITION) {
            printf("     The composition of the Bulk Phases will be calculated\n");
        } else if (m_bulkFunc == BULK_ETCH) {
            printf("     Bulk Phases have fixed compositions\n");
        } else {
            throw CanteraError("solveSP::print_header",
                               "Unknown bulkFunc flag = " + int2str(m_bulkFunc));
        }

        if (damping) {
            printf("     Damping is ON   \n");
        } else {
            printf("     Damping is OFF  \n");
        }

        printf("     Reltol = %9.3e, Abstol = %9.3e\n", reltol, abstol);
    }

    if (ioflag == 1) {
        printf("\n\n\t Iter    Time       Del_t      Damp      DelX   "
               "     Resid    Name-Time    Name-Damp\n");
        printf("\t -----------------------------------------------"
               "------------------------------------\n");
    }
}

void solveSP::printIteration(int ioflag, doublereal damp, int label_d,
                             int label_t, doublereal inv_t, doublereal t_real,
                             size_t iter, doublereal update_norm,
                             doublereal resid_norm, bool do_time, bool final)
{
    size_t i, k;
    string nm;
    if (ioflag == 1) {
        if (final) {
            printf("\tFIN%3s ", int2str(iter).c_str());
        } else {
            printf("\t%6s ", int2str(iter).c_str());
        }
        if (do_time) {
            printf("%9.4e %9.4e ", t_real, 1.0/inv_t);
        } else
            for (i = 0; i < 22; i++) {
                printf(" ");
            }
        if (damp < 1.0) {
            printf("%9.4e ", damp);
        } else
            for (i = 0; i < 11; i++) {
                printf(" ");
            }
        printf("%9.4e %9.4e", update_norm, resid_norm);
        if (do_time) {
            k = m_kinSpecIndex[label_t];
            size_t isp = m_kinObjIndex[label_t];
            InterfaceKinetics* m_kin = m_objects[isp];
            nm = m_kin->kineticsSpeciesName(k);
            printf(" %-16s", nm.c_str());
        } else {
            for (i = 0; i < 16; i++) {
                printf(" ");
            }
        }
        if (label_d >= 0) {
            k = m_kinSpecIndex[label_d];
            size_t isp = m_kinObjIndex[label_d];
            InterfaceKinetics* m_kin = m_objects[isp];
            nm = m_kin->kineticsSpeciesName(k);
            printf(" %-16s", nm.c_str());
        }
        if (final) {
            printf(" -- success");
        }
        printf("\n");
    }
} /* printIteration */

}
