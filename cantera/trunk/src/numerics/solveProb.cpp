/**
 * @file: solveProb.cpp Implicit solver for nonlinear problems
 */

/*
 * Copyright 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "cantera/numerics/solveProb.h"
#include "cantera/base/clockWC.h"
#include "cantera/numerics/ctlapack.h"
#include "cantera/base/stringUtils.h"

using namespace std;
namespace Cantera
{

/***************************************************************************
 *       STATIC ROUTINES DEFINED IN THIS FILE
 ***************************************************************************/

static doublereal calcWeightedNorm(const doublereal [], const doublereal dx[], size_t);

solveProb::solveProb(ResidEval* resid) :
    m_residFunc(resid),
    m_neq(0),
    m_atol(0),
    m_rtol(1.0E-4),
    m_maxstep(1000),
    m_ioflag(0)
{
    m_neq =   m_residFunc->nEquations();

    // Dimension solution vector
    size_t dim1 = std::max<size_t>(1, m_neq);

    m_atol.resize(dim1, 1.0E-9);
    m_netProductionRatesSave.resize(dim1, 0.0);
    m_numEqn1.resize(dim1, 0.0);
    m_numEqn2.resize(dim1, 0.0);
    m_CSolnSave.resize(dim1, 0.0);
    m_CSolnSP.resize(dim1, 0.0);
    m_CSolnSPInit.resize(dim1, 0.0);
    m_CSolnSPOld.resize(dim1, 0.0);
    m_wtResid.resize(dim1, 0.0);
    m_wtSpecies.resize(dim1, 0.0);
    m_resid.resize(dim1, 0.0);
    m_ipiv.resize(dim1, 0);
    m_topBounds.resize(dim1, 1.0);
    m_botBounds.resize(dim1, 0.0);

    m_Jac.resize(dim1, dim1, 0.0);
    m_JacCol.resize(dim1, 0);
    for (size_t k = 0; k < dim1; k++) {
        m_JacCol[k] = m_Jac.ptrColumn(k);
    }

}

solveProb::~solveProb()
{
}

int solveProb::solve(int ifunc, doublereal time_scale,
                     doublereal reltol)
{
    /*
     * The following calculation is a Newton's method to get the surface fractions
     * of the surface and bulk species by requiring that the surface species
     * production rate = 0 and that the bulk fractions are proportional to their
     * production rates.
     */
    doublereal EXTRA_ACCURACY = 0.001;
    if (ifunc == SOLVEPROB_JACOBIAN) {
        EXTRA_ACCURACY *= 0.001;
    }
    int info = 0;
    size_t label_t = npos; /* Species IDs for time control */
    size_t label_d; /* Species IDs for damping control */
    size_t label_t_old = npos;
    doublereal label_factor = 1.0;
    int iter=0; // iteration number on numlinear solver
    int iter_max=1000; // maximum number of nonlinear iterations
    int nrhs=1;
    doublereal deltaT = 1.0E-10; // Delta time step
    doublereal damp=1.0, tmp;
    //  Weighted L2 norm of the residual.  Currently, this is only
    //  used for IO purposes. It doesn't control convergence.
    //  Therefore, it is turned off when DEBUG_SOLVEPROB isn't defined.
    doublereal  resid_norm;
    doublereal inv_t = 0.0;
    doublereal t_real = 0.0, update_norm = 1.0E6;

    bool do_time = false, not_converged = true;

#ifdef DEBUG_SOLVEPROB
#ifdef DEBUG_SOLVEPROB_TIME
    doublereal         t1;
#endif
#else
    if (m_ioflag > 1) {
        m_ioflag = 1;
    }
#endif

#ifdef DEBUG_SOLVEPROB
#ifdef DEBUG_SOLVEPROB_TIME
    Cantera::clockWC wc;
    if (m_ioflag) {
        t1 = wc.secondsWC();
    }
#endif
#endif

    /*
     *       Set the initial value of the do_time parameter
     */
    if (ifunc == SOLVEPROB_INITIALIZE || ifunc == SOLVEPROB_TRANSIENT) {
        do_time = true;
    }

    /*
     *  upload the initial conditions
     */
    m_residFunc->getInitialConditions(t_real, DATA_PTR(m_CSolnSP), DATA_PTR(m_numEqn1));

    /*
     * Store the initial guess in the soln vector,
     *  CSolnSP, and in an separate vector CSolnSPInit.
     */
    std::copy(m_CSolnSP.begin(), m_CSolnSP.end(), m_CSolnSPInit.begin());



    if (m_ioflag) {
        print_header(m_ioflag, ifunc, time_scale, reltol,
                     DATA_PTR(m_netProductionRatesSave));
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
        //  if (iter%5 == 4) {
        //    evalSurfLarge(DATA_PTR(m_CSolnSP));
        // }

        /*
         *    Calculate the value of the time step
         *       - heuristics to stop large oscillations in deltaT
         */
        if (do_time) {
            /* don't hurry increase in time step at the same time as damping */
            if (damp < 1.0) {
                label_factor = 1.0;
            }
            tmp = calc_t(DATA_PTR(m_netProductionRatesSave), DATA_PTR(m_CSolnSP),
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

            if (ifunc == SOLVEPROB_TRANSIENT) {
                tmp = t_real + 1.0/inv_t;
                if (tmp > time_scale) {
                    inv_t = 1.0/(time_scale - t_real);
                }
            }
        } else {
            /* make steady state calc a step of 1 million seconds to
               prevent singular jacobians for some pathological cases */
            inv_t = 1.0e-6;
        }
        deltaT = 1.0/inv_t;

        /*
         * Call the routine to numerically evaluation the jacobian
         * and residual for the current iteration.
         */
        resjac_eval(m_JacCol, DATA_PTR(m_resid), DATA_PTR(m_CSolnSP),
                    DATA_PTR(m_CSolnSPOld), do_time, deltaT);

        /*
         * Calculate the weights. Make sure the calculation is carried
         * out on the first iteration.
         */
        if (iter%4 == 1) {
            calcWeights(DATA_PTR(m_wtSpecies), DATA_PTR(m_wtResid),
                        DATA_PTR(m_CSolnSP));
        }

        /*
         *    Find the weighted norm of the residual
         */
        resid_norm = calcWeightedNorm(DATA_PTR(m_wtResid), DATA_PTR(m_resid), m_neq);

#ifdef DEBUG_SOLVEPROB
        if (m_ioflag > 1) {
            printIterationHeader(m_ioflag, damp, inv_t, t_real, iter, do_time);
            /*
             *    Print out the residual and jacobian
             */
            printResJac(m_ioflag, m_neq, m_Jac, DATA_PTR(m_resid),
                        DATA_PTR(m_wtResid), resid_norm);
        }
#endif

        /*
         *  Solve Linear system (with LAPACK).  The solution is in resid[]
         */

        ct_dgetrf(m_neq, m_neq, m_JacCol[0], m_neq, DATA_PTR(m_ipiv), info);
        if (info==0) {
            ct_dgetrs(ctlapack::NoTranspose, m_neq, nrhs, m_JacCol[0],
                      m_neq, DATA_PTR(m_ipiv), DATA_PTR(m_resid), m_neq,
                      info);
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

                printf("-----\n");
            }
            if (do_time) {
                t_real += time_scale;
            }
#ifdef DEBUG_SOLVEPROB
            if (m_ioflag) {
                printf("\nResidual is small, forcing convergence!\n");
            }
#endif
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


        if (do_time) {
            t_real += damp/inv_t;
        }

        if (m_ioflag) {
            printIteration(m_ioflag, damp, label_d, label_t,  inv_t, t_real, iter,
                           update_norm, resid_norm,
                           DATA_PTR(m_netProductionRatesSave),
                           DATA_PTR(m_CSolnSP), DATA_PTR(m_resid),
                           DATA_PTR(m_wtSpecies), m_neq, do_time);
        }

        if (ifunc == SOLVEPROB_TRANSIENT) {
            not_converged = (t_real < time_scale);
        } else {
            if (do_time) {
                if (t_real > time_scale ||
                        (resid_norm < 1.0e-7 &&
                         update_norm*time_scale/t_real < EXTRA_ACCURACY))  {
                    do_time = false;
#ifdef DEBUG_SOLVEPROB
                    if (m_ioflag > 1) {
                        printf("\t\tSwitching to steady solve.\n");
                    }
#endif
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
            printf("#$#$#$# Error in solveProb $#$#$#$ \n");
            printf("Newton iter on surface species did not converge, "
                   "update_norm = %e \n", update_norm);
            printf("Continuing anyway\n");
        }
    }
#ifdef DEBUG_SOLVEPROB
#ifdef DEBUG_SOLVEPROB_TIME
    if (m_ioflag) {
        printf("\nEnd of solve, time used: %e\n", wc.secondsWC()-t1);
    }
#endif
#endif

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
        printFinal(m_ioflag, damp, label_d, label_t,  inv_t, t_real, iter,
                   update_norm, resid_norm, DATA_PTR(m_netProductionRatesSave),
                   DATA_PTR(m_CSolnSP), DATA_PTR(m_resid),
                   DATA_PTR(m_wtSpecies),
                   DATA_PTR(m_wtResid), m_neq, do_time);
    }

    /*
     *        Return with the appropriate flag
     */
    if (update_norm > 1.0) {
        return -1;
    }
    return 0;
}

void solveProb::reportState(doublereal* const CSolnSP) const
{
    std::copy(m_CSolnSP.begin(), m_CSolnSP.end(), CSolnSP);
}

void solveProb::fun_eval(doublereal* const resid, const doublereal* const CSoln,
                         const doublereal* const CSolnOld,  const bool do_time,
                         const doublereal deltaT)
{
    /*
     *   This routine uses the m_numEqn1 and m_netProductionRatesSave vectors
     *   as temporary internal storage.
     */
    if (do_time) {
        m_residFunc->evalSimpleTD(0.0, CSoln, CSolnOld, deltaT, resid);
    } else {
        m_residFunc->evalSS(0.0, CSoln, resid);
    }
}

void solveProb::resjac_eval(std::vector<doublereal*> &JacCol,
                            doublereal resid[], doublereal CSoln[],
                            const doublereal CSolnOld[], const bool do_time,
                            const doublereal deltaT)
{
    doublereal dc, cSave, sd;
    doublereal* col_j;
    /*
     * Calculate the residual
     */
    fun_eval(resid, CSoln, CSolnOld, do_time, deltaT);
    /*
     * Now we will look over the columns perturbing each unknown.
     */

    for (size_t kCol = 0; kCol < m_neq; kCol++) {
        cSave = CSoln[kCol];
        sd = fabs(cSave) + fabs(CSoln[kCol]) + m_atol[kCol] * 1.0E6;
        if (sd < 1.0E-200) {
            sd = 1.0E-4;
        }
        dc = std::max(1.0E-11 * sd, fabs(cSave) * 1.0E-6);
        CSoln[kCol] += dc;
        // Use the m_numEqn2 vector as temporary internal storage.
        fun_eval(DATA_PTR(m_numEqn2), CSoln, CSolnOld, do_time, deltaT);
        col_j = JacCol[kCol];
        for (size_t i = 0; i < m_neq; i++) {
            col_j[i] = (m_numEqn2[i] - resid[i])/dc;
        }
        CSoln[kCol] = cSave;
    }

}

doublereal solveProb::calc_damping(doublereal x[], doublereal dxneg[], size_t dim, size_t* label)
{
    const doublereal APPROACH = 0.50;
    doublereal  damp = 1.0, xnew, xtop, xbot;
    static doublereal damp_old = 1.0;
    *label = npos;

    for (size_t i = 0; i < dim; i++) {
        doublereal topBounds = m_topBounds[i];
        doublereal botBounds = m_botBounds[i];
        /*
         * Calculate the new suggested new value of x[i]
         */
        double delta_x = - dxneg[i];
        xnew = x[i] - damp * dxneg[i];

        /*
         *  Calculate the allowed maximum and minimum values of x[i]
         *   - Only going to allow x[i] to converge to  the top and bottom bounds by a
         *     single order of magnitude at one time
         */
        bool canCrossOrigin = false;
        if (topBounds > 0.0 && botBounds < 0.0) {
            canCrossOrigin = true;
        }

        xtop = topBounds - 0.1 * fabs(topBounds - x[i]);

        xbot = botBounds + 0.1 * fabs(x[i] - botBounds);

        if (xnew > xtop) {
            damp = - APPROACH * (xtop - x[i]) / dxneg[i];
            *label = i;
        } else if (xnew < xbot) {
            damp = APPROACH * (x[i] - xbot) / dxneg[i];
            *label = i;
        }
        // else  if (fabs(xnew) > 2.0*MAX(fabs(x[i]), 1.0E-10)) {
        //    damp = 0.5 * MAX(fabs(x[i]), 1.0E-9)/ fabs(xnew);
        //    *label = i;
        //     }
        double denom = fabs(x[i]) + 1.0E5 * m_atol[i];
        if ((fabs(delta_x) / denom) > 0.3) {
            double newdamp = 0.3 * denom / fabs(delta_x);
            if (canCrossOrigin) {
                if (xnew * x[i] < 0.0) {
                    if (fabs(x[i]) < 1.0E8 * m_atol[i]) {
                        newdamp = 2.0 * fabs(x[i]) / fabs(delta_x);
                    }
                }
            }
            damp = std::min(damp, newdamp);
        }

    }

    /*
     * Only allow the damping parameter to increase by a factor of three each
     * iteration. Heuristic to avoid oscillations in the value of damp
     */
    if (damp > damp_old*3) {
        damp = damp_old*3;
        *label = npos;
    }

    /*
     *      Save old value of the damping parameter for use
     *      in subsequent calls.
     */
    damp_old = damp;
    return damp;

}

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

void solveProb::calcWeights(doublereal wtSpecies[], doublereal wtResid[],
                            const doublereal CSoln[])
{
    /*
     * First calculate the weighting factor
     */

    for (size_t k = 0; k < m_neq; k++) {
        wtSpecies[k] = m_atol[k] + m_rtol * fabs(CSoln[k]);
    }
    /*
     * Now do the residual Weights. Since we have the Jacobian, we
     * will use it to generate a number based on the what a significant
     * change in a solution variable does to each residual.
     * This is a row sum scale operation.
     */
    for (size_t k = 0; k < m_neq; k++) {
        wtResid[k] = 0.0;
        for (size_t jcol = 0; jcol < m_neq; jcol++) {
            wtResid[k] += fabs(m_Jac(k,jcol) * wtSpecies[jcol]);
        }
    }
}

doublereal solveProb::
calc_t(doublereal netProdRateSolnSP[], doublereal Csoln[],
       size_t* label, size_t* label_old, doublereal* label_factor, int ioflag)
{
    doublereal tmp, inv_timeScale=0.0;
    for (size_t k = 0; k < m_neq; k++) {
        if (Csoln[k] <= 1.0E-10) {
            tmp = 1.0E-10;
        } else {
            tmp = Csoln[k];
        }
        tmp = fabs(netProdRateSolnSP[k]/ tmp);


        if (netProdRateSolnSP[k]> 0.0) {
            tmp /= 100.;
        }
        if (tmp > inv_timeScale) {
            inv_timeScale = tmp;
            *label = k;
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
    inv_timeScale = inv_timeScale / *label_factor;
#ifdef DEBUG_SOLVEPROB
    if (ioflag > 1) {
        if (*label_factor > 1.0) {
            printf("Delta_t increase due to repeated controlling species = %e\n",
                   *label_factor);
        }
        int kkin = m_kinSpecIndex[*label];

        string sn = "  "
                    printf("calc_t: spec=%d(%s)  sf=%e  pr=%e  dt=%e\n",
                           *label, sn.c_str(), XMolSolnSP[*label],
                           netProdRateSolnSP[*label], 1.0/inv_timeScale);
    }
#endif

    return inv_timeScale;

}

void solveProb::setBounds(const doublereal botBounds[], const doublereal topBounds[])
{
    for (size_t k = 0; k < m_neq; k++) {
        m_botBounds[k] = botBounds[k];
        m_topBounds[k] = topBounds[k];
    }
}

#ifdef DEBUG_SOLVEPROB
void solveProb::printResJac(int ioflag, int neq, const Array2D& Jac,
                            doublereal resid[], doublereal wtRes[],
                            doublereal norm)
{

}
#endif

void solveProb::print_header(int ioflag, int ifunc, doublereal time_scale,
                             doublereal reltol,
                             doublereal netProdRate[])
{
    int damping = 1;
    if (ioflag) {
        printf("\n================================ SOLVEPROB CALL SETUP "
               "========================================\n");
        if (ifunc == SOLVEPROB_INITIALIZE) {
            printf("\n  SOLVEPROB Called with Initialization turned on\n");
            printf("     Time scale input = %9.3e\n", time_scale);
        } else if (ifunc == SOLVEPROB_RESIDUAL) {
            printf("\n   SOLVEPROB Called to calculate steady state residual\n");
            printf("           from a good initial guess\n");
        } else if (ifunc == SOLVEPROB_JACOBIAN)  {
            printf("\n   SOLVEPROB Called to calculate steady state jacobian\n");
            printf("           from a good initial guess\n");
        } else if (ifunc == SOLVEPROB_TRANSIENT) {
            printf("\n   SOLVEPROB Called to integrate surface in time\n");
            printf("           for a total of %9.3e sec\n", time_scale);
        } else {
            fprintf(stderr,"Unknown ifunc flag = %d\n", ifunc);
            exit(EXIT_FAILURE);
        }



        if (damping) {
            printf("     Damping is ON   \n");
        } else {
            printf("     Damping is OFF  \n");
        }

        printf("     Reltol = %9.3e, Abstol = %9.3e\n", reltol, m_atol[0]);
    }

    /*
     *   Print out the initial guess
     */
#ifdef DEBUG_SOLVEPROB
    if (ioflag > 1) {
        printf("\n================================ INITIAL GUESS "
               "========================================\n");
        int kindexSP = 0;
        for (int isp = 0; isp < m_numSurfPhases; isp++) {
            InterfaceKinetics* m_kin = m_objects[isp];
            int surfIndex = m_kin->surfacePhaseIndex();
            int nPhases = m_kin->nPhases();
            m_kin->getNetProductionRates(netProdRate);
            updateMFKinSpecies(XMolKinSpecies, isp);

            printf("\n IntefaceKinetics Object # %d\n\n", isp);

            printf("\t  Number of Phases = %d\n", nPhases);
            printf("\t  Phase:SpecName      Prod_Rate  MoleFraction   kindexSP\n");
            printf("\t  -------------------------------------------------------"
                   "----------\n");

            int kspindex = 0;
            bool inSurfacePhase = false;
            for (int ip = 0; ip < nPhases; ip++) {
                if (ip == surfIndex) {
                    inSurfacePhase = true;
                } else {
                    inSurfacePhase = false;
                }
                ThermoPhase& THref = m_kin->thermo(ip);
                int nsp = THref.nSpecies();
                string pname = THref.id();
                for (int k = 0; k < nsp; k++) {
                    string sname = THref.speciesName(k);
                    string cname = pname + ":" + sname;
                    if (inSurfacePhase) {
                        printf("\t  %-24s %10.3e %10.3e      %d\n", cname.c_str(),
                               netProdRate[kspindex], XMolKinSpecies[kspindex],
                               kindexSP);
                        kindexSP++;
                    } else {
                        printf("\t  %-24s %10.3e %10.3e\n", cname.c_str(),
                               netProdRate[kspindex], XMolKinSpecies[kspindex]);
                    }
                    kspindex++;
                }
            }
            printf("=========================================================="
                   "=================================\n");
        }
    }
#endif
    if (ioflag == 1) {
        printf("\n\n\t Iter    Time       Del_t      Damp      DelX   "
               "     Resid    Name-Time    Name-Damp\n");
        printf("\t -----------------------------------------------"
               "------------------------------------\n");
    }
}
//================================================================================================
void solveProb::printIteration(int ioflag, doublereal damp, size_t label_d,
                               size_t label_t,
                               doublereal inv_t, doublereal t_real, int iter,
                               doublereal update_norm, doublereal resid_norm,
                               doublereal netProdRate[], doublereal CSolnSP[],
                               doublereal resid[],
                               doublereal wtSpecies[], size_t dim, bool do_time)
{
    size_t i, k;
    string nm;
    if (ioflag == 1) {

        printf("\t%6d ", iter);
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
            k = label_t;
            printf(" %s", int2str(k).c_str());
        } else {
            for (i = 0; i < 16; i++) {
                printf(" ");
            }
        }
        if (label_d != npos) {
            k = label_d;
            printf(" %s", int2str(k).c_str());
        }
        printf("\n");
    }
#ifdef DEBUG_SOLVEPROB
    else if (ioflag > 1) {

        updateMFSolnSP(XMolSolnSP);
        printf("\n\t      Weighted norm of update = %10.4e\n", update_norm);

        printf("\t  Name            Prod_Rate        XMol       Conc   "
               "  Conc_Old     wtConc");
        if (damp < 1.0) {
            printf(" UnDamped_Conc");
        }
        printf("\n");
        printf("\t---------------------------------------------------------"
               "-----------------------------\n");
        int kindexSP = 0;
        for (int isp = 0; isp < m_numSurfPhases; isp++) {
            int nsp = m_nSpeciesSurfPhase[isp];
            InterfaceKinetics* m_kin = m_objects[isp];
            //int surfPhaseIndex = m_kinObjPhaseIDSurfPhase[isp];
            m_kin->getNetProductionRates(DATA_PTR(m_numEqn1));
            for (int k = 0; k < nsp; k++, kindexSP++) {
                int kspIndex = m_kinSpecIndex[kindexSP];
                nm = m_kin->kineticsSpeciesName(kspIndex);
                printf("\t%-16s  %10.3e   %10.3e  %10.3e  %10.3e %10.3e ",
                       nm.c_str(),
                       m_numEqn1[kspIndex],
                       XMolSolnSP[kindexSP],
                       CSolnSP[kindexSP], CSolnSP[kindexSP]+damp*resid[kindexSP],
                       wtSpecies[kindexSP]);
                if (damp < 1.0) {
                    printf("%10.4e ", CSolnSP[kindexSP]+(damp-1.0)*resid[kindexSP]);
                    if (label_d == kindexSP) {
                        printf(" Damp ");
                    }
                }
                if (label_t == kindexSP) {
                    printf(" Tctrl");
                }
                printf("\n");
            }

        }

        printf("\t--------------------------------------------------------"
               "------------------------------\n");
    }
#endif
}

void solveProb::printFinal(int ioflag, doublereal damp, size_t label_d, size_t label_t,
                           doublereal inv_t, doublereal t_real, int iter,
                           doublereal update_norm, doublereal resid_norm,
                           doublereal netProdRateKinSpecies[], const doublereal CSolnSP[],
                           const doublereal resid[],
                           const doublereal wtSpecies[], const doublereal wtRes[],
                           size_t dim, bool do_time)
{
    size_t i, k;
    string nm;
    if (ioflag == 1) {

        printf("\tFIN%3d ", iter);
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
            k = label_t;
            printf(" %s", int2str(k).c_str());
        } else {
            for (i = 0; i < 16; i++) {
                printf(" ");
            }
        }
        if (label_d != npos) {
            k = label_d;

            printf(" %s", int2str(k).c_str());
        }
        printf(" -- success\n");
    }
#ifdef DEBUG_SOLVEPROB
    else if (ioflag > 1) {


        printf("\n================================== FINAL RESULT ========="
               "==================================================\n");

        printf("\n        Weighted norm of solution update = %10.4e\n", update_norm);
        printf("        Weighted norm of residual update = %10.4e\n\n", resid_norm);

        printf("  Name            Prod_Rate        XMol       Conc   "
               "  wtConc      Resid   Resid/wtResid   wtResid");
        if (damp < 1.0) {
            printf(" UnDamped_Conc");
        }
        printf("\n");
        printf("---------------------------------------------------------------"
               "---------------------------------------------\n");

        for (int k = 0; k < m_neq; k++, k++) {
            printf("%-16s  %10.3e   %10.3e  %10.3e  %10.3e %10.3e  %10.3e %10.3e",
                   nm.c_str(),
                   m_numEqn1[k],
                   XMolSolnSP[k],
                   CSolnSP[k],
                   wtSpecies[k],
                   resid[k],
                   resid[k]/wtRes[k], wtRes[k]);
            if (damp < 1.0) {
                printf("%10.4e ", CSolnSP[k]+(damp-1.0)*resid[k]);
                if (label_d == k) {
                    printf(" Damp ");
                }
            }
            if (label_t == k) {
                printf(" Tctrl");
            }
            printf("\n");
        }

        printf("\n");
        printf("==============================================================="
               "============================================\n\n");
    }
#endif
}

#ifdef DEBUG_SOLVEPROB
void solveProb::
printIterationHeader(int ioflag, doublereal damp,doublereal inv_t, doublereal t_real,
                     int iter, bool do_time)
{
    if (ioflag > 1) {
        printf("\n===============================Iteration %5d "
               "=================================\n", iter);
        if (do_time) {
            printf("   Transient step with: Real Time_n-1 = %10.4e sec,", t_real);
            printf(" Time_n = %10.4e sec\n", t_real + 1.0/inv_t);
            printf("                        Delta t = %10.4e sec", 1.0/inv_t);
        } else {
            printf("   Steady Solve ");
        }
        if (damp < 1.0) {
            printf(", Damping value =  %10.4e\n", damp);
        } else {
            printf("\n");
        }
    }
}
#endif

void solveProb::setAtol(const doublereal atol[])
{
    for (size_t k = 0; k < m_neq; k++, k++) {
        m_atol[k] = atol[k];
    }
}

void solveProb::setAtolConst(const doublereal atolconst)
{
    for (size_t k = 0; k < m_neq; k++, k++) {
        m_atol[k] = atolconst;
    }
}

}
