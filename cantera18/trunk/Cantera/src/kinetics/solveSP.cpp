/*
 * @file: solveSP.cpp Implicit surface site concentration solver
 */
/*
 * $Id: solveSP.cpp,v 1.4 2008/12/17 17:09:37 hkmoffa Exp $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "solveSP.h"
#include "clockWC.h"
#include "ctlapack.h"

/* Standard include files */

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <vector>

using namespace std;
namespace Cantera {

    /***************************************************************************
     *       STATIC ROUTINES DEFINED IN THIS FILE
     ***************************************************************************/

    static double calc_damping(double *x, double *dx, int dim, int *);
    static double calcWeightedNorm(const double [], const double dx[], int);

    /***************************************************************************
     *                    LAPACK PROTOTYPES
     ***************************************************************************/
    //#define FSUB_TYPE void
    // extern "C" {
    //  extern FSUB_TYPE dgetrf_(int *, int *, double *, int *, int [], int *);
    // extern FSUB_TYPE dgetrs_(char *, int *, int *, double *, int *, int [],
    //			     double [], int *, int *, unsigned int);
    // }
    /*****************************************************************************
     *   PROTOTYPES and PREPROC DIRECTIVES FOR MISC. ROUTINES 
     *****************************************************************************/

#ifndef MAX
#  define MAX(x,y) (( (x) > (y) ) ? (x) : (y))     /* max function */
#endif

#ifndef DAMPING
#  define DAMPING true
#endif

    /***************************************************************************
     *  solveSP Class Definitinos
     ***************************************************************************/

    // Main constructor
    solveSP::solveSP(ImplicitSurfChem *surfChemPtr, int bulkFunc) :
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
        int numPossibleSurfPhases = m_objects.size();
        for (int n = 0; n < numPossibleSurfPhases; n++) {
            InterfaceKinetics *m_kin = m_objects[n];
            int surfPhaseIndex = m_kin->surfacePhaseIndex();
            if (surfPhaseIndex >= 0) {
                m_numSurfPhases++;
                m_indexKinObjSurfPhase.push_back(n);
                m_kinObjPhaseIDSurfPhase.push_back(surfPhaseIndex);
            } else {
                throw CanteraError("solveSP", 
                    "InterfaceKinetics object has no surface phase");
            }
            ThermoPhase *tp = &(m_kin->thermo(surfPhaseIndex));
            SurfPhase *sp = dynamic_cast<SurfPhase *>(tp);
            if (!sp) {
                throw CanteraError("solveSP", 
                    "Inconsistent ThermoPhase object within "
                    "InterfaceKinetics object");
            }

            m_ptrsSurfPhase.push_back(sp);
            int nsp = sp->nSpecies();
            m_nSpeciesSurfPhase.push_back(nsp);
            m_numTotSurfSpecies += nsp;

        }
        /*
         * We rely on ordering to figure things out
         */
        if (1) {
            //m_numBulkPhases = m_kin0->nPhases() - 1 - m_numSurfPhases;
            // Disable the capability until we figure out what is going on
            m_numBulkPhasesSS = 0;
            //if (m_numBulkPhasesSS > 0) {
            //m_numBulkSpecies.resize(m_numBulkPhasesSS, 0);
            //m_bulkPhasePtrs.resize(m_numBulkPhasesSS, 0);
            //m_bulkIndex = 1;
            //if (m_bulkIndex == surfPhaseIndex) {
            // m_bulkIndex += m_numSurfPhases;
            //}
	
            //for (i = 0; i < m_numBulkPhasesSS; i++) {
            //  m_bulkPhasePtrs[i] = &(m_kin0->thermo(m_bulkIndex + i));
            // m_numBulkSpecies[i] =  m_bulkPhasePtrs[i]->nSpecies();
            // m_numTotBulkSpeciesSS += m_numBulkSpecies[i];
            //}
            //}
        }

        if (bulkFunc == BULK_DEPOSITION) {
            m_neq = m_numTotSurfSpecies + m_numTotBulkSpeciesSS;
        } else {
            m_neq = m_numTotSurfSpecies;
        }
    
        m_maxTotSpecies = 0;
        for (int n = 0; n < m_numSurfPhases; n++) {
            int tsp =  m_objects[n]->nTotalSpecies();
            m_maxTotSpecies = MAX(m_maxTotSpecies, tsp);
        }
        m_maxTotSpecies = MAX(m_maxTotSpecies, m_neq);

 
        m_netProductionRatesSave.resize(m_maxTotSpecies, 0.0);
        m_numEqn1.resize(m_maxTotSpecies, 0.0);
        m_numEqn2.resize(m_maxTotSpecies, 0.0);
        m_XMolKinSpecies.resize(m_maxTotSpecies, 0.0);
        m_CSolnSave.resize(m_neq, 0.0);
    
        m_spSurfLarge.resize(m_numSurfPhases, 0);
   
        m_kinSpecIndex.resize(m_numTotSurfSpecies + m_numTotBulkSpeciesSS, 0);
        m_kinObjIndex.resize(m_numTotSurfSpecies + m_numTotBulkSpeciesSS, 0);
        m_eqnIndexStartSolnPhase.resize(m_numSurfPhases + m_numBulkPhasesSS, 0);
  
        int kindexSP = 0;
        int isp, k, nsp, kstart;
        for (isp = 0; isp < m_numSurfPhases; isp++) {
            int iKinObject = m_indexKinObjSurfPhase[isp];
            InterfaceKinetics *m_kin = m_objects[iKinObject];
            int surfPhaseIndex = m_kinObjPhaseIDSurfPhase[isp];
            kstart = m_kin->kineticsSpeciesIndex(0, surfPhaseIndex);
            nsp =  m_nSpeciesSurfPhase[isp];
            m_eqnIndexStartSolnPhase[isp] = kindexSP;
            for (k = 0; k < nsp; k++, kindexSP++) {
                m_kinSpecIndex[kindexSP] = kstart + k;
                m_kinObjIndex[kindexSP] = isp;
            }
        }
        if (0) {
            //for (isp = 0; isp < m_numBulkPhasesSS; isp++) {
            //nt iKinObject = m_bulkKinObjID[isp];
            //InterfaceKinetics *m_kin = m_objects[iKinObject];
            //int bulkIndex = m_bulkKinObjPhaseID[isp];
            //kstart = m_kin->kineticsSpeciesIndex(0, bulkIndex);
            //	nsp =  m_numBulkSpecies[isp];
            //m_eqnIndexStartSolnPhase[isp] = kindexSP;
            //for (k = 0; k < nsp; k++, kindexSP++) {
            //  m_kinSpecIndex[kindexSP] = kstart + k;
            //  m_kinObjIndex[kindexSP] = m_numSurfPhases + isp;
            //}
            //}
        }

        // Dimension solution vector
        int dim1 = MAX(1, m_neq);
        m_CSolnSP.resize(dim1, 0.0);
        m_CSolnSPInit.resize(dim1, 0.0);
        m_CSolnSPOld.resize(dim1, 0.0);
        m_wtResid.resize(dim1, 0.0);
        m_wtSpecies.resize(dim1, 0.0);
        m_resid.resize(dim1, 0.0);
        m_ipiv.resize(dim1, 0);
    
        m_Jac.resize(dim1, dim1, 0.0);
        m_JacCol.resize(dim1, 0);
        for (int k = 0; k < dim1; k++) {
            m_JacCol[k] = m_Jac.ptrColumn(k);
        }
    }

    // Empty destructor
    solveSP::~solveSP() {
    }

    /*
     * The following calculation is a Newton's method to
     * get the surface fractions of the surface and bulk species by 
     * requiring that the
     * surface species production rate = 0 and that the bulk fractions are
     * proportional to their production rates. 
     */
    int solveSP::solveSurfProb(int ifunc, double time_scale, double TKelvin, 
        double PGas, double reltol, double abstol)
    {
        double EXTRA_ACCURACY = 0.001;
        if (ifunc == SFLUX_JACOBIAN) {
            EXTRA_ACCURACY *= 0.001;
        }
        int k, irow;
        int  jcol, info = 0;
        int label_t=-1; /* Species IDs for time control */
        int label_d; /* Species IDs for damping control */
        int        label_t_old=-1;
        double     label_factor = 1.0;
        int iter=0; // iteration number on numlinear solver
        int iter_max=1000; // maximum number of nonlinear iterations
        int nrhs=1;
        double deltaT = 1.0E-10; // Delta time step
        double damp=1.0, tmp;
        //  Weighted L2 norm of the residual.  Currently, this is only
        //  used for IO purposes. It doesn't control convergence.
        //  Therefore, it is turned off when DEBUG_SOLVESP isn't defined.
        double  resid_norm;
        double inv_t = 0.0;
        double t_real = 0.0, update_norm = 1.0E6;

        bool do_time = false, not_converged = true;

#ifdef DEBUG_SOLVESP
#ifdef DEBUG_SOLVESP_TIME
        double         t1;
#endif
#else
        if (m_ioflag > 1) {
            m_ioflag = 1;
        }
#endif

#ifdef DEBUG_SOLVESP
#ifdef DEBUG_SOLVESP_TIME
        Cantera::clockWC wc;
        if (m_ioflag)  t1 = wc.secondsWC();
#endif
#endif

        /*
         *       Set the initial value of the do_time parameter
         */
        if (ifunc == SFLUX_INITIALIZE || ifunc == SFLUX_TRANSIENT) do_time = true;
 
        /*
         *    Store the initial guess for the surface problem in the soln vector,
         *  CSoln, and in an separate vector CSolnInit.
         */
        int loc = 0;
        for (int n = 0; n < m_numSurfPhases; n++) {
            SurfPhase *sf_ptr =  m_ptrsSurfPhase[n];
            sf_ptr->getConcentrations(DATA_PTR(m_numEqn1));
            int nsp = m_nSpeciesSurfPhase[n];
            for (k = 0; k <nsp; k++) {
                m_CSolnSP[loc] = m_numEqn1[k];
                loc++;
            }
        }
  
        if (m_bulkFunc == BULK_DEPOSITION) {
            //for (isp = 0; isp < m_numBulkPhasesSS; isp++) {
            //ThermoPhase *bf_ptr =  m_bulkPhasePtrs[isp];
            //bf_ptr->getConcentrations(DATA_PTR(m_numEqn1));
            //int nsp = m_numBulkSpecies[isp];
            //for (k = 0; k < nsp; k++, kindex++) {
            //  m_CSolnSP[loc] = m_numEqn1[k];
            // loc++;
            //}
            //}
        }
        std::copy(m_CSolnSP.begin(), m_CSolnSP.end(), m_CSolnSPInit.begin());

        // Calculate the largest species in each phase
        evalSurfLarge(DATA_PTR(m_CSolnSP));
        /*
         * Get the net production rate of all species in the kinetics manager.
         */
        // m_kin->getNetProductionRates(DATA_PTR(m_netProductionRatesSave));
  
        if (m_ioflag) {
            print_header(m_ioflag, ifunc, time_scale, DAMPING, reltol, abstol, 
                TKelvin, PGas, DATA_PTR(m_netProductionRatesSave),
                DATA_PTR(m_XMolKinSpecies));
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
                if (damp < 1.0) label_factor = 1.0;
                tmp = calc_t(DATA_PTR(m_netProductionRatesSave), 
                    DATA_PTR(m_XMolKinSpecies),
                    &label_t, &label_t_old,  &label_factor, m_ioflag);
                if (iter < 10)
                    inv_t = tmp;
                else if (tmp > 2.0*inv_t)
                    inv_t =  2.0*inv_t;
                else {
                    inv_t = tmp;
                }

                /*
                 *   Check end condition
                 */

                if (ifunc == SFLUX_TRANSIENT) {
                    tmp = t_real + 1.0/inv_t;
                    if (tmp > time_scale) inv_t = 1.0/(time_scale - t_real);
                }
            }
            else {
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
		    m_Jac, DATA_PTR(m_CSolnSP), abstol, reltol);
            }
 
            /*
             *    Find the weighted norm of the residual
             */
            resid_norm = calcWeightedNorm(DATA_PTR(m_wtResid), 
                DATA_PTR(m_resid), m_neq);

#ifdef DEBUG_SOLVESP
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
            // (void) dgetrf_(&m_neq, &m_neq, m_JacCol[0], &m_neq, 
            //	     DATA_PTR(m_ipiv), &info);
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
                for (jcol = 0; jcol < m_neq; jcol++) m_resid[jcol] = 0.0;

                /* print out some helpful info */
                if (m_ioflag > 1) {
                    printf("-----\n");
                    printf("solveSurfProb: iter %d t_real %g delta_t %g\n\n",
                        iter,t_real, 1.0/inv_t);
                    printf("solveSurfProb: init guess, current concentration,"
                        "and prod rate:\n");
                    for (jcol = 0; jcol < m_neq; jcol++) {
                        printf("\t%d  %g %g %g\n", jcol, m_CSolnSPInit[jcol], m_CSolnSP[jcol], 
                            m_netProductionRatesSave[m_kinSpecIndex[jcol]]);
                    }
                    printf("-----\n");
                }
                if (do_time) t_real += time_scale;
#ifdef DEBUG_SOLVESP
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

#ifdef DAMPING
            damp = calc_damping( DATA_PTR(m_CSolnSP), DATA_PTR(m_resid), m_neq, &label_d);
#endif

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
            for (irow = 0; irow < m_neq; irow++) m_CSolnSP[irow] -= damp * m_resid[irow];
            for (irow = 0; irow < m_neq; irow++) {
                m_CSolnSP[irow] = MAX(0.0, m_CSolnSP[irow]);
            }
            updateState( DATA_PTR(m_CSolnSP));
  
            if (do_time) t_real += damp/inv_t;

            if (m_ioflag) {
                printIteration(m_ioflag, damp, label_d, label_t,  inv_t, t_real, iter,
                    update_norm, resid_norm,
                    DATA_PTR(m_netProductionRatesSave),
                    DATA_PTR(m_CSolnSP), DATA_PTR(m_resid), 
                    DATA_PTR(m_XMolKinSpecies), DATA_PTR(m_wtSpecies),
                    m_neq, do_time);
            }

            if (ifunc == SFLUX_TRANSIENT)
                not_converged = (t_real < time_scale);
            else {
                if (do_time) {
                    if (t_real > time_scale ||
                        (resid_norm < 1.0e-7 && 
                            update_norm*time_scale/t_real < EXTRA_ACCURACY) )  {
                        do_time = false;
#ifdef DEBUG_SOLVESP
                        if (m_ioflag > 1) {
                            printf("\t\tSwitching to steady solve.\n");
                        }
#endif
                    }
                }
                else {
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
#ifdef DEBUG_SOLVESP
#ifdef DEBUG_SOLVESP_TIME
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
                DATA_PTR(m_XMolKinSpecies), DATA_PTR(m_wtSpecies),
                DATA_PTR(m_wtResid), m_neq, do_time,
                TKelvin, PGas);
        }

        /*
         *        Return with the appropriate flag
         */
        if (update_norm > 1.0) {
            return -1;
        }
        return 1;
    }

#undef DAMPING


    /*
     * Update the surface states of the surface phases.
     */
    void solveSP::updateState(const double *CSolnSP) {
        int loc = 0;
        for (int n = 0; n < m_numSurfPhases; n++) {
            m_ptrsSurfPhase[n]->setConcentrations(CSolnSP + loc);
            loc += m_nSpeciesSurfPhase[n];
        }
        //if (m_bulkFunc == BULK_DEPOSITION) {
        //  for (int n = 0; n < m_numBulkPhasesSS; n++) {
        //	m_bulkPhasePtrs[n]->setConcentrations(CSolnSP + loc);
        //	loc += m_numBulkSpecies[n];
        //   }
        //}
    }

    /*
     * Update the mole fractions for phases which are part of the equation set
     */
    void solveSP::updateMFSolnSP(double *XMolSolnSP) {
        for (int isp = 0; isp < m_numSurfPhases; isp++) {
            int keqnStart = m_eqnIndexStartSolnPhase[isp];
            m_ptrsSurfPhase[isp]->getMoleFractions(XMolSolnSP + keqnStart);
        }
        //if (m_bulkFunc == BULK_DEPOSITION) {
        // for (int isp = 0; isp < m_numBulkPhasesSS; isp++) {
        //	int keqnStart = m_eqnIndexStartSolnPhase[isp + m_numSurfPhases];
        //	m_bulkPhasePtrs[isp]->getMoleFractions(XMolSolnSP + keqnStart);
        // }
        //}
    }

    /*
     * Update the mole fractions for phases which are part of a single
     * interfacial kinetics object
     */
    void solveSP::updateMFKinSpecies(double *XMolKinSpecies, int isp) {
        InterfaceKinetics *m_kin = m_objects[isp];
        int nph = m_kin->nPhases();
        for (int iph = 0; iph < nph; iph++) {
            int ksi = m_kin->kineticsSpeciesIndex(0, iph);
            ThermoPhase &thref = m_kin->thermo(iph);
            thref.getMoleFractions(XMolKinSpecies + ksi);
        }
    }
 
    /*
     * Update the vector that keeps track of the largest species in each
     * surface phase.
     */
    void solveSP::evalSurfLarge(const double *CSolnSP) {
        int kindexSP = 0;
        for (int isp = 0; isp < m_numSurfPhases; isp++) {
            int nsp = m_nSpeciesSurfPhase[isp];
            double Clarge = CSolnSP[kindexSP];
            m_spSurfLarge[isp] = 0;
            kindexSP++;
            for (int k = 1; k < nsp; k++, kindexSP++) {
                if (CSolnSP[kindexSP] > Clarge) {
                    Clarge = CSolnSP[kindexSP];
                    m_spSurfLarge[isp] = k;
                }
            }
        }
    }
 
    /*
     * This calculates the net production rates of all species
     *
     * This calculates the function eval.
     *      (should switch to special_species formulation for sum condition)
     *
     * @internal
     *   This routine uses the m_numEqn1 and m_netProductionRatesSave vectors
     *   as temporary internal storage.
     */
    void solveSP::fun_eval(double* resid, const double *CSoln, 
        const double *CSolnOld,  const bool do_time,
        const double deltaT)
    {
        int isp, nsp, kstart, k, kindexSP, kins, kspecial;
        double lenScale = 1.0E-9;
        double sd = 0.0;
        double grRate;
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
                    InterfaceKinetics *kinPtr = m_objects[isp];
                    int surfIndex = kinPtr->surfacePhaseIndex();
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
                    InterfaceKinetics *kinPtr = m_objects[isp];
                    int surfIndex = kinPtr->surfacePhaseIndex();
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
                    double *XBlk = DATA_PTR(m_numEqn1);
                    //ThermoPhase *THptr = m_bulkPhasePtrs[isp];
                    //THptr->getMoleFractions(XBlk);
                    nsp = m_nSpeciesSurfPhase[isp];
                    int surfPhaseIndex = m_indexKinObjSurfPhase[isp];
                    InterfaceKinetics *m_kin = m_objects[isp];
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

    /*
     * Calculate the Jacobian and residual
     *
     * @internal
     *   This routine uses the m_numEqn2 vector
     *   as temporary internal storage.
     */
    void solveSP::resjac_eval(std::vector<double*> &JacCol,
        double resid[], double CSoln[], 
        const double CSolnOld[], const bool do_time, 
        const double deltaT)
    {
        int kColIndex = 0, nsp, jsp, i, kCol;
        double dc, cSave, sd;
        double *col_j;
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
                dc = fmaxx(1.0E-10 * sd, fabs(cSave) * 1.0E-7);
                CSoln[kColIndex] += dc;
                fun_eval(DATA_PTR(m_numEqn2), CSoln, CSolnOld, do_time, deltaT);
                col_j = JacCol[kColIndex];
                for (i = 0; i < m_neq; i++) {
                    col_j[i] = (m_numEqn2[i] - resid[i])/dc;
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
                    dc = fmaxx(1.0E-10 * sd, fabs(cSave) * 1.0E-7);
                    CSoln[kColIndex] += dc;
                    fun_eval(DATA_PTR(m_numEqn2), CSoln, CSolnOld, do_time, deltaT);
                    col_j = JacCol[kColIndex];
                    for (i = 0; i < m_neq; i++) {
                        col_j[i] = (m_numEqn2[i] - resid[i])/dc;
                    }
                    CSoln[kColIndex] = cSave;
                    kColIndex++;
                }
            }
        }
    }


#define APPROACH 0.80

    static double calc_damping(double x[], double dxneg[], int dim, int *label)

        /* This function calculates a damping factor for the Newton iteration update
         * vector, dxneg, to insure that all site and bulk fractions, x, remain
         * bounded between zero and one.
         *
         *      dxneg[] = negative of the update vector.
         *
         * The constant "APPROACH" sets the fraction of the distance to the boundary
         * that the step can take.  If the full step would not force any fraction
         * outside of 0-1, then Newton's method is allowed to operate normally.
         */

    {
        int       i;
        double    damp = 1.0, xnew, xtop, xbot;
        static double damp_old = 1.0;

        *label = -1;

        for (i = 0; i < dim; i++) {

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
            if (xnew > xtop )  {
                damp = - APPROACH * (1.0 - x[i]) / dxneg[i];
                *label = i;
            }
            else if (xnew < xbot) {
                damp = APPROACH * x[i] / dxneg[i];
                *label = i;
            } else  if (xnew > 3.0*MAX(x[i], 1.0E-10)) {
                damp = - 2.0 * MAX(x[i], 1.0E-10) / dxneg[i];
                *label = i;
            }
        }

        if (damp < 1.0e-2) damp = 1.0e-2;
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
#undef APPROACH

    /*
     *    This function calculates the norm  of an update, dx[], 
     *    based on the weighted values of x.
     */
    static double calcWeightedNorm(const double wtX[], const double dx[], int dim) {
        double norm = 0.0;
        double tmp;
        if (dim == 0) return 0.0;
        for (int i = 0; i < dim; i++) {
            tmp = dx[i] / wtX[i];
            norm += tmp * tmp;
        }
        return (sqrt(norm/dim));
    } 
 
    /*
     * Calculate the weighting factors for norms wrt both the species
     * concentration unknowns and the residual unknowns.
     *
     */
    void solveSP::calcWeights(double wtSpecies[], double wtResid[],
        const Array2D &Jac, const double CSoln[],
        const double abstol, const double reltol)
    {
        int k, jcol, kindex, isp, nsp;
        double sd;
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
 
    /*
     *    This routine calculates a pretty conservative 1/del_t based
     *    on  MAX_i(sdot_i/(X_i*SDen0)).  This probably guarantees
     *    diagonal dominance.
     *
     *     Small surface fractions are allowed to intervene in the del_t
     *     determination, no matter how small.  This may be changed.
     *     Now minimum changed to 1.0e-12,
     *
     *     Maximum time step set to time_scale.
     */
    double solveSP::
    calc_t(double netProdRateSolnSP[], double XMolSolnSP[],
        int *label, int *label_old, double *label_factor, int ioflag)
    {
        int  k, isp, nsp, kstart;
        double   inv_timeScale = 1.0E-10;
        double sden, tmp;
        int kindexSP = 0;
        *label = 0;
        int ispSpecial = 0;
        int kspSpecial = 0;
        updateMFSolnSP(XMolSolnSP);
        for (isp = 0; isp < m_numSurfPhases; isp++) {
            nsp = m_nSpeciesSurfPhase[isp];
    
            // Get the interface kinetics associated with this surface
            InterfaceKinetics *m_kin = m_objects[isp];
      
            // Calcuate the start of the species index for surfaces within
            // the InterfaceKinetics object
            int surfIndex = m_kin->surfacePhaseIndex();
            kstart = m_kin->kineticsSpeciesIndex(0, surfIndex);
            ThermoPhase& THref = m_kin->thermo(surfIndex);
 
            m_kin->getNetProductionRates(DATA_PTR(m_numEqn1));

            sden = THref.molarDensity();
            for (k = 0; k < nsp; k++, kindexSP++) {
                int kspindex = kstart + k;
                netProdRateSolnSP[kindexSP] = m_numEqn1[kspindex];
                if (XMolSolnSP[kindexSP] <= 1.0E-10) {
                    tmp = 1.0E-10;
                } else { 
                    tmp = XMolSolnSP[kindexSP];
                }
                tmp *= sden;
                tmp = fabs(netProdRateSolnSP[kindexSP]/ tmp);
                if (netProdRateSolnSP[kindexSP]> 0.0) tmp /= 100.;
                if (tmp > inv_timeScale) {
                    inv_timeScale = tmp;
                    *label = kindexSP;
                    ispSpecial = isp;
                    kspSpecial = k;
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
        inv_timeScale = inv_timeScale / *label_factor;
#ifdef DEBUG_SOLVESP
        if (ioflag > 1) {
            if (*label_factor > 1.0) {
                printf("Delta_t increase due to repeated controlling species = %e\n",
                    *label_factor);
            }
            int kkin = m_kinSpecIndex[*label];
            InterfaceKinetics *m_kin = m_objects[ispSpecial];
            string sn = m_kin->kineticsSpeciesName(kkin);
            printf("calc_t: spec=%d(%s)  sf=%e  pr=%e  dt=%e\n",
                *label, sn.c_str(), XMolSolnSP[*label], 
                netProdRateSolnSP[*label], 1.0/inv_timeScale);
        }
#endif
    
        return (inv_timeScale);
    
    } /* calc_t */


    /**
     *  printResJac(): prints out the residual and Jacobian.
     *
     */
#ifdef DEBUG_SOLVESP
    void solveSP::printResJac(int ioflag, int neq, const Array2D &Jac,
        double resid[], double wtRes[],
        double norm)
    {
        int i, j, isp, nsp, irowKSI, irowISP;
        int kstartKSI;
        int kindexSP = 0;
        string sname, pname, cname;
        if (ioflag > 1) {
            printf("     Printout of residual and jacobian\n");
            printf("\t  Residual: weighted norm = %10.4e\n", norm);
            printf("\t  Index       Species_Name      Residual     "
                "Resid/wtRes      wtRes\n");
            for (isp = 0; isp < m_numSurfPhases; isp++) {
                nsp = m_nSpeciesSurfPhase[isp];
                InterfaceKinetics *m_kin = m_objects[isp];
                int surfPhaseIndex = m_kinObjPhaseIDSurfPhase[isp];
                m_kin->getNetProductionRates(DATA_PTR(m_numEqn1));
                kstartKSI = m_kin->kineticsSpeciesIndex(0, surfPhaseIndex);
                SurfPhase *sp_ptr = m_ptrsSurfPhase[isp];
                pname = sp_ptr->id();

                for (int k = 0; k < nsp; k++, kindexSP++) {
                    sname = sp_ptr->speciesName(k);
                    cname = pname + ":" + sname;
                    printf("\t  %d: %-24s: %11.3e  %11.3e  %11.3e\n", kindexSP, 
                        cname.c_str(), resid[kindexSP],
                        resid[kindexSP]/wtRes[kindexSP], wtRes[kindexSP]);
                }
            }
            if (m_bulkFunc == BULK_DEPOSITION) {
                for (isp = 0; isp < m_numBulkPhasesSS; isp++) {
                    // fill in
                }
            }
            if (ioflag > 2) {
                printf("\t  Jacobian:\n");
                for (i = 0; i < m_neq; i++) {
                    irowISP = m_kinObjIndex[i];
                    InterfaceKinetics *m_kin = m_objects[irowISP];
                    irowKSI = m_kinSpecIndex[i];
                    ThermoPhase& THref = m_kin->speciesPhase(irowKSI);
                    int phaseIndex = m_kin->speciesPhaseIndex(irowKSI);
                    kstartKSI = m_kin->kineticsSpeciesIndex(0, phaseIndex);
                    int klocal = i - m_eqnIndexStartSolnPhase[irowISP];
                    sname = THref.speciesName(klocal);
                    printf("\t   Row %d:%-16s:\n", i, sname.c_str());
                    printf("\t     ");
                    for (j = 0; j < m_neq; j++) {
                        printf("%10.4e ", Jac(i,j));
                    }
                    printf("\n");
                }
            }
        }
    } /* printResJac */
#endif
 
    /*
     * Optional printing at the start of the solveSP problem
     */
    void solveSP::print_header(int ioflag, int ifunc, double time_scale, 
        int damping, double reltol, double abstol,  
        double TKelvin,
        double PGas, double netProdRate[],
        double XMolKinSpecies[]) {
        if (ioflag) {
            printf("\n================================ SOLVESP CALL SETUP "
                "========================================\n");
            if (ifunc == SFLUX_INITIALIZE) {
                printf("\n  SOLVESP Called with Initialization turned on\n");
                printf("     Time scale input = %9.3e\n", time_scale);
            }
            else if (ifunc == SFLUX_RESIDUAL) {
                printf("\n   SOLVESP Called to calculate steady state residual\n");
                printf( "           from a good initial guess\n");
            }
            else if (ifunc == SFLUX_JACOBIAN)  {
                printf("\n   SOLVESP Called to calculate steady state jacobian\n");
                printf( "           from a good initial guess\n");
            }
            else if (ifunc == SFLUX_TRANSIENT) {
                printf("\n   SOLVESP Called to integrate surface in time\n");
                printf( "           for a total of %9.3e sec\n", time_scale);
            }
            else {
                fprintf(stderr,"Unknown ifunc flag = %d\n", ifunc);
                exit(EXIT_FAILURE);
            }

            if (m_bulkFunc == BULK_DEPOSITION)
                printf("     The composition of the Bulk Phases will be calculated\n");
            else if (m_bulkFunc == BULK_ETCH)
                printf("     Bulk Phases have fixed compositions\n");
            else {
                fprintf(stderr,"Unknown bulkFunc flag = %d\n", m_bulkFunc);
                exit(EXIT_FAILURE);
            }

            if (damping)
                printf("     Damping is ON   \n");
            else
                printf("     Damping is OFF  \n");

            printf("     Reltol = %9.3e, Abstol = %9.3e\n", reltol, abstol);
        }

        /*
         *   Print out the initial guess
         */
#ifdef DEBUG_SOLVESP
        if (ioflag > 1) {
            printf("\n================================ INITIAL GUESS "
                "========================================\n");
            int kindexSP = 0;
            for (int isp = 0; isp < m_numSurfPhases; isp++) {
                InterfaceKinetics *m_kin = m_objects[isp];
                int surfIndex = m_kin->surfacePhaseIndex();
                int nPhases = m_kin->nPhases();
                m_kin->getNetProductionRates(netProdRate);
                updateMFKinSpecies(XMolKinSpecies, isp);

                printf("\n IntefaceKinetics Object # %d\n\n", isp);

                printf("\t  Number of Phases = %d\n", nPhases);
                printf("\t  Temperature = %10.3e Kelvin\n", TKelvin);
                printf("\t  Pressure    = %10.3g Pa\n\n", PGas);
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
                    ThermoPhase &THref = m_kin->thermo(ip);
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
            printf(    "\t -----------------------------------------------"
                "------------------------------------\n");
        }
    } 

    void solveSP::printIteration(int ioflag, double damp, int label_d,
        int label_t,
        double inv_t, double t_real, int iter,
        double update_norm, double resid_norm,
        double netProdRate[], double CSolnSP[],
        double resid[], double XMolSolnSP[],
        double wtSpecies[], int dim, bool do_time)
    {
        int i, k;
        string nm;
        if (ioflag == 1) {

            printf("\t%6d ", iter);
            if (do_time)
                printf("%9.4e %9.4e ", t_real, 1.0/inv_t);
            else
                for (i = 0; i < 22; i++) printf(" ");
            if (damp < 1.0)
                printf("%9.4e ", damp);
            else
                for (i = 0; i < 11; i++) printf(" ");
            printf("%9.4e %9.4e", update_norm, resid_norm);
            if (do_time) {
                k = m_kinSpecIndex[label_t];
                int isp = m_kinObjIndex[label_t];
                InterfaceKinetics *m_kin = m_objects[isp];
                nm = m_kin->kineticsSpeciesName(k);
                printf(" %-16s", nm.c_str());
            } else {
                for (i = 0; i < 16; i++) printf(" ");
            }
            if (label_d >= 0) {
                k = m_kinSpecIndex[label_d];
                int isp = m_kinObjIndex[label_d];
                InterfaceKinetics *m_kin = m_objects[isp];
                nm = m_kin->kineticsSpeciesName(k);
                printf(" %-16s", nm.c_str());
            }
            printf("\n");
        }
#ifdef DEBUG_SOLVESP
        else if (ioflag > 1) {

            updateMFSolnSP(XMolSolnSP);
            printf("\n\t      Weighted norm of update = %10.4e\n", update_norm);

            printf("\t  Name            Prod_Rate        XMol       Conc   "
                "  Conc_Old     wtConc");
            if (damp < 1.0) printf(" UnDamped_Conc");
            printf("\n");
            printf("\t---------------------------------------------------------"
                "-----------------------------\n");
            int kindexSP = 0;
            for (int isp = 0; isp < m_numSurfPhases; isp++) {
                int nsp = m_nSpeciesSurfPhase[isp];
                InterfaceKinetics *m_kin = m_objects[isp];
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
                        if (label_d == kindexSP) printf(" Damp ");
                    }
                    if (label_t == kindexSP) printf(" Tctrl");
                    printf("\n");
                }

            }
  
            printf("\t--------------------------------------------------------"
                "------------------------------\n");
        }
#endif
    } /* printIteration */


    void solveSP::printFinal(int ioflag, double damp, int label_d, int label_t,
        double inv_t, double t_real, int iter,
        double update_norm, double resid_norm,
        double netProdRateKinSpecies[], const double CSolnSP[],
        const double resid[], double XMolSolnSP[],
        const double wtSpecies[], const double wtRes[],
        int dim, bool do_time,
        double TKelvin, double PGas)
    {
        int i, k;
        string nm;
        if (ioflag == 1) {

            printf("\tFIN%3d ", iter);
            if (do_time)
                printf("%9.4e %9.4e ", t_real, 1.0/inv_t);
            else
                for (i = 0; i < 22; i++) printf(" ");
            if (damp < 1.0)
                printf("%9.4e ", damp);
            else
                for (i = 0; i < 11; i++) printf(" ");
            printf("%9.4e %9.4e", update_norm, resid_norm);
            if (do_time) {
                k = m_kinSpecIndex[label_t];
                int isp = m_kinObjIndex[label_t];
                InterfaceKinetics *m_kin = m_objects[isp];
                nm = m_kin->kineticsSpeciesName(k);
                printf(" %-16s", nm.c_str());
            } else {
                for (i = 0; i < 16; i++) printf(" ");
            }
            if (label_d >= 0) {
                k = m_kinSpecIndex[label_d];
                int isp = m_kinObjIndex[label_d];
                InterfaceKinetics *m_kin = m_objects[isp];
                nm = m_kin->kineticsSpeciesName(k);
                printf(" %-16s", nm.c_str());
            }
            printf(" -- success\n");
        }
#ifdef DEBUG_SOLVESP
        else if (ioflag > 1) {

   
            printf("\n================================== FINAL RESULT ========="
                "==================================================\n");
            updateMFSolnSP(XMolSolnSP);
            printf("\n        Weighted norm of solution update = %10.4e\n", update_norm);
            printf("        Weighted norm of residual update = %10.4e\n\n", resid_norm);

            printf("  Name            Prod_Rate        XMol       Conc   "
                "  wtConc      Resid   Resid/wtResid   wtResid");
            if (damp < 1.0) printf(" UnDamped_Conc");
            printf("\n");
            printf("---------------------------------------------------------------"
                "---------------------------------------------\n");
            int kindexSP = 0;
            for (int isp = 0; isp < m_numSurfPhases; isp++) {
                int nsp = m_nSpeciesSurfPhase[isp];
                InterfaceKinetics *m_kin = m_objects[isp];
                //int surfPhaseIndex = m_kinObjPhaseIDSurfPhase[isp];
                m_kin->getNetProductionRates(DATA_PTR(m_numEqn1));
                for (int k = 0; k < nsp; k++, kindexSP++) {
                    int kspIndex = m_kinSpecIndex[kindexSP];
                    nm = m_kin->kineticsSpeciesName(kspIndex);
                    printf("%-16s  %10.3e   %10.3e  %10.3e  %10.3e %10.3e  %10.3e %10.3e",
                        nm.c_str(),
                        m_numEqn1[kspIndex], 
                        XMolSolnSP[kindexSP], 
                        CSolnSP[kindexSP],
                        wtSpecies[kindexSP],
                        resid[kindexSP],
                        resid[kindexSP]/wtRes[kindexSP], wtRes[kindexSP]);
                    if (damp < 1.0) {
                        printf("%10.4e ", CSolnSP[kindexSP]+(damp-1.0)*resid[kindexSP]);
                        if (label_d == kindexSP) printf(" Damp ");
                    }
                    if (label_t == kindexSP) printf(" Tctrl");
                    printf("\n");
                }

            }
            printf("---------------------------------------------------------------"
                "---------------------------------------------\n");
            double *XMolKinSpecies = DATA_PTR(m_numEqn2);
            kindexSP = 0;
            for (int isp = 0; isp < m_numSurfPhases; isp++) {
                InterfaceKinetics *m_kin = m_objects[isp];
                int surfIndex = m_kin->surfacePhaseIndex();
                int nPhases = m_kin->nPhases();
                m_kin->getNetProductionRates(netProdRateKinSpecies);

                updateMFKinSpecies(XMolKinSpecies, isp);

                printf("\n    IntefaceKinetics Object # %d\n\n", isp);

                printf("\t   Number of Phases = %d\n", nPhases);
                printf("\t   Temperature = %10.3e Kelvin\n", TKelvin);
                printf("\t   Pressure    = %10.3g Pa\n\n", PGas);
                printf("\t   Phase:SpecName      Prod_Rate  MoleFraction   kindexSP\n");
                printf("\t--------------------------------------------------------------"
                    "---\n");
      
                int kspindex = 0;
                bool inSurfacePhase = false;
                for (int ip = 0; ip < nPhases; ip++) {
                    if (ip == surfIndex) {
                        inSurfacePhase = true;
                    } else {
                        inSurfacePhase = false;
                    }
                    ThermoPhase &THref = m_kin->thermo(ip);
                    int nsp = THref.nSpecies();
                    string pname = THref.id();
                    for (k = 0; k < nsp; k++) {
                        string sname = THref.speciesName(k);
                        string cname = pname + ":" + sname;
                        if (inSurfacePhase) {
                            printf("\t%-24s %10.3e %10.3e      %d\n", cname.c_str(), 
                                netProdRateKinSpecies[kspindex], XMolKinSpecies[kspindex], kindexSP);
                            kindexSP++;
                        } else {
                            printf("\t%-24s %10.3e %10.3e\n", cname.c_str(), 
                                netProdRateKinSpecies[kspindex], XMolKinSpecies[kspindex]);
                        }
                        kspindex++;
                    }
                }
            }
            printf("\n");
            printf("==============================================================="
                "============================================\n\n");
        }
#endif
    }

#ifdef DEBUG_SOLVESP
    void solveSP::
    printIterationHeader(int ioflag, double damp,double inv_t, double t_real,
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
 
}
