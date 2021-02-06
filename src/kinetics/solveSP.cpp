/*
 * @file: solveSP.cpp Implicit surface site concentration solver
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/solveSP.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/kinetics/ImplicitSurfChem.h"

using namespace std;
namespace Cantera
{

// STATIC ROUTINES DEFINED IN THIS FILE

static doublereal calc_damping(doublereal* x, doublereal* dx, size_t dim, int*);
static doublereal calcWeightedNorm(const doublereal [], const doublereal dx[], size_t);

// solveSP Class Definitions

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
    for (size_t n = 0; n < m_objects.size(); n++) {
        InterfaceKinetics* kin = m_objects[n];
        size_t surfPhaseIndex = kin->surfacePhaseIndex();
        if (surfPhaseIndex != npos) {
            m_numSurfPhases++;
            m_indexKinObjSurfPhase.push_back(n);
            m_kinObjPhaseIDSurfPhase.push_back(surfPhaseIndex);
        } else {
            throw CanteraError("solveSP::solveSP",
                               "InterfaceKinetics object has no surface phase");
        }
        SurfPhase* sp = dynamic_cast<SurfPhase*>(&kin->thermo(surfPhaseIndex));
        if (!sp) {
            throw CanteraError("solveSP::solveSP",
                               "Inconsistent ThermoPhase object within "
                               "InterfaceKinetics object");
        }

        m_ptrsSurfPhase.push_back(sp);
        size_t nsp = sp->nSpecies();
        m_nSpeciesSurfPhase.push_back(nsp);
        m_numTotSurfSpecies += nsp;
    }
    // We rely on ordering to figure things out
    m_numBulkPhasesSS = 0;

    if (bulkFunc == BULK_DEPOSITION) {
        m_neq = m_numTotSurfSpecies + m_numTotBulkSpeciesSS;
    } else {
        m_neq = m_numTotSurfSpecies;
    }

    m_maxTotSpecies = 0;
    for (size_t n = 0; n < m_numSurfPhases; n++) {
        size_t tsp = m_objects[n]->nTotalSpecies();
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
    for (size_t isp = 0; isp < m_numSurfPhases; isp++) {
        size_t iKinObject = m_indexKinObjSurfPhase[isp];
        InterfaceKinetics* kin = m_objects[iKinObject];
        size_t surfPhaseIndex = m_kinObjPhaseIDSurfPhase[isp];
        size_t kstart = kin->kineticsSpeciesIndex(0, surfPhaseIndex);
        size_t nsp = m_nSpeciesSurfPhase[isp];
        m_eqnIndexStartSolnPhase[isp] = kindexSP;
        for (size_t k = 0; k < nsp; k++, kindexSP++) {
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
    int label_t=-1; // Species IDs for time control
    int label_d = -1; // Species IDs for damping control
    int label_t_old=-1;
    doublereal label_factor = 1.0;
    int iter=0; // iteration number on nonlinear solver
    int iter_max=1000; // maximum number of nonlinear iterations
    doublereal deltaT = 1.0E-10; // Delta time step
    doublereal damp=1.0;
    doublereal inv_t = 0.0;
    doublereal t_real = 0.0, update_norm = 1.0E6;
    bool do_time = false, not_converged = true;
    m_ioflag = std::min(m_ioflag, 1);

    // Set the initial value of the do_time parameter
    if (ifunc == SFLUX_INITIALIZE || ifunc == SFLUX_TRANSIENT) {
        do_time = true;
    }

    // Store the initial guess for the surface problem in the soln vector,
    // CSoln, and in an separate vector CSolnInit.
    size_t loc = 0;
    for (size_t n = 0; n < m_numSurfPhases; n++) {
        m_ptrsSurfPhase[n]->getConcentrations(m_numEqn1.data());
        for (size_t k = 0; k < m_nSpeciesSurfPhase[n]; k++) {
            m_CSolnSP[loc] = m_numEqn1[k];
            loc++;
        }
    }

    m_CSolnSPInit = m_CSolnSP;

    // Calculate the largest species in each phase
    evalSurfLarge(m_CSolnSP.data());

    if (m_ioflag) {
        print_header(m_ioflag, ifunc, time_scale, true, reltol, abstol);
    }

    // Quick return when there isn't a surface problem to solve
    if (m_neq == 0) {
        not_converged = false;
        update_norm = 0.0;
    }

    // Start of Newton's method
    while (not_converged && iter < iter_max) {
        iter++;
        // Store previous iteration's solution in the old solution vector
        m_CSolnSPOld = m_CSolnSP;

        // Evaluate the largest surface species for each surface phase every
        // 5 iterations.
        if (iter%5 == 4) {
            evalSurfLarge(m_CSolnSP.data());
        }

        // Calculate the value of the time step
        // - heuristics to stop large oscillations in deltaT
        if (do_time) {
            // don't hurry increase in time step at the same time as damping
            if (damp < 1.0) {
                label_factor = 1.0;
            }
            double tmp = calc_t(m_netProductionRatesSave.data(),
                         m_XMolKinSpecies.data(),
                         &label_t, &label_t_old, &label_factor, m_ioflag);
            if (iter < 10) {
                inv_t = tmp;
            } else if (tmp > 2.0*inv_t) {
                inv_t = 2.0*inv_t;
            } else {
                inv_t = tmp;
            }

            // Check end condition
            if (ifunc == SFLUX_TRANSIENT) {
                tmp = t_real + 1.0/inv_t;
                if (tmp > time_scale) {
                    inv_t = 1.0/(time_scale - t_real);
                }
            }
        } else {
            // make steady state calc a step of 1 million seconds to prevent
            // singular Jacobians for some pathological cases
            inv_t = 1.0e-6;
        }
        deltaT = 1.0/inv_t;

        // Call the routine to numerically evaluation the Jacobian and residual
        // for the current iteration.
        resjac_eval(m_Jac, m_resid.data(), m_CSolnSP.data(),
                    m_CSolnSPOld.data(), do_time, deltaT);

        // Calculate the weights. Make sure the calculation is carried out on
        // the first iteration.
        if (iter%4 == 1) {
            calcWeights(m_wtSpecies.data(), m_wtResid.data(),
                        m_Jac, m_CSolnSP.data(), abstol, reltol);
        }

        // Find the weighted norm of the residual
        double resid_norm = calcWeightedNorm(m_wtResid.data(), m_resid.data(), m_neq);

        // Solve Linear system.  The solution is in m_resid
        solve(m_Jac, m_resid.data());

        // Calculate the Damping factor needed to keep all unknowns between 0
        // and 1, and not allow too large a change (factor of 2) in any unknown.
        damp = calc_damping(m_CSolnSP.data(), m_resid.data(), m_neq, &label_d);

        // Calculate the weighted norm of the update vector Here, resid is the
        // delta of the solution, in concentration units.
        update_norm = calcWeightedNorm(m_wtSpecies.data(),
                                       m_resid.data(), m_neq);

        // Update the solution vector and real time Crop the concentrations to
        // zero.
        for (size_t irow = 0; irow < m_neq; irow++) {
            m_CSolnSP[irow] -= damp * m_resid[irow];
        }
        for (size_t irow = 0; irow < m_neq; irow++) {
            m_CSolnSP[irow] = std::max(0.0, m_CSolnSP[irow]);
        }
        updateState(m_CSolnSP.data());

        if (do_time) {
            t_real += damp/inv_t;
        }

        if (m_ioflag) {
            printIteration(m_ioflag, damp, label_d, label_t, inv_t, t_real, iter,
                           update_norm, resid_norm, do_time);
        }

        if (ifunc == SFLUX_TRANSIENT) {
            not_converged = (t_real < time_scale);
        } else {
            if (do_time) {
                if (t_real > time_scale ||
                        (resid_norm < 1.0e-7 &&
                         update_norm*time_scale/t_real < EXTRA_ACCURACY)) {
                    do_time = false;
                }
            } else {
                not_converged = ((update_norm > EXTRA_ACCURACY) ||
                                 (resid_norm  > EXTRA_ACCURACY));
            }
        }
    } // End of Newton's Method while statement

    // End Newton's method. If not converged, print error message and
    // recalculate sdot's at equal site fractions.
    if (not_converged && m_ioflag) {
        writelog("#$#$#$# Error in solveSP $#$#$#$ \n");
        writelogf("Newton iter on surface species did not converge, "
                  "update_norm = %e \n", update_norm);
        writelog("Continuing anyway\n");
    }

    // Decide on what to return in the solution vector. Right now, will always
    // return the last solution no matter how bad
    if (m_ioflag) {
        fun_eval(m_resid.data(), m_CSolnSP.data(), m_CSolnSPOld.data(),
                 false, deltaT);
        double resid_norm = calcWeightedNorm(m_wtResid.data(), m_resid.data(), m_neq);
        printIteration(m_ioflag, damp, label_d, label_t, inv_t, t_real, iter,
                       update_norm, resid_norm, do_time, true);
    }

    // Return with the appropriate flag
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
    InterfaceKinetics* kin = m_objects[isp];
    for (size_t iph = 0; iph < kin->nPhases(); iph++) {
        size_t ksi = kin->kineticsSpeciesIndex(0, iph);
        kin->thermo(iph).getMoleFractions(XMolKinSpecies + ksi);
    }
}

void solveSP::evalSurfLarge(const doublereal* CSolnSP)
{
    size_t kindexSP = 0;
    for (size_t isp = 0; isp < m_numSurfPhases; isp++) {
        doublereal Clarge = CSolnSP[kindexSP];
        m_spSurfLarge[isp] = 0;
        kindexSP++;
        for (size_t k = 1; k <  m_nSpeciesSurfPhase[isp]; k++, kindexSP++) {
            if (CSolnSP[kindexSP] > Clarge) {
                Clarge = CSolnSP[kindexSP];
                m_spSurfLarge[isp] = k;
            }
        }
    }
}

void solveSP::fun_eval(doublereal* resid, const doublereal* CSoln,
                       const doublereal* CSolnOld, const bool do_time,
                       const doublereal deltaT)
{
    size_t k;
    doublereal lenScale = 1.0E-9;
    if (m_numSurfPhases > 0) {
        // update the surface concentrations with the input surface
        // concentration vector
        updateState(CSoln);

        // Get the net production rates of all of the species in the
        // surface kinetics mechanism
        //
        // HKM Should do it here for all kinetics objects so that
        //     bulk will eventually work.
        if (do_time) {
            size_t kindexSP = 0;
            for (size_t isp = 0; isp < m_numSurfPhases; isp++) {
                size_t nsp = m_nSpeciesSurfPhase[isp];
                InterfaceKinetics* kinPtr = m_objects[isp];
                size_t surfIndex = kinPtr->surfacePhaseIndex();
                size_t kstart = kinPtr->kineticsSpeciesIndex(0, surfIndex);
                size_t kins = kindexSP;
                kinPtr->getNetProductionRates(m_netProductionRatesSave.data());
                for (k = 0; k < nsp; k++, kindexSP++) {
                    resid[kindexSP] =
                        (CSoln[kindexSP] - CSolnOld[kindexSP]) / deltaT
                        - m_netProductionRatesSave[kstart + k];
                }

                size_t kspecial = kins + m_spSurfLarge[isp];
                double sd = m_ptrsSurfPhase[isp]->siteDensity();
                resid[kspecial] = sd;
                for (k = 0; k < nsp; k++) {
                    resid[kspecial] -= CSoln[kins + k];
                }
            }
        } else {
            size_t kindexSP = 0;
            for (size_t isp = 0; isp < m_numSurfPhases; isp++) {
                size_t nsp = m_nSpeciesSurfPhase[isp];
                InterfaceKinetics* kinPtr = m_objects[isp];
                size_t surfIndex = kinPtr->surfacePhaseIndex();
                size_t kstart = kinPtr->kineticsSpeciesIndex(0, surfIndex);
                size_t kins = kindexSP;
                kinPtr->getNetProductionRates(m_netProductionRatesSave.data());
                for (k = 0; k < nsp; k++, kindexSP++) {
                    resid[kindexSP] = - m_netProductionRatesSave[kstart + k];
                }
                size_t kspecial = kins + m_spSurfLarge[isp];
                double sd = m_ptrsSurfPhase[isp]->siteDensity();
                resid[kspecial] = sd;
                for (k = 0; k < nsp; k++) {
                    resid[kspecial] -= CSoln[kins + k];
                }
            }
        }

        if (m_bulkFunc == BULK_DEPOSITION) {
            size_t kindexSP = m_numTotSurfSpecies;
            for (size_t isp = 0; isp < m_numBulkPhasesSS; isp++) {
                doublereal* XBlk = m_numEqn1.data();
                size_t nsp = m_nSpeciesSurfPhase[isp];
                size_t surfPhaseIndex = m_indexKinObjSurfPhase[isp];
                InterfaceKinetics* kin = m_objects[isp];
                double grRate = 0.0;
                size_t kstart = kin->kineticsSpeciesIndex(0, surfPhaseIndex);
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
                    //! @todo the appearance of k in this formula is suspicious
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

void solveSP::resjac_eval(DenseMatrix& jac,
                          doublereal resid[], doublereal CSoln[],
                          const doublereal CSolnOld[], const bool do_time,
                          const doublereal deltaT)
{
    size_t kColIndex = 0;
    // Calculate the residual
    fun_eval(resid, CSoln, CSolnOld, do_time, deltaT);
    // Now we will look over the columns perturbing each unknown.
    for (size_t jsp = 0; jsp < m_numSurfPhases; jsp++) {
        size_t nsp = m_nSpeciesSurfPhase[jsp];
        double sd = m_ptrsSurfPhase[jsp]->siteDensity();
        for (size_t kCol = 0; kCol < nsp; kCol++) {
            double cSave = CSoln[kColIndex];
            double dc = std::max(1.0E-10 * sd, fabs(cSave) * 1.0E-7);
            CSoln[kColIndex] += dc;
            fun_eval(m_numEqn2.data(), CSoln, CSolnOld, do_time, deltaT);
            for (size_t i = 0; i < m_neq; i++) {
                jac(i, kColIndex) = (m_numEqn2[i] - resid[i])/dc;
            }
            CSoln[kColIndex] = cSave;
            kColIndex++;
        }
    }

    if (m_bulkFunc == BULK_DEPOSITION) {
        for (size_t jsp = 0; jsp < m_numBulkPhasesSS; jsp++) {
            size_t nsp = m_numBulkSpecies[jsp];
            double sd = m_bulkPhasePtrs[jsp]->molarDensity();
            for (size_t kCol = 0; kCol < nsp; kCol++) {
                double cSave = CSoln[kColIndex];
                double dc = std::max(1.0E-10 * sd, fabs(cSave) * 1.0E-7);
                CSoln[kColIndex] += dc;
                fun_eval(m_numEqn2.data(), CSoln, CSolnOld, do_time, deltaT);
                for (size_t i = 0; i < m_neq; i++) {
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
    doublereal damp = 1.0;
    static doublereal damp_old = 1.0; //! @todo this variable breaks thread safety
    *label = -1;

    for (size_t i = 0; i < dim; i++) {
        // Calculate the new suggested new value of x[i]
        double xnew = x[i] - damp * dxneg[i];

        // Calculate the allowed maximum and minimum values of x[i]
        //  - Only going to allow x[i] to converge to zero by a
        //    single order of magnitude at a time
        double xtop = 1.0 - 0.1*fabs(1.0-x[i]);
        double xbot = fabs(x[i]*0.1) - 1.0e-16;
        if (xnew > xtop)  {
            damp = - APPROACH * (1.0 - x[i]) / dxneg[i];
            *label = int(i);
        } else if (xnew < xbot) {
            damp = APPROACH * x[i] / dxneg[i];
            *label = int(i);
        } else if (xnew > 3.0*std::max(x[i], 1.0E-10)) {
            damp = - 2.0 * std::max(x[i], 1.0E-10) / dxneg[i];
            *label = int(i);
        }
    }
    damp = std::max(damp, 1e-2);

    // Only allow the damping parameter to increase by a factor of three each
    // iteration. Heuristic to avoid oscillations in the value of damp
    if (damp > damp_old*3) {
        damp = damp_old*3;
        *label = -1;
    }

    // Save old value of the damping parameter for use in subsequent calls.
    damp_old = damp;
    return damp;

} /* calc_damping */

/*
 *    This function calculates the norm of an update, dx[], based on the
 *    weighted values of x.
 */
static doublereal calcWeightedNorm(const doublereal wtX[], const doublereal dx[], size_t dim)
{
    doublereal norm = 0.0;
    if (dim == 0) {
        return 0.0;
    }
    for (size_t i = 0; i < dim; i++) {
        norm += pow(dx[i] / wtX[i], 2);
    }
    return sqrt(norm/dim);
}

void solveSP::calcWeights(doublereal wtSpecies[], doublereal wtResid[],
                          const Array2D& Jac, const doublereal CSoln[],
                          const doublereal abstol, const doublereal reltol)
{
    // First calculate the weighting factor for the concentrations of the
    // surface species and bulk species.
    size_t kindex = 0;
    for (size_t isp = 0; isp < m_numSurfPhases; isp++) {
        double sd = m_ptrsSurfPhase[isp]->siteDensity();
        for (size_t k = 0; k < m_nSpeciesSurfPhase[isp]; k++, kindex++) {
            wtSpecies[kindex] = abstol * sd + reltol * fabs(CSoln[kindex]);
        }
    }
    if (m_bulkFunc == BULK_DEPOSITION) {
        for (size_t isp = 0; isp < m_numBulkPhasesSS; isp++) {
            double sd = m_bulkPhasePtrs[isp]->molarDensity();
            for (size_t k = 0; k < m_numBulkSpecies[isp]; k++, kindex++) {
                wtSpecies[kindex] = abstol * sd + reltol * fabs(CSoln[kindex]);
            }
        }
    }

    // Now do the residual Weights. Since we have the Jacobian, we will use it
    // to generate a number based on the what a significant change in a solution
    // variable does to each residual. This is a row sum scale operation.
    for (size_t k = 0; k < m_neq; k++) {
        wtResid[k] = 0.0;
        for (size_t jcol = 0; jcol < m_neq; jcol++) {
            wtResid[k] += fabs(Jac(k,jcol) * wtSpecies[jcol]);
        }
    }
}

doublereal solveSP::calc_t(doublereal netProdRateSolnSP[],
                          doublereal XMolSolnSP[],
                          int* label, int* label_old,
                          doublereal* label_factor, int ioflag)
{
    doublereal inv_timeScale = 1.0E-10;
    size_t kindexSP = 0;
    *label = 0;
    updateMFSolnSP(XMolSolnSP);
    for (size_t isp = 0; isp < m_numSurfPhases; isp++) {
        // Get the interface kinetics associated with this surface
        InterfaceKinetics* kin = m_objects[isp];

        // Calculate the start of the species index for surfaces within
        // the InterfaceKinetics object
        size_t surfIndex = kin->surfacePhaseIndex();
        size_t kstart = kin->kineticsSpeciesIndex(0, surfIndex);
        kin->getNetProductionRates(m_numEqn1.data());
        double sden = kin->thermo(surfIndex).molarDensity();
        for (size_t k = 0; k < m_nSpeciesSurfPhase[isp]; k++, kindexSP++) {
            size_t kspindex = kstart + k;
            netProdRateSolnSP[kindexSP] = m_numEqn1[kspindex];
            double tmp = std::max(XMolSolnSP[kindexSP], 1.0e-10);
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

    // Increase time step exponentially as same species repeatedly controls time
    // step
    if (*label == *label_old) {
        *label_factor *= 1.5;
    } else {
        *label_old = *label;
        *label_factor = 1.0;
    }
    return inv_timeScale / *label_factor;
} // calc_t

void solveSP::print_header(int ioflag, int ifunc, doublereal time_scale,
                           int damping, doublereal reltol, doublereal abstol)
{
    if (ioflag) {
        writelog("\n================================ SOLVESP CALL SETUP "
                 "========================================\n");
        if (ifunc == SFLUX_INITIALIZE) {
            writelog("\n  SOLVESP Called with Initialization turned on\n");
            writelogf("     Time scale input = %9.3e\n", time_scale);
        } else if (ifunc == SFLUX_RESIDUAL) {
            writelog("\n   SOLVESP Called to calculate steady state residual\n");
            writelog("           from a good initial guess\n");
        } else if (ifunc == SFLUX_JACOBIAN)  {
            writelog("\n   SOLVESP Called to calculate steady state Jacobian\n");
            writelog("           from a good initial guess\n");
        } else if (ifunc == SFLUX_TRANSIENT) {
            writelog("\n   SOLVESP Called to integrate surface in time\n");
            writelogf("           for a total of %9.3e sec\n", time_scale);
        } else {
            throw CanteraError("solveSP::print_header",
                               "Unknown ifunc flag = {}", ifunc);
        }

        if (m_bulkFunc == BULK_DEPOSITION) {
            writelog("     The composition of the Bulk Phases will be calculated\n");
        } else if (m_bulkFunc == BULK_ETCH) {
            writelog("     Bulk Phases have fixed compositions\n");
        } else {
            throw CanteraError("solveSP::print_header",
                               "Unknown bulkFunc flag = {}", m_bulkFunc);
        }

        if (damping) {
            writelog("     Damping is ON   \n");
        } else {
            writelog("     Damping is OFF  \n");
        }

        writelogf("     Reltol = %9.3e, Abstol = %9.3e\n", reltol, abstol);
    }

    if (ioflag == 1) {
        writelog("\n\n\t Iter    Time       Del_t      Damp      DelX   "
                 "     Resid    Name-Time    Name-Damp\n");
        writelog("\t -----------------------------------------------"
                 "------------------------------------\n");
    }
}

void solveSP::printIteration(int ioflag, doublereal damp, int label_d,
                             int label_t, doublereal inv_t, doublereal t_real,
                             size_t iter, doublereal update_norm,
                             doublereal resid_norm, bool do_time, bool final)
{
    if (ioflag == 1) {
        if (final) {
            writelogf("\tFIN%3d ", iter);
        } else {
            writelogf("\t%6d ", iter);
        }
        if (do_time) {
            writelogf("%9.4e %9.4e ", t_real, 1.0/inv_t);
        } else {
            writeline(' ', 22, false);
        }
        if (damp < 1.0) {
            writelogf("%9.4e ", damp);
        } else {
            writeline(' ', 11, false);
        }
        writelogf("%9.4e %9.4e", update_norm, resid_norm);
        if (do_time) {
            size_t k = m_kinSpecIndex[label_t];
            size_t isp = m_kinObjIndex[label_t];
            writelog(" %-16s", m_objects[isp]->kineticsSpeciesName(k));
        } else {
            writeline(' ', 16, false);
        }
        if (label_d >= 0) {
            size_t k = m_kinSpecIndex[label_d];
            size_t isp = m_kinObjIndex[label_d];
            writelogf(" %-16s", m_objects[isp]->kineticsSpeciesName(k));
        }
        if (final) {
            writelog(" -- success");
        }
        writelog("\n");
    }
} // printIteration

}
