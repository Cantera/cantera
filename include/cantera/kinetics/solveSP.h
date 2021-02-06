/**
 * @file solveSP.h Header file for implicit surface problem solver (see \ref
 *       chemkinetics and class \link Cantera::solveSP solveSP\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef SOLVESP_H
#define SOLVESP_H

#include "cantera/kinetics/InterfaceKinetics.h"
#include "cantera/numerics/DenseMatrix.h"

//! @defgroup solvesp_methods Surface Problem Solver Methods
//! @{

//! This assumes that the initial guess supplied to the routine is far from
//! the correct one. Substantial work plus transient time-stepping is to be
//! expected to find a solution.
const int SFLUX_INITIALIZE = 1;

//! Need to solve the surface problem in order to calculate the surface fluxes
//! of gas-phase species. (Can expect a moderate change in the solution
//! vector; try to solve the system by direct methods with no damping first,
//! then try time-stepping if the first method fails). A "time_scale" supplied
//! here is used in the algorithm to determine when to shut off time-stepping.
const int SFLUX_RESIDUAL = 2;

//! Calculation of the surface problem is due to the need for a numerical
//! Jacobian for the gas-problem. The solution is expected to be very close to
//! the initial guess, and accuracy is needed because solution variables have
//! been perturbed from nominal values to create Jacobian entries.
const int SFLUX_JACOBIAN = 3;

//! The transient calculation is performed here for an amount of time
//! specified by "time_scale".  It is not guaranteed to be time-accurate -
//! just stable and fairly fast. The solution after del_t time is returned,
//! whether it's converged to a steady state or not. This is a poor man's time
//! stepping algorithm.
const int SFLUX_TRANSIENT = 4;
// @}

//! @defgroup solvesp_bulkFunc Surface Problem Bulk Phase Mode
//! Functionality expected from the bulk phase. This changes the equations
//! that will be used to solve for the bulk mole fractions.
//! @{

//! Deposition of a bulk phase is to be expected. Bulk mole fractions are
//! determined from ratios of growth rates of bulk species.
const int BULK_DEPOSITION = 1;

//! Etching of a bulk phase is to be expected. Bulk mole fractions are assumed
//! constant, and given by the initial conditions. This is also used whenever
//! the condensed phase is part of the larger solution.
const int BULK_ETCH = 2;
// @}

namespace Cantera
{

//! Method to solve a pseudo steady state surface problem
/*!
 *  The following class handles solving the surface problem. The calculation
 *  uses Newton's method to obtain the surface fractions of the surface and
 *  bulk species by requiring that the surface species production rate = 0 and
 *  that the either the bulk fractions are proportional to their production
 *  rates or they are constants.
 *
 *  Currently, the bulk mole fractions are treated as constants. Implementation
 *  of their being added to the unknown solution vector is delayed.
 *
 *  Lets introduce the unknown vector for the "surface problem". The surface
 *  problem is defined as the evaluation of the surface site fractions for
 *  multiple surface phases. The unknown vector will consist of the vector of
 *  surface concentrations for each species in each surface vector. Species
 *  are grouped first by their surface phases
 *
 *  - C_i_j = Concentration of the ith species in the jth surface phase
 *  - Nj = number of surface species in the jth surface phase
 *
 *  The unknown solution vector is defined as follows:
 *
 *  C_i_j     | kindexSP
 *  --------- | ----------
 *  C_0_0     |   0
 *  C_1_0     |   1
 *  C_2_0     |   2
 *   . . .    |  ...
 *  C_N0-1_0  | N0-1
 *  C_0_1     | N0
 *  C_1_1     | N0+1
 *  C_2_1     | N0+2
 *   . . .    | ...
 *  C_N1-1_1  | NO+N1-1
 *
 *  Note there are a couple of different types of species indices floating
 *  around in the formulation of this object.
 *
 *  kindexSP: This is the species index in the contiguous vector of unknowns
 *            for the surface problem.
 *
 *  Note, in the future, BULK_DEPOSITION systems will be added, and the
 *  solveSP unknown vector will get more complicated. It will include the mole
 *  fraction and growth rates of specified bulk phases
 *
 *  Indices which relate to individual kinetics objects use the suffix KSI
 *  (kinetics species index).
 *
 *  ## Solution Method
 *
 *  This routine is typically used within a residual calculation in a large code.
 *  It's typically invoked millions of times for large calculations, and it must
 *  work every time. Therefore, requirements demand that it be robust but also
 *  efficient.
 *
 *  The solution methodology is largely determined by the `ifunc` parameter,
 *  that is input to the solution object. This parameter may have one of the
 *  values defined in @ref solvesp_methods.
 *
 *  ### Pseudo time stepping algorithm:
 *  The time step is determined from sdot[], so that the time step
 *  doesn't ever change the value of a variable by more than 100%.
 *
 *  This algorithm does use a damped Newton's method to relax the equations.
 *  Damping is based on a "delta damping" technique. The solution unknowns
 *  are not allowed to vary too much between iterations.
 *
 *  `EXTRA_ACCURACY`: A constant that is the ratio of the required update norm
 *  in this Newton iteration compared to that in the nonlinear solver. A value
 *  of 0.1 is used so surface species are safely overconverged.
 *
 *  Functions called:
 *  - `ct_dgetrf` -- First half of LAPACK direct solve of a full Matrix
 *  - `ct_dgetrs` -- Second half of LAPACK direct solve of a full matrix.
 *    Returns solution vector in the right-hand-side vector, resid.
 */
class solveSP
{
public:
    //! Constructor for the object
    /*!
     *  @param surfChemPtr  Pointer to the ImplicitSurfChem object that
     *                      defines the surface problem to be solved.
     *  @param bulkFunc     Integer representing how the bulk phases should be
     *                      handled. See @ref solvesp_bulkFunc. Currently,
     *                      only the default value of BULK_ETCH is supported.
     */
    solveSP(ImplicitSurfChem* surfChemPtr, int bulkFunc = BULK_ETCH);

    //! Destructor. Deletes the integrator.
    ~solveSP() {}

private:
    //! Unimplemented private copy constructor
    solveSP(const solveSP& right);

    //! Unimplemented private assignment operator
    solveSP& operator=(const solveSP& right);

public:
    //! Main routine that actually calculates the pseudo steady state
    //! of the surface problem
    /*!
     * The actual converged solution is returned as part of the internal state
     * of the InterfaceKinetics objects.
     *
     * Uses Newton's method to get the surface fractions of the surface and
     * bulk species by requiring that the surface species production rate = 0
     * and that the bulk fractions are proportional to their production rates.
     *
     * @param ifunc Determines the type of solution algorithm to be used. See
     *                  @ref solvesp_methods for possible values.
     * @param time_scale  Time over which to integrate the surface equations,
     *                    where applicable
     * @param TKelvin     Temperature (kelvin)
     * @param PGas        Pressure (pascals)
     * @param reltol      Relative tolerance to use
     * @param abstol      absolute tolerance.
     * @return  1 if the surface problem is successfully solved.
     *          -1 if the surface problem wasn't solved successfully.
     *          Note the actual converged solution is returned as part of the
     *          internal state of the InterfaceKinetics objects.
     */
    int solveSurfProb(int ifunc, doublereal time_scale, doublereal TKelvin,
                      doublereal PGas, doublereal reltol, doublereal abstol);

private:
    //! Printing routine that optionally gets called at the start of every
    //! invocation
    void print_header(int ioflag, int ifunc, doublereal time_scale,
                      int damping, doublereal reltol, doublereal abstol);

    //! Printing routine that gets called after every iteration
    void printIteration(int ioflag, doublereal damp, int label_d, int label_t,
                        doublereal inv_t, doublereal t_real, size_t iter,
                        doublereal update_norm, doublereal resid_norm,
                        bool do_time, bool final=false);

    //! Calculate a conservative delta T to use in a pseudo-steady state
    //! algorithm
    /*!
     *  This routine calculates a pretty conservative 1/del_t based on
     *  MAX_i(sdot_i/(X_i*SDen0)). This probably guarantees diagonal dominance.
     *
     *  Small surface fractions are allowed to intervene in the del_t
     *  determination, no matter how small.  This may be changed.
     *  Now minimum changed to 1.0e-12,
     *
     *  Maximum time step set to time_scale.
     *
     *  @param netProdRateSolnSP  Output variable. Net production rate of all
     *      of the species in the solution vector.
     *  @param XMolSolnSP output variable. Mole fraction of all of the species
     *      in the solution vector
     *  @param label Output variable. Pointer to the value of the species
     *      index (kindexSP) that is controlling the time step
     *  @param label_old Output variable. Pointer to the value of the species
     *      index (kindexSP) that controlled the time step at the previous
     *      iteration
     *  @param label_factor Output variable. Pointer to the current factor
     *      that is used to indicate the same species is controlling the time
     *      step.
     * @param ioflag Level of the output requested.
     * @returns the 1. /  delta T to be used on the next step
     */
    doublereal calc_t(doublereal netProdRateSolnSP[], doublereal XMolSolnSP[],
                      int* label, int* label_old,
                      doublereal* label_factor, int ioflag);

    //! Calculate the solution and residual weights
    /*!
     *  @param wtSpecies Weights to use for the soln unknowns. These are in
     *      concentration units
     *  @param wtResid    Weights to sue for the residual unknowns.
     *  @param Jac        Jacobian. Row sum scaling is used for the Jacobian
     *  @param CSolnSP    Solution vector for the surface problem
     *  @param abstol     Absolute error tolerance
     *  @param reltol     Relative error tolerance
     */
    void calcWeights(doublereal wtSpecies[], doublereal wtResid[],
                     const Array2D& Jac, const doublereal CSolnSP[],
                     const doublereal abstol, const doublereal reltol);

    /**
     * Update the surface states of the surface phases.
     */
    void updateState(const doublereal* cSurfSpec);

    //! Update mole fraction vector consisting of unknowns in surface problem
    /*!
     * @param XMolSolnSP  Vector of mole fractions for the unknowns in the
     *                    surface problem.
     */
    void updateMFSolnSP(doublereal* XMolSolnSP);

    //! Update the mole fraction vector for a specific kinetic species vector
    //! corresponding to one InterfaceKinetics object
    /*!
     * @param XMolKinSp Mole fraction vector corresponding to a particular
     *                  kinetic species for a single InterfaceKinetics Object
     *                  This is a vector over all the species in all of the
     *                  phases in the InterfaceKinetics object
     * @param isp       ID of the InterfaceKinetics Object.
     */
    void updateMFKinSpecies(doublereal* XMolKinSp, int isp);

    //! Update the vector that keeps track of the largest species in each
    //! surface phase.
    /*!
     * @param CSolnSP Vector of the current values of the surface concentrations
     *                in all of the surface species.
     */
    void evalSurfLarge(const doublereal* CSolnSP);

    //! Main Function evaluation
    /*!
     *  @param resid output Vector of residuals, length = m_neq
     *  @param CSolnSP  Vector of species concentrations, unknowns in the
     *                  problem, length = m_neq
     *  @param CSolnOldSP Old Vector of species concentrations, unknowns in the
     *                  problem, length = m_neq
     *  @param do_time Calculate a time dependent residual
     *  @param deltaT  Delta time for time dependent problem.
     */
    void fun_eval(doublereal* resid, const doublereal* CSolnSP,
                  const doublereal* CSolnOldSP, const bool do_time, const doublereal deltaT);

    //! Main routine that calculates the current residual and Jacobian
    /*!
     *  @param jac     Jacobian to be evaluated.
     *  @param resid   output Vector of residuals, length = m_neq
     *  @param CSolnSP  Vector of species concentrations, unknowns in the
     *                  problem, length = m_neq. These are tweaked in order
     *                  to derive the columns of the Jacobian.
     *  @param CSolnSPOld Old Vector of species concentrations, unknowns in the
     *                  problem, length = m_neq
     *  @param do_time Calculate a time dependent residual
     *  @param deltaT  Delta time for time dependent problem.
     */
    void resjac_eval(DenseMatrix& jac, doublereal* resid,
                     doublereal* CSolnSP,
                     const doublereal* CSolnSPOld, const bool do_time,
                     const doublereal deltaT);

    //! Pointer to the manager of the implicit surface chemistry problem
    /*!
     *  This object actually calls the current object. Thus, we are providing a
     *  loop-back functionality here.
     */
    ImplicitSurfChem* m_SurfChemPtr;

    //! Vector of interface kinetics objects
    /*!
     * Each of these is associated with one and only one surface phase.
     */
    std::vector<InterfaceKinetics*> &m_objects;

    //! Total number of equations to solve in the implicit problem.
    /*!
     * Note, this can be zero, and frequently is
     */
    size_t m_neq;

    //! This variable determines how the bulk phases are to be handled
    /*!
     *  Possible values are given in @ref solvesp_bulkFunc.
     */
    int m_bulkFunc;

    //! Number of surface phases in the surface problem
    /*!
     * This number is equal to the number of InterfaceKinetics objects
     * in the problem. (until further noted)
     */
    size_t m_numSurfPhases;

    //! Total number of surface species in all surface phases.
    /*!
     * This is also the number of equations to solve for m_mode=0 system
     * It's equal to the sum of the number of species in each of the
     * m_numSurfPhases.
     */
    size_t m_numTotSurfSpecies;

    //! Mapping between the surface phases and the InterfaceKinetics objects
    /*!
     *  Currently this is defined to be a 1-1 mapping (and probably assumed
     *  in some places)
     *  m_surfKinObjID[i] = i
     */
    std::vector<size_t> m_indexKinObjSurfPhase;

    //! Vector of length number of surface phases containing
    //! the number of surface species in each phase
    /*!
     *  Length is equal to the number of surface phases, m_numSurfPhases
     */
    std::vector<size_t> m_nSpeciesSurfPhase;

    //! Vector of surface phase pointers
    /*!
     *  This is created during the constructor
     *  Length is equal to the number of surface phases, m_numSurfPhases
     */
    std::vector<SurfPhase*> m_ptrsSurfPhase;

    //! Index of the start of the unknowns for each solution phase
    /*!
     *        i_eqn = m_eqnIndexStartPhase[isp]
     *
     *  isp is the phase id in the list of phases solved by the
     *  surface problem.
     *
     *  i_eqn is the equation number of the first unknown in the
     *  solution vector corresponding to isp'th phase.
     */
    std::vector<size_t> m_eqnIndexStartSolnPhase;

    //! Phase ID in the InterfaceKinetics object of the surface phase
    /*!
     *  For each surface phase, this lists the PhaseId of the
     *  surface phase in the corresponding InterfaceKinetics object
     *
     * Length is equal to m_numSurfPhases
     */
    std::vector<size_t> m_kinObjPhaseIDSurfPhase;

    //! Total number of volumetric condensed phases included in the steady state
    //! problem handled by this routine.
    /*!
     * This is equal to or less than the total number of volumetric phases in
     * all of the InterfaceKinetics objects. We usually do not include bulk
     * phases. Bulk phases are only included in the calculation when their
     * domain isn't included in the underlying continuum model conservation
     * equation system.
     *
     * This is equal to 0, for the time being
     */
    size_t m_numBulkPhasesSS;

    //! Vector of number of species in the m_numBulkPhases phases.
    /*!
     * Length is number of bulk phases
     */
    std::vector<size_t> m_numBulkSpecies;

    //! Total number of species in all bulk phases.
    /*!
     *  This is also the number of bulk equations to solve when bulk equation
     *  solving is turned on.
     */
    size_t m_numTotBulkSpeciesSS;

    //! Vector of bulk phase pointers, length is equal to m_numBulkPhases.
    std::vector<ThermoPhase*> m_bulkPhasePtrs;

    //! Index between the equation index and the position in the kinetic
    //! species array for the appropriate kinetics operator
    /*!
     *  Length = m_neq.
     *
     *  ksp = m_kinSpecIndex[ieq]
     *  ksp is the kinetic species index for the ieq'th equation.
     */
    std::vector<size_t> m_kinSpecIndex;

    //! Index between the equation index and the index of the
    //! InterfaceKinetics object
    /*!
     *   Length m_neq
     */
    std::vector<size_t> m_kinObjIndex;

    //! Vector containing the indices of the largest species
    //! in each surface phase
    /*!
     *  `k = m_spSurfLarge[i]` where `k` is the local species index, i.e., it
     *  varies from 0 to (num species in phase - 1) and `i` is the surface
     *  phase index in the problem. Length is equal to #m_numSurfPhases.
     */
    std::vector<size_t> m_spSurfLarge;

    //! The absolute tolerance in real units. units are (kmol/m2)
    doublereal m_atol;

    //! The relative error tolerance.
    doublereal m_rtol;

    //! maximum value of the time step. units = seconds
    doublereal m_maxstep;

    //! Maximum number of species in any single kinetics operator
    //! -> also maxed wrt the total # of solution species
    size_t m_maxTotSpecies;

    //! Temporary vector with length equal to max m_maxTotSpecies
    vector_fp m_netProductionRatesSave;

    //! Temporary vector with length equal to max m_maxTotSpecies
    vector_fp m_numEqn1;

    //! Temporary vector with length equal to max m_maxTotSpecies
    vector_fp m_numEqn2;

    //! Temporary vector with length equal to max m_maxTotSpecies
    vector_fp m_CSolnSave;

    //! Solution vector. length MAX(1, m_neq)
    vector_fp m_CSolnSP;

    //! Saved initial solution vector. length MAX(1, m_neq)
    vector_fp m_CSolnSPInit;

    //! Saved  solution vector at the old time step. length MAX(1, m_neq)
    vector_fp m_CSolnSPOld;

    //! Weights for the residual norm calculation. length MAX(1, m_neq)
    vector_fp m_wtResid;

    //! Weights for the species concentrations norm calculation
    /*!
     * length MAX(1, m_neq)
     */
    vector_fp m_wtSpecies;

    //!  Residual for the surface problem
    /*!
     *  The residual vector of length "dim" that, that has the value of "sdot"
     *  for surface species.  The residuals for the bulk species are a function
     *  of the sdots for all species in the bulk phase. The last residual of
     *  each phase enforces {Sum(fractions) = 1}. After linear solve (dgetrf_ &
     *  dgetrs_), resid holds the update vector.
     *
     * length MAX(1, m_neq)
     */
    vector_fp m_resid;

    //! Vector of mole fractions. length m_maxTotSpecies
    vector_fp m_XMolKinSpecies;

    //! Jacobian. m_neq by m_neq computed Jacobian matrix for the local
    //! Newton's method.
    DenseMatrix m_Jac;

public:
    int m_ioflag;
};
}
#endif
