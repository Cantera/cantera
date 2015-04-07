/**
 * @file solveProb.h Header file for implicit nonlinear solver with the option
 *       of a pseudotransient (see \ref numerics and class \link
 *       Cantera::solveProb solveProb\endlink).
 */

/*
 * Copyright 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#ifndef SOLVEPROB_H
#define SOLVEPROB_H
/**
 * @defgroup solverGroup Solvers for Equation Systems
 */

#include "cantera/base/Array.h"
#include "ResidEval.h"

//! Solution Methods
/*!
 * Flag to specify the solution method
 *
 *  1: SOLVEPROB_INITIALIZE   = This assumes that the initial guess supplied to the
 *                          routine is far from the correct one. Substantial
 *                          work plus transient time-stepping is to be expected
 *                          to find a solution.
 *  2:  SOLVEPROB_RESIDUAL    = Need to solve the surface problem in order to
 *                          calculate the surface fluxes of gas-phase species.
 *                          (Can expect a moderate change in the solution
 *                           vector -> try to solve the system by direct
 *                            methods
 *                           with no damping first -> then, try time-stepping
 *                           if the first method fails)
 *                          A "time_scale" supplied here is used in the
 *                          algorithm to determine when to shut off
 *                          time-stepping.
 *  3:  SOLVEPROB_JACOBIAN    = Calculation of the surface problem is due to the
 *                          need for a numerical jacobian for the gas-problem.
 *                          The solution is expected to be very close to the
 *                          initial guess, and accuracy is needed.
 *  4:  SOLVEPROB_TRANSIENT   = The transient calculation is performed here for an
 *                          amount of time specified by "time_scale".  It is
 *                          not guaranteed to be time-accurate - just stable
 *                          and fairly fast. The solution after del_t time is
 *                          returned, whether it's converged to a steady
 *                          state or not.
 */
const int SOLVEPROB_INITIALIZE = 1;
const int SOLVEPROB_RESIDUAL   = 2;
const int SOLVEPROB_JACOBIAN   = 3;
const int SOLVEPROB_TRANSIENT  = 4;

namespace Cantera
{
//! Method to solve a pseudo steady state of a nonlinear problem
/*!
 *   The following class handles the solution of a nonlinear problem.
 *
 *      Res_ss(C) =   - Res(C) = 0
 *
 *   Optionally a pseudo transient algorithm may be used to relax the residual if
 *   it is available.
 *
 *      Res_td(C) =  dC/dt  - Res(C) = 0;
 *
 *    Res_ss(C) is the steady state residual to be solved. Res_td(C) is the
 *    time dependent residual which leads to the steady state residual.
 *
 *  Solution Method
 *
 *  This routine is typically used within a residual calculation in a large code.
 *  It's typically invoked millions of times for large calculations, and it must
 *  work every time. Therefore, requirements demand that it be robust but also
 *  efficient.
 *
 *  The solution methodology is largely determined by the <TT>ifunc</TT> parameter,
 *  that is input to the solution object. This parameter may have the following
 *  4 values:
 *
 *  1: SOLVEPROB_INITIALIZE   = This assumes that the initial guess supplied to the
 *                          routine is far from the correct one. Substantial
 *                          work plus transient time-stepping is to be expected
 *                          to find a solution.
 *
 *  2:  SOLVEPROB_RESIDUAL    = Need to solve the nonlinear problem in order to
 *                          calculate quantities for a residual calculation
 *                          (Can expect a moderate change in the solution
 *                           vector -> try to solve the system by direct methods
 *                           with no damping first -> then, try time-stepping
 *                           if the first method fails)
 *                          A "time_scale" supplied here is used in the
 *                          algorithm to determine when to shut off
 *                          time-stepping.
 *
 *  3:  SOLVEPROB_JACOBIAN    = Calculation of the surface problem is due to the
 *                          need for a numerical jacobian for the gas-problem.
 *                          The solution is expected to be very close to the
 *                          initial guess, and extra accuracy is needed because
 *                          solution variables have been delta'd from
 *                          nominal values to create jacobian entries.
 *
 *  4: SOLVEPROB_TRANSIENT   = The transient calculation is performed here for an
 *                          amount of time specified by "time_scale".  It is
 *                          not guaranteed to be time-accurate - just stable
 *                          and fairly fast. The solution after del_t time is
 *                          returned, whether it's converged to a steady
 *                          state or not. This is a poor man's time stepping
 *                          algorithm.
 *
 * Pseudo time stepping algorithm:
 *  The time step is determined from sdot[],  so that the time step
 *   doesn't ever change the value of a variable by more than 100%.
 *
 *  This algorithm does use a damped Newton's method to relax the equations.
 *  Damping is based on a "delta damping" technique. The solution unknowns
 *  are not allowed to vary too much between iterations.
 *
 *   EXTRA_ACCURACY:A constant that is the ratio of the required update norm in
 *    this Newton iteration compared to that in the nonlinear solver.
 *     A value of 0.1 is used so surface species are safely  overconverged.
 *
 *  Functions called:
 *----------------------------------------------------------------------------
 *
 * ct_dgetrf    -- First half of LAPACK direct solve of a full Matrix
 *
 * ct_dgetrs    -- Second half of LAPACK direct solve of a full matrix. Returns
 *                 solution vector in the right-hand-side vector, resid.
 *
 *----------------------------------------------------------------------------
 *
 *  @ingroup solverGroup
 */
class solveProb
{
public:

    //! Constructor for the object
    solveProb(ResidEval* resid);

    virtual ~solveProb();

private:

    //! Unimplemented private copy constructor
    solveProb(const solveProb& right);

    //! Unimplemented private assignment operator
    solveProb& operator=(const solveProb& right);

public:
    //! Main routine that actually calculates the pseudo steady state
    //! of the surface problem
    /*!
     *   The actual converged solution is returned as part of the
     *   internal state of the InterfaceKinetics objects.
     *
     * @param ifunc Determines the type of solution algorithm to be
     *                  used.  Possible values are  SOLVEPROB_INITIALIZE  ,
     *                  SOLVEPROB_RESIDUAL SOLVEPROB_JACOBIAN  SOLVEPROB_TRANSIENT   .
     *
     * @param time_scale  Time over which to integrate the surface equations,
     *                    where applicable
     *
     * @param reltol      Relative tolerance to use
     *
     * @return  Returns 1 if the surface problem is successfully solved.
     *          Returns -1 if the surface problem wasn't solved successfully.
     *          Note the actual converged solution is returned as part of the
     *          internal state of the InterfaceKinetics objects.
     */
    int solve(int ifunc, doublereal time_scale,  doublereal reltol);

    //! Report the current state of the solution
    /*!
     *  @param[out] CSoln solution vector for the nonlinear problem
     */
    virtual void reportState(doublereal* const CSoln) const;

    //! Set the bottom and top bounds on the solution vector
    /*!
     *  The default is for the bottom is 0.0, while the default for the top is 1.0
     *
     *   @param botBounds Vector of bottom bounds
     *   @param topBounds vector of top bounds
     */
    virtual void setBounds(const doublereal botBounds[], const doublereal topBounds[]);

    void setAtol(const doublereal atol[]);
    void setAtolConst(const doublereal atolconst);

private:
    //! Printing routine that gets called at the start of every invocation
    virtual void print_header(int ioflag, int ifunc, doublereal time_scale,
                              doublereal reltol,
                              doublereal netProdRate[]);

#ifdef DEBUG_SOLVEPROB
    //! Prints out the residual and Jacobian
    virtual void printResJac(int ioflag, int neq, const Array2D& Jac,
                             doublereal resid[], doublereal wtResid[], doublereal norm);
#endif

    //! Printing routine that gets called after every iteration
    virtual void printIteration(int ioflag, doublereal damp, size_t label_d, size_t label_t,
                                doublereal inv_t, doublereal t_real, int iter,
                                doublereal update_norm, doublereal resid_norm,
                                doublereal netProdRate[], doublereal CSolnSP[],
                                doublereal resid[],
                                doublereal wtSpecies[], size_t dim, bool do_time);

    //! Print a summary of the solution
    virtual void printFinal(int ioflag, doublereal damp, size_t label_d, size_t label_t,
                            doublereal inv_t, doublereal t_real, int iter,
                            doublereal update_norm, doublereal resid_norm,
                            doublereal netProdRateKinSpecies[], const doublereal CSolnSP[],
                            const doublereal resid[],
                            const doublereal wtSpecies[], const doublereal wtRes[],
                            size_t dim, bool do_time);

    //! Calculate a conservative delta T to use in a pseudo-steady state
    //! algorithm
    /*!
     *    This routine calculates a pretty conservative 1/del_t based
     *    on  MAX_i(sdot_i/(X_i*SDen0)).  This probably guarantees
     *    diagonal dominance.
     *
     *     Small surface fractions are allowed to intervene in the del_t
     *     determination, no matter how small.  This may be changed.
     *     Now minimum changed to 1.0e-12,
     *
     *     Maximum time step set to time_scale.
     *
     *    @param netProdRateSolnSP  Output variable. Net production rate
     *             of all of the species in the solution vector.
     *    @param Csoln output variable.
     *            Mole fraction of all of the species in the  solution vector
     *    @param label Output variable. Pointer to the value of the
     *                 species index (kindexSP) that is controlling
     *                 the time step
     *    @param label_old Output variable. Pointer to the value of the
     *                 species index (kindexSP) that controlled
     *                 the time step at the previous iteration
     *    @param label_factor Output variable. Pointer to the current
     *                 factor that is used to indicate the same species
     *                 is controlling the time step.
     *
     *    @param ioflag Level of the output requested.
     *
     *    @return  Returns the 1. /  delta T to be used on the next step
     */
    virtual doublereal calc_t(doublereal netProdRateSolnSP[], doublereal Csoln[],
                              size_t* label, size_t* label_old,
                              doublereal* label_factor, int ioflag);

    //! Calculate the solution and residual weights
    /*!
     * Calculate the weighting factors for norms wrt both the species
     * concentration unknowns and the residual unknowns.
     *  @param wtSpecies Weights to use for the soln unknowns. These
     *                    are in concentration units
     *  @param wtResid    Weights to sue for the residual unknowns.
     *
     *  @param CSolnSP    Solution vector for the surface problem
     */
    virtual void calcWeights(doublereal wtSpecies[], doublereal wtResid[],
                             const doublereal CSolnSP[]);

#ifdef DEBUG_SOLVEPROB
    //! Utility routine to print a header for high lvls of debugging
    /*!
     *  @param ioflag Lvl of debugging
     *  @param damp   lvl of damping
     *  @param inv_t  Inverse of the value of delta T
     *  @param t_real Value of the time
     *  @param iter   Iteration number
     *  @param do_time boolean indicating whether time stepping is taking
     *                 place
     */
    virtual void printIterationHeader(int ioflag, doublereal damp,
                                      doublereal inv_t, doublereal t_real, int iter,
                                      bool do_time);
#endif

    //! Main Function evaluation
    /*!
     * This calculates the net production rates of all species
     *
     *  @param resid output Vector of residuals, length = m_neq
     *  @param CSolnSP  Vector of species concentrations, unknowns in the
     *                  problem, length = m_neq
     *  @param CSolnSPOld Old Vector of species concentrations, unknowns in the
     *                  problem, length = m_neq
     *  @param do_time Calculate a time dependent residual
     *  @param deltaT  Delta time for time dependent problem.
     */
    virtual void  fun_eval(doublereal* const resid, const doublereal* const CSolnSP,
                           const doublereal* const CSolnSPOld, const bool do_time, const doublereal deltaT);

    //! Main routine that calculates the current residual and Jacobian
    /*!
     *  @param JacCol  Vector of pointers to the tops of columns of the
     *                 Jacobian to be evaluated.
     *  @param resid   output Vector of residuals, length = m_neq
     *  @param CSolnSP  Vector of species concentrations, unknowns in the
     *                  problem, length = m_neq. These are tweaked in order
     *                  to derive the columns of the jacobian.
     *  @param CSolnSPOld Old Vector of species concentrations, unknowns in the
     *                  problem, length = m_neq
     *  @param do_time Calculate a time dependent residual
     *  @param deltaT  Delta time for time dependent problem.
     */
    virtual void resjac_eval(std::vector<doublereal*>& JacCol, doublereal* resid,
                             doublereal* CSolnSP,
                             const doublereal* CSolnSPOld,  const bool do_time,
                             const doublereal deltaT);

    //!  This function calculates a damping factor for the Newton iteration update
    //!  vector, dxneg, to insure that all solution components stay within prescribed bounds
    /*!
     *  The default for this class is that all solution components are bounded between zero and one.
     *  this is because the original unknowns were mole fractions and surface site fractions.
     *
     *      dxneg[] = negative of the update vector.
     *
     * The constant "APPROACH" sets the fraction of the distance to the boundary
     * that the step can take.  If the full step would not force any fraction
     * outside of the bounds, then Newton's method is mostly allowed to operate normally.
     * There is also some solution damping employed.
     *
     *  @param x       Vector of the current solution components
     *  @param dxneg   Vector of the negative of the full solution update vector.
     *  @param dim     Size of the solution vector
     *  @param label   return int, stating which solution component caused the most damping.
     */
    virtual doublereal calc_damping(doublereal x[], doublereal dxneg[], size_t dim, size_t* label);

    //! residual function pointer to be solved.
    ResidEval* m_residFunc;

    //! Total number of equations to solve in the implicit problem.
    /*!
     * Note, this can be zero, and frequently is
     */
    size_t m_neq;

    //! m_atol is the absolute tolerance in real units.
    vector_fp m_atol;

    //! m_rtol is the relative error tolerance.
    doublereal m_rtol;

    //! maximum value of the time step
    /*!
     * units = seconds
     */
    doublereal m_maxstep;

    //! Temporary vector with length MAX(1, m_neq)
    vector_fp m_netProductionRatesSave;

    //! Temporary vector with length MAX(1, m_neq)
    vector_fp m_numEqn1;

    //! Temporary vector with  length MAX(1, m_neq)
    vector_fp m_numEqn2;

    //! Temporary vector with length MAX(1, m_neq)
    vector_fp m_CSolnSave;

    //! Solution vector
    /*!
     * length MAX(1, m_neq)
     */
    vector_fp m_CSolnSP;

    //! Saved initial solution vector
    /*!
     * length MAX(1, m_neq)
     */
    vector_fp m_CSolnSPInit;

    //! Saved  solution vector at the old time step
    /*!
     * length MAX(1, m_neq)
     */
    vector_fp m_CSolnSPOld;

    //!  Weights for the residual norm calculation
    /*!
     * length MAX(1, m_neq)
     */
    vector_fp m_wtResid;

    //!  Weights for the species concentrations norm calculation
    /*!
     * length MAX(1, m_neq)
     */
    vector_fp m_wtSpecies;

    //!  Residual for the surface problem
    /*!
     *  The residual vector of length "dim" that, that has the value
     *  of "sdot" for surface species.  The residuals for the bulk
     *  species are a function of the sdots for all species in the bulk
     *  phase. The last residual of each phase enforces {Sum(fractions)
     *  = 1}. After linear solve (dgetrf_ & dgetrs_), resid holds the
     *  update vector.
     *
     * length MAX(1, m_neq)
     */
    vector_fp m_resid;

    //!  pivots
    /*!
     * length MAX(1, m_neq)
     */
    vector_int m_ipiv;

    //! Vector of pointers to the top of the columns of the jacobians
    /*!
     *   The "dim" by "dim" computed Jacobian matrix for the
     *   local Newton's method.
     */
    std::vector<doublereal*> m_JacCol;

    //! Jacobian
    /*!
     *  m_neq by m_neq computed Jacobian matrix for the local Newton's method.
     */
    Array2D m_Jac;

    //! Top bounds for the solution vector
    /*!
     *  This defaults to 1.0
     */
    vector_fp m_topBounds;

    //! Bottom bounds for the solution vector
    /*!
     *  This defaults to 0.0
     */
    vector_fp m_botBounds;

public:
    int m_ioflag;
};
}
#endif
