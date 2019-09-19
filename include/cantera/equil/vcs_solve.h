/**
 * @file vcs_solve.h Header file for the internal object that holds the vcs
 *    equilibrium problem (see Class \link Cantera::VCS_SOLVE
 *    VCS_SOLVE\endlink and \ref equilfunctions ).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef _VCS_SOLVE_H
#define _VCS_SOLVE_H

/*
* Index of Symbols
* -------------------
*     irxn -> refers to the species or rxn between the species and
*             the components in the problem
*     k    -> refers to the species
*     j    -> refers to the element or component
*
*     ###  -> to be eliminated
*/

#include "cantera/base/ct_defs.h"
#include "cantera/equil/vcs_defs.h"
#include "cantera/equil/vcs_internal.h"
#include "cantera/base/Array.h"

namespace Cantera
{

class vcs_VolPhase;
class VCS_SPECIES_THERMO;
class VCS_COUNTERS;
class MultiPhase;

//! This is the main structure used to hold the internal data
//! used in vcs_solve_TP(), and to solve TP systems.
/*!
 * The indices of information in this structure may change when the species
 * basis changes or when phases pop in and out of existence. Both of these
 * operations change the species ordering.
 */
class VCS_SOLVE
{
public:
    //! Initialize the sizes within the VCS_SOLVE object
    /*!
     * This resizes all of the internal arrays within the object.
     */
    VCS_SOLVE(MultiPhase* mphase, int printLvl=0);

    ~VCS_SOLVE();

    //! Solve an equilibrium problem
    /*!
     * This is the main interface routine to the equilibrium solver
     *
     * @param ipr Printing of results
     *     ipr = 1 -> Print problem statement and final results to
     *                standard output
     *           0 -> don't report on anything
     * @param ip1 Printing of intermediate results
     *     IP1 = 1 -> Print intermediate results.
     * @param maxit  Maximum number of iterations for the algorithm
     * @return nonzero value: failure to solve the problem at hand. zero :
     *       success
     */
    int vcs(int ipr, int ip1, int maxit);

    //! Main routine that solves for equilibrium at constant T and P using a
    //! variant of the VCS method
    /*!
     * This is the main routine that solves for equilibrium at constant T and P
     * using a variant of the VCS method. Nonideal phases can be accommodated as
     * well.
     *
     * Any number of single-species phases and multi-species phases can be
     * handled by the present version.
     *
     * @param print_lvl     1 -> Print results to standard output;
     *                      0 -> don't report on anything
     * @param printDetails  1 -> Print intermediate results.
     * @param maxit         Maximum number of iterations for the algorithm
     * @return
     *     * 0 = Equilibrium Achieved
     *     * 1 = Range space error encountered. The element abundance criteria
     *       are only partially satisfied. Specifically, the first NC= (number
     *       of components) conditions are satisfied. However, the full NE
     *       (number of elements) conditions are not satisfied. The
     *       equilibrium condition is returned.
     *     * -1 = Maximum number of iterations is exceeded. Convergence was
     *       not found.
     */
    int vcs_solve_TP(int print_lvl, int printDetails, int maxit);

    /*!
     * We make decisions on the initial mole number, and major-minor status
     * here. We also fix up the total moles in a phase.
     *
     * irxn = id of the noncomponent species formation reaction for the
     *        species to be added in.
     *
     * The algorithm proceeds to implement these decisions in the previous
     * position of the species. Then, vcs_switch_pos is called to move the
     * species into the last active species slot, incrementing the number of
     * active species at the same time.
     *
     * This routine is responsible for the global data manipulation only.
     */
    void vcs_reinsert_deleted(size_t kspec);

    //! Choose the optimum species basis for the calculations
    /*!
     * This is done by choosing the species with the largest mole fraction not
     * currently a linear combination of the previous components. Then,
     * calculate the stoichiometric coefficient matrix for that basis.
     *
     * Rearranges the solution data to put the component data at the
     * front of the species list.
     *
     * Then, calculates `m_stoichCoeffRxnMatrix(jcomp,irxn)` the formation
     * reactions for all noncomponent species in the mechanism. Also
     * calculates DNG(I) and DNL(I), the net mole change for each formation
     * reaction. Also, initializes IR(I) to the default state.
     *
     * @param[in] doJustComponents  If true, the #m_stoichCoeffRxnMatrix and
     *                              #m_deltaMolNumPhase are not calculated.
     * @param[in] aw     Vector of mole fractions which will be used to
     *                   construct an optimal basis from.
     * @param[in] sa     Gram-Schmidt orthog work space (nc in length) sa[j]
     * @param[in] ss     Gram-Schmidt orthog work space (nc in length) ss[j]
     * @param[in] sm     QR matrix work space (nc*ne in length)         sm[i+j*ne]
     * @param[in] test   This is a small negative number dependent upon whether
     *                   an estimate is supplied or not.
     * @param[out] usedZeroedSpecies  If true, then a species with a zero
     *                                concentration was used as a component.
     *                                The problem may be converged. Or, the
     *                                problem may have a range space error and
     *                                may not have a proper solution.
     * @returns VCS_SUCCESS if everything went ok. Returns
     *     VCS_FAILED_CONVERGENCE if there is a problem.
     *
     * ### Internal Variables calculated by this routine:
     *
     * - #m_numComponents:  Number of component species. This routine
     *   calculates the #m_numComponents species. It switches their positions
     *   in the species vector so that they occupy the first #m_numComponents
     *   spots in the species vector.
     * - `m_stoichCoeffRxnMatrix(jcomp,irxn)` Stoichiometric coefficient
     *   matrix for the reaction mechanism expressed in Reduced Canonical
     *   Form. jcomp refers to the component number, and irxn refers to the
     *   irxn_th non-component species.
     * - `m_deltaMolNumPhase(iphase,irxn)`: Change in the number of moles in
     *   phase, iphase, due to the noncomponent formation reaction, irxn.
     * - `m_phaseParticipation(iphase,irxn)`: This is 1 if the phase, iphase,
     *   participates in the formation reaction, irxn, and zero otherwise.
     */
    int vcs_basopt(const bool doJustComponents, double aw[], double sa[], double sm[],
                   double ss[], double test, bool* const usedZeroedSpecies);

    //!  Choose a species to test for the next component
    /*!
     * We make the choice based on testing (molNum[i] * spSize[i]) for its
     * maximum value. Preference for single species phases is also made.
     *
     * @param molNum  Mole number vector
     * @param j       index into molNum[] that indicates where the search will
     *     start from Previous successful components are swapped into the front
     *     of molNum[].
     * @param n       Length of molNum[]
     */
    size_t vcs_basisOptMax(const double* const molNum, const size_t j, const size_t n);

    //! Evaluate the species category for the indicated species
    /*!
     * All evaluations are done using the "old" version of the solution.
     *
     * @param kspec   Species to be evaluated
     * @returns the calculated species type
     */
    int vcs_species_type(const size_t kspec) const;

    //!  This routine evaluates the species type for all species
    /*!
     *  - #VCS_SPECIES_MAJOR: Major species
     *  - #VCS_SPECIES_MINOR: Minor species
     *  - #VCS_SPECIES_SMALLMS: The species lies in a multicomponent phase
     *    that exists. Its concentration is currently very low, necessitating
     *    a different method of calculation.
     *  - #VCS_SPECIES_ZEROEDMS: The species lies in a multicomponent phase
     *    which currently doesn't exist. Its concentration is currently zero.
     *  - #VCS_SPECIES_ZEROEDSS: Species lies in a single-species phase which
     *       is currently zeroed out.
     *  - #VCS_SPECIES_DELETED: Species has such a small mole fraction it is
     *    deleted even though its phase may possibly exist. The species is
     *    believed to have such a small mole fraction that it best to throw
     *    the calculation of it out.  It will be added back in at the end of
     *    the calculation.
     *  - #VCS_SPECIES_INTERFACIALVOLTAGE: Species refers to an electron in
     *    the metal The unknown is equal to the interfacial voltage drop
     *    across the interface on the SHE (standard hydrogen electrode) scale
     *    (volts).
     *  - #VCS_SPECIES_ZEROEDPHASE: Species lies in a multicomponent phase
     *    that is zeroed atm and will stay deleted due to a choice from a
     *    higher level. These species will formally always have zero mole
     *    numbers in the solution vector.
     *  - #VCS_SPECIES_ACTIVEBUTZERO: The species lies in a multicomponent
     *    phase which currently does exist.  Its concentration is currently
     *    identically zero, though the phase exists. Note, this is a temporary
     *    condition that exists at the start of an equilibrium problem. The
     *    species is soon "birthed" or "deleted".
     *  - #VCS_SPECIES_STOICHZERO: The species lies in a multicomponent phase
     *    which currently does exist.  Its concentration is currently
     *    identically zero, though the phase exists. This is a permanent
     *    condition due to stoich constraints
     */
    bool vcs_evaluate_speciesType();

    //! Calculate the dimensionless chemical potentials of all species or
    //! of certain groups of species, at a fixed temperature and pressure.
    /*!
     * We calculate the dimensionless chemical potentials of all species or
     * certain groups of species here, at a fixed temperature and pressure, for
     * the input mole vector z[] in the parameter list. Nondimensionalization is
     * achieved by division by RT.
     *
     * Note, for multispecies phases which are currently zeroed out, the
     * chemical potential is filled out with the standard chemical potential.
     *
     * For species in multispecies phases whose concentration is zero, we need
     * to set the mole fraction to a very low value. Its chemical potential is
     * then calculated using the VCS_DELETE_MINORSPECIES_CUTOFF concentration
     * to keep numbers positive.
     *
     * For formulas, see vcs_chemPotPhase().
     *
     *  Handling of Small Species:
     * ------------------------------
     * As per the discussion above, for small species where the mole fraction
     *
     *     z(i) < VCS_DELETE_MINORSPECIES_CUTOFF
     *
     * The chemical potential is calculated as:
     *
     *      m_feSpecies(I)(I) = m_SSfeSpecies(I) + ln(ActCoeff[i](VCS_DELETE_MINORSPECIES_CUTOFF))
     *
     * Species in the following categories are treated as "small species"
     *   - #VCS_SPECIES_DELETED
     *   - #VCS_SPECIES_ACTIVEBUTZERO
     *
     * For species in multispecies phases which are currently not active, the
     * treatment is different. These species are in the following species
     * categories:
     *   - #VCS_SPECIES_ZEROEDMS
     *   - #VCS_SPECIES_ZEROEDPHASE
     *
     * For these species, the `ln( ActCoeff[I] X[I])` term is dropped
     * altogether. The following equation is used:
     *
     *     m_feSpecies(I) = m_SSfeSpecies(I)
     *                    + Charge[I] * Faraday_dim * phasePhi[iphase];
     *
     * Handling of "Species" Representing Interfacial Voltages
     * ---------------------------------------------------------
     *
     * These species have species types of
     * #VCS_SPECIES_TYPE_INTERFACIALVOLTAGE The chemical potentials for these
     * "species" refer to electrons in metal electrodes. They have the
     * following formula
     *
     *     m_feSpecies(I) = m_SSfeSpecies(I) - F z[I] / RT
     *
     * - `F` is Faraday's constant.
     * - `R` = gas constant
     * - `T` = temperature
     * - `V` = potential of the interface = phi_electrode - phi_solution
     *
     * For these species, the solution vector unknown, z[I], is V, the phase
     * voltage, in volts.
     *
     * @param ll     Determine which group of species gets updated
     *     - `ll = 0`: Calculate for all species
     *     - `ll < 0`: calculate for components and for major non-components
     *     - `ll = 1`: calculate for components and for minor non-components
     *
     * @param lbot    Restricts the calculation of the chemical potential
     *                to the species between LBOT <= i < LTOP. Usually
     *                 LBOT and LTOP will be equal to 0 and MR, respectively.
     * @param ltop    Top value of the loops
     *
     * @param stateCalc   Determines whether z is old or new or tentative:
     *            - 1: Use the tentative values for the total number of
     *                 moles in the phases, i.e., use TG1 instead of TG etc.
     *            - 0: Use the base values of the total number of
     *                 moles in each system.
     *
     *  Also needed:
     *     ff     : standard state chemical potentials. These are the
     *              chemical potentials of the standard states at
     *              the same T and P as the solution.
     *     tg     : Total Number of moles in the phase.
     */
    void vcs_dfe(const int stateCalc, const int ll, const size_t lbot, const size_t ltop);

    //! This routine uploads the state of the system into all of the
    //! vcs_VolumePhase objects in the current problem.
    /*!
     * @param stateCalc Determines where to get the mole numbers from.
     *               -  VCS_STATECALC_OLD -> from m_molNumSpecies_old
     *               -  VCS_STATECALC_NEW -> from m_molNumSpecies_new
     */
    void vcs_updateVP(const int stateCalc);

    //! Utility function that evaluates whether a phase can be popped
    //! into existence
    /*!
     * A phase can be popped iff the stoichiometric coefficients for the
     * component species, whose concentrations will be lowered during the
     * process, are positive by at least a small degree.
     *
     * If one of the phase species is a zeroed component, then the phase can
     * be popped if the component increases in mole number as the phase moles
     * are increased.
     *
     * @param iphasePop  id of the phase, which is currently zeroed,
     * @returns true if the phase can come into existence and false otherwise.
     */
    bool vcs_popPhasePossible(const size_t iphasePop) const;

    //! Decision as to whether a phase pops back into existence
    /*!
     * @param  phasePopPhaseIDs Vector containing the phase ids of the phases
     *         that will be popped this step.
     * @returns the phase id of the phase that pops back into existence. Returns
     *         -1 if there are no phases
     */
    size_t vcs_popPhaseID(std::vector<size_t> &phasePopPhaseIDs);

    //! Calculates the deltas of the reactions due to phases popping
    //! into existence
    /*!
     * Updates #m_deltaMolNumSpecies[irxn] : reaction adjustments, where irxn
     * refers to the irxn'th species formation reaction. This adjustment is
     * for species irxn + M, where M is the number of components.
     *
     * @param iphasePop  Phase id of the phase that will come into existence
     * @returns an int representing the status of the step
     *            -  0 : normal return
     *            -  1 : A single species phase species has been zeroed out
     *                   in this routine. The species is a noncomponent
     *            -  2 : Same as one but, the zeroed species is a component.
     *            -  3 : Nothing was done because the phase couldn't be birthed
     *                   because a needed component is zero.
     */
    int vcs_popPhaseRxnStepSizes(const size_t iphasePop);

    //! Calculates formation reaction step sizes.
    /*!
     * This is equation 6.4-16, p. 143 in Smith and Missen.
     *
     * Output
     * -------
     * m_deltaMolNumSpecies(irxn) : reaction adjustments, where irxn refers to
     *                              the irxn'th species formation reaction.
     *                              This adjustment is for species irxn + M,
     *                              where M is the number of components.
     *
     * Special branching occurs sometimes. This causes the component basis
     * to be reevaluated
     *
     * @param forceComponentCalc  integer flagging whether a component
     *                            recalculation needs to be carried out.
     * @param kSpecial            species number of phase being zeroed.
     * @returns an int representing which phase may need to be zeroed
     */
    size_t vcs_RxnStepSizes(int& forceComponentCalc, size_t& kSpecial);

    //! Calculates the total number of moles of species in all phases.
    /*!
     * Also updates the variable m_totalMolNum and Reconciles Phase existence
     * flags with total moles in each phase.
     */
    double vcs_tmoles();

    void check_tmoles() const;

    //! This subroutine calculates reaction free energy changes for
    //! all noncomponent formation reactions.
    /*!
     * Formation reactions are
     * reactions which create each noncomponent species from the component
     * species. `m_stoichCoeffRxnMatrix(jcomp,irxn)` are the stoichiometric
     * coefficients for these reactions. A stoichiometric coefficient of
     * one is assumed for species irxn in this reaction.
     *
     * @param L
     *    - `L < 0`: Calculate reactions corresponding to major noncomponent
     *      and zeroed species only
     *    - `L = 0`: Do all noncomponent reactions, i, between
     *      0 <= i < irxnl
     *    - `L > 0`: Calculate reactions corresponding to minor noncomponent
     *      and zeroed species only
     *
     * @param doDeleted   Do deleted species
     * @param vcsState    Calculate deltaG corresponding to either old or new
     *                     free energies
     * @param alterZeroedPhases boolean indicating whether we should
     *                          add in a special section for zeroed phases.
     *
     * Note we special case one important issue. If the component has zero
     * moles, then we do not allow deltaG < 0.0 for formation reactions which
     * would lead to the loss of more of that same component. This dG < 0.0
     * condition feeds back into the algorithm in several places, and leads
     * to a infinite loop in at least one case.
     */
    void vcs_deltag(const int L, const bool doDeleted, const int vcsState,
                    const bool alterZeroedPhases = true);

    void vcs_printDeltaG(const int stateCalc);

    //!  Swaps the indices for all of the global data for two species, k1
    //!  and k2.
    /*!
     * @param  ifunc  If true, switch the species data and the noncomponent
     *     reaction data. This must be called for a non-component species only.
     *     If false, switch the species data only. Typically, we use this option
     *     when determining the component species and at the end of the
     *     calculation, when we want to return unscrambled results. All rxn data
     *     will be out-of-date.
     *  @param k1     First species index
     *  @param k2     Second species index
     */
    void vcs_switch_pos(const bool ifunc, const size_t k1, const size_t k2);

    //! Main program to test whether a deleted phase should be brought
    //! back into existence
    /*!
     * @param iph Phase id of the deleted phase
     */
    double vcs_phaseStabilityTest(const size_t iph);

    //! Solve an equilibrium problem at a particular fixed temperature
    //! and pressure
    /*!
     * The actual problem statement is assumed to be in the structure already.
     * This is a wrapper around the solve_TP() function. In this wrapper, we
     * calculate the standard state Gibbs free energies of the species
     * and we decide whether to we need to use the initial guess algorithm.
     *
     * @param ipr = 1 -> Print results to standard output;
     *              0 -> don't report on anything
     * @param ip1 = 1 -> Print intermediate results;
     *              0 -> Dont print any intermediate results
     * @param maxit  Maximum number of iterations for the algorithm
     * @param T    Value of the Temperature (Kelvin)
     * @param pres Value of the Pressure
     * @returns an integer representing the success of the algorithm
     * * 0 = Equilibrium Achieved
     * * 1 = Range space error encountered. The element abundance criteria are
     *   only partially satisfied. Specifically, the first NC= (number of
     *   components) conditions are satisfied. However, the full NE (number of
     *   elements) conditions are not satisfied. The equilibrium condition is
     *   returned.
     * * -1 = Maximum number of iterations is exceeded. Convergence was not
     *   found.
     */
    int vcs_TP(int ipr, int ip1, int maxit, double T, double pres);

    /*!
     * Evaluate the standard state free energies at the current temperature
     * and pressure. Ideal gas pressure contribution is added in here.
     *
     * @param  ipr   1 -> Print results to standard output;  0 -> don't report
     *     on anything
     * @param  ip1   1 -> Print intermediate results; 0 -> don't.
     * @param  Temp  Temperature (Kelvin)
     * @param  pres  Pressure (Pascal)
     */
    int vcs_evalSS_TP(int ipr, int ip1, double Temp, double pres);

    //! Initialize the chemical potential of single species phases
    /*!
     * For single species phases, initialize the chemical potential with the
     * value of the standard state chemical potential. This value doesn't
     * change during the calculation
     */
    void vcs_fePrep_TP();

    //! Calculation of the total volume and the partial molar volumes
    /*!
     * This function calculates the partial molar volume for all species, kspec,
     * in the thermo problem at the temperature TKelvin and pressure, Pres, pres
     * is in atm. And, it calculates the total volume of the combined system.
     *
     * @param[in] tkelvin   Temperature in kelvin()
     * @param[in] pres      Pressure in Pascal
     * @param[in] w         w[] is the vector containing the current mole
     *                      numbers in units of kmol.
     * @param[out] volPM[]  For species in all phase, the entries are the
     *                      partial molar volumes units of M**3 / kmol.
     * @returns the total volume of the entire system in units of m**3.
     */
    double vcs_VolTotal(const double tkelvin, const double pres,
                        const double w[], double volPM[]);

    //! This routine is mostly concerned with changing the private data to be
    //! consistent with what's needed for solution. It is called one time for
    //! each new problem structure definition.
    /*!
     * The problem structure refers to:
     *
     *  - the number and identity of the species.
     *  - the formula matrix and thus the number of components.
     *  - the number and identity of the phases.
     *  - the equation of state
     *  - the method and parameters for determining the standard state
     *  - The method and parameters for determining the activity coefficients.
     *
     * Tasks:
     *    1. Fill in the SSPhase[] array.
     *    2. Check to see if any multispecies phases actually have only one
     *       species in that phase. If true, reassign that phase and species
     *       to be a single-species phase.
     *    3. Determine the number of components in the problem if not already
     *       done so. During this process the order of the species is changed
     *       in the private data structure. All references to the species
     *       properties must employ the ind[] index vector.
     *    4. Initialization of arrays to zero.
     *    5. Check to see if the problem is well posed (If all the element 
     *       abundances are zero, the algorithm will fail)
     *
     * @param printLvl Print level of the routine
     * @return
     *     VCS_SUCCESS = everything went OK;
     *     VCS_PUB_BAD = There is an irreconcilable difference in the
     *                   public data structure from when the problem was
     *                   initially set up.
     */
    int vcs_prep(int printLvl);

    //! Rearrange the constraint equations represented by the Formula
    //! Matrix so that the operational ones are in the front
    /*!
     * This subroutine handles the rearrangement of the constraint equations
     * represented by the Formula Matrix. Rearrangement is only necessary when
     * the number of components is less than the number of elements. For this
     * case, some constraints can never be satisfied exactly, because the
     * range space represented by the Formula Matrix of the components can't
     * span the extra space. These constraints, which are out of the range
     * space of the component Formula matrix entries, are migrated to the back
     * of the Formula matrix.
     *
     * A prototypical example is an extra element column in FormulaMatrix[],
     * which is identically zero. For example, let's say that argon is has an
     * element column in FormulaMatrix[], but no species in the mechanism
     * actually contains argon. Then, nc < ne. Also, without perturbation of
     * FormulaMatrix[] vcs_basopt[] would produce a zero pivot because the
     * matrix would be singular (unless the argon element column was already
     * the last column of FormulaMatrix[].
     *
     * This routine borrows heavily from vcs_basopt's algorithm. It finds nc
     * constraints which span the range space of the Component Formula matrix,
     * and assigns them as the first nc components in the formula matrix. This
     * guarantees that vcs_basopt[] has a nonsingular matrix to invert.
     *
     * Other Variables
     *  @param aw   Mole fraction work space        (ne in length)
     *  @param sa   Gram-Schmidt orthog work space (ne in length)
     *  @param sm   QR matrix work space (ne*ne in length)
     *  @param ss   Gram-Schmidt orthog work space (ne in length)
     */
    int vcs_elem_rearrange(double* const aw, double* const sa,
                           double* const sm, double* const ss);

    //! Swaps the indices for all of the global data for two elements, ipos
    //! and jpos.
    /*!
     * This function knows all of the element information with VCS_SOLVE, and
     * can therefore switch element positions
     *
     * @param ipos  first global element index
     * @param jpos  second global element index
     */
    void vcs_switch_elem_pos(size_t ipos, size_t jpos);

    //! Calculates the diagonal contribution to the Hessian due to
    //! the dependence of the activity coefficients on the mole numbers.
    /*!
     * (See framemaker notes, Eqn. 20 - VCS Equations document)
     *
     * We allow the diagonal to be increased positively to any degree.
     * We allow the diagonal to be decreased to 1/3 of the ideal solution
     * value, but no more -> it must remain positive.
     */
    double vcs_Hessian_diag_adj(size_t irxn, double hessianDiag_Ideal);

    //! Calculates the diagonal contribution to the Hessian due to
    //! the dependence of the activity coefficients on the mole numbers.
    /*!
     * (See framemaker notes, Eqn. 20 - VCS Equations document)
     */
    double vcs_Hessian_actCoeff_diag(size_t irxn);

    //! Recalculate all of the activity coefficients in all of the phases
    //! based on input mole numbers
    /*
     * @param moleSpeciesVCS kmol of species to be used in the update.
     *
     * NOTE: This routine needs to be regulated.
     */
    void vcs_CalcLnActCoeffJac(const double* const moleSpeciesVCS);

    //! Print out a report on the state of the equilibrium problem to
    //! standard output.
    /*!
     * @param iconv Indicator of convergence, to be printed out in the report:
     *    -   0 converged
     *    -   1 range space error
     *    -  -1 not converged
     */
    int vcs_report(int iconv);

    //! Computes the current elemental abundances vector
    /*!
     * Computes the elemental abundances vector, m_elemAbundances[], and stores
     * it back into the global structure
     */
    void vcs_elab();

    /*!
     * Checks to see if the element abundances are in compliance. If they are,
     * then TRUE is returned. If not, FALSE is returned. Note the number of
     * constraints checked is usually equal to the number of components in the
     * problem. This routine can check satisfaction of all of the constraints
     * in the problem, which is equal to ne. However, the solver can't fix
     * breakage of constraints above nc, because that nc is the range space by
     * definition. Satisfaction of extra constraints would have had to occur
     * in the problem specification.
     *
     * The constraints should be broken up into 2 sections. If a constraint
     * involves a formula matrix with positive and negative signs, and eaSet =
     * 0.0, then you can't expect that the sum will be zero. There may be
     * roundoff that inhibits this. However, if the formula matrix is all of
     * one sign, then this requires that all species with nonzero entries in
     * the formula matrix be identically zero. We put this into the logic
     * below.
     *
     * @param ibound  1: Checks constraints up to the number of elements;
     *                0: Checks constraints up to the number of components.
     */
    bool vcs_elabcheck(int ibound);

    /*!
     * Computes the elemental abundances vector for a single phase,
     * elemAbundPhase[], and returns it through the argument list. The mole
     * numbers of species are taken from the current value in
     * m_molNumSpecies_old[].
     */
    void vcs_elabPhase(size_t iphase, double* const elemAbundPhase);

    /*!
     * This subroutine corrects for element abundances. At the end of the
     * subroutine, the total moles in all phases are recalculated again,
     * because we have changed the number of moles in this routine.
     *
     * @param aa  temporary work vector, length ne*ne
     * @param x   temporary work vector, length ne
     *
     * @returns
     * - 0 = Nothing of significance happened, Element abundances were and
     *   still are good.
     * - 1 = The solution changed significantly; The element abundances are
     *   now good.
     * - 2 = The solution changed significantly, The element abundances are
     *   still bad.
     * - 3 = The solution changed significantly, The element abundances are
     *   still bad and a component species got zeroed out.
     *
     *  Internal data to be worked on::
     *  - ga    Current element abundances
     *  - m_elemAbundancesGoal   Required elemental abundances
     *  - m_molNumSpecies_old     Current mole number of species.
     *  - m_formulaMatrix[][]  Formula matrix of the species
     *  - ne    Number of elements
     *  - nc    Number of components.
     *
     * NOTES:
     *  This routine is turning out to be very problematic. There are
     *  lots of special cases and problems with zeroing out species.
     *
     *  Still need to check out when we do loops over nc vs. ne.
     */
    int vcs_elcorr(double aa[], double x[]);

    //! Create an initial estimate of the solution to the thermodynamic
    //! equilibrium problem.
    /*!
     * @return  value indicates success:
     *   -    0: successful initial guess
     *   -   -1: Unsuccessful initial guess; the elemental abundances aren't
     *           satisfied.
     */
    int vcs_inest_TP();

    //! Estimate the initial mole numbers by constrained linear programming
    /*!
     * This is done by running each reaction as far forward or backward as
     * possible, subject to the constraint that all mole numbers remain non-
     * negative. Reactions for which \f$ \Delta \mu^0 \f$ are positive are run
     * in reverse, and ones for which it is negative are run in the forward
     * direction. The end result is equivalent to solving the linear
     * programming problem of minimizing the linear Gibbs function subject to
     * the element and non-negativity constraints.
     */
    int vcs_setMolesLinProg();

    //! Calculate the total dimensionless Gibbs free energy
    /*!
     * Inert species are handled as if they had a standard free energy of
     * zero. Note, for this algorithm this function should be monotonically
     * decreasing.
     */
    double vcs_Total_Gibbs(double* w, double* fe, double* tPhMoles);

    //! Calculate the total dimensionless Gibbs free energy of a single phase
    /*!
     * Inert species are handled as if they had a standard free energy of
     * zero and if they obeyed ideal solution/gas theory.
     *
     * @param iphase   ID of the phase
     * @param w        Species mole number vector for all species
     * @param fe       vector of partial molar free energies of all of the
     *                 species
     */
    double vcs_GibbsPhase(size_t iphase, const double* const w,
                          const double* const fe);

    //! Transfer the results of the equilibrium calculation back from VCS_SOLVE
    void vcs_prob_update();

    //! Fully specify the problem to be solved
    void vcs_prob_specifyFully();

private:
    //! Zero out the concentration of a species.
    /*!
     * Make sure to conserve elements and keep track of the total moles in all
     * phases.
     * - w[]
     * - m_tPhaseMoles_old[]
     *
     * @param kspec  Species index
     *  @return:
     *      1: succeeded
     *      0: failed.
     */
    int vcs_zero_species(const size_t kspec);

    //! Change a single species from active to inactive status
    /*!
     * Rearrange data when species is added or removed. The Lth species is
     * moved to the back of the species vector. The back of the species
     * vector is indicated by the value of MR, the current number of
     * active species in the mechanism.
     *
     * @param kspec   Species Index
     * @return
     *     Returns 0 unless.
     *     The return is 1 when the current number of
     *     noncomponent species is equal to zero. A recheck of deleted species
     *     is carried out in the main code.
     */
    int vcs_delete_species(const size_t kspec);

    //! This routine handles the bookkeeping involved with the deletion of
    //! multiphase phases from the problem.
    /*!
     * When they are deleted, all of their species become active species, even
     * though their mole numbers are set to zero. The routine does not make the
     * decision to eliminate multiphases.
     *
     * Note, species in phases with zero mole numbers are still considered
     * active. Whether the phase pops back into existence or not is checked as
     * part of the main iteration loop.
     *
     * @param iph Phase to be deleted
     * @returns whether the operation was successful or not
     */
    bool vcs_delete_multiphase(const size_t iph);

    //! Change the concentration of a species by delta moles.
    /*!
     * Make sure to conserve elements and keep track of the total kmoles in
     * all phases.
     *
     * @param kspec The species index
     * @param delta_ptr   pointer to the delta for the species. This may
     *                     change during the calculation
     * @return
     *     1: succeeded without change of dx
     *     0: Had to adjust dx, perhaps to zero, in order to do the delta.
     */
    int delta_species(const size_t kspec, double* const delta_ptr);

    //! Provide an estimate for the deleted species in phases that are not
    //! zeroed out
    /*!
     * Try to add back in all deleted species. An estimate of the kmol numbers
     * are obtained and the species is added back into the equation system, into
     * the old state vector.
     *
     * This routine is called at the end of the calculation, just before
     * returning to the user.
     */
    size_t vcs_add_all_deleted();

    //! Recheck deleted species in multispecies phases.
    /*!
     * We are checking the equation:
     *
     *         sum_u = sum_j_comp [ sigma_i_j * u_j ]
     *               = u_i_O + log((AC_i * W_i)/m_tPhaseMoles_old)
     *
     * by first evaluating:
     *
     *          DG_i_O = u_i_O - sum_u.
     *
     * Then, if TL is zero, the phase pops into existence if DG_i_O < 0. Also,
     * if the phase exists, then we check to see if the species can have a mole
     * number larger than VCS_DELETE_SPECIES_CUTOFF (default value = 1.0E-32).
     */
    int vcs_recheck_deleted();

    //! Minor species alternative calculation
    /*!
     * This is based upon the following approximation: The mole fraction
     * changes due to these reactions don't affect the mole numbers of the
     * component species. Therefore the following approximation is valid for
     * a small component of an ideal phase:
     *
     *       0 = m_deltaGRxn_old(I) + log(molNum_new(I)/molNum_old(I))
     *
     * `m_deltaGRxn_old` contains the contribution from
     *
     *       m_feSpecies_old(I) =
     *              m_SSfeSpecies(I) +
     *              log(ActCoeff[i] * molNum_old(I) / m_tPhaseMoles_old(iph))
     * Thus,
     *
     *       molNum_new(I)= molNum_old(I) * EXP(-m_deltaGRxn_old(I))
     *
     * Most of this section is mainly restricting the update to reasonable
     * values. We restrict the update a factor of 1.0E10 up and 1.0E-10 down
     * because we run into trouble with the addition operator due to roundoff if
     * we go larger than ~1.0E15. Roundoff will then sometimes produce zero mole
     * fractions.
     *
     * Note: This routine was generalized to incorporate nonideal phases and
     * phases on the molality basis
     *
     * @param[in] kspec   The current species and corresponding formation
     *                    reaction number.
     * @param[in] irxn    The current species and corresponding formation
     *                    reaction number.
     * @param[out] do_delete:  BOOLEAN which if true on return, then we
     *                         branch to the section that deletes a species
     *                         from the current set of active species.
     */
    double vcs_minor_alt_calc(size_t kspec, size_t irxn, bool* do_delete,
                              char* ANOTE=0) const;

    //! This routine optimizes the minimization of the total Gibbs free energy
    //! by making sure the slope of the following functional stays negative:
    /*!
     * The slope of the following functional is equivalent to the slope of the
     * total Gibbs free energy of the system:
     *
     *      d_Gibbs/ds = sum_k( m_deltaGRxn * m_deltaMolNumSpecies[k] )
     *
     * along the current direction m_deltaMolNumSpecies[], by choosing a value, al: (0<al<1)
     * such that the a parabola approximation to Gibbs(al) fit to the
     * end points al = 0 and al = 1 is minimized.
     *  - s1 = slope of Gibbs function at al = 0, which is the previous
     *    solution = d(Gibbs)/d(al).
     *  - s2 = slope of Gibbs function at al = 1, which is the current
     *    solution = d(Gibbs)/d(al).
     *
     * Only if there has been an inflection point (i.e., s1 < 0 and s2 > 0),
     * does this code section kick in. It finds the point on the parabola
     * where the slope is equal to zero.
     */
    bool vcs_globStepDamp();

    //! Calculate the norm of a deltaGibbs free energy vector
    /*!
     * Positive DG for species which don't exist are ignored.
     *
     * @param dg Vector of local delta G's.
     */
    double l2normdg(double dg[]) const;

    void checkDelta1(double* const ds, double* const delTPhMoles, size_t kspec);

    //! Estimate equilibrium compositions
    /*!
     * Algorithm covered in a section of Smith and Missen's Book.
     *
     * Linear programming module is based on using dbolm.
     *
     * @param aw   aw[i[  Mole fraction work space        (ne in length)
     * @param sa   sa[j] = Gram-Schmidt orthog work space (ne in length)
     * @param sm   sm[i+j*ne] = QR matrix work space (ne*ne in length)
     * @param ss   ss[j] = Gram-Schmidt orthog work space (ne in length)
     * @param test This is a small negative number.
     */
    void vcs_inest(double* const aw, double* const sa, double* const sm,
                   double* const ss, double test);

    //! Calculate the status of single species phases.
    void vcs_SSPhase();

    //! Delete memory that isn't just resizable STL containers
    /*!
     * This gets called by the destructor or by InitSizes().
     */
    void vcs_delete_memory();

    //! Initialize the internal counters
    /*!
     * Initialize the internal counters containing the subroutine call
     * values and times spent in the subroutines.
     *
     *  ifunc = 0     Initialize only those counters appropriate for the top of
     *                vcs_solve_TP().
     *        = 1     Initialize all counters.
     */
    void vcs_counters_init(int ifunc);

    //! Create a report on the plog file containing timing and its information
    /*!
     * @param timing_print_lvl If 0, just report the iteration count. If larger
     *     than zero, report the timing information
     */
    void vcs_TCounters_report(int timing_print_lvl = 1);

    void vcs_setFlagsVolPhases(const bool upToDate, const int stateCalc);

    void vcs_setFlagsVolPhase(const size_t iph, const bool upToDate, const int stateCalc);

    //! Update all underlying vcs_VolPhase objects
    /*!
     * Update the mole numbers and the phase voltages of all phases in the vcs
     * problem
     *
     * @param stateCalc Location of the update (either VCS_STATECALC_NEW or
     *        VCS_STATECALC_OLD).
     */
    void vcs_updateMolNumVolPhases(const int stateCalc);

    // Helper functions used internally by vcs_solve_TP
    int solve_tp_component_calc(bool& allMinorZeroedSpecies);
    void solve_tp_inner(size_t& iti, size_t& it1, bool& uptodate_minors,
                        bool& allMinorZeroedSpecies, int& forceComponentCalc,
                        int& stage, bool printDetails, char* ANOTE);
    void solve_tp_equilib_check(bool& allMinorZeroedSpecies, bool& uptodate_minors,
                                bool& giveUpOnElemAbund, int& solveFail,
                                size_t& iti, size_t& it1, int maxit,
                                int& stage, bool& lec);
    void solve_tp_elem_abund_check(size_t& iti, int& stage, bool& lec,
                                   bool& giveUpOnElemAbund,
                                   int& finalElemAbundAttempts,
                                   int& rangeErrorFound);

    // data used by vcs_solve_TP and it's helper functions
    vector_fp m_sm;
    vector_fp m_ss;
    vector_fp m_sa;
    vector_fp m_aw;
    vector_fp m_wx;

public:
    //! Print level for print routines
    int m_printLvl;

    MultiPhase* m_mix;

    //! Print out the  problem specification in all generality as it currently
    //! exists in the VCS_SOLVE object
    /*!
     *  @param print_lvl Parameter lvl for printing
     *      * 0 - no printing
     *      * 1 - all printing
     */
    void prob_report(int print_lvl);

    //! Add elements to the local element list
    /*!
     * This routine sorts through the elements defined in the vcs_VolPhase
     * object. It then adds the new elements to the VCS_SOLVE object, and creates
     * a global map, which is stored in the vcs_VolPhase object. Id and matching
     * of elements is done strictly via the element name, with case not
     * mattering.
     *
     * The routine also fills in the position of the element in the vcs_VolPhase
     * object's ElGlobalIndex field.
     *
     * @param volPhase  Object containing the phase to be added. The elements in
     *     this phase are parsed for addition to the global element list
     */
    void addPhaseElements(vcs_VolPhase* volPhase);

    //! This routines adds entries for the formula matrix for one species
    /*!
     * This routines adds entries for the formula matrix for this object for one
     * species
     *
     * This object also fills in the index filed, IndSpecies, within the
     * volPhase object.
     *
     * @param volPhase object containing the species
     * @param k        Species number within the volPhase k
     * @param kT       global Species number within this object
     *
     */
    size_t addOnePhaseSpecies(vcs_VolPhase* volPhase, size_t k, size_t kT);

    //! This routine resizes the number of elements in the VCS_SOLVE object by
    //! adding a new element to the end of the element list
    /*!
     * The element name is added. Formula vector entries ang element abundances
     * for the new element are set to zero.
     *
     * @param elNameNew New name of the element
     * @param elType    Type of the element
     * @param elactive  boolean indicating whether the element is active
     * @returns the index number of the new element
     */
    size_t addElement(const char* elNameNew, int elType, int elactive);

    void reportCSV(const std::string& reportFile);

    //!  Total number of species in the problems
    size_t m_nsp;

    //! Number of element constraints in the problem
    size_t m_nelem;

    //! Number of components calculated for the problem
    size_t m_numComponents;

    //! Total number of non-component species in the problem
    size_t m_numRxnTot;

    //! Current number of species in the problems. Species can be deleted if
    //! they aren't stable under the current conditions
    size_t m_numSpeciesRdc;

    //! Current number of non-component species in the problem. Species can be
    //! deleted if they aren't stable under the current conditions
    size_t m_numRxnRdc;

    //! Number of active species which are currently either treated as
    //! minor species
    size_t m_numRxnMinorZeroed;

    //! Number of Phases in the problem
    size_t m_numPhases;

    //! Formula matrix for the problem
    /*!
     *  FormulaMatrix(kspec,j) =  Number of elements, j, in the kspec species
     *
     *  Both element and species indices are swapped.
     */
    Array2D m_formulaMatrix;

    //! Stoichiometric coefficient matrix for the reaction mechanism expressed
    //! in Reduced Canonical Form.
    /*!
     * This is the stoichiometric coefficient matrix for the reaction which
     * forms species kspec from the component species. A stoichiometric
     * coefficient of one is assumed for the species kspec in this mechanism.
     *
     * NOTE: kspec = irxn + m_numComponents
     *
     * `m_stoichCoeffRxnMatrix(j,irxn)` : j refers to the component number, and
     *     irxn refers to the irxn_th non-component species. The stoichiometric
     *     coefficients multiplied by the Formula coefficients of the component
     *     species add up to the negative value of the number of elements in the
     *     species kspec.
     *
     * size = nelements0 x nspecies0
     */
    Array2D m_stoichCoeffRxnMatrix;

    //! Absolute size of the stoichiometric coefficients
    /*!
     * scSize[irxn] = abs(Size) of the stoichiometric coefficients. These are
     *                used to determine whether a given species should be
     *                handled by the alt_min treatment or should be handled as a
     *                major species.
     */
    vector_fp m_scSize;

    //! total size of the species
    /*!
     *  This is used as a multiplier to the mole number in figuring out which
     *  species should be components.
     */
    vector_fp m_spSize;

    //! Standard state chemical potentials for species K at the current
    //! temperature and pressure.
    /*!
     *  The first NC entries are for components. The following NR entries are
     *  for the current non-component species in the mechanism.
     */
    vector_fp m_SSfeSpecies;

    //! Free energy vector from the start of the current iteration
    /*!
     *  The free energies are saved at the start of the current iteration.
     *  Length = number of species
     */
    vector_fp m_feSpecies_old;

    //! Dimensionless new free energy for all the species in the mechanism
    //! at the new tentative T, P, and mole numbers.
    /*!
     *   The first NC entries are for components. The following
     *   NR entries are for the current non-component species in the mechanism.
     *  Length = number of species
     */
    vector_fp m_feSpecies_new;

    //! Setting for whether to do an initial estimate
    /*!
     *  Initial estimate:
     *       0 Do not estimate the solution at all. Use the
     *         supplied mole numbers as is.
     *       1 Only do an estimate if the element abundances aren't satisfied.
     *      -1 Force an estimate of the soln. Throw out the input mole numbers.
     */
    int m_doEstimateEquil;

    //! Total moles of the species
    /*!
     *  Total number of moles of the kth species.
     *  Length = Total number of species = m
     */
    vector_fp m_molNumSpecies_old;

    //! Specifies the species unknown type
    /*!
     *  There are two types. One is the straightforward species, with the mole
     *  number w[k], as the unknown. The second is the an interfacial voltage
     *  where w[k] refers to the interfacial voltage in volts.
     *
     *  These species types correspond to metallic electrons corresponding to
     *  electrodes. The voltage and other interfacial conditions sets up an
     *  interfacial current, which is set to zero in this initial treatment.
     *  Later we may have non-zero interfacial currents.
     */
    vector_int m_speciesUnknownType;

    //! Change in the number of moles of phase, iphase, due to the
    //! noncomponent formation reaction, irxn, for species, k:
    /*!
     *       m_deltaMolNumPhase(iphase,irxn) = k = nc + irxn
     */
    Array2D m_deltaMolNumPhase;

    //!  This is 1 if the phase, iphase, participates in the formation reaction
    //!  irxn, and zero otherwise.  PhaseParticipation(iphase,irxn)
    Array2D m_phaseParticipation;

    //! electric potential of the iph phase
    vector_fp m_phasePhi;

    //! Tentative value of the mole number vector. It's also used to store the
    //! mole fraction vector.
    vector_fp m_molNumSpecies_new;

    //! Delta G(irxn) for the noncomponent species in the mechanism.
    /*!
     * Computed by the subroutine deltaG. m_deltaGRxn is the free energy change
     * for the reaction which forms species K from the component species. This
     * vector has length equal to the number of noncomponent species in the
     * mechanism. It starts with the first current noncomponent species in the
     * mechanism.
     */
    vector_fp m_deltaGRxn_new;

    //!  Last deltag[irxn] from the previous step
    vector_fp m_deltaGRxn_old;

    //! Last deltag[irxn] from the previous step with additions for
    //! possible births of zeroed phases.
    vector_fp m_deltaGRxn_Deficient;

    //! Temporary vector of Rxn DeltaG's
    /*!
     *  This is used from time to time, for printing purposes
     */
    vector_fp m_deltaGRxn_tmp;

    //! Reaction Adjustments for each species during the current step
    /*!
     *  delta Moles for each species during the current step.
     *  Length = number of species
     */
    vector_fp m_deltaMolNumSpecies;

    //! Element abundances vector
    /*!
     *  Vector of moles of each element actually in the solution vector. Except
     *  for certain parts of the algorithm, this is a constant. Note other
     *  constraint conditions are added to this vector. This is input from the
     *  input file and is considered a constant from thereon. units = kmoles
     */
    vector_fp m_elemAbundances;

    //! Element abundances vector Goals
    /*!
     * Vector of moles of each element that are the goals of the simulation.
     * This is a constant in the problem. Note other constraint conditions are
     * added to this vector. This is input from the input file and is considered
     * a constant from thereon. units = kmoles
     */
    vector_fp m_elemAbundancesGoal;

    //! Total number of kmoles in all phases
    /*!
     * This number includes the inerts.
     *            -> Don't use this except for scaling purposes
     */
    double m_totalMolNum;

    //! Total kmols of species in each phase
    /*!
     * This contains the total number of moles of species in each phase
     *
     * Length = number of phases
     */
    vector_fp m_tPhaseMoles_old;

    //! total kmols of species in each phase in the tentative soln vector
    /*!
     * This contains the total number of moles of species in each phase
     * in the tentative solution vector
     *
     * Length = number of phases
     */
    vector_fp m_tPhaseMoles_new;

    //! Temporary vector of length NPhase
    mutable vector_fp m_TmpPhase;

    //! Temporary vector of length NPhase
    mutable vector_fp m_TmpPhase2;

    //! Change in the total moles in each phase
    /*!
     *  Length number of phases.
     */
    vector_fp m_deltaPhaseMoles;

    //! Temperature (Kelvin)
    double m_temperature;

    //! Pressure
    double m_pressurePA;

    //! Total kmoles of inert to add to each phase
    /*!
     * TPhInertMoles[iph] = Total kmoles of inert to add to each phase
     * length = number of phases
     */
    vector_fp TPhInertMoles;

    //! Tolerance requirement for major species
    double m_tolmaj;

    //! Tolerance requirements for minor species
    double m_tolmin;

    //! Below this, major species aren't refined any more
    double m_tolmaj2;

    //! Below this, minor species aren't refined any more
    double m_tolmin2;

    //! Index vector that keeps track of the species vector rearrangement
    /*!
     * At the end of each run, the species vector and associated data gets put
     * back in the original order.
     *
     * Example
     *
     *     k = m_speciesMapIndex[kspec]
     *
     *     kspec = current order in the vcs_solve object
     *     k     = original order in the MultiPhase object
     */
    std::vector<size_t> m_speciesMapIndex;

    //! Index that keeps track of the index of the species within the local
    //! phase
    /*!
     * This returns the local index of the species within the phase. Its
     * argument is the global species index within the VCS problem.
     *
     * k = m_speciesLocalPhaseIndex[kspec]
     *
     * k varies between 0 and the nSpecies in the phase
     *
     * Length = number of species
     */
    std::vector<size_t> m_speciesLocalPhaseIndex;

    //! Index vector that keeps track of the rearrangement of the elements
    /*!
     * At the end of each run, the element vector and associated data gets
     * put back in the original order.
     *
     * Example
     *
     *     e    = m_elementMapIndex[eNum]
     *     eNum  = current order in the vcs_solve object
     *     e     = original order in the MultiPhase object
     */
    std::vector<size_t> m_elementMapIndex;

    //! Mapping between the species index for noncomponent species and the
    //! full species index.
    /*!
     * ir[irxn]   = Mapping between the reaction index for noncomponent
     *              formation reaction of a species and the full species
     *              index.
     *
     * Initially set to a value of K = NC + I This vector has length equal to
     * number of noncomponent species in the mechanism. It starts with the
     * first current noncomponent species in the mechanism. kspec = ir[irxn]
     */
    std::vector<size_t> m_indexRxnToSpecies;

    //! Major -Minor status vector for the species in the problem
    /*!
     * The index for this is species. The reaction that this is referring to
     * is `kspec = irxn + m_numComponents`. For possible values and their
     * meanings, see vcs_evaluate_speciesType().
     */
    vector_int m_speciesStatus;

    //! Mapping from the species number to the phase number
    std::vector<size_t> m_phaseID;

    //! Boolean indicating whether a species belongs to a single-species phase
    // vector<bool> can't be used here because it doesn't work with std::swap
    std::vector<char> m_SSPhase;

    //! Species string name for the kth species
    std::vector<std::string> m_speciesName;

    //! Vector of strings containing the element names
    /*!
     *   ElName[j]  = String containing element names
     */
    std::vector<std::string> m_elementName;

    //! Type of the element constraint
    /*!
     *  * 0 - #VCS_ELEM_TYPE_ABSPOS Normal element that is positive or zero in
     *    all species.
     *  * 1 - #VCS_ELEM_TYPE_ELECTRONCHARGE element dof that corresponds to
     *    the electronic charge DOF.
     *  * 2 - #VCS_ELEM_TYPE_CHARGENEUTRALITY element dof that corresponds to
     *    a required charge neutrality constraint on the phase. The element
     *    abundance is always exactly zero.
     *  * 3 - #VCS_ELEM_TYPE_OTHERCONSTRAINT Other constraint which may mean
     *    that a species has neg 0 or pos value of that constraint (other than
     *    charge)
     */
    vector_int m_elType;

    //! Specifies whether an element constraint is active
    /*!
     * The default is true. Length = nelements
     */
    vector_int m_elementActive;

    //! Array of Phase Structures. Length = number of phases.
    std::vector<std::unique_ptr<vcs_VolPhase>> m_VolPhaseList;

    //! specifies the activity convention of the phase containing the species
    /*!
     *  * 0 = molar based
     *  * 1 = molality based
     *
     * length = number of species
     */
    vector_int m_actConventionSpecies;

    //! specifies the activity convention of the phase.
    /*!
     *  * 0 = molar based
     *  * 1 = molality based
     *
     * length = number of phases
     */
    vector_int m_phaseActConvention;

    //! specifies the ln(Mnaught) used to calculate the chemical potentials
    /*!
     * For molar based activity conventions this will be equal to 0.0.
     * length = number of species.
     */
    vector_fp m_lnMnaughtSpecies;

    //! Molar-based Activity Coefficients for Species.
    //! Length = number of species
    vector_fp m_actCoeffSpecies_new;

    //! Molar-based Activity Coefficients for Species based on old mole numbers
    /*!
     * These activity coefficients are based on the m_molNumSpecies_old
     * values Molar based activity coefficients. Length = number of species
     */
    vector_fp m_actCoeffSpecies_old;

    //! Change in the log of the activity coefficient with respect to the mole
    //! number multiplied by the phase mole number
    /*!
     * size = nspecies x nspecies
     *
     * This is a temporary array that gets regenerated every time it's needed.
     * It is not swapped wrt species.
     */
    Array2D m_np_dLnActCoeffdMolNum;

    //! Molecular weight of each species
    /*!
     * units = kg/kmol. length = number of species.
     *
     * note: this is a candidate for removal. I don't think we use it.
     */
    vector_fp m_wtSpecies;

    //! Charge of each species. Length = number of species.
    vector_fp m_chargeSpecies;

    std::vector<std::vector<size_t> > phasePopProblemLists_;

    //! Vector of pointers to thermo structures which identify the model
    //! and parameters for evaluating the thermodynamic functions for that
    //! particular species.
    /*!
     * SpeciesThermo[k] pointer to the thermo information for the kth species
     */
    std::vector<std::unique_ptr<VCS_SPECIES_THERMO>> m_speciesThermoList;

    //! Choice of Hessians
    /*!
     * If this is true, then we will use a better approximation to the Hessian
     * based on Jacobian of the ln(ActCoeff) with respect to mole numbers
     */
    int m_useActCoeffJac;

    //! Total volume of all phases. Units are m^3
    double m_totalVol;

    //! Partial molar volumes of the species
    /*!
     *  units = mks (m^3/kmol)
     *  Length = number of species
     */
    vector_fp m_PMVolumeSpecies;

    //! dimensionless value of Faraday's constant, F / RT  (1/volt)
    double m_Faraday_dim;

    //! Timing and iteration counters for the vcs object
    VCS_COUNTERS* m_VCount;

    //! Debug printing lvl
    /*!
     *  Levels correspond to the following guidelines
     *     * 0  No printing at all
     *     * 1  Serious warnings or fatal errors get one line
     *     * 2  one line per each successful vcs package call
     *     * 3  one line per every successful solve_TP calculation
     *     * 4  one line for every successful operation -> solve_TP gets a summary report
     *     * 5  each iteration in solve_TP gets a report with one line per species
     *     * 6  Each decision in solve_TP gets a line per species in addition to 4
     *     * 10 Additionally Hessian matrix is printed out
     */
    int m_debug_print_lvl;

    //! printing level of timing information
    /*!
     *  * 1 allowing printing of timing
     *  * 0 do not allow printing of timing -> everything is printed as a NA.
     */
    int m_timing_print_lvl;

    //! Disable printing of timing information. Used to generate consistent
    //! output for tests.
    static void disableTiming();

    friend class vcs_phaseStabilitySolve;
};

}
#endif
