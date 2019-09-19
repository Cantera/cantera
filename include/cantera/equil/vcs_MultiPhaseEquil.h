/**
 *  @file  vcs_MultiPhaseEquil.h
 *  Interface class for the vcsnonlinear solver
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef VCS_MULTIPHASEEQUIL_H
#define VCS_MULTIPHASEEQUIL_H

#include "MultiPhase.h"
#include "vcs_solve.h"

namespace Cantera
{

//! %Cantera's Interface to the Multiphase chemical equilibrium solver.
/*!
 * Class vcs_MultiPhaseEquil is designed to be used to set a mixture containing
 * one or more phases to a state of chemical equilibrium.
 *
 * Note, as currently constructed, the underlying ThermoPhase objects are shared
 * between the MultiPhase object and this object. Therefore, mix is not a const
 * argument, and the return parameters are contained in underlying ThermoPhase
 * objects.
 *
 * @ingroup equilfunctions
 */
class vcs_MultiPhaseEquil
{
public:
    //! Constructor for the multiphase equilibrium solver
    /*!
     * This constructor will initialize the object with a MultiPhase object,
     * setting up the internal equilibration problem. Note, as currently
     * constructed, the underlying ThermoPhase objects are shared between the
     * MultiPhase object and this object. Therefore, mix is not a const
     * argument, and the return parameters are contained in underlying
     * ThermoPhase objects.
     *
     * @param mix Object containing the MultiPhase object
     * @param printLvl Determines the amount of printing to stdout
     *            that occurs for each call:
     *        - 0: No printing
     *        - 1: Only printing to the .csv file
     *        - 2: print the soln only
     *        - 3: Print the setup and then the soln only
     *        - 4: Print a table for each iteration
     *        - 5: Print more than a table for each iteration
     */
    vcs_MultiPhaseEquil(MultiPhase* mix, int printLvl);

    virtual ~vcs_MultiPhaseEquil() {}

    //! return the number of iterations
    int iterations() const {
        return m_iter;
    }

    //! Equilibrate the solution using the current element abundances
    //! stored in the MultiPhase object
    /*!
     * Use the vcs algorithm to equilibrate the current multiphase mixture.
     *
     * @param XY       Integer representing what two thermo quantities are
     *                 held constant during the equilibration
     * @param estimateEquil integer indicating whether the solver should
     *     estimate its own initial condition.
     *     - If 0, the initial mole fraction vector in the ThermoPhase object is
     *       used as the initial condition.
     *     - If 1, the initial mole fraction vector is used if the element
     *       abundances are satisfied.
     *     - if -1, the initial mole fraction vector is thrown out, and an
     *       estimate is formulated.
     * @param printLvl Determines the amount of printing that gets sent to
     *     stdout from the vcs package (Note, you may have to compile with debug
     *     flags to get some printing).
     * @param err      Internal error level
     * @param maxsteps max steps allowed.
     * @param loglevel Determines the amount of printing to the output file.
     */
    int equilibrate(int XY, int estimateEquil = 0,
                    int printLvl= 0, doublereal err = 1.0e-6,
                    int maxsteps = VCS_MAXSTEPS, int loglevel=-99);

    //! Equilibrate the solution using the current element abundances
    //! stored in the MultiPhase object using constant T and P
    /*!
     * Use the vcs algorithm to equilibrate the current multiphase mixture.
     *
     * @param estimateEquil integer indicating whether the solver should
     *     estimate its own initial condition.
     *     - If 0, the initial mole fraction vector in the ThermoPhase object is
     *       used as the initial condition.
     *     - If 1, the initial mole fraction vector is used if the element
     *       abundances are satisfied.
     *     - if -1, the initial mole fraction vector is thrown out, and an
     *       estimate is formulated.
     * @param printLvl  Determines the amount of printing that gets sent to
     *     stdout from the vcs package (Note, you may have to compile with debug
     *     flags to get some printing).
     * @param err     Internal error level
     * @param maxsteps max steps allowed.
     * @param loglevel Determines the amount of printing to the output file.
     */
    int equilibrate_TP(int estimateEquil = 0,
                       int printLvl= 0, doublereal err = 1.0e-6,
                       int maxsteps = VCS_MAXSTEPS, int loglevel=-99);

    //! Equilibrate the solution using the current element abundances
    //! stored in the MultiPhase object using either constant H and P
    //! or constant U and P.
    /*!
     * Use the vcs algorithm to equilibrate the current multiphase mixture. The
     * pressure of the calculation is taken from the current pressure stored
     * with the MultiPhase object.
     *
     * @param Htarget Value of the total mixture enthalpy or total internal
     *     energy that will be kept constant. Note, this is and must be an
     *     extensive quantity.  units = Joules
     * @param XY      Integer flag indicating what is held constant. Must be
     *     either HP or UP.
     * @param Tlow    Lower limit of the temperature. It's an error condition
     *     if the temperature falls below Tlow.
     * @param Thigh   Upper limit of the temperature. It's an error condition
     *     if the temperature goes higher than Thigh.
     * @param estimateEquil integer indicating whether the solver
     *     should estimate its own initial condition.
     *     - If 0, the initial mole fraction vector in the ThermoPhase object is
     *       used as the initial condition.
     *     - If 1, the initial mole fraction vector is used if the element
     *       abundances are satisfied.
     *     - if -1, the initial mole fraction vector is thrown out, and an
     *       estimate is formulated.
     * @param printLvl  Determines the amount of printing that gets sent to
     *                  stdout from the vcs package (Note, you may have to
     *                  compile with debug flags to get some printing). See main
     *                  constructor call for meaning of the levels.
     * @param err     Internal error level
     * @param maxsteps max steps allowed.
     * @param loglevel Determines the amount of printing to the output file.
     */
    int equilibrate_HP(doublereal Htarget, int XY, double Tlow, double Thigh,
                       int estimateEquil = 0,
                       int printLvl = 0, doublereal err = 1.0E-6,
                       int maxsteps = VCS_MAXSTEPS, int loglevel=-99);

    //! Equilibrate the solution using the current element abundances stored in
    //! the MultiPhase object using constant S and P.
    /*!
     * Use the vcs algorithm to equilibrate the current multiphase mixture. The
     * pressure of the calculation is taken from the current pressure stored
     * with the MultiPhase object.
     *
     * @param Starget Value of the total mixture entropy that will be kept
     *     constant. Note, this is and must be an extensive quantity.
     *     units = Joules/K
     * @param Tlow    Lower limit of the temperature. It's an error condition if
     *     the temperature falls below Tlow.
     * @param Thigh   Upper limit of the temperature. It's an error condition if
     *     the temperature goes higher than Thigh.
     * @param estimateEquil integer indicating whether the solver should
     *     estimate its own initial condition.
     *     - If 0, the initial mole fraction vector in the ThermoPhase object is
     *       used as the initial condition.
     *     - If 1, the initial mole fraction vector is used if the element
     *       abundances are satisfied.
     *     - If -1, the initial mole fraction vector is thrown out, and an
     *       estimate is formulated.
     * @param printLvl  Determines the amount of printing that gets sent to
     *                  stdout from the vcs package (Note, you may have to
     *                  compile with debug flags to get some printing). See main
     *                  constructor call for meaning of the levels.
     * @param err     Internal error level
     * @param maxsteps max steps allowed.
     * @param loglevel Determines the amount of printing to the output file.
     */
    int equilibrate_SP(doublereal Starget, double Tlow, double Thigh,
                       int estimateEquil = 0,
                       int printLvl = 0, doublereal err = 1.0E-6,
                       int maxsteps = VCS_MAXSTEPS, int loglevel=-99);

    //! Equilibrate the solution using the current element abundances stored
    //! in the MultiPhase object using constant V and constant T, H, U or S.
    /*!
     * Use the vcs algorithm to equilibrate the current multiphase mixture. The
     * pressure of the calculation is taken from the current pressure stored
     * with the MultiPhase object.
     *
     * @param XY      Integer flag indicating what is held constant.
     *     Must be either TV, HV, UV, or SV.
     * @param xtarget Value of the total thermodynamic parameter to be held
     *     constant in addition to V. Note, except for T, this must be an
     *     extensive quantity.  units = Joules/K or Joules
     * @param estimateEquil integer indicating whether the solver should
     *     estimate its own initial condition.
     *     - If 0, the initial mole fraction vector in the ThermoPhase object is
     *       used as the initial condition.
     *     - If 1, the initial mole fraction vector is used if the element
     *       abundances are satisfied.
     *     - if -1, the initial mole fraction vector is thrown out, and an
     *       estimate is formulated.
     * @param printLvl  Determines the amount of printing that gets sent to
     *     stdout from the vcs package (Note, you may have to compile with debug
     *     flags to get some printing). See main constructor call for meaning of
     *     the levels.
     * @param err      Internal error level
     * @param maxsteps max steps allowed.
     * @param logLevel Determines the amount of printing to the output file.
     */
    int equilibrate_TV(int XY, doublereal xtarget,
                       int estimateEquil = 0,
                       int printLvl = 0, doublereal err = 1.0E-6,
                       int maxsteps = VCS_MAXSTEPS, int logLevel = -99);

    //! Report the equilibrium answer in a comma separated table format
    /*!
     * This routine is used for in the test suite.
     *
     * @param reportFile Base name of the file to get the report. File name is
     *     incremented by 1 for each report.
     */
    void reportCSV(const std::string& reportFile);

protected:
    //! Vector that takes into account of the current sorting of the species
    /*!
     * The index of m_order is the original k value of the species in the
     * multiphase.  The value of m_order, k_sorted, is the current value of the
     * species index.
     *
     * `m_order[korig] = k_sorted`
     */
    vector_int m_order;

    //! Pointer to the MultiPhase mixture that will be equilibrated.
    /*!
     *  Equilibrium solutions will be returned via this variable.
     */
    MultiPhase* m_mix;

    //! Print level from the VCSnonlinear package
    /*!
     * (Note, you may have to compile with debug flags to get some printing).
     *
     * - 0: No IO from the routine whatsoever
     * - 1: file IO from reportCSV() carried out. One line print statements
     *   from equilibrate_XY() functions
     * - 2: Problem statement information from vcs_Cantera_update_vprob();
     *   Final state of the system from vcs_solve_TP()
     * - 3: Several more setup tables; Problem initialization routine
     * - 4: One table for each iteration within vcs_solve_Tp()
     * - 5: Multiple tables for each iteration within vcs_solve_TP()
     */
    int m_printLvl;

    //! Stoichiometric matrix
    DenseMatrix m_N;

    //! Iteration Count
    int m_iter;

    //! Vector of indices for species that are included in the calculation. This
    //! is used to exclude pure-phase species with invalid thermo data
    vector_int m_species;

    //! The object that contains the problem statement and does all of the equilibration work
    /*!
     * The problem statement may contain some subtleties. For example, the
     * element constraints may be different than just an element conservation
     * contraint equations. There may be kinetically frozen degrees of freedom.
     * There may be multiple electrolyte phases with zero charge constraints.
     * All of these make the problem statement different than the simple element
     * conservation statement.
     *
     * VCS_SOLVE will have different ordering for species and element constraints
     * than this object.
     */
    VCS_SOLVE m_vsolve;
};

}

#endif
