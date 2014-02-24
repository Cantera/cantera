/**
 *  @file  vcs_MultiPhaseEquil.h
 *  Interface class for the vcsnonlinear solver
 */
#ifndef VCS_MULTIPHASEEQUIL_H
#define VCS_MULTIPHASEEQUIL_H

#include "cantera/base/ct_defs.h"
#include "MultiPhase.h"
#include "vcs_defs.h"

namespace Cantera
{
//!  Set a single-phase chemical solution to chemical equilibrium.
/*!
 *  The function uses the element abundance vector that is currently
 *  consistent with the composition within the phase itself. Two other
 *  thermodynamic quantities, determined by the XY string,  are held constant
 *  during the equilibration. This is a convenience function that uses one or
 *  the other of the two chemical equilibrium solvers.
 *
 *  @param s The object to set to an equilibrium state
 *  @param XY An integer specifying the two properties to be held constant.
 *  @param estimateEquil integer indicating whether the solver
 *                   should estimate its own initial condition.
 *                   - If 0, the initial mole fraction vector in the
 *                     ThermoPhase object is used as the initial condition.
 *                   - If 1, the initial mole fraction vector is used if the
 *                     element abundances are satisfied.
 *                   - if -1, the initial mole fraction vector is thrown out,
 *                     and an estimate is formulated.
 *  @param printLvl Determines the amount of printing that gets sent to stdout
 *                  from the vcs package (Note, you may have to compile with
 *                  debug flags to get some printing).
 *  @param solver The equilibrium solver to use. If solver = 0, the ChemEquil
 *                solver will be used, and if solver = 1, the
 *                vcs_MultiPhaseEquil solver will be used (slower than
 *                ChemEquil, but more stable). If solver < 0 (default, then
 *                ChemEquil will be tried first, and if it fails
 *                vcs_MultiPhaseEquil will be tried.
 *  @param rtol  Relative tolerance of the solve. Defaults to 1.0E-9.
 *  @param maxsteps The maximum number of steps to take to find the solution.
 *  @param maxiter For the MultiPhaseEquil solver only, this is the maximum
 *                 number of outer temperature or pressure iterations to take
 *                 when T and/or P is not held fixed.
 *  @param loglevel Controls amount of diagnostic output. loglevel
 *                  = 0 suppresses diagnostics, and increasingly-verbose
 *                  messages are written as loglevel increases. The
 *                  messages are written to a file in HTML format for viewing
 *                  in a web browser. @see HTML_logs
 *
 *  @ingroup equilfunctions
 */
int vcs_equilibrate(thermo_t& s, const char* XY,
                    int estimateEquil = 0, int printLvl = 0,
                    int solver = -1, doublereal rtol = 1.0e-9,
                    int maxsteps = VCS_MAXSTEPS,
                    int maxiter = 100, int loglevel = -99);

//!  Set a multi-phase chemical solution to chemical equilibrium.
/*!
 *  This function uses the vcs_MultiPhaseEquil interface to the vcs solver.
 *  The function uses the element abundance vector that is currently
 *  consistent with the composition within the phases themselves. Two other
 *  thermodynamic quantities, determined by the XY string,  are held constant
 *  during the equilibration.
 *
 *  @param s The object to set to an equilibrium state
 *  @param XY A character string representing the unknowns to be held constant
 *  @param estimateEquil integer indicating whether the solver
 *                   should estimate its own initial condition.
 *                   - If 0, the initial mole fraction vector in the
 *                     ThermoPhase object is used as the initial condition.
 *                   - If 1, the initial mole fraction vector is used if the
 *                     element abundances are satisfied.
 *                   - If -1, the initial mole fraction vector is thrown out,
 *                     and an estimate is formulated.
 *  @param printLvl Determines the amount of printing that gets sent to stdout
 *                  from the vcs package (Note, you may have to compile with
 *                  debug flags to get some printing).
 *  @param solver   Determines which solver is used.
 *                 - 1 MultiPhaseEquil solver
 *                 - 2 VCSnonideal Solver (default)
 *  @param rtol  Relative tolerance of the solve. Defaults to 1.0E-9.
 *  @param maxsteps The maximum number of steps to take to find the solution.
 *  @param maxiter For the MultiPhaseEquil solver only, this is the maximum
 *                 number of outer temperature or pressure iterations to take
 *                 when T and/or P is not held fixed.
 *  @param loglevel Controls amount of diagnostic output. loglevel
 *                  = 0 suppresses diagnostics, and increasingly-verbose
 *                  messages are written as loglevel increases. The
 *                  messages are written to a file in HTML format for viewing
 *                  in a web browser. @see HTML_logs
 *
 *  @ingroup equilfunctions
 */
int vcs_equilibrate(MultiPhase& s, const char* XY,
                    int estimateEquil = 0, int printLvl = 0,
                    int solver = 2,
                    doublereal rtol = 1.0e-9, int maxsteps = VCS_MAXSTEPS,
                    int maxiter = 100, int loglevel = -99);

//!  Set a multi-phase chemical solution to chemical equilibrium.
/*!
 *  This function uses the vcs_MultiPhaseEquil interface to the vcs solver.
 *  The function uses the element abundance vector that is currently
 *  consistent with the composition within the phases themselves. Two other
 *  thermodynamic quantities, determined by the XY string,  are held constant
 *  during the equilibration.
 *
 *  @param s The MultiPhase object to be set to an equilibrium state
 *  @param ixy An integer specifying the two properties to be held constant.
 *  @param estimateEquil integer indicating whether the solver
 *                   should estimate its own initial condition.
 *                   - If 0, the initial mole fraction vector in the
 *                     ThermoPhase object is used as the initial condition.
 *                   - If 1, the initial mole fraction vector is used if the
 *                     element abundances are satisfied.
 *                   - if -1, the initial mole fraction vector is thrown out,
 *                     and an estimate is formulated.
 *  @param printLvl Determines the amount of printing that gets sent to stdout
 *                  from the vcs package (Note, you may have to compile with
 *                  debug flags to get some printing).
 *  @param solver   Determines which solver is used.
 *                 - 1 MultiPhaseEquil solver
 *                 - 2 VCSnonideal Solver (default)
 *  @param rtol  Relative tolerance of the solve. Defaults to 1.0E-9.
 *  @param maxsteps The maximum number of steps to take to find the solution.
 *  @param maxiter For the MultiPhaseEquil solver only, this is
 *                 the maximum number of outer temperature or
 *                 pressure iterations to take when T and/or P is
 *                 not held fixed.
 *  @param loglevel Controls amount of diagnostic output. loglevel
 *                  = 0 suppresses diagnostics, and increasingly-verbose
 *                  messages are written as loglevel increases. The
 *                  messages are written to a file in HTML format for viewing
 *                  in a web browser. @see HTML_logs
 *
 *  @ingroup equilfunctions
 */
int vcs_equilibrate_1(MultiPhase& s, int ixy,
                      int estimateEquil = 0, int printLvl = 0,
                      int solver = 2,
                      doublereal rtol = 1.0e-9, int maxsteps = VCS_MAXSTEPS,
                      int maxiter = 100, int loglevel = -99);

//! Determine the phase stability of a single phase given the current conditions
//! in a MultiPhase object
/*!
 *  @param s         The MultiPhase object to be set to an equilibrium state
 *  @param iphase    Phase index within the multiphase object to be
 *                   tested for stability.
 *  @param funcStab  Function value that tests equilibrium. > 0 indicates stable
 *                   < 0 indicates unstable
 *  @param printLvl   Determines the amount of printing that gets sent to
 *                  stdout from the vcs package (Note, you may have to compile
 *                  with debug flags to get some printing).
 *  @param loglevel Controls amount of diagnostic output. loglevel
 *                  = 0 suppresses diagnostics, and increasingly-verbose
 *                  messages are written as loglevel increases. The
 *                  messages are written to a file in HTML format for viewing
 *                  in a web browser. @see HTML_logs
 */
int vcs_determine_PhaseStability(MultiPhase& s, int iphase,
                                 double& funcStab, int printLvl, int loglevel);

}

//! Contains classes and functions implementing the VCS multi-phase
//! equilibrium solver.
namespace VCSnonideal
{
class VCS_PROB;
class VCS_SOLVE;

//! Translate a MultiPhase object into a VCS_PROB problem definition object
/*!
 *  @param mphase MultiPhase object that is the source for all of the information
 *  @param vprob  VCS_PROB problem definition that gets all of the information
 *
 *  Note, both objects share the underlying Thermophase objects. So, neither
 *  can be const objects.
 */
int vcs_Cantera_to_vprob(Cantera::MultiPhase* mphase,
                         VCSnonideal::VCS_PROB* vprob);

//! Translate a MultiPhase information into a VCS_PROB problem definition object
/*!
 *  This version updates the problem statement information only. All species and
 *  phase definitions remain the same.
 *
 *  @param mphase MultiPhase object that is the source for all of the information
 *  @param vprob  VCS_PROB problem definition that gets all of the information
 */
int vcs_Cantera_update_vprob(Cantera::MultiPhase* mphase,
                             VCSnonideal::VCS_PROB* vprob);

//! %Cantera's Interface to the Multiphase chemical equilibrium solver.
/*!
 *  Class vcs_MultiPhaseEquil is designed to be used to set a mixture
 *  containing one or more phases to a state of chemical equilibrium.
 *
 * Note, as currently constructed, the underlying ThermoPhase objects are
 * shared between the MultiPhase object and this object. Therefore, mix is not
 * a const argument, and the return parameters are contained in underlying
 * ThermoPhase objects.
 *
 * @ingroup equilfunctions
 */
class vcs_MultiPhaseEquil
{
public:
    vcs_MultiPhaseEquil();

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
    vcs_MultiPhaseEquil(Cantera::MultiPhase* mix, int printLvl);

    virtual ~vcs_MultiPhaseEquil();

    //! Return the index of the ith component
    /*!
     *  Returns the index of the ith component in the equilibrium
     *  calculation. The index refers to the ordering of the species
     *  in the MultiPhase object.
     *
     *  @param m Index of the component. Must be between 0 and the
     *           number of components, which can be obtained from the
     *           numComponents() command.
     */
    size_t component(size_t m) const ;

    //! Get the stoichiometric reaction coefficients for a single
    //! reaction index
    /*!
     * This returns a stoichiometric reaction vector for a single
     * formation reaction for a noncomponent species. There are
     * (nSpecies() - nComponents) formation reactions. Each
     * formation reaction will have a value of 1.0 for the species
     * that is being formed, and the other non-zero coefficients will
     * all involve the components of the mixture.
     *
     * @param rxn Reaction number.
     * @param nu  Vector of coefficients for the formation reaction. Length is
     *            equal to the number of species in the MultiPhase object.
     */
    void getStoichVector(size_t rxn, Cantera::vector_fp& nu);

    //! return the number of iterations
    int iterations() const {
        return m_iter;
    }

    //! Equilibrate the solution using the current element abundances
    //! stored in the MultiPhase object
    /*!
     *  Use the vcs algorithm to equilibrate the current multiphase mixture.
     *
     *  @param XY       Integer representing what two thermo quantities are
     *                  held constant during the equilibration
     *  @param estimateEquil integer indicating whether the solver should
     *                  estimate its own initial condition.
     *                  - If 0, the initial mole fraction vector in the
     *                    ThermoPhase object is used as the initial condition.
     *                  - If 1, the initial mole fraction vector is used if
     *                    the element abundances are satisfied.
     *                  - if -1, the initial mole fraction vector is thrown
     *                    out, and an estimate is formulated.
     *  @param printLvl Determines the amount of printing that gets sent to
     *                  stdout from the vcs package (Note, you may have to
     *                  compile with debug flags to get some printing).
     *  @param err      Internal error level
     *  @param maxsteps max steps allowed.
     *  @param loglevel Determines the amount of printing to the HTML
     *                  output file.
     */
    int equilibrate(int XY, int estimateEquil = 0,
                    int printLvl= 0, doublereal err = 1.0e-6,
                    int maxsteps = VCS_MAXSTEPS, int loglevel=-99);

    //! Equilibrate the solution using the current element abundances
    //! stored in the MultiPhase object using constant T and P
    /*!
     *  Use the vcs algorithm to equilibrate the current multiphase mixture.
     *
     *  @param estimateEquil integer indicating whether the solver should
     *                 estimate its own initial condition.
     *                 - If 0, the initial mole fraction vector in the
     *                   ThermoPhase object is used as the initial condition.
     *                 - If 1, the initial mole fraction vector is used if the
     *                   element abundances are satisfied.
     *                 - if -1, the initial mole fraction vector is thrown
     *                   out, and an estimate is formulated.
     *  @param printLvl  Determines the amount of printing that gets sent to
     *                  stdout from the vcs package (Note, you may have to
     *                  compile with debug flags to get some printing).
     *  @param err     Internal error level
     *  @param maxsteps max steps allowed.
     *  @param loglevel Determines the amount of printing to the HTML
     *                  output file.
     */
    int equilibrate_TP(int estimateEquil = 0,
                       int printLvl= 0, doublereal err = 1.0e-6,
                       int maxsteps = VCS_MAXSTEPS, int loglevel=-99);

    //! Equilibrate the solution using the current element abundances
    //! stored in the MultiPhase object using either constant H and P
    //! or constant U and P.
    /*!
     *  Use the vcs algorithm to equilibrate the current multiphase
     *  mixture. The pressure of the calculation is taken from
     *  the current pressure stored with the MultiPhase object.
     *
     *  @param Htarget Value of the total mixture enthalpy or total internal
     *                 energy that will be kept constant. Note, this is and
     *                 must be an extensive quantity.  units = Joules
     *  @param XY      Integer flag indicating what is held constant.
     *                 Must be either HP or UP.
     *  @param Tlow    Lower limit of the temperature. It's an error condition
     *                 if the temperature falls below Tlow.
     *  @param Thigh   Upper limit of the temperature. It's an error condition
     *                 if the temperature goes higher than Thigh.
     *  @param estimateEquil integer indicating whether the solver
     *                 should estimate its own initial condition.
     *                 - If 0, the initial mole fraction vector in the
     *                   ThermoPhase object is used as the initial condition.
     *                 - If 1, the initial mole fraction vector is used if the
     *                   element abundances are satisfied.
     *                 - if -1, the initial mole fraction vector is thrown
     *                   out, and an estimate is formulated.
     *  @param printLvl  Determines the amount of printing that
     *                  gets sent to stdout from the vcs package
     *                  (Note, you may have to compile with debug
     *                   flags to get some printing). See main
     *                   constructor call for meaning of the levels.
     *  @param err     Internal error level
     *  @param maxsteps max steps allowed.
     *  @param loglevel Determines the amount of printing to the HTML
     *                  output file.
     */
    int equilibrate_HP(doublereal Htarget, int XY, double Tlow, double Thigh,
                       int estimateEquil = 0,
                       int printLvl = 0, doublereal err = 1.0E-6,
                       int maxsteps = VCS_MAXSTEPS, int loglevel=-99);

    //! Equilibrate the solution using the current element abundances
    //! stored in the MultiPhase object using constant S and P.
    /*!
     *  Use the vcs algorithm to equilibrate the current multiphase
     *  mixture. The pressure of the calculation is taken from
     *  the current pressure stored with the MultiPhase object.
     *
     *  @param Starget Value of the total mixture entropy that will be kept
     *                 constant. Note, this is and must be an extensive
     *                 quantity.  units = Joules/K
     *  @param Tlow    Lower limit of the temperature. It's an error condition
     *                 if the temperature falls below Tlow.
     *  @param Thigh   Upper limit of the temperature. It's an error condition
     *                 if the temperature goes higher than Thigh.
     *  @param estimateEquil integer indicating whether the solver should
     *                 estimate its own initial condition.
     *                 - If 0, the initial mole fraction vector in the
     *                   %ThermoPhase object is used as the initial condition.
     *                 - If 1, the initial mole fraction vector is used if the
     *                   element abundances are satisfied.
     *                 - If -1, the initial mole fraction vector is thrown
     *                   out, and an estimate is formulated.
     *  @param printLvl  Determines the amount of printing that
     *                  gets sent to stdout from the vcs package
     *                  (Note, you may have to compile with debug
     *                   flags to get some printing). See main
     *                   constructor call for meaning of the levels.
     *  @param err     Internal error level
     *  @param maxsteps max steps allowed.
     *  @param loglevel Determines the amount of printing to the HTML
     *                  output file.
     */
    int equilibrate_SP(doublereal Starget, double Tlow, double Thigh,
                       int estimateEquil = 0,
                       int printLvl = 0, doublereal err = 1.0E-6,
                       int maxsteps = VCS_MAXSTEPS, int loglevel=-99);

    //! Equilibrate the solution using the current element abundances stored
    //! in the MultiPhase object using constant V and constant T, H, U or S.
    /*!
     *  Use the vcs algorithm to equilibrate the current multiphase
     *  mixture. The pressure of the calculation is taken from
     *  the current pressure stored with the MultiPhase object.
     *
     *  @param XY      Integer flag indicating what is held constant.
     *                 Must be either TV, HV, UV, or SV.
     *  @param xtarget Value of the total thermodynamic parameter to
     *                 be held constant in addition to V.
     *                 Note, except for T, this must be an extensive
     *                 quantity.  units = Joules/K or Joules
     *  @param estimateEquil integer indicating whether the solver should
     *                 estimate its own initial condition.
     *                 - If 0, the initial mole fraction vector in the
     *                   ThermoPhase object is used as the initial condition.
     *                 - If 1, the initial mole fraction vector is used if the
     *                   element abundances are satisfied.
     *                 - if -1, the initial mole fraction vector is thrown
     *                   out, and an estimate is formulated.
     *  @param printLvl  Determines the amount of printing that gets sent to
     *                  stdout from the vcs package (Note, you may have to
     *                  compile with debug flags to get some printing). See
     *                  main constructor call for meaning of the levels.
     *  @param err     Internal error level
     *  @param maxsteps max steps allowed.
     *  @param logLevel Determines the amount of printing to the HTML
     *                  output file.
     */
    int equilibrate_TV(int XY, doublereal xtarget,
                       int  estimateEquil = 0,
                       int printLvl = 0, doublereal err = 1.0E-6,
                       int maxsteps = VCS_MAXSTEPS, int logLevel = -99);

    //! Determine the phase stability of a phase at the current conditions
    /*!
     * Equilibration of the solution is not done before the determination is made.
     *
     *  @param iph       Phase number to determine the equilibrium. If the phase
     *                   has a non-zero mole number....
     *  @param funcStab  Value of the phase pop function
     *  @param printLvl  Determines the amount of printing that gets sent to
     *                   stdout from the vcs package (Note, you may have to
     *                   compile with debug flags to get some printing).
     *  @param logLevel Determines the amount of printing to the HTML output file.
     */
    int determine_PhaseStability(int iph, double& funcStab,  int printLvl= 0, int logLevel = -99);

    //! Report the equilibrium answer in a comma separated table format
    /*!
     *  This routine is used for in the test suite.
     *
     *  @param reportFile Base name of the file to get the report.
     *         File name is incremented by 1 for each report.
     */
    void reportCSV(const std::string& reportFile);

    //! reports the number of components in the equilibration problem
    /*!
     *  @return returns the number of components. If an equilibrium
     *          problem hasn't been solved yet, it returns -1.
     */
    size_t numComponents() const;

    //! Reports the number of element constraints in the equilibration problem
    /*!
     *  @return returns the number of element constraints. If an equilibrium
     *          problem hasn't been solved yet, it returns -1.
     */
    size_t numElemConstraints() const;

    // Friend functions
    friend int vcs_Cantera_to_vprob(Cantera::MultiPhase* mphase,
                                    VCSnonideal::VCS_PROB* vprob);
    friend int  vcs_Cantera_update_vprob(Cantera::MultiPhase* mphase,
                                         VCSnonideal::VCS_PROB* vprob);

protected:
    //! Vector that takes into account of the current sorting of the species
    /*!
     *  The index of m_order is the original k value of the species in the
     *  multiphase.  The value of m_order, k_sorted, is the current value of
     *  the species index.
     *
     *  `m_order[korig] = k_sorted`
     */
    Cantera::vector_int m_order;

    //! Object which contains the problem statement
    /*!
     *  The problem statement may contain some subtleties. For example, the
     *  element constraints may be different than just an element conservation
     *  contraint equations. There may be kinetically frozen degrees of
     *  freedom. There may be multiple electrolyte phases with zero charge
     *  constraints. All of these make the problem statement different than
     *  the simple element conservation statement.
     */
    VCSnonideal::VCS_PROB* m_vprob;

    //! Pointer to the MultiPhase mixture that will be equilibrated.
    /*!
     *  Equilibrium solutions will be returned via this variable.
     */
    Cantera::MultiPhase* m_mix;

    //! Print level from the VCSnonlinear package
    /*!
     *  (Note, you may have to compile with debug flags to get some printing).
     *
     *  - 0: No IO from the routine whatsoever
     *  - 1: file IO from reportCSV() carried out. One line print statements
     *    from equilibrate_XY() functions
     *  - 2: Problem statement information from vcs_Cantera_update_vprob();
     *    Final state of the system from vcs_solve_TP()
     *  - 3: Several more setup tables; Problem initialization routine
     *  - 4: One table for each iteration within vcs_solve_Tp()
     *  - 5: Multiple tables for each iteration within vcs_solve_TP()
     */
    int m_printLvl;

    //! Stoichiometric matrix
    Cantera::DenseMatrix m_N;

    //! Iteration Count
    int m_iter;

    //! Vector of indices for species that are included in the calculation.
    /*!
     *   This is used to exclude pure-phase species with invalid thermo data
     */
    Cantera::vector_int m_species;

    //! Pointer to the object that does all of the equilibration work.
    /*!
     * VCS_SOLVE will have different ordering for species and element constraints
     * than this object or the VCS_PROB object. This object owns the pointer.
     */
    VCSnonideal::VCS_SOLVE* m_vsolvePtr;
};

//! Global hook for turning on and off time printing.
/*!
 * Default is to allow printing. But, you can assign this to zero globally to
 * turn off all time printing. This is helpful for test suite purposes where
 * you are interested in differences in text files.
 */
extern int vcs_timing_print_lvl;

}
#endif
