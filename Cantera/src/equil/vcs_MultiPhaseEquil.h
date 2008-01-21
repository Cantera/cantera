/**
 *  @file  vcs_MultiPhase.h
 *  Interface class for the vcsnonlinear solver
 */

/*
 *  $Author$
 *  $Revision$
 *  $Date$
 */


#ifndef VCS_MULTIPHASEEQUIL_H
#define VCS_MULTIPHASEEQUIL_H

#include "ct_defs.h"
#include "MultiPhase.h"


namespace Cantera {

   
  //!  Set a single-phase chemical solution to chemical equilibrium.
  /*!
   *  The function uses the element abundance vector that is 
   *  currently consistent with the composition within the phase
   *  itself. Two other thermodynamic quantities, determined by the
   *  XY string,  are held constant during the equilibration.
   *  This is a convenience function that uses one or the other of
   *  the two chemical equilibrium solvers.
   *
   *  @param s The object to set to an equilibrium state
   *
   *  @param XY An integer specifying the two properties to be held
   *            constant.
   *
   *  @param estimateEquil Boolean indicating whether the solver
   *                   should estimate its own initial condition.
   *                   If false, the initial mole fraction vector
   *                   in the %ThermoPhase object is used as the 
   *                   initial condition.
   *
   *  @param printLvl Determines the amount of printing that
   *                  gets sent to stdout from the vcs package
   *                  (Note, you may have to compile with debug
   *                   flags to get some printing).
   *
   *  @param solver The equilibrium solver to use. If solver = 0,
   *                the ChemEquil solver will be used, and if
   *                solver = 1, the vcs_MultiPhaseEquil solver will
   *                be used (slower than ChemEquil,
   *                but more stable). If solver < 0 (default, then 
   *                ChemEquil will be tried first, and if it fails 
   *                vcs_MultiPhaseEquil will be tried.
   *
   *  @param maxsteps The maximum number of steps to take to find
   *                  the solution.
   *
   *  @param maxiter For the MultiPhaseEquil solver only, this is
   *                 the maximum number of outer temperature or 
   *                 pressure iterations to take when T and/or P is 
   *                 not held fixed.
   *
   *  @param loglevel Controls amount of diagnostic output. loglevel
   *                  = 0 suppresses diagnostics, and increasingly-verbose
   *                  messages are written as loglevel increases. The 
   *                  messages are written to a file in HTML format for viewing 
   *                  in a web browser. @see HTML_logs
   *
   *  @ingroup equilfunctions
   */
  int vcs_equilibrate(thermo_t& s, const char* XY,
		      bool estimateEquil = false, int printLvl = 0,
		      int solver = -1, doublereal rtol = 1.0e-9, 
		      int maxsteps = 1000,
		      int maxiter = 100, int loglevel = -99);


  //!  Set a multi-phase chemical solution to chemical equilibrium.
  /*!
   *  This function uses the vcs_MultiPhaseEquil interface to the
   *  vcs solver.
   *  The function uses the element abundance vector that is 
   *  currently consistent with the composition within the phases
   *  themselves. Two other thermodynamic quantities, determined by the
   *  XY string,  are held constant during the equilibration.
   *
   *  @param s The object to set to an equilibrium state
   *
   *  @param XY A character string representing the unknowns
   *              to be held constant
   *
   *  @param estimateEquil Boolean indicating whether the solver
   *                   should estimate its own initial condition.
   *                   If false, the initial mole fraction vector
   *                   in the %ThermoPhase object is used as the 
   *                   initial condition.
   *
   *  @param printLvl Determines the amount of printing that
   *                  gets sent to stdout from the vcs package
   *                  (Note, you may have to compile with debug
   *                   flags to get some printing).
   *
   *  @param solver   Determines which solver is used. 
   *                 - 1 MultiPhaseEquil solver
   *                 - 2 VCSnonideal Solver (default)
   *
   *  @param maxsteps The maximum number of steps to take to find
   *                  the solution.
   *
   *  @param maxiter For the MultiPhaseEquil solver only, this is
   *                 the maximum number of outer temperature or 
   *                 pressure iterations to take when T and/or P is 
   *                 not held fixed.
   *
   *  @param loglevel Controls amount of diagnostic output. loglevel
   *                  = 0 suppresses diagnostics, and increasingly-verbose
   *                  messages are written as loglevel increases. The 
   *                  messages are written to a file in HTML format for viewing 
   *                  in a web browser. @see HTML_logs
   *
   *  @ingroup equilfunctions
   */
  int vcs_equilibrate(MultiPhase& s, const char* XY, 
		      bool estimateEquil = false, int printLvl = 0,
		      int solver = 2,
		      doublereal rtol = 1.0e-9, int maxsteps = 1000, 
		      int maxiter = 100, int loglevel = -99);

  //!  Set a multi-phase chemical solution to chemical equilibrium.
  /*!
   *  This function uses the vcs_MultiPhaseEquil interface to the
   *  vcs solver.
   *  The function uses the element abundance vector that is 
   *  currently consistent with the composition within the phases
   *  themselves. Two other thermodynamic quantities, determined by the
   *  XY string,  are held constant during the equilibration.
   *
   *  @param s The object to set to an equilibrium state
   *
   *  @param XY An integer specifying the two properties to be held
   *            constant.
   *
   *  @param estimateEquil Boolean indicating whether the solver
   *                   should estimate its own initial condition.
   *                   If false, the initial mole fraction vector
   *                   in the %ThermoPhase object is used as the 
   *                   initial condition.
   *
   *  @param printLvl Determines the amount of printing that
   *                  gets sent to stdout from the vcs package
   *                  (Note, you may have to compile with debug
   *                   flags to get some printing).
   *
   *  @param solver   Determines which solver is used. 
   *                 - 1 MultiPhaseEquil solver
   *                 - 2 VCSnonideal Solver (default)
   *
   *  @param maxsteps The maximum number of steps to take to find
   *                  the solution.
   *
   *  @param maxiter For the MultiPhaseEquil solver only, this is
   *                 the maximum number of outer temperature or 
   *                 pressure iterations to take when T and/or P is 
   *                 not held fixed.
   *
   *  @param loglevel Controls amount of diagnostic output. loglevel
   *                  = 0 suppresses diagnostics, and increasingly-verbose
   *                  messages are written as loglevel increases. The 
   *                  messages are written to a file in HTML format for viewing 
   *                  in a web browser. @see HTML_logs
   *
   *  @ingroup equilfunctions
   */
  int vcs_equilibrate_1(MultiPhase& s, int ixy, 
			bool estimateEquil = false, int printLvl = 0,
			int solver = 2,
			doublereal rtol = 1.0e-9, int maxsteps = 1000, 
			int maxiter = 100, int loglevel = -99);

}

namespace VCSnonideal {

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
  int vcs_Cantera_to_vprob(Cantera::MultiPhase *mphase, 
			   VCSnonideal::VCS_PROB *vprob);

  //! Translate a MultiPhase information into a VCS_PROB problem definition object
  /*!
   *  This version updates the problem statement information only. All species and
   *  phase definitions remain the same.
   *
   *  @param mphase MultiPhase object that is the source for all of the information
   *  @param vprob  VCS_PROB problem definition that gets all of the information
   *
   */
  int vcs_Cantera_update_vprob(Cantera::MultiPhase *mphase, 
			       VCSnonideal::VCS_PROB *vprob);


  //! Cantera's Interface to the Multiphase chemical equilibrium solver.
  /*!
   *  Class MultiPhaseEquil is designed to be used to set a mixture
   *  containing one or more phases to a state of chemical equilibrium. 
   * 
   * @ingroup equilfunctions
   */
  class vcs_MultiPhaseEquil {
  public:
    //! Shorthand for the MultiPhase mixture object used by Cantera
    //! to store information about multiple phases
    typedef Cantera::MultiPhase mix_t;

    //! Typedef for an index variable
    typedef size_t index_t;

    //! Typedef for a dense 2d matrix.
    typedef Cantera::DenseMatrix matrix_t;

    //! Default empty constructor
    vcs_MultiPhaseEquil();

    
    vcs_MultiPhaseEquil(mix_t* mix, int printLvl, bool start=true);

    //! Destructor for the class
    virtual ~vcs_MultiPhaseEquil();

    int constituent(index_t m) { 
      if (m < m_nel) return m_order[m]; 
      else return -1;
    }

    //! Get the stoichiometric matrix for a single reaction index
    /*!
     * This returns a stoichiometric reaction matrix for a single
     * formation reaction.
     */
    void getStoichVector(index_t rxn, Cantera::vector_fp& nu);

    //! return the number of iterations
    int iterations() const { return m_iter; }

    //! Equilibrate the solution using the current element abundances
    //! storred in the MultiPhase object
    /*!
     *  Use the vcs algorithm to equilibrate the current multiphase
     *  mixture.
     *
     *  @param XY  Integer representing what two thermo quantities
     *             are held constant during the equilibration
     *
     *  @param estimateEquil Boolean indicating whether the solver
     *                   should estimate its own initial condition.
     *                   If false, the initial mole fraction vector
     *                   in the %ThermoPhase object is used as the 
     *                   initial condition.
     *
     *  @param printLvl  Determines the amount of printing that
     *                  gets sent to stdout from the vcs package
     *                  (Note, you may have to compile with debug
     *                   flags to get some printing).
     *  @param err     Internal error level
     *  @param maxsteps max steps allowed.
     *  @param  loglevel for 
     */
    int equilibrate(int XY,  bool estimateEquil = false,
		    int printLvl= 0, doublereal err = 1.0e-6, 
		    int maxsteps = 1000, int loglevel=-99);

    //! Equilibrate the solution using the current element abundances
    //! storred in the MultiPhase object using constant T and P
    /*!
     *  Use the vcs algorithm to equilibrate the current multiphase
     *  mixture.
     *
     *  @param estimateEquil Boolean indicating whether the solver
     *                   should estimate its own initial condition.
     *                   If false, the initial mole fraction vector
     *                   in the %ThermoPhase object is used as the 
     *                   initial condition.
     *
     *  @param printLvl  Determines the amount of printing that
     *                  gets sent to stdout from the vcs package
     *                  (Note, you may have to compile with debug
     *                   flags to get some printing).
     *  @param err     Internal error level
     *  @param maxsteps max steps allowed.
     *  @param  loglevel for 
     */
    int equilibrate_TP(bool estimateEquil = false,
		       int printLvl= 0, doublereal err = 1.0e-6, 
		       int maxsteps = 1000, int loglevel=-99);

    //! Equilibrate the solution using the current element abundances
    //! storred in the MultiPhase object using constant H and P
    /*!
     *  Use the vcs algorithm to equilibrate the current multiphase
     *  mixture.
     *
     *  @param estimateEquil Boolean indicating whether the solver
     *                   should estimate its own initial condition.
     *                   If false, the initial mole fraction vector
     *                   in the %ThermoPhase object is used as the 
     *                   initial condition.
     *
     *  @param printLvl  Determines the amount of printing that
     *                  gets sent to stdout from the vcs package
     *                  (Note, you may have to compile with debug
     *                   flags to get some printing).
     *  @param err     Internal error level
     *  @param maxsteps max steps allowed.
     *  @param  loglevel for 
     */
    int equilibrate_HP(doublereal Htarget, int XY, double Tlow, double Thigh,
		       bool estimateEquil = false,
		       int printLvl = 0, doublereal err = 1.0E-6, 
		       int maxsteps = 1000, int loglevel=-99);


    int equilibrate_SP(doublereal Starget, double Tlow, double Thigh,
		       bool estimateEquil = false,
		       int printLvl = 0, doublereal err = 1.0E-6, 
		       int maxsteps = 1000, int loglevel=-99);

    int equilibrate_TV(int XY, doublereal xtarget,
		       bool estimateEquil,
		       int printLvl, doublereal err, 
		       int maxsteps, int loglevel);

    //! Report the equilibrium answer in a comma separated table format
    /*!
     *  This routine is used for in the test suite.
     * 
     *  @param reportFile Base name of the file to get the report.
     *         File name is incremented by 1 for each report.
     */
    void reportCSV(const std::string &reportFile);

    //! reports the number of components in the equilibration problem
    /*!
     *  @return returns the number of components. If an equilibrium
     *          problem hasn't been solved yet, it returns -1.
     */
    int numComponents() const;

    //! Reports the number of element contraints in the equilibration problem
    /*!
     *  @return returns the number of element constraints. If an equilibrium
     *          problem hasn't been solved yet, it returns -1.
     */
    int numElemConstraints() const;
  
    //index_t componentIndex(index_t n) { return m_species[m_order[n]]; }

    friend int vcs_Cantera_to_vprob(Cantera::MultiPhase *mphase, 
				    VCSnonideal::VCS_PROB *vprob);
    friend int  vcs_Cantera_update_vprob(Cantera::MultiPhase *mphase, 
					 VCSnonideal::VCS_PROB *vprob);

  protected:


    //!  Number of elements in the combined element object describing all of the
    //!  phases.
    index_t m_nel;

    //! Number of species in the combined multiphase object
    index_t m_nsp;

    //!  Vector that takes into account of the current sorting of the species
    /*!
     *   The index of m_order is the original k value of the species in the
     *   multiphase.  The value of m_order, k_sorted, is the current value of the
     *   species index.
     *  
     *       m_order[korig] = k_sorted
     */
    Cantera::vector_int m_order;

    //! Object which contains the problem statement
    /*!
     *  The problem statement may contain some subtleties. For example,
     *  the element constraints may be different than just an element 
     *  conservation contraint equations. 
     *  There may be kinetically frozen degrees of freedom. 
     *  There may be multiple electrolyte phases with zero charge constraints.
     *  All of these make the problem statement different than the
     *  simple element conservation statement.
     */
    VCSnonideal::VCS_PROB *m_vprob;

    //! Pointer to the MultiPhase mixture that will be equilibrated.
    /*!
     *  Equilibrium solutions will be returned via this variable.
     */
    mix_t *m_mix;

    //! Print level from the VCSnonlinear package
    /*!
     *  (Note, you may have to compile with debug
     *                   flags to get some printing).
     *
     *    - 0 No IO from the routine whatsoever
     *    - 1 file IO from reportCSV() carried out.
     *        One line print statements from equilibrate_XY() functions
     *    - 2 Problem statement information from vcs_Cantera_update_vprob()
     *        - Final state of the system from vcs_solve_TP()
     *    - 3 Several more setup tables
     *        - Problem initialization routine
     *    - 4 One table for each iteration within vcs_solve_Tp()
     *    - 5 Multiple tables for each iteration within vcs_solve_TP()
     *        - full discussion of decisions made for each variable.
     */
    int m_printLvl;


    //! Stoichiometric matrix
    /*!
     *
     */
    matrix_t m_N;

    int m_iter;

    //! Vector of indices for species that are included in the
    //! calculation. 
    /*!
     *   This is used to exclude pure-phase species
     *   with invalid thermo data
     */
    Cantera::vector_int m_species;

    //! Pointer to the object that does all of the equilibration work.
    /*!
     * VCS_SOLVE will have different ordering for species and element constraints
     * than this object or the VCS_PROB object.
     * This object owns the pointer.
     */
    VCSnonideal::VCS_SOLVE *m_vsolvePtr;
  };

}


#endif
