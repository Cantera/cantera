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
/*
 * VCS_PROB is outside of Cantera namespace
 */
namespace VCSnonideal {
  class VCS_PROB;
  class VCS_SOLVE;
}

namespace Cantera {

  int vcs_Cantera_to_vprob(MultiPhase *mphase, VCSnonideal::VCS_PROB *vprob);

  int vcs_Cantera_update_vprob(MultiPhase *mphase, 
			       VCSnonideal::VCS_PROB *vprob);
   
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
    typedef MultiPhase       mix_t;
    typedef size_t           index_t;
    typedef DenseMatrix      matrix_t;

    vcs_MultiPhaseEquil();
    vcs_MultiPhaseEquil(mix_t* mix, int printLvl, bool start=true);

    virtual ~vcs_MultiPhaseEquil();

    int constituent(index_t m) { 
      if (m < m_nel) return m_order[m]; 
      else return -1;
    }

    void getStoichVector(index_t rxn, vector_fp& nu) {
      index_t k;
      nu.resize(m_nsp, 0.0);
      if (rxn > m_nsp - m_nel) return;
      for (k = 0; k < m_nsp; k++) {
	nu[m_order[k]] = m_N(k, rxn);
      }
    }

    int iterations() { return m_iter; }

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

    void reportCSV(const std::string &reportFile);
  
    index_t componentIndex(index_t n) { return m_species[m_order[n]]; }

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
    vector_int m_order;

    VCSnonideal::VCS_PROB *m_vprob;

    //! Pointer to the MultiPhase mixture that will be equilibrated.
    /*!
     *  Solutions will be returned in this variable.
     */
    mix_t *m_mix;
    int m_printLvl;

    matrix_t m_N;

    int m_iter;

    // Vector of indices for species that are included in the
    // calculation.  This is used to exclude pure-phase species
    // with invalid thermo data
    vector_int m_species;

    //! Pointer to the object that does all of the equilibration work.
    /*!
     * This object owns the pointer.
     */
    VCSnonideal::VCS_SOLVE *m_vsolvePtr;
  };

}


#endif
