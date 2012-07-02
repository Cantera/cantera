/**
 * @file vcs_phasePopProblem 
 *   Header for the object representing each phase within vcs
 */
/*
 * $Id: vcs_VolPhase.h 657 2010-12-18 18:26:16Z hkmoffa $ 
 */
/*
 * Copywrite (2005) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#ifndef VCS_PHASEPOPPROBLEM_H
#define VCS_PHASEPOPPROBLEM_H



#include <vector>
#include <string>

/*
 * Forward references
 */
// Forward reference for ThermoPhase object within the Cantera namespace
namespace Cantera {
  class ThermoPhase;
}

namespace VCSnonideal {

  struct VCS_SPECIES;
  class vcs_SpeciesProperties;
  class VCS_SOLVE;


  //!  Phase information and Phase calculations for vcs.
  /*!
   *
   */
  class vcs_phasePopProblem {
  public:

    /*************************************************************************
     *              FUNCTIONS                                                *
     ************************************************************************/

    //! Base constructor for the class
    vcs_phasePopProblem(VCS_SOLVE * owningSolverObject = 0);

    //! Copy constructor
    /*!
     *  @param b object to be copied
     */
    vcs_phasePopProblem(const vcs_phasePopProblem& b);

    //! Assignment operator
    /*!
     *  @param b object to be copied
     */
    vcs_phasePopProblem& operator=(const vcs_phasePopProblem& b);

    //! Destructor
    ~vcs_phasePopProblem();

    //! Number of phases that are popped or destroyed at one time
    int numPopPhases_;

    std::vector<int> phasePopIDs_;
    
    std::vector<doublereal> xMFVector_;

    std::vector<doublereal> molVolPhases_;
  
    std::vector<doublereal> xMFVectorFixed_;

    std::vector<doublereal> molVolPhasesFixed_;

    std::vector<int> kIndexStart_;

    //! Source of the mole fraction  vectors and phase volumes
    /*!
     *  0   From the current solution vector
     *  1   Optimization problem, whose current solution is storred in xMFVector_;
     *  2   VolPhase object
     *  3   Value of xMFVectorFixed_
     */
    int MFsource;
    

      
  
  };

}

#endif
