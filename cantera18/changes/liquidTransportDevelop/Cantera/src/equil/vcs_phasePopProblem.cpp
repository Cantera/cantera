/**
 * @file vcs_VolPhase.cpp
 */

/* $Id: vcs_VolPhase.cpp 714 2011-04-14 18:25:24Z hkmoffa $ */

/*
 * Copywrite (2005) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
#include "vcs_VolPhase.h"
#include "vcs_internal.h"
#include "vcs_SpeciesProperties.h"
#include "vcs_species_thermo.h"
#include "vcs_solve.h"

#include "ThermoPhase.h"
#include "mix_defs.h"

#include <string>
#include <cstdio>
#include <cstdlib>

namespace VCSnonideal {

  //! Base constructor for the class
  vcs_phasePopProblem::vcs_phasePopProblem(VCS_SOLVE * owningSolverObject = 0)
  {

  }

    //! Copy constructor
    /*!
     *  @param b object to be copied
     */
  vcs_phasePopProblem::vcs_phasePopProblem(const vcs_phasePopProblem& b)
  {

  }

    //! Assignment operator
    /*!
     *  @param b object to be copied
     */
  vcs_phasePopProblem& vcs_phasePopProblem::operator=(const vcs_phasePopProblem& b)
  {

  }

    //! Destructor
  vcs_phasePopProblem::~vcs_phasePopProblem()
  {

  }


}

