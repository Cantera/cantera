/**
 * @file equilibrium.h
 * cxx layer - Header file providing support for chemical equilibrium calculations
 *    (see \ref equilfunctions)
 *    @deprecated Equilibrium solvers are directly available through class
 *        Cantera::ThermoPhase and class Cantera::MultiPhase
 */
#ifndef CT_EQUIL_INCL
#define CT_EQUIL_INCL
#pragma message("cantera/equilibrium.h is deprecated")
#include "equil/equil.h"
#include "equil/ChemEquil.h"
#include "equil/MultiPhaseEquil.h"
#include "equil/vcs_MultiPhaseEquil.h"
#endif


