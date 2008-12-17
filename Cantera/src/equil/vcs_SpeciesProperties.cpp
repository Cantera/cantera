/**
 * @file vcs_SpeciesProperties.cpp
 */
/*
 * $Id $
 */

#include "vcs_defs.h"
#include "vcs_SpeciesProperties.h"
#include "vcs_VolPhase.h"
#include "vcs_species_thermo.h"
#include "vcs_internal.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>

using namespace std;

namespace VCSnonideal {

/*****************************************************************************
 *
 * constructor():
 */
vcs_SpeciesProperties::vcs_SpeciesProperties(int indexPhase, 
					     int indexSpeciesPhase,
					     vcs_VolPhase *owning) :
  IndexPhase(indexPhase),
  IndexSpeciesPhase(indexSpeciesPhase),
  OwningPhase(owning),
  SpeciesThermo(0),
  WtSpecies(0.0),
  Charge(0.0),
  SurfaceSpecies(0),
  VolPM(0.0),
  ReferenceMoleFraction(1.0E-6)
{
}

/******************************************************************************
 *
 * destructor
 */
vcs_SpeciesProperties::~vcs_SpeciesProperties() 
{
}

/*****************************************************************************
 *
 * Copy Constructor vcs_SpeciesProperties
 */
vcs_SpeciesProperties::vcs_SpeciesProperties(const vcs_SpeciesProperties& b) :
  IndexPhase(b.IndexPhase),
  IndexSpeciesPhase(b.IndexSpeciesPhase),
  OwningPhase(b.OwningPhase),
  NumElements(b.NumElements),
  SpeciesThermo(b.SpeciesThermo),
  WtSpecies(b.WtSpecies),
  Charge(b.Charge),
  SurfaceSpecies(b.SurfaceSpecies),
  VolPM(b.VolPM),
  ReferenceMoleFraction(b.ReferenceMoleFraction)
{
  SpName = b.SpName;
  FormulaMatrixCol = b.FormulaMatrixCol;
}

/*****************************************************************************
 *
 * Assignment operator for vcs_SpeciesProperties
 */
vcs_SpeciesProperties& 
vcs_SpeciesProperties::operator=(const vcs_SpeciesProperties& b)
{
  if (&b != this) {
    IndexPhase              = b.IndexPhase;
    IndexSpeciesPhase       = b.IndexSpeciesPhase;
    OwningPhase             = b.OwningPhase;
    NumElements             = b.NumElements;
    SpName                  = b.SpName;
    WtSpecies             = b.WtSpecies;
    FormulaMatrixCol      = b.FormulaMatrixCol;
    Charge                = b.Charge;
    SurfaceSpecies        = b.SurfaceSpecies;
    VolPM                 = b.VolPM;
    ReferenceMoleFraction = b.ReferenceMoleFraction;
  }
  return *this;
}

/*****************************************************************************/
}

