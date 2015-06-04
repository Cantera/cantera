/**
 * @file vcs_SpeciesProperties.cpp
 */
#include "cantera/equil/vcs_SpeciesProperties.h"

using namespace std;

namespace Cantera
{

vcs_SpeciesProperties::vcs_SpeciesProperties(size_t indexPhase,
        size_t indexSpeciesPhase,
        vcs_VolPhase* owning) :
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

}
