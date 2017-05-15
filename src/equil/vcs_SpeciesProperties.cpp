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

}
