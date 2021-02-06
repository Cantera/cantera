/**
 * @file vcs_SpeciesProperties.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

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
