//! @file vcs_SpeciesProperties.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef VCS_SPECIES_PROPERTIES_H
#define VCS_SPECIES_PROPERTIES_H

#include "cantera/base/ct_defs.h"

namespace Cantera
{

class VCS_SPECIES_THERMO;
class vcs_VolPhase;

//! Properties of a single species.
class vcs_SpeciesProperties
{
public:
    size_t IndexPhase;
    size_t IndexSpeciesPhase;
    vcs_VolPhase* OwningPhase;
    size_t NumElements;

    //! Name of the species
    std::string SpName;

    //! Pointer to the thermo structure for this species
    VCS_SPECIES_THERMO* SpeciesThermo;

    //! Molecular Weight of the species (gm/mol)
    double WtSpecies;

    //! Column of the formula matrix, comprising the
    //! element composition of the species
    vector_fp FormulaMatrixCol;

    //! Charge state of the species -> This may be duplication of what's in the
    //! FormulaMatrixCol entries. However, it's prudent to separate it out.
    double Charge;

    //! True if this species belongs to a surface phase
    int SurfaceSpecies;

    /*
     * Various Calculated Quantities that are appropriate to keep copies of at
     * this level.
     */

    //! Partial molar volume of the species
    double VolPM;

    //! Representative value of the mole fraction of this species in a phase.
    //! This value is used for convergence issues and for calculation of
    //! numerical derivatives
    double ReferenceMoleFraction;

    vcs_SpeciesProperties(size_t indexPhase, size_t indexSpeciesPhase,
                          vcs_VolPhase* owning);
    virtual ~vcs_SpeciesProperties() {}
};

}

#endif
