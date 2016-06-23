/**
 *  @file SpeciesThermo.h
 *  @deprecated To be removed after Cantera 2.3.
 */

#ifndef CT_SPECIESTHERMO_H
#define CT_SPECIESTHERMO_H

#pragma message "Deprecated. SpeciesThermo.h will be removed after Cantera 2.3. Use MultiSpeciesThermo.h instead."

#include "cantera/base/ct_defs.h"
#include "MultiSpeciesThermo.h"

namespace Cantera
{

//! @deprecated To be removed after Cantera 2.3. Use class MultiSpeciesThermo
//!     instead.
class SpeciesThermo : public MultiSpeciesThermo
{
    SpeciesThermo() {
        warn_deprecated("class SpeciesThermo", "To be removed after Cantera 2.3. "
            "Use class MultiSpeciesThermo instead.");
    }
};

}

#endif
