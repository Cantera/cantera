/**
 *  @file GeneralSpeciesThermo.h
 *  @deprecated To be removed after Cantera 2.3.
 */

#ifndef CT_GENERALSPECIESTHERMO_H
#define CT_GENERALSPECIESTHERMO_H

#pragma message "Deprecated. GeneralSpeciesThermo.h will be removed after Cantera 2.3. Use MultiSpeciesThermo.h instead."

#include "cantera/base/ct_defs.h"
#include "MultiSpeciesThermo.h"

namespace Cantera
{

//! @deprecated To be removed after Cantera 2.3. Use class MultiSpeciesThermo
//!     instead.
class GeneralSpeciesThermo : public MultiSpeciesThermo
{
    GeneralSpeciesThermo() {
        warn_deprecated("class GeneralSpeciesThermo", "To be removed after"
            " Cantera 2.3. Use class MultiSpeciesThermo instead.");
    }
};

}

#endif
