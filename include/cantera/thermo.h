/**
 * @file thermo.h
 *
 * Support for thermo property calculation from C++ application programs.
 * This header file includes several headers needed to create and use objects
 * which evaluate thermo properties.
 */

#ifndef CT_THERMO_INCL
#define CT_THERMO_INCL

#pragma message("warning: thermo.h is deprecated and will be removed after " \
                "Cantera 3.2. Use core.h and/or thermo/ThermoFactory.h instead.")

#include "thermo/ThermoPhase.h"
#include "thermo/Species.h"
#include "thermo/ThermoFactory.h"
#include "thermo/SpeciesThermoInterpType.h"
#include "thermo/SpeciesThermoFactory.h"

#endif
