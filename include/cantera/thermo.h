/**
 * @file thermo.h
 *
 * Support for thermo property calculation from C++ application programs.
 * This header file includes several headers from the Cantera kernel needed
 * to evaluate thermo properties.
 */

#ifndef CT_THERMO_INCL
#define CT_THERMO_INCL

#include "thermo/ThermoFactory.h"
#include "thermo/SurfPhase.h"
#include "thermo/EdgePhase.h"

#include "thermo/GibbsExcessVPSSTP.h"
#include "thermo/MargulesVPSSTP.h"

#include "electrolyteThermo.h"

#include "thermo/LatticePhase.h"
#include "thermo/LatticeSolidPhase.h"

#endif
