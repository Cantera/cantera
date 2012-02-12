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


#ifdef WITH_IDEAL_SOLUTIONS

#include "thermo/GibbsExcessVPSSTP.h"
#include "thermo/MargulesVPSSTP.h"

#endif

#ifdef WITH_ELECTROLYTES

#include "electrolyteThermo.h"

#endif


#ifdef WITH_LATTICE_SOLID

#include "thermo/LatticePhase.h"
#include "thermo/LatticeSolidPhase.h"

#endif

#ifdef WITH_PURE_FLUIDS


#endif



#endif
