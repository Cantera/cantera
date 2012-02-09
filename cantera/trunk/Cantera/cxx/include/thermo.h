/**
 * @file thermo.h
 *
 * Support for thermo property calculation from C++ application programs.
 * This header file includes several headers from the Cantera kernel needed
 * to evaluate thermo properties.
 */

#ifndef CT_THERMO_INCL
#define CT_THERMO_INCL

#include "kernel/ThermoFactory.h"
#include "importPhase.h"
#include "kernel/SurfPhase.h"
#include "kernel/EdgePhase.h"


#ifdef WITH_IDEAL_SOLUTIONS

#include "kernel/GibbsExcessVPSSTP.h"  
#include "kernel/MargulesVPSSTP.h"

#endif

#ifdef WITH_ELECTROLYTES

#include "electrolyteThermo.h"

#endif


#ifdef WITH_LATTICE_SOLID

#include "kernel/LatticePhase.h"
#include "kernel/LatticeSolidPhase.h"

#endif

#ifdef WITH_PURE_FLUIDS


#endif



#endif
