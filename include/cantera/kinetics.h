/**
 * @file kinetics.h
 *
 * Support for chemical kinetics calculation from C++ application programs.
 * This header file includes headers needed to create and use objects for
 * evaluating chemical kinetic mechanisms.
 */

#ifndef CXX_INCL_KINETICS
#define CXX_INCL_KINETICS

#pragma message("warning: kinetics.h is deprecated and will be removed after " \
                "Cantera 3.2. Use core.h and/or kinetics/KineticsFactory.h instead.")

#include "kinetics/Kinetics.h"
#include "kinetics/Reaction.h"
#include "kinetics/KineticsFactory.h"

#endif
