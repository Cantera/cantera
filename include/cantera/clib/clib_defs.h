/**
 * @file clib_defs.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CTC_DEFS_H
#define CTC_DEFS_H

#include "cantera/base/config.h"
#include <stdlib.h>

#ifdef _WIN32
// Windows (MSVC or MinGW)
# ifdef CANTERA_USE_INTERNAL
#  define CANTERA_CAPI extern __declspec(dllexport)
# else
#  define CANTERA_CAPI extern __declspec(dllimport)
# endif
#else
// Non-Windows platform
# ifdef CANTERA_USE_INTERNAL
#  define CANTERA_CAPI extern
# else
#  define CANTERA_CAPI
# endif
#endif

// Values returned for error conditions
#ifndef ERR
# define ERR -999
#endif

#ifndef DERR
# define DERR -999.999
#endif

#endif
