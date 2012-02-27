/**
 * @file clib_defs.h
 */
#ifndef CTC_DEFS_H
#define CTC_DEFS_H

#include "cantera/base/ct_defs.h"

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
