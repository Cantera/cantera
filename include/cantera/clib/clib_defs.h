/**
 * @file clib_defs.h
 *
 * @deprecated  Deprecated in %Cantera 3.2 and to be removed thereafter.
 *      The legacy CLib is superseded by the generated CLib.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CLIB_DEFS_H
#define CLIB_DEFS_H  // use same guards as generated CLib

#include "cantera/base/config.h"
#include <stdlib.h>

// Legacy attribute applied to clib functions. Currently, used only to identify
// functions that should be considered by the 'sourcegen' parser for inclusion in the
// C# interface.
#define CANTERA_CAPI

// Values returned for error conditions
#ifndef ERR
# define ERR -999
#endif

#ifndef DERR
# define DERR -999.999
#endif

// Used by external logger
enum LogLevel { INFO, WARN , ERROR };

//! Represents a callback that is invoked to produce log output.
//! TODO: Only needed in the main CLib library. Should be moved once the
//! legacy CLib is removed.
typedef void
    (*LogCallback)(enum LogLevel logLevel, const char* category, const char* message);

#endif
