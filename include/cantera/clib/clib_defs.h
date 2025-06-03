/**
 * @file clib_defs.h
 *
 * @warning  This module is an experimental part of the %Cantera API and
 *      may be changed or removed without notice.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef __CLIB_DEFS_H__
#define __CLIB_DEFS_H__  // use same guards as clib_generated

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
