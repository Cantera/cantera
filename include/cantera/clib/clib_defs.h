/**
 * @file clib_defs.h
 *
 * @warning  This module is an experimental part of the %Cantera API and
 *      may be changed or removed without notice.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CTC_DEFS_H
#define CTC_DEFS_H

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
typedef void
    (*LogCallback)(enum LogLevel logLevel, const char* category, const char* message);

#endif
