/**
 * CLIB3_DEFS - Common definitions of experimental CLib %Cantera interface.
 *
 * @file clib3_defs.h
 *
 * @warning  This module is an experimental part of the %Cantera API and
 *      may be changed or removed without notice.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef __CLIB3_DEFS_H__
#define __CLIB3_DEFS_H__

#include "cantera/base/config.h"
#include <stdlib.h>

// Values returned for error conditions
#ifndef ERR
#define ERR -999
#endif

#ifndef DERR
#define DERR -999.999
#endif

// Used by external logger
enum LogLevel { INFO, WARN, ERROR };

//! Represents a callback that is invoked to produce log output.
typedef void
    (*LogCallback)(enum LogLevel logLevel, const char* category, const char* message);

#endif // __CLIB3_DEFS_H__
