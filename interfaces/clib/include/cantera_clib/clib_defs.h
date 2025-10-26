/**
 *  @file clib_defs.h
 *
 *  @warning  The generated CLib API is an experimental part of %Cantera and
 *      may be changed without notice.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

/**
 *  @defgroup CAPIindex Generated C Libraries
 *  @brief %Cantera's generated C interface (CLib).
 *  @warning  The generated CLib API is an experimental part of %Cantera and
 *      may be changed without notice.
 *
 *  %Cantera classes and generated CLib libraries have a one-to-one relationship.
 *  This section provides a list of all currently available libraries.
 */

#ifndef CLIB_DEFS_H
#define CLIB_DEFS_H

#include <stdlib.h>

// Values returned for error conditions
#ifndef ERR
    #define ERR -999
#endif

#ifndef DERR
    #define DERR -999.999
#endif

// Used by external logger
enum LogLevel { INFO, WARN , ERROR };

//! Represents a callback that is invoked to produce log output.
typedef void
    (*LogCallback)(enum LogLevel logLevel, const char* category, const char* message);

#endif
