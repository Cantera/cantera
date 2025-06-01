/**
 * @file clib3_defs.h
 *
 * @warning  This module is an experimental part of the %Cantera API and
 *      may be changed or removed without notice.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

/**
 * @defgroup clibGroup Experimental CLib API
 * @brief %Cantera's experimental CLib interface.
 */

#ifndef __CLIB3_DEFS_H__
#define __CLIB3_DEFS_H__

#include "cantera/base/config.h"
#include <stdlib.h>

// Values returned for error conditions
#ifndef ERR
# define ERR -999
#endif

#ifndef DERR
# define DERR -999.999
#endif

// // Used by external logger
// enum LogLevel { INFO, WARN , ERROR };

// //! Represents a callback that is invoked to produce log output.
// //! TODO: Only needed in the main CLib library. Should be moved once the
// //! traditional CLib is removed.
// typedef void
//     (*LogCallback)(enum LogLevel logLevel, const char* category, const char* message);

#endif
