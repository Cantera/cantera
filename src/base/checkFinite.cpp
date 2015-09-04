/**
 *   @file checkFinite.cpp Declarations for routines that check for the
 *       presence of NaNs in the code.
 */
/*
 * Copyright 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "cantera/base/ct_defs.h"

#include <stdexcept>
#include <cstdio>
#include "cantera/base/stringUtils.h"
#include "cantera/base/ctexceptions.h"

// We expect that there will be special casing based on the computer
// system here

#ifdef SOLARIS
#include <ieeefp.h>
#include <sunmath.h>
#endif

// Compiler-dependent names for 'isnan' and 'finite'
#if defined(USE_UNDERSCORE_ISNAN)
    // Windows
    #include <float.h>
    #define isnan(x) _isnan(x)
    #define finite(x) _finite(x)
#elif defined(USE_GLOBAL_ISNAN)
    // From C99
    using ::isnan;
    using ::finite;
#elif defined(USE_STD_ISNAN)
    // From C++11
    using std::isnan;
    #define finite(x) std::isfinite(x)
#endif

namespace Cantera {

void checkFinite(const double tmp)
{
    if (!finite(tmp)) {
        if (isnan(tmp)) {
            printf("checkFinite() ERROR: we have encountered a nan!\n");
        } else if (tmp > 0) {
            printf("checkFinite() ERROR: we have encountered a pos inf!\n");
        } else {
            printf("checkFinite() ERROR: we have encountered a neg inf!\n");
        }
        throw std::range_error("checkFinite()");
    }
}

void checkFinite(const std::string& name, double* values, size_t N)
{
    for (size_t i = 0; i < N; i++) {
        if (!finite(values[i])) {
            std::string message = name + " contains non-finite elements:\n\n";
            for (size_t j = 0; j < N; j++) {
                if (!finite(values[j])) {
                    message += name + "[" + int2str(j) + "] = " +
                               fp2str(values[j]) + "\n";
                }
            }
            throw CanteraError("checkFinite", message);
        }
    }
}

}
