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

// We expect that there will be special casing based on the computer
// system here

#ifdef SOLARIS
#include <ieeefp.h>
#include <sunmath.h>
#endif

#ifdef _WIN32
#include <float.h>
#endif

using namespace std;

namespace Cantera {

void checkFinite(const double tmp)
{
#if defined _WIN32
    if (_finite(tmp)) {
        if (_isnan(tmp)) {
            printf("checkFinite() ERROR: we have encountered a nan!\n");
        } else if (_fpclass(tmp) == _FPCLASS_PINF) {
            printf("checkFinite() ERROR: we have encountered a pos inf!\n");
        } else {
            printf("checkFinite() ERROR: we have encountered a neg inf!\n");
        }
        throw std::range_error("checkFinite()");
    }
#elif defined __CYGWIN__
    if (!finite(tmp)) {
        if (isnan(tmp)) {
            printf("checkFinite() ERROR: we have encountered a nan!\n");
        } else if (isinf(tmp) == 1) {
            printf("checkFinite() ERROR: we have encountered a pos inf!\n");
        } else {
            printf("checkFinite() ERROR: we have encountered a neg inf!\n");
        }
        throw std::range_error("checkFinite()");
    }
#else
    if (!::finite(tmp)) {
        if (::isnan(tmp)) {
            printf("checkFinite() ERROR: we have encountered a nan!\n");
        } else if (::isinf(tmp) == 1) {
            printf("checkFinite() ERROR: we have encountered a pos inf!\n");
        } else {
            printf("checkFinite() ERROR: we have encountered a neg inf!\n");
        }
        throw std::range_error("checkFinite()");
    }
#endif
}

}
