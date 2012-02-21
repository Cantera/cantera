/**
 *   @file checkFinite.cpp
 *   Declarations for Multi Dimensional Pointer (mdp) routines that
 *   check for the presence of NaNs in the code.
 */
/*
 * Copyright 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "cantera/base/ct_defs.h"

#include <stdexcept>
#include <string>

#include <cmath>
#include <cstdlib>
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

namespace mdp
{

// Utility routine to check to see that a number is finite.
/*
 *  @param tmp number to be checked
 */
#ifdef _WIN32
void checkFinite(const double tmp)
{
    if (_finite(tmp)) {
        if (_isnan(tmp)) {
            printf("checkFinite() ERROR: we have encountered a nan!\n");
        } else if (_fpclass(tmp) == _FPCLASS_PINF) {
            printf("checkFinite() ERROR: we have encountered a pos inf!\n");
        } else {
            printf("checkFinite() ERROR: we have encountered a neg inf!\n");
        }
        const std::string s = "checkFinite()";
        throw std::range_error(s);
    }
}
#else
void checkFinite(const double tmp)
{
    if (! finite(tmp)) {
        if (isnan(tmp)) {
            printf("checkFinite() ERROR: we have encountered a nan!\n");
        } else if (isinf(tmp) == 1) {
            printf("checkFinite() ERROR: we have encountered a pos inf!\n");
        } else {
            printf("checkFinite() ERROR: we have encountered a neg inf!\n");
        }
        const std::string s = "checkFinite()";
        throw std::range_error(s);
    }
}
#endif


// Utility routine to link checkFinte() to fortran program
/*
 *  This routine is accessible from fortran, usually
 *
 * @param tmp Pointer to the number to check
 *
 * @todo link it into the usual way Cantera handles Fortran calls
 */
extern "C" void checkfinite_(double* tmp)
{
    checkFinite(*tmp);
}


// Utility routine to check that a double stays bounded
/*
 *   This routine checks to see if a number stays bounded. The absolute
 *   value of the number is required to stay below the trigger.
 *
 * @param tmp     Number to be checked
 * @param trigger bounds on the number. Defaults to 1.0E20
 */
void checkMagnitude(const double tmp, const double trigger)
{
    checkFinite(tmp);
    if (fabs(tmp) >= trigger) {
        char sbuf[64];
        sprintf(sbuf, "checkMagnitude() ERROR: Trigger %g exceeded: %g\n", trigger,
                tmp);
        throw std::range_error(sbuf);
    }
}

// Utility routine to check to see that a number is neither zero
// nor indefinite.
/*
 *  This check can be used before using the number in a denominator.
 *
 *  @param tmp number to be checked
 */
#ifdef _WIN32
void checkZeroFinite(const double tmp)
{
    if ((tmp == 0.0) || (! _finite(tmp))) {
        if (tmp == 0.0) {
            printf("checkZeroFinite() ERROR: we have encountered a zero!\n");
        } else if (_isnan(tmp)) {
            printf("checkZeroFinite() ERROR: we have encountered a nan!\n");
        } else if (_fpclass(tmp) == _FPCLASS_PINF) {
            printf("checkZeroFinite() ERROR: we have encountered a pos inf!\n");
        } else {
            printf("checkZeroFinite() ERROR: we have encountered a neg inf!\n");
        }
        char sbuf[64];
        sprintf(sbuf, "checkZeroFinite() ERROR: zero or indef exceeded: %g\n",
                tmp);
        throw std::range_error(sbuf);
    }
}
#else
void checkZeroFinite(const double tmp)
{
    if ((tmp == 0.0) || (! finite(tmp))) {
        if (tmp == 0.0) {
            printf("checkZeroFinite() ERROR: we have encountered a zero!\n");
        } else if (isnan(tmp)) {
            printf("checkZeroFinite() ERROR: we have encountered a nan!\n");
        } else if (isinf(tmp) == 1) {
            printf("checkZeroFinite() ERROR: we have encountered a pos inf!\n");
        } else {
            printf("checkZeroFinite() ERROR: we have encountered a neg inf!\n");
        }
        char sbuf[64];
        sprintf(sbuf, "checkZeroFinite() ERROR: zero or indef exceeded: %g\n",
                tmp);
        throw std::range_error(sbuf);
    }
}
#endif
}
