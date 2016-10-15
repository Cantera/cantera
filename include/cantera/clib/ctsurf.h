/**
 * @file ctsurf.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CTC_SURF_H
#define CTC_SURF_H

#include "clib_defs.h"

#ifdef __cplusplus
extern "C" {
#endif

    CANTERA_CAPI int surf_setcoverages(int i, const double* c, int norm);
    CANTERA_CAPI int surf_getcoverages(int i, double* c);
    CANTERA_CAPI int surf_setconcentrations(int i, const double* c);
    CANTERA_CAPI int surf_getconcentrations(int i, double* c);
    CANTERA_CAPI int surf_setsitedensity(int i, double s0);
    CANTERA_CAPI double surf_sitedensity(int i);
    CANTERA_CAPI int surf_setcoveragesbyname(int i, const char* c);

#ifdef __cplusplus
}
#endif

#endif
