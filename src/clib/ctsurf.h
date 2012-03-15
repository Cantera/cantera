/**
 * @file ctsurf.h
 */
#ifndef CTC_SURF_H
#define CTC_SURF_H

#include "clib_defs.h"
#include "cantera/base/config.h"

extern "C" {
    CANTERA_CAPI int surf_setcoverages(int i, double* c);
    CANTERA_CAPI int surf_getcoverages(int i, double* c);
    CANTERA_CAPI int surf_setconcentrations(int i, double* c);
    CANTERA_CAPI int surf_getconcentrations(int i, double* c);
    CANTERA_CAPI int surf_setsitedensity(int i, double s0);
    CANTERA_CAPI double surf_sitedensity(int i);
    CANTERA_CAPI int surf_setcoveragesbyname(int i, char* c);
}

#endif
