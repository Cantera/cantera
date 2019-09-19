/**
 * @file ctfunc.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CTC_FUNC1_H
#define CTC_FUNC1_H

#include "clib_defs.h"

#ifdef __cplusplus
extern "C" {
#endif

    CANTERA_CAPI int func_new(int type, size_t n, size_t lenp, const double* p);
    CANTERA_CAPI int func_del(int i);
    CANTERA_CAPI double func_value(int i, double t);
    CANTERA_CAPI int func_derivative(int i);
    CANTERA_CAPI int func_duplicate(int i);
    CANTERA_CAPI int func_write(int i, size_t lennm, const char* arg, char* nm);
    CANTERA_CAPI int ct_clearFunc();

#ifdef __cplusplus
}
#endif

#endif
