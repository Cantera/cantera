/**
 * @file ctfunc.h
 */
#ifndef CTC_FUNC1_H
#define CTC_FUNC1_H

#include "clib_defs.h"

extern "C" {
    CANTERA_CAPI int func_new(int type, size_t n, size_t lenp, double* p);
    CANTERA_CAPI int func_del(int i);
    CANTERA_CAPI int func_copy(int i);
    CANTERA_CAPI int func_assign(int i, int j);
    CANTERA_CAPI double func_value(int i, double t);
    CANTERA_CAPI int func_derivative(int i);
    CANTERA_CAPI int func_duplicate(int i);
    CANTERA_CAPI int func_write(int i, size_t lennm, const char* arg, char* nm);
}

#endif
