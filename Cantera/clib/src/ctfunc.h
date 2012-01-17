/**
 * @file ctfunc.h
 */
#ifndef CTC_FUNC1_H
#define CTC_FUNC1_H

#include "clib_defs.h"

extern "C" {
    EEXXTT int DLL_CPREFIX func_new(int type, size_t n, size_t lenp, double* p);
    EEXXTT int DLL_CPREFIX func_del(int i);
    EEXXTT int DLL_CPREFIX func_copy(int i);
    EEXXTT int DLL_CPREFIX func_assign(int i, int j);
    EEXXTT double DLL_CPREFIX func_value(int i, double t);
    EEXXTT int DLL_CPREFIX func_derivative(int i);
    EEXXTT int DLL_CPREFIX func_duplicate(int i);
    EEXXTT int DLL_CPREFIX func_write(int i, size_t lennm, const char* arg, char* nm);
}

#endif
