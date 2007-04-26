#ifndef CTC_FUNC1_H
#define CTC_FUNC1_H

#include "clib_defs.h"

extern "C" {
    int DLL_IMPORT func_new(int type, int n, int lenp, double* p);
    int DLL_IMPORT func_del(int i);
    int DLL_IMPORT func_copy(int i);
    int DLL_IMPORT func_assign(int i, int j);
    double DLL_IMPORT func_value(int i, double t);
    int DLL_IMPORT func_derivative(int i);
    int DLL_IMPORT func_duplicate(int i);
    int DLL_EXPORT func_write(int i, int lennm, const char* arg, char* nm);
}

#endif
