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

    const int FourierFuncType = 1;
    const int PolyFuncType = 2;
    const int ArrheniusFuncType = 3;
    const int GaussianFuncType = 4;
    const int SumFuncType = 20;
    const int DiffFuncType = 25;
    const int ProdFuncType = 30;
    const int RatioFuncType = 40;
    const int PeriodicFuncType = 50;
    const int CompositeFuncType = 60;
    const int TimesConstantFuncType = 70;
    const int PlusConstantFuncType = 80;
    const int SinFuncType = 100;
    const int CosFuncType = 102;
    const int ExpFuncType = 104;
    const int PowFuncType = 106;
    const int ConstFuncType = 110;
    const int TabulatedFuncType = 120;

    CANTERA_CAPI int func_new(int type, size_t n, size_t lenp, const double* p);
    CANTERA_CAPI int func_new_basic(const char* type, double c);
    CANTERA_CAPI int func_new_advanced(const char* type, size_t lenp, const double* p);
    CANTERA_CAPI int func_new_compound(const char* type, int a, int b);
    CANTERA_CAPI int func_new_modified(const char* type, int a, double c);
    CANTERA_CAPI int func_del(int i);
    CANTERA_CAPI int func_type(int i, size_t lennm, char* nm);
    CANTERA_CAPI double func_value(int i, double t);
    CANTERA_CAPI int func_derivative(int i);
    CANTERA_CAPI int func_duplicate(int i);
    CANTERA_CAPI int func_write(int i, size_t lennm, const char* arg, char* nm);
    CANTERA_CAPI int func_write3(int i, const char* arg, size_t lennm, char* nm);
    CANTERA_CAPI int ct_clearFunc();

#ifdef __cplusplus
}
#endif

#endif
