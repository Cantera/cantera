/**
 * @file ctfunc.h
 */
/*
 *      $Id: ctfunc.h,v 1.4 2009/07/11 17:16:09 hkmoffa Exp $
 */

#ifndef CTC_FUNC1_H
#define CTC_FUNC1_H

#include "clib_defs.h"

extern "C" {
    EEXXTT int DLL_CPREFIX func_new(int type, int n, int lenp, double* p);
    EEXXTT int DLL_CPREFIX func_del(int i);
    EEXXTT int DLL_CPREFIX func_copy(int i);
    EEXXTT int DLL_CPREFIX func_assign(int i, int j);
    EEXXTT double DLL_CPREFIX func_value(int i, double t);
    EEXXTT int DLL_CPREFIX func_derivative(int i);
    EEXXTT int DLL_CPREFIX func_duplicate(int i);
    EEXXTT int DLL_CPREFIX func_write(int i, int lennm, const char* arg, char* nm);
}

#endif
