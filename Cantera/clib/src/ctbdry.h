#ifndef CTC_BDRY_H
#define CTC_BDRY_H

#include "clib_defs.h"

extern "C" {  

    int DLL_IMPORT bndry_new(int itype);
    int DLL_IMPORT bndry_del(int i);
    double DLL_IMPORT bndry_temperature(int i);
    int DLL_IMPORT bndry_settemperature(int i, double t);
    int DLL_IMPORT bndry_setmdot(int i, double mdot);
    double DLL_IMPORT bndry_mdot(int i);
    int DLL_IMPORT bndry_setxin(int i, double* xin);
    int DLL_IMPORT bndry_setxinbyname(int i, char* xin);
}
#endif
