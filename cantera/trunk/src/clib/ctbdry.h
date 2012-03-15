/**
 * @file ctbdry.h
 */
#ifndef CTC_BDRY_H
#define CTC_BDRY_H

#include "clib_defs.h"

extern "C" {
    CANTERA_CAPI int bndry_new(int itype);
    CANTERA_CAPI int bndry_del(int i);
    CANTERA_CAPI double bndry_temperature(int i);
    CANTERA_CAPI int bndry_settemperature(int i, double t);
    CANTERA_CAPI double bndry_spreadrate(int i);
    CANTERA_CAPI int bndry_setSpreadRate(int i, double v);
    CANTERA_CAPI int bndry_setmdot(int i, double mdot);
    CANTERA_CAPI double bndry_mdot(int i);
    CANTERA_CAPI int bndry_setxin(int i, double* xin);
    CANTERA_CAPI int bndry_setxinbyname(int i, char* xin);
}
#endif
