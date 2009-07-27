/**
 * @file ctbdry.h
 */
/*
 *      $Id: ctbdry.h,v 1.4 2009/07/11 17:16:09 hkmoffa Exp $
 */

#ifndef CTC_BDRY_H
#define CTC_BDRY_H

#include "clib_defs.h"

extern "C" {

    EEXXTT int DLL_CPREFIX bndry_new(int itype);
    EEXXTT int DLL_CPREFIX bndry_del(int i);
    EEXXTT double DLL_CPREFIX bndry_temperature(int i);
    EEXXTT int DLL_CPREFIX bndry_settemperature(int i, double t);
    EEXXTT double DLL_CPREFIX bndry_spreadrate(int i);
    EEXXTT int DLL_CPREFIX bndry_setSpreadRate(int i, double v);
    EEXXTT int DLL_CPREFIX bndry_setmdot(int i, double mdot);
    EEXXTT double DLL_CPREFIX bndry_mdot(int i);
    EEXXTT int DLL_CPREFIX bndry_setxin(int i, double* xin);
    EEXXTT int DLL_CPREFIX bndry_setxinbyname(int i, char* xin);
    EEXXTT int DLL_CPREFIX bndry_setkinetics(int i, int j);
}
#endif
