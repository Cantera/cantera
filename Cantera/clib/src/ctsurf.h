/**
 * @file ctsurf.h
 */
/*
 *      $Id: ctsurf.h,v 1.3 2009/07/11 17:16:09 hkmoffa Exp $
 */

#ifndef CTC_SURF_H
#define CTC_SURF_H


#include "clib_defs.h"

#ifdef CANTERA_USE_INTERNAL
#include "config.h"
#else
#include "cantera/config.h"
#endif

extern "C" {

    EEXXTT int DLL_CPREFIX surface_new(int ikin);
    EEXXTT int DLL_CPREFIX surface_del(int i);

    EEXXTT int DLL_CPREFIX surfphase_new();
    EEXXTT int DLL_CPREFIX surf_del(int i);
    EEXXTT int DLL_CPREFIX surf_copy(int i);
    EEXXTT int DLL_CPREFIX surf_assign(int i, int j);
    EEXXTT int DLL_CPREFIX surf_addspecies(int i, char* nm, double sz);
    EEXXTT int DLL_CPREFIX surf_nspecies(int i);
    EEXXTT int DLL_CPREFIX surf_freezespecies(int i);
    EEXXTT int DLL_CPREFIX surf_setcoverages(int i, double* c);
    EEXXTT int DLL_CPREFIX surf_getcoverages(int i, double* c);
    EEXXTT int DLL_CPREFIX surf_setconcentrations(int i, double* c);
    EEXXTT int DLL_CPREFIX surf_getconcentrations(int i, double* c);
    EEXXTT int DLL_CPREFIX surf_setsitedensity(int i, double s0);
    EEXXTT double DLL_CPREFIX surf_sitedensity(int i);
    EEXXTT int DLL_CPREFIX surf_doc(int i, char* key, char* value);

    EEXXTT int DLL_CPREFIX surfkin_new(int iph, int ith1, int ith2);
    EEXXTT int DLL_CPREFIX surfkin_del(int i);
    EEXXTT int DLL_CPREFIX surfkin_addreaction(int i, int nr, int* r, 
        int* rst, int* ro, int np, int* p, int* pst,  
        int nrate, double* rateParams);
    EEXXTT int DLL_CPREFIX surfkin_nreactions(int i);
    EEXXTT int DLL_CPREFIX surfkin_getratesofprogress(int i, double* rop);
    EEXXTT int DLL_CPREFIX surfkin_getnetproductionrates(int i, double* sdot);
    EEXXTT int DLL_CPREFIX surfkin_integrate(int i, double dt);
    EEXXTT int DLL_CPREFIX surfkin_save(int i, char* fname, char* id, char* comment);

    EEXXTT int DLL_CPREFIX surface_settolerances(int i, int nr, 
        double* rtol, int na, double* atol);
    EEXXTT int DLL_CPREFIX surface_fixspecies(int i, int k, double c=-1.0);
    EEXXTT int DLL_CPREFIX surface_solvespecies(int i, int k);
    EEXXTT int DLL_CPREFIX surface_setmultiplier(int i, int k, double f);
    EEXXTT double DLL_CPREFIX surface_multiplier(int i, int k);
    EEXXTT double DLL_CPREFIX surface_temperature(int i);
    EEXXTT int DLL_CPREFIX surface_settemperature(int i, double t);
    EEXXTT int DLL_CPREFIX surface_setcoverages(int i, double* c);
    EEXXTT int DLL_EXPORT surf_setcoveragesbyname(int i, char* c);
}

#endif
