/**
 * @file ctsurf.h
 */
#ifndef CTC_SURF_H
#define CTC_SURF_H


#include "clib_defs.h"

#ifdef CANTERA_USE_INTERNAL
#include "cantera/base/config.h"
#else
#include "cantera/base/config.h"
#endif

extern "C" {

    CANTERA_CAPI int surface_new(int ikin);
    CANTERA_CAPI int surface_del(int i);

    CANTERA_CAPI int surfphase_new();
    CANTERA_CAPI int surf_del(int i);
    CANTERA_CAPI int surf_copy(int i);
    CANTERA_CAPI int surf_assign(int i, int j);
    CANTERA_CAPI int surf_addspecies(int i, char* nm, double sz);
    CANTERA_CAPI int surf_nspecies(int i);
    CANTERA_CAPI int surf_freezespecies(int i);
    CANTERA_CAPI int surf_setcoverages(int i, double* c);
    CANTERA_CAPI int surf_getcoverages(int i, double* c);
    CANTERA_CAPI int surf_setconcentrations(int i, double* c);
    CANTERA_CAPI int surf_getconcentrations(int i, double* c);
    CANTERA_CAPI int surf_setsitedensity(int i, double s0);
    CANTERA_CAPI double surf_sitedensity(int i);
    CANTERA_CAPI int surf_doc(int i, char* key, char* value);

    CANTERA_CAPI int surfkin_new(int iph, int ith1, int ith2);
    CANTERA_CAPI int surfkin_del(int i);
    CANTERA_CAPI int surfkin_addreaction(int i, int nr, int* r,
            int* rst, int* ro, int np, int* p, int* pst,
            int nrate, double* rateParams);
    CANTERA_CAPI int surfkin_nreactions(int i);
    CANTERA_CAPI int surfkin_getratesofprogress(int i, double* rop);
    CANTERA_CAPI int surfkin_getnetproductionrates(int i, double* sdot);
    CANTERA_CAPI int surfkin_integrate(int i, double dt);
    CANTERA_CAPI int surfkin_save(int i, char* fname, char* id, char* comment);

    CANTERA_CAPI int surface_settolerances(int i, int nr,
            double* rtol, int na, double* atol);
    CANTERA_CAPI int surface_fixspecies(int i, int k, double c=-1.0);
    CANTERA_CAPI int surface_solvespecies(int i, int k);
    CANTERA_CAPI int surface_setmultiplier(int i, int k, double f);
    CANTERA_CAPI double surface_multiplier(int i, int k);
    CANTERA_CAPI double surface_temperature(int i);
    CANTERA_CAPI int surface_settemperature(int i, double t);
    CANTERA_CAPI int surface_setcoverages(int i, double* c);
    CANTERA_CAPI int surf_setcoveragesbyname(int i, char* c);
}

#endif
