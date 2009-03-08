#ifndef CTC_SURF_H
#define CTC_SURF_H

// Cantera includes
//#include "surface.h"

//#include "Cabinet.h"
//#include "Storage.h"
#include "clib_defs.h"

extern "C" {

    int DLL_IMPORT surface_new(int ikin);
    int DLL_IMPORT surface_del(int i);

    int DLL_IMPORT surfphase_new();
    int DLL_IMPORT surf_del(int i);
    int DLL_IMPORT surf_copy(int i);
    int DLL_IMPORT surf_assign(int i, int j);
    int DLL_IMPORT surf_addspecies(int i, char* nm, double sz);
    int DLL_IMPORT surf_nspecies(int i);
    int DLL_IMPORT surf_freezespecies(int i);
    int DLL_IMPORT surf_setcoverages(int i, double* c);
    int DLL_IMPORT surf_getcoverages(int i, double* c);
    int DLL_IMPORT surf_setconcentrations(int i, double* c);
    int DLL_IMPORT surf_getconcentrations(int i, double* c);
    int DLL_IMPORT surf_setsitedensity(int i, double s0);
    double DLL_IMPORT surf_sitedensity(int i);
    int DLL_IMPORT surf_doc(int i, char* key, char* value);

    int DLL_IMPORT surfkin_new(int iph, int ith1, int ith2);
    int DLL_IMPORT surfkin_del(int i);
    int DLL_IMPORT surfkin_addreaction(int i, int nr, int* r, 
        int* rst, int* ro, int np, int* p, int* pst,  
        int nrate, double* rateParams);
    int DLL_IMPORT surfkin_nreactions(int i);
    int DLL_IMPORT surfkin_getratesofprogress(int i, double* rop);
    int DLL_IMPORT surfkin_getnetproductionrates(int i, double* sdot);
    int DLL_IMPORT surfkin_integrate(int i, double dt);
    int DLL_IMPORT surfkin_save(int i, char* fname, char* id, char* comment);

    int DLL_IMPORT surface_settolerances(int i, int nr, 
        double* rtol, int na, double* atol);
    int DLL_IMPORT surface_fixspecies(int i, int k, double c=-1.0);
    int DLL_IMPORT surface_solvespecies(int i, int k);
    int DLL_IMPORT surface_setmultiplier(int i, int k, double f);
    double DLL_IMPORT surface_multiplier(int i, int k);
    double DLL_IMPORT surface_temperature(int i);
    int DLL_IMPORT surface_settemperature(int i, double t);
    int DLL_IMPORT surface_setcoverages(int i, double* c);
    int DLL_EXPORT surf_setcoveragesbyname(int i, char* c);
}

#endif
