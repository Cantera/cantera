/**
 * @file ctsurf.cpp
 */
/*
 *      $Id: ctsurf.cpp,v 1.6 2009/07/11 17:16:09 hkmoffa Exp $
 */

// turn off warnings under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

// clib header information
#define CANTERA_USE_INTERNAL
#include "ctsurf.h"


// Cantera includes
#include "SurfPhase.h"
#include "InterfaceKinetics.h"
#include "ImplicitSurfChem.h"

#include "Cabinet.h"
#include "Storage.h"

using namespace std;
using namespace Cantera;



//Cabinet<Surf1D>*         Cabinet<Surf1D>::__storage = 0;

//inline Surf1D* _surface(int i) {
//   return Cabinet<Surf1D>::cabinet()->item(i);
//}

inline SurfPhase* _surfphase(int n) {
    return (SurfPhase*)Storage::__storage->__thtable[n];
}

inline InterfaceKinetics* _surfkin(int n) {
    return (InterfaceKinetics*)Storage::__storage->__ktable[n];
}


extern "C" {  

//     int DLL_EXPORT surface_new(int ikin) {
//         InterfaceKinetics* sk = 0;
//         if (ikin > 0) sk = _surfkin(ikin);
//         Surf1D* s = new Surf1D(sk);
//         return Cabinet<Surf1D>::cabinet()->add(s);
//     }

//     int DLL_EXPORT surface_del(int i) {
//         Cabinet<Surf1D>::cabinet()->del(i);
//         return 0;
//     }


    int DLL_EXPORT surf_setsitedensity(int i, double s0) {
        _surfphase(i)->setSiteDensity(s0);
        return 0;
    }

    double DLL_EXPORT surf_sitedensity(int i) {
        return _surfphase(i)->siteDensity();
    }

    int DLL_EXPORT surf_setcoverages(int i, double* c) {
        _surfphase(i)->setCoverages(c);
        return 0;
    }

    int DLL_EXPORT surf_setcoveragesbyname(int i, char* c) {
        _surfphase(i)->setCoveragesByName(string(c));
        return 0;
    }

    int DLL_EXPORT surf_getcoverages(int i, double* c) {
        _surfphase(i)->getCoverages(c);
        return 0;
    }

    int DLL_EXPORT surf_setconcentrations(int i, double* c) {
        _surfphase(i)->setConcentrations(c);
        return 0;
    }

    int DLL_EXPORT surf_getconcentrations(int i, double* c) {
        _surfphase(i)->getConcentrations(c);
        return 0;
    }

//     int DLL_EXPORT surface_setcoverages(int i, double* c) {
//         _surface(i)->setCoverages(c);
//         return 0;
//     }
}
