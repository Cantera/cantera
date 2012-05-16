/**
 * @file ctbdry.cpp
 */
/*
 *      $Id$
 */


#define CANTERA_USE_INTERNAL
#include "ctbdry.h"


#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif


// Cantera includes
#include "OneDim.h"
#include "Inlet1D.h"
#include "InterfaceKinetics.h"
#include "Cabinet.h"
#include "Storage.h"


using namespace std;
using namespace Cantera;

typedef Cabinet<Bdry1D> BoundaryCabinet;
template<> BoundaryCabinet* BoundaryCabinet::__storage = 0;

//Cabinet<Bdry1D>*  Cabinet<Bdry1D>::__storage = 0;

inline Bdry1D* _bndry(int i) {
    return Cabinet<Bdry1D>::cabinet()->item(i);
}

//inline Phase* _phase(int n) {
//    return Storage::__storage->__phasetable[n];
//}

inline ThermoPhase* _thermo(int n) {
    return Storage::__storage->__thtable[n];
}

inline Kinetics* _kin(int n) {
    return Storage::__storage->__ktable[n];
}

extern "C" {  

    int DLL_EXPORT bndry_new(int itype) {
        Bdry1D* s;
        switch (itype) {
        case 1:
            s = new Inlet1D(); break;
        case 2: 
            s = new Symm1D(); break;
        case 3:
            s = new Surf1D(); break;
        case 4:
            s = new ReactingSurf1D(); break;
        default:
            return -2;
        }
        int i = Cabinet<Bdry1D>::cabinet()->add(s);
        return i;
    }

    int DLL_EXPORT bndry_del(int i) {
        Cabinet<Bdry1D>::cabinet()->del(i);
        return 0;
    }

    double DLL_EXPORT bndry_temperature(int i) {
        return _bndry(i)->temperature();
    }

    int DLL_EXPORT bndry_settemperature(int i, double t) {
        try {
            _bndry(i)->setTemperature(t);
        }
        catch (CanteraError) {return -1;}
        return 0;
    }

    double DLL_EXPORT bndry_spreadrate(int i) {
        try {
            return ((Inlet1D*)_bndry(i))->spreadRate();
        }
        catch (CanteraError) {return -1;}
        return 0;
    }

    int DLL_EXPORT bndry_setSpreadRate(int i, double v) {
        try {
            ((Inlet1D*)_bndry(i))->setSpreadRate(v);
        }
        catch (CanteraError) {return -1;}
        return 0;
    }

    int DLL_EXPORT bndry_setmdot(int i, double mdot) {
        try {
            _bndry(i)->setMdot(mdot);
        }
        catch (CanteraError) {return -1;}
        return 0;
    }


    double DLL_EXPORT bndry_mdot(int i) {
        return _bndry(i)->mdot();
        return 0;
    }

    int DLL_EXPORT bndry_setxin(int i, double* xin) {
        try {
            _bndry(i)->setMoleFractions(xin);
        }
        catch (CanteraError) {return -1;}
        return 0;
    }

    int DLL_EXPORT bndry_setxinbyname(int i, char* xin) {
        try {
            _bndry(i)->setMoleFractions(string(xin));
        }
        catch (CanteraError) {return -1;}
        return 0;
    }

    int DLL_EXPORT surf_setkinetics(int i, int j) {
        try {
            ReactingSurf1D* srf = (ReactingSurf1D*)_bndry(i);
            InterfaceKinetics* k = (InterfaceKinetics*)_kin(j);
            srf->setKineticsMgr(k);
        }
        catch (CanteraError) {return -1;}
        return 0;
    }
}
