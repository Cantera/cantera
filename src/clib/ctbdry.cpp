/**
 * @file ctbdry.cpp
 */
#define CANTERA_USE_INTERNAL
#include "ctbdry.h"

// Cantera includes
#include "cantera/oneD/OneDim.h"
#include "cantera/oneD/Inlet1D.h"
#include "cantera/kinetics/InterfaceKinetics.h"
#include "Cabinet.h"

using namespace std;
using namespace Cantera;

typedef Cabinet<Bdry1D> BoundaryCabinet;
template<> BoundaryCabinet* BoundaryCabinet::__storage = 0;

extern "C" {

    int DLL_EXPORT bndry_new(int itype)
    {
        Bdry1D* s;
        switch (itype) {
        case 1:
            s = new Inlet1D();
            break;
        case 2:
            s = new Symm1D();
            break;
        case 3:
            s = new Surf1D();
            break;
        case 4:
            s = new ReactingSurf1D();
            break;
        default:
            return -2;
        }
        int i = BoundaryCabinet::add(s);
        return i;
    }

    int DLL_EXPORT bndry_del(int i)
    {
        BoundaryCabinet::del(i);
        return 0;
    }

    double DLL_EXPORT bndry_temperature(int i)
    {
        return BoundaryCabinet::item(i).temperature();
    }

    int DLL_EXPORT bndry_settemperature(int i, double t)
    {
        try {
            BoundaryCabinet::item(i).setTemperature(t);
        } catch (CanteraError) {
            return -1;
        }
        return 0;
    }

    double DLL_EXPORT bndry_spreadrate(int i)
    {
        try {
            return dynamic_cast<Inlet1D*>(&BoundaryCabinet::item(i))->spreadRate();
        } catch (CanteraError) {
            return -1;
        }
        return 0;
    }

    int DLL_EXPORT bndry_setSpreadRate(int i, double v)
    {
        try {
            dynamic_cast<Inlet1D*>(&BoundaryCabinet::item(i))->setSpreadRate(v);
        } catch (CanteraError) {
            return -1;
        }
        return 0;
    }

    int DLL_EXPORT bndry_setmdot(int i, double mdot)
    {
        try {
            BoundaryCabinet::item(i).setMdot(mdot);
        } catch (CanteraError) {
            return -1;
        }
        return 0;
    }


    double DLL_EXPORT bndry_mdot(int i)
    {
        return BoundaryCabinet::item(i).mdot();
    }

    int DLL_EXPORT bndry_setxin(int i, double* xin)
    {
        try {
            BoundaryCabinet::item(i).setMoleFractions(xin);
        } catch (CanteraError) {
            return -1;
        }
        return 0;
    }

    int DLL_EXPORT bndry_setxinbyname(int i, char* xin)
    {
        try {
            BoundaryCabinet::item(i).setMoleFractions(string(xin));
        } catch (CanteraError) {
            return -1;
        }
        return 0;
    }

    int DLL_EXPORT surf_setkinetics(int i, int j)
    {
        try {
            ReactingSurf1D* srf =
                dynamic_cast<ReactingSurf1D*>(&BoundaryCabinet::item(i));
            InterfaceKinetics* k =
                dynamic_cast<InterfaceKinetics*>(&Cabinet<Kinetics>::item(j));
            srf->setKineticsMgr(k);
        } catch (CanteraError) {
            return -1;
        }
        return 0;
    }
}
