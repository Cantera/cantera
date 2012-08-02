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

    int bndry_new(int itype)
    {
        try {
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
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int bndry_del(int i)
    {
        try {
            BoundaryCabinet::del(i);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    double bndry_temperature(int i)
    {
        try {
            return BoundaryCabinet::item(i).temperature();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    int bndry_settemperature(int i, double t)
    {
        try {
            BoundaryCabinet::item(i).setTemperature(t);
        } catch (...) {
            return Cantera::handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    double bndry_spreadrate(int i)
    {
        try {
            return BoundaryCabinet::get<Inlet1D>(i).spreadRate();
        } catch (...) {
            return Cantera::handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    int bndry_setSpreadRate(int i, double v)
    {
        try {
            BoundaryCabinet::get<Inlet1D>(i).setSpreadRate(v);
        } catch (...) {
            return Cantera::handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    int bndry_setmdot(int i, double mdot)
    {
        try {
            BoundaryCabinet::item(i).setMdot(mdot);
        } catch (...) {
            return Cantera::handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    double bndry_mdot(int i)
    {
        try {
            return BoundaryCabinet::item(i).mdot();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    int bndry_setxin(int i, double* xin)
    {
        try {
            BoundaryCabinet::item(i).setMoleFractions(xin);
        } catch (...) {
            return Cantera::handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    int bndry_setxinbyname(int i, char* xin)
    {
        try {
            BoundaryCabinet::item(i).setMoleFractions(xin);
        } catch (...) {
            return Cantera::handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    int surf_setkinetics(int i, int j)
    {
        try {
            InterfaceKinetics& k = Cabinet<Kinetics>::get<InterfaceKinetics>(j);
            BoundaryCabinet::get<ReactingSurf1D>(i).setKineticsMgr(&k);
        } catch (...) {
            return Cantera::handleAllExceptions(-1, ERR);
        }
        return 0;
    }
}
