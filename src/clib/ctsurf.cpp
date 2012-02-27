/**
 * @file ctsurf.cpp
 */
// clib header information
#define CANTERA_USE_INTERNAL
#include "ctsurf.h"

// Cantera includes
#include "cantera/thermo/SurfPhase.h"
#include "cantera/kinetics/InterfaceKinetics.h"
#include "kinetics/ImplicitSurfChem.h"
#include "Cabinet.h"

using namespace std;
using namespace Cantera;

inline SurfPhase* _surfphase(int n)
{
    return dynamic_cast<SurfPhase*>(&Cabinet<ThermoPhase>::item(n));
}

inline InterfaceKinetics* _surfkin(int n)
{
    return dynamic_cast<InterfaceKinetics*>(&Cabinet<Kinetics>::item(n));
}

extern "C" {

    int surf_setsitedensity(int i, double s0)
    {
        _surfphase(i)->setSiteDensity(s0);
        return 0;
    }

    double surf_sitedensity(int i)
    {
        return _surfphase(i)->siteDensity();
    }

    int surf_setcoverages(int i, double* c)
    {
        _surfphase(i)->setCoverages(c);
        return 0;
    }

    int surf_setcoveragesbyname(int i, char* c)
    {
        _surfphase(i)->setCoveragesByName(string(c));
        return 0;
    }

    int surf_getcoverages(int i, double* c)
    {
        _surfphase(i)->getCoverages(c);
        return 0;
    }

    int surf_setconcentrations(int i, double* c)
    {
        _surfphase(i)->setConcentrations(c);
        return 0;
    }

    int surf_getconcentrations(int i, double* c)
    {
        _surfphase(i)->getConcentrations(c);
        return 0;
    }

    //     int surface_setcoverages(int i, double* c) {
    //         _surface(i)->setCoverages(c);
    //         return 0;
    //     }
}
