/**
 * @file ctsurf.cpp
 */
// clib header information
#define CANTERA_USE_INTERNAL
#include "ctsurf.h"

// Cantera includes
#include "cantera/thermo/SurfPhase.h"
#include "cantera/kinetics/InterfaceKinetics.h"
#include "cantera/kinetics/ImplicitSurfChem.h"
#include "Cabinet.h"

using namespace std;
using namespace Cantera;

typedef Cabinet<ThermoPhase> ThermoCabinet;

extern "C" {

    int surf_setsitedensity(int i, double s0)
    {
        try {
            ThermoCabinet::get<SurfPhase>(i).setSiteDensity(s0);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    double surf_sitedensity(int i)
    {
        try {
            return ThermoCabinet::get<SurfPhase>(i).siteDensity();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    int surf_setcoverages(int i, double* c)
    {
        try {
            ThermoCabinet::get<SurfPhase>(i).setCoverages(c);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int surf_setcoveragesbyname(int i, char* c)
    {
        try {
            ThermoCabinet::get<SurfPhase>(i).setCoveragesByName(c);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int surf_getcoverages(int i, double* c)
    {
        try {
            ThermoCabinet::get<SurfPhase>(i).getCoverages(c);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int surf_setconcentrations(int i, double* c)
    {
        try {
            ThermoCabinet::get<SurfPhase>(i).setConcentrations(c);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int surf_getconcentrations(int i, double* c)
    {
        try {
            ThermoCabinet::get<SurfPhase>(i).getConcentrations(c);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }
}
