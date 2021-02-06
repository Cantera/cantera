/**
 * @file ctsurf.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

// clib header information
#define CANTERA_USE_INTERNAL
#include "cantera/clib/ctsurf.h"

// Cantera includes
#include "cantera/thermo/SurfPhase.h"
#include "Cabinet.h"

using namespace std;
using namespace Cantera;

typedef Cabinet<ThermoPhase> ThermoCabinet;

extern "C" {

    int surf_setSiteDensity(int i, double s0)
    {
        try {
            ThermoCabinet::get<SurfPhase>(i).setSiteDensity(s0);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    double surf_siteDensity(int i)
    {
        try {
            return ThermoCabinet::get<SurfPhase>(i).siteDensity();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    int surf_setCoverages(int i, const double* c, int norm)
    {
        try {
            if(norm){
                ThermoCabinet::get<SurfPhase>(i).setCoverages(c);
            } else {
                ThermoCabinet::get<SurfPhase>(i).setCoveragesNoNorm(c);
            }
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int surf_setCoveragesByName(int i, const char* c)
    {
        try {
            ThermoCabinet::get<SurfPhase>(i).setCoveragesByName(c);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int surf_getCoverages(int i, double* c)
    {
        try {
            ThermoCabinet::get<SurfPhase>(i).getCoverages(c);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int surf_setConcentrations(int i, const double* c)
    {
        try {
            ThermoCabinet::get<SurfPhase>(i).setConcentrations(c);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int surf_getConcentrations(int i, double* c)
    {
        try {
            ThermoCabinet::get<SurfPhase>(i).getConcentrations(c);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }
}
