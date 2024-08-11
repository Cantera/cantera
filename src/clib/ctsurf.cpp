/**
 * @file ctsurf.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

// clib header information
#include "cantera/clib/ctsurf.h"

// Cantera includes
#include "cantera/thermo/SurfPhase.h"
#include "clib_utils.h"

using namespace Cantera;

typedef SharedCabinet<ThermoPhase> ThermoCabinet;
template<> ThermoCabinet* ThermoCabinet::s_storage; // defined in ct.cpp


extern "C" {

    int surf_setSiteDensity(int i, double s0)
    {
        try {
            ThermoCabinet::as<SurfPhase>(i)->setSiteDensity(s0);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    double surf_siteDensity(int i)
    {
        try {
            return ThermoCabinet::as<SurfPhase>(i)->siteDensity();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    int surf_setCoverages(int i, const double* c, int norm)
    {
        try {
            if(norm){
                ThermoCabinet::as<SurfPhase>(i)->setCoverages(c);
            } else {
                ThermoCabinet::as<SurfPhase>(i)->setCoveragesNoNorm(c);
            }
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int surf_setCoveragesByName(int i, const char* c)
    {
        try {
            ThermoCabinet::as<SurfPhase>(i)->setCoveragesByName(c);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int surf_getCoverages(int i, double* c)
    {
        try {
            ThermoCabinet::as<SurfPhase>(i)->getCoverages(c);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int surf_setConcentrations(int i, const double* c)
    {
        try {
            ThermoCabinet::as<SurfPhase>(i)->setConcentrations(c);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int surf_getConcentrations(int i, double* c)
    {
        try {
            ThermoCabinet::as<SurfPhase>(i)->getConcentrations(c);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }
}
