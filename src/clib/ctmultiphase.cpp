/**
 * @file ctmultiphase.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/clib/ctmultiphase.h"

// Cantera includes
#include "cantera/equil/MultiPhase.h"
#include "cantera/thermo/ThermoPhase.h"
#include "clib_utils.h"

using namespace Cantera;

typedef Cabinet<MultiPhase> mixCabinet;
typedef Cabinet<ThermoPhase> ThermoCabinet;

template<> mixCabinet* mixCabinet::s_storage = 0;
template<> ThermoCabinet* ThermoCabinet::s_storage; // defined in ct.cpp

extern "C" {

    int mix_new()
    {
        try {
            return mixCabinet::add(make_shared<MultiPhase>());
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int mix_del(int i)
    {
        try {
            mixCabinet::del(i);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int ct_clearMix()
    {
        try {
            mixCabinet::clear();
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int mix_addPhase(int i, int j, double moles)
    {
        try {
            mixCabinet::at(i)->addPhase(ThermoCabinet::at(j).get(), moles);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int mix_init(int i)
    {
        try {
            mixCabinet::at(i)->init();
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int mix_updatePhases(int i)
    {
        try {
            mixCabinet::at(i)->updatePhases();
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    size_t mix_nElements(int i)
    {
        try {
            return mixCabinet::at(i)->nElements();
        } catch (...) {
            return handleAllExceptions(npos, npos);
        }
    }

    size_t mix_elementIndex(int i, const char* name)
    {
        try {
            return mixCabinet::at(i)->elementIndex(name);
        } catch (...) {
            return handleAllExceptions(npos, npos);
        }
    }

    size_t mix_nSpecies(int i)
    {
        try {
            return mixCabinet::at(i)->nSpecies();
        } catch (...) {
            return handleAllExceptions(npos, npos);
        }
    }

    size_t mix_speciesIndex(int i, int k, int p)
    {
        try {
            auto& mix = mixCabinet::at(i);
            mix->checkPhaseIndex(p);
            mix->checkSpeciesIndex(k);
            return mix->speciesIndex(k, p);
        } catch (...) {
            return handleAllExceptions(npos, npos);
        }
    }

    double mix_nAtoms(int i, int k, int m)
    {
        try {
            auto& mix = mixCabinet::at(i);
            mix->checkSpeciesIndex(k);
            mix->checkElementIndex(m);
            return mixCabinet::at(i)->nAtoms(k,m);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    size_t mix_nPhases(int i)
    {
        try {
            return mixCabinet::at(i)->nPhases();
        } catch (...) {
            return handleAllExceptions(npos, npos);
        }
    }

    double mix_phaseMoles(int i, int n)
    {
        try {
            auto& mix = mixCabinet::at(i);
            mix->checkPhaseIndex(n);
            return mix->phaseMoles(n);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    int mix_setPhaseMoles(int i, int n, double v)
    {
        try {
            auto& mix = mixCabinet::at(i);
            mix->checkPhaseIndex(n);
            if (v < 0.0) {
                throw CanteraError("mix_setPhaseMoles",
                                   "Mole number must be non-negative.");
            }
            mix->setPhaseMoles(n, v);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int mix_setMoles(int i, size_t nlen, const double* n)
    {
        try {
            auto& mix = mixCabinet::at(i);
            mix->checkSpeciesArraySize(nlen);
            mix->setMoles(n);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }


    int mix_setMolesByName(int i, const char* n)
    {
        try {
            mixCabinet::at(i)->setMolesByName(n);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int mix_setTemperature(int i, double t)
    {
        try {
            if (t < 0.0) {
                throw CanteraError("mix_setTemperature",
                                   "Temperature must be positive.");
            }
            mixCabinet::at(i)->setTemperature(t);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    double mix_temperature(int i)
    {
        try {
            return mixCabinet::at(i)->temperature();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double mix_minTemp(int i)
    {
        try {
            return mixCabinet::at(i)->minTemp();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double mix_maxTemp(int i)
    {
        try {
            return mixCabinet::at(i)->maxTemp();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double mix_charge(int i)
    {
        try {
            return mixCabinet::at(i)->charge();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double mix_phaseCharge(int i, int p)
    {
        try {
            auto& mix = mixCabinet::at(i);
            mix->checkPhaseIndex(p);
            return mix->phaseCharge(p);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    int mix_setPressure(int i, double p)
    {
        try {
            if (p < 0.0) {
                throw CanteraError("mix_setPressure",
                                   "Pressure must be positive.");
            }
            mixCabinet::at(i)->setPressure(p);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    double mix_pressure(int i)
    {
        try {
            return mixCabinet::at(i)->pressure();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double mix_speciesMoles(int i, int k)
    {
        try {
            auto& mix = mixCabinet::at(i);
            mix->checkSpeciesIndex(k);
            return mix->speciesMoles(k);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double mix_elementMoles(int i, int m)
    {
        try {
            auto& mix = mixCabinet::at(i);
            mix->checkElementIndex(m);
            return mix->elementMoles(m);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double mix_equilibrate(int i, const char* XY, double rtol,
                           int maxsteps, int maxiter, int loglevel)
    {
        try {
            mixCabinet::at(i)->equilibrate(XY, "auto", rtol, maxsteps, maxiter,
                                            0, loglevel);
            return 0;
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    int mix_getChemPotentials(int i, size_t lenmu, double* mu)
    {
        try {
            auto& mix = mixCabinet::at(i);
            mix->checkSpeciesArraySize(lenmu);
            mix->getChemPotentials(mu);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    double mix_enthalpy(int i)
    {
        try {
            return mixCabinet::at(i)->enthalpy();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double mix_entropy(int i)
    {
        try {
            return mixCabinet::at(i)->entropy();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double mix_gibbs(int i)
    {
        try {
            return mixCabinet::at(i)->gibbs();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double mix_cp(int i)
    {
        try {
            return mixCabinet::at(i)->cp();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double mix_volume(int i)
    {
        try {
            return mixCabinet::at(i)->volume();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    size_t mix_speciesPhaseIndex(int i, int k)
    {
        try {
            auto& mix = mixCabinet::at(i);
            mix->checkSpeciesIndex(k);
            return mix->speciesPhaseIndex(k);
        } catch (...) {
            return handleAllExceptions(npos, npos);
        }
    }

    double mix_moleFraction(int i, int k)
    {
        try {
            auto& mix = mixCabinet::at(i);
            mix->checkSpeciesIndex(k);
            return mix->moleFraction(k);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }
}
