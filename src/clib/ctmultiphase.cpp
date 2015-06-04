/**
 * @file ctmultiphase.cpp
 */
#define CANTERA_USE_INTERNAL
#include "ctmultiphase.h"

// Cantera includes
#include "cantera/equil/MultiPhase.h"
#include "Cabinet.h"

using namespace std;
using namespace Cantera;

typedef Cabinet<MultiPhase> mixCabinet;
template<> mixCabinet* mixCabinet::s_storage = 0;

extern "C" {

    int mix_new()
    {
        try {
            MultiPhase* m = new MultiPhase;
            return mixCabinet::add(m);
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

    int mix_clear()
    {
        try {
            mixCabinet::clear();
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int mix_copy(int i)
    {
        try {
            return mixCabinet::newCopy(i);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int mix_assign(int i, int j)
    {
        try {
            return mixCabinet::assign(i,j);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int mix_addPhase(int i, int j, double moles)
    {
        try {
            mixCabinet::item(i).addPhase(&Cabinet<ThermoPhase>::item(j), moles);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int mix_init(int i)
    {
        try {
            mixCabinet::item(i).init();
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    size_t mix_nElements(int i)
    {
        try {
            return mixCabinet::item(i).nElements();
        } catch (...) {
            return handleAllExceptions(npos, npos);
        }
    }

    size_t mix_elementIndex(int i, char* name)
    {
        try {
            return mixCabinet::item(i).elementIndex(name);
        } catch (...) {
            return handleAllExceptions(npos, npos);
        }
    }

    size_t mix_nSpecies(int i)
    {
        try {
            return mixCabinet::item(i).nSpecies();
        } catch (...) {
            return handleAllExceptions(npos, npos);
        }
    }

    size_t mix_speciesIndex(int i, int k, int p)
    {
        try {
            MultiPhase& mix = mixCabinet::item(i);
            mix.checkPhaseIndex(p);
            mix.checkSpeciesIndex(k);
            return mix.speciesIndex(k, p);
        } catch (...) {
            return handleAllExceptions(npos, npos);
        }
    }

    doublereal mix_nAtoms(int i, int k, int m)
    {
        try {
            MultiPhase& mix = mixCabinet::item(i);
            mix.checkSpeciesIndex(k);
            mix.checkElementIndex(m);
            return mixCabinet::item(i).nAtoms(k,m);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    size_t mix_nPhases(int i)
    {
        try {
            return mixCabinet::item(i).nPhases();
        } catch (...) {
            return handleAllExceptions(npos, npos);
        }
    }

    doublereal mix_phaseMoles(int i, int n)
    {
        try {
            MultiPhase& mix = mixCabinet::item(i);
            mix.checkPhaseIndex(n);
            return mix.phaseMoles(n);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    int mix_setPhaseMoles(int i, int n, double v)
    {
        try {
            MultiPhase& mix = mixCabinet::item(i);
            mix.checkPhaseIndex(n);
            if (v < 0.0) {
                throw CanteraError("setPhaseMoles",
                                   "Mole number must be non-negative.");
            }
            mix.setPhaseMoles(n, v);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int mix_setMoles(int i, size_t nlen, double* n)
    {
        try {
            MultiPhase& mix = mixCabinet::item(i);
            mix.checkSpeciesArraySize(nlen);
            mix.setMoles(n);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }


    int mix_setMolesByName(int i, char* n)
    {
        try {
            mixCabinet::item(i).setMolesByName(n);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int mix_setTemperature(int i, double t)
    {
        try {
            if (t < 0.0) {
                throw CanteraError("setTemperature",
                                   "Temperature must be positive.");
            }
            mixCabinet::item(i).setTemperature(t);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    doublereal mix_temperature(int i)
    {
        try {
            return mixCabinet::item(i).temperature();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    doublereal mix_minTemp(int i)
    {
        try {
            return mixCabinet::item(i).minTemp();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    doublereal mix_maxTemp(int i)
    {
        try {
            return mixCabinet::item(i).maxTemp();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    doublereal mix_charge(int i)
    {
        try {
            return mixCabinet::item(i).charge();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    doublereal mix_phaseCharge(int i, int p)
    {
        try {
            MultiPhase& mix = mixCabinet::item(i);
            mix.checkPhaseIndex(p);
            return mix.phaseCharge(p);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    int mix_setPressure(int i, double p)
    {
        try {
            if (p < 0.0) {
                throw CanteraError("setPressure",
                                   "Pressure must be positive.");
            }
            mixCabinet::item(i).setPressure(p);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    doublereal mix_pressure(int i)
    {
        try {
            return mixCabinet::item(i).pressure();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    doublereal mix_speciesMoles(int i, int k)
    {
        try {
            MultiPhase& mix = mixCabinet::item(i);
            mix.checkSpeciesIndex(k);
            return mix.speciesMoles(k);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    doublereal mix_elementMoles(int i, int m)
    {
        try {
            MultiPhase& mix = mixCabinet::item(i);
            mix.checkElementIndex(m);
            return mix.elementMoles(m);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    doublereal mix_equilibrate(int i, char* XY, doublereal rtol, int maxsteps,
                               int maxiter, int loglevel)
    {
        try {
            mixCabinet::item(i).equilibrate(XY, "auto", rtol, maxsteps, maxiter,
                                            0, loglevel);
            return 0;
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    doublereal mix_vcs_equilibrate(int i, char* XY, int estimateEquil,
                                   int printLvl, int solver, doublereal rtol,
                                   int maxsteps, int maxiter, int loglevel)
    {
        try {
            string ssolver;
            if (solver < 0) {
                ssolver = "auto";
            } else if (solver == 1) {
                ssolver = "gibbs";
            } else if (solver == 2) {
                ssolver = "vcs";
            } else {
                throw CanteraError("mix_vcs_equilibrate",
                    "Invalid equilibrium solver specified.");
            }
            mixCabinet::item(i).equilibrate(XY, ssolver, rtol, maxsteps,
                                            maxiter, estimateEquil, loglevel);
            return 0;
        } catch (...) {
            return handleAllExceptions(ERR, ERR);
        }
    }

    int mix_getChemPotentials(int i, size_t lenmu, double* mu)
    {
        try {
            MultiPhase& mix = mixCabinet::item(i);
            mix.checkSpeciesArraySize(lenmu);
            mix.getChemPotentials(mu);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int mix_getValidChemPotentials(int i, double bad_mu,
                                   int standard, size_t lenmu, double* mu)
    {
        try {
            bool st = (standard == 1);
            MultiPhase& mix = mixCabinet::item(i);
            mix.checkSpeciesArraySize(lenmu);
            mixCabinet::item(i).getValidChemPotentials(bad_mu, mu, st);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    double mix_enthalpy(int i)
    {
        try {
            return mixCabinet::item(i).enthalpy();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double mix_entropy(int i)
    {
        try {
            return mixCabinet::item(i).entropy();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double mix_gibbs(int i)
    {
        try {
            return mixCabinet::item(i).gibbs();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double mix_cp(int i)
    {
        try {
            return mixCabinet::item(i).cp();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double mix_volume(int i)
    {
        try {
            return mixCabinet::item(i).volume();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    size_t mix_speciesPhaseIndex(int i, int k)
    {
        try {
            MultiPhase& mix = mixCabinet::item(i);
            mix.checkSpeciesIndex(k);
            return mix.speciesPhaseIndex(k);
        } catch (...) {
            return handleAllExceptions(npos, npos);
        }
    }

    double mix_moleFraction(int i, int k)
    {
        try {
            MultiPhase& mix = mixCabinet::item(i);
            mix.checkSpeciesIndex(k);
            return mix.moleFraction(k);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }
}
