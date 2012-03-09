/**
 * @file ctmultiphase.cpp
 */
#define CANTERA_USE_INTERNAL
#include "ctmultiphase.h"

// Cantera includes
#include "cantera/equil/equil.h"
#include "cantera/equil/MultiPhase.h"
#include "cantera/equil/MultiPhaseEquil.h"
#include "cantera/equil/vcs_MultiPhaseEquil.h"
#include "Cabinet.h"

using namespace std;
using namespace Cantera;

typedef Cabinet<MultiPhase> mixCabinet;
template<> mixCabinet* mixCabinet::__storage = 0;

static bool checkSpecies(int i, size_t k)
{
    try {
        if (k >= mixCabinet::item(i).nSpecies())
            throw CanteraError("checkSpecies",
                               "illegal species index ("+int2str(k)+") ");
        return true;
    } catch (...) {
        return Cantera::handleAllExceptions(false, false);
    }
}

static bool checkElement(int i, size_t m)
{
    try {
        if (m >= mixCabinet::item(i).nElements())
            throw CanteraError("checkElement",
                               "illegal element index ("+int2str(m)+") ");
        return true;
    } catch (...) {
        return Cantera::handleAllExceptions(false, false);
    }
}

static bool checkPhase(int i, int n)
{
    try {
        if (n < 0 || n >= int(mixCabinet::item(i).nPhases()))
            throw CanteraError("checkPhase",
                               "illegal phase index ("+int2str(n)+") ");
        return true;
    } catch (...) {
        return Cantera::handleAllExceptions(false, false);
    }
}

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
            return mixCabinet::item(i).elementIndex(string(name));
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
            return mixCabinet::item(i).speciesIndex(k, p);
        } catch (...) {
            return handleAllExceptions(npos, npos);
        }
    }

    doublereal mix_nAtoms(int i, int k, int m)
    {
        try {
            bool ok = (checkSpecies(i,k) && checkElement(i,m));
            if (ok) {
                return mixCabinet::item(i).nAtoms(k,m);
            } else {
                return DERR;
            }
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
            if (!checkPhase(i, n)) {
                return DERR;
            }
            return mixCabinet::item(i).phaseMoles(n);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    int mix_setPhaseMoles(int i, int n, double v)
    {
        try {
            if (!checkPhase(i, n)) {
                return ERR;
            }
            if (v < 0.0) {
                return -1;
            }
            mixCabinet::item(i).setPhaseMoles(n, v);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int mix_setMoles(int i, size_t nlen, double* n)
    {
        try {
            if (nlen < mixCabinet::item(i).nSpecies()) {
                throw CanteraError("setMoles","array size too small.");
            }
            mixCabinet::item(i).setMoles(n);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }


    int mix_setMolesByName(int i, char* n)
    {
        try {
            mixCabinet::item(i).setMolesByName(string(n));
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int mix_setTemperature(int i, double t)
    {
        try {
            if (t < 0.0) {
                return -1;
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
            if (!checkPhase(i,p)) {
                return DERR;
            }
            return mixCabinet::item(i).phaseCharge(p);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    int mix_setPressure(int i, double p)
    {
        try {
            if (p < 0.0) {
                return -1;
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
            if (!checkSpecies(i,k)) {
                return DERR;
            }
            return mixCabinet::item(i).speciesMoles(k);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    doublereal mix_elementMoles(int i, int m)
    {
        try {
            if (!checkElement(i,m)) {
                return DERR;
            }
            return mixCabinet::item(i).elementMoles(m);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }


    doublereal mix_equilibrate(int i, char* XY,
                                          doublereal rtol, int maxsteps,
                                          int maxiter, int loglevel)
    {
        try {
            return equilibrate(mixCabinet::item(i), XY,
                               rtol, maxsteps, maxiter, loglevel);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }


    doublereal mix_vcs_equilibrate(int i, char* XY, int estimateEquil,
            int printLvl, int solver,
            doublereal rtol, int maxsteps,
            int maxiter, int loglevel)
    {
        try {
#ifdef WITH_VCSNONIDEAL
            int retn = vcs_equilibrate(mixCabinet::item(i), XY, estimateEquil, printLvl, solver,
                                       rtol, maxsteps, maxiter, loglevel);
#else
            int retn = -1;
            throw CanteraError("mix_vcs_equilibrate",
                               "The VCS NonIdeal equilibrium solver isn't compiled in\n"
                               " To use this feature add export WITH_VCS_NONIDEAL='y' to the preconfig file");
#endif
            return (double) retn;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int mix_getChemPotentials(int i, size_t lenmu, double* mu)
    {
        try {
            if (lenmu < mixCabinet::item(i).nSpecies()) {
                throw CanteraError("getChemPotentials","array too small");
            }
            mixCabinet::item(i).getChemPotentials(mu);
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
            if (lenmu < mixCabinet::item(i).nSpecies()) {
                throw CanteraError("getChemPotentials","array too small");
            }
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
            return mixCabinet::item(i).speciesPhaseIndex(k);
        } catch (...) {
            return handleAllExceptions(npos, npos);
        }
    }

    double mix_moleFraction(int i, int k)
    {
        try {
            return mixCabinet::item(i).moleFraction(k);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }
}
