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
    } catch (CanteraError& err) {
        err.save();
        return false;
    }
}

static bool checkElement(int i, size_t m)
{
    try {
        if (m >= mixCabinet::item(i).nElements())
            throw CanteraError("checkElement",
                               "illegal element index ("+int2str(m)+") ");
        return true;
    } catch (CanteraError& err) {
        err.save();
        return false;
    }
}

static bool checkPhase(int i, int n)
{
    try {
        if (n < 0 || n >= int(mixCabinet::item(i).nPhases()))
            throw CanteraError("checkPhase",
                               "illegal phase index ("+int2str(n)+") ");
        return true;
    } catch (CanteraError& err) {
        err.save();
        return false;
    }
}

namespace Cantera
{
int _equilflag(const char* xy);
}

extern "C" {

    int mix_new()
    {
        MultiPhase* m = new MultiPhase;
        return mixCabinet::add(m);
    }

    int mix_del(int i)
    {
        mixCabinet::del(i);
        return 0;
    }

    int mix_copy(int i)
    {
        return mixCabinet::newCopy(i);
    }

    int mix_assign(int i, int j)
    {
        return mixCabinet::assign(i,j);
    }

    int mix_addPhase(int i, int j, double moles)
    {
        mixCabinet::item(i).addPhase(&Cabinet<ThermoPhase>::item(j), moles);
        return 0;
    }

    int mix_init(int i)
    {
        mixCabinet::item(i).init();
        return 0;
    }

    size_t mix_nElements(int i)
    {
        return mixCabinet::item(i).nElements();
    }

    size_t mix_elementIndex(int i, char* name)
    {
        return mixCabinet::item(i).elementIndex(string(name));
    }

    size_t mix_nSpecies(int i)
    {
        return mixCabinet::item(i).nSpecies();
    }

    size_t mix_speciesIndex(int i, int k, int p)
    {
        return mixCabinet::item(i).speciesIndex(k, p);
    }

    doublereal mix_nAtoms(int i, int k, int m)
    {
        bool ok = (checkSpecies(i,k) && checkElement(i,m));
        if (ok) {
            return mixCabinet::item(i).nAtoms(k,m);
        } else {
            return DERR;
        }
    }

    size_t mix_nPhases(int i)
    {
        return mixCabinet::item(i).nPhases();
    }

    doublereal mix_phaseMoles(int i, int n)
    {
        if (!checkPhase(i, n)) {
            return DERR;
        }
        return mixCabinet::item(i).phaseMoles(n);
    }

    int mix_setPhaseMoles(int i, int n, double v)
    {
        if (!checkPhase(i, n)) {
            return ERR;
        }
        if (v < 0.0) {
            return -1;
        }
        mixCabinet::item(i).setPhaseMoles(n, v);
        return 0;
    }

    int mix_setMoles(int i, size_t nlen, double* n)
    {
        try {
            if (nlen < mixCabinet::item(i).nSpecies()) {
                throw CanteraError("setMoles","array size too small.");
            }
            mixCabinet::item(i).setMoles(n);
            return 0;
        } catch (CanteraError& err) {
            err.save();
            return ERR;
        }
    }


    int mix_setMolesByName(int i, char* n)
    {
        try {
            mixCabinet::item(i).setMolesByName(string(n));
            return 0;
        } catch (CanteraError& err) {
            err.save();
            return -1;
        }
    }

    int mix_setTemperature(int i, double t)
    {
        if (t < 0.0) {
            return -1;
        }
        mixCabinet::item(i).setTemperature(t);
        return 0;
    }

    doublereal mix_temperature(int i)
    {
        return mixCabinet::item(i).temperature();
    }

    doublereal mix_minTemp(int i)
    {
        return mixCabinet::item(i).minTemp();
    }

    doublereal mix_maxTemp(int i)
    {
        return mixCabinet::item(i).maxTemp();
    }

    doublereal mix_charge(int i)
    {
        return mixCabinet::item(i).charge();
    }

    doublereal mix_phaseCharge(int i, int p)
    {
        if (!checkPhase(i,p)) {
            return DERR;
        }
        return mixCabinet::item(i).phaseCharge(p);
    }

    int mix_setPressure(int i, double p)
    {
        if (p < 0.0) {
            return -1;
        }
        mixCabinet::item(i).setPressure(p);
        return 0;
    }

    doublereal mix_pressure(int i)
    {
        return mixCabinet::item(i).pressure();
    }

    doublereal mix_speciesMoles(int i, int k)
    {
        if (!checkSpecies(i,k)) {
            return DERR;
        }
        return mixCabinet::item(i).speciesMoles(k);
    }

    doublereal mix_elementMoles(int i, int m)
    {
        if (!checkElement(i,m)) {
            return DERR;
        }
        return mixCabinet::item(i).elementMoles(m);
    }


    doublereal mix_equilibrate(int i, char* XY,
                                          doublereal rtol, int maxsteps,
                                          int maxiter, int loglevel)
    {
        try {
            return equilibrate(mixCabinet::item(i), XY,
                               rtol, maxsteps, maxiter, loglevel);
        } catch (CanteraError& err) {
            err.save();
            return DERR;
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
        } catch (CanteraError& err) {
            err.save();
            return DERR;
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
        } catch (CanteraError& err) {
            err.save();
            return -1;
        }
    }

    int mix_getValidChemPotentials(int i, double bad_mu,
            int standard, size_t lenmu, double* mu)
    {
        bool st = (standard == 1);
        try {
            if (lenmu < mixCabinet::item(i).nSpecies()) {
                throw CanteraError("getChemPotentials","array too small");
            }
            mixCabinet::item(i).getValidChemPotentials(bad_mu, mu, st);
            return 0;
        } catch (CanteraError& err) {
            err.save();
            return -1;
        }
    }

    double mix_enthalpy(int i)
    {
        return mixCabinet::item(i).enthalpy();
    }

    double mix_entropy(int i)
    {
        return mixCabinet::item(i).entropy();
    }

    double mix_gibbs(int i)
    {
        return mixCabinet::item(i).gibbs();
    }

    double mix_cp(int i)
    {
        return mixCabinet::item(i).cp();
    }

    double mix_volume(int i)
    {
        return mixCabinet::item(i).volume();
    }

    size_t mix_speciesPhaseIndex(int i, int k)
    {
        return mixCabinet::item(i).speciesPhaseIndex(k);
    }

    double mix_moleFraction(int i, int k)
    {
        return mixCabinet::item(i).moleFraction(k);
    }
}
