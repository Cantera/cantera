/**
 * @file ctmultiphase.cpp
 */
/*
 *      $Id: ctmultiphase.cpp,v 1.15 2009/07/11 17:16:09 hkmoffa Exp $
 */


#define CANTERA_USE_INTERNAL
#include "ctmultiphase.h"

// Cantera includes
#include "equil.h"
#include "MultiPhase.h"
#include "MultiPhaseEquil.h"
#include "vcs_MultiPhaseEquil.h"

#include "Cabinet.h"
#include "Storage.h"

using namespace std;
using namespace Cantera;

typedef MultiPhase  mix_t;

template<> Cabinet<mix_t>*    Cabinet<mix_t>::__storage = 0;

inline mix_t* _mix(int i) {
    return Cabinet<mix_t>::cabinet()->item(i);
}

inline ThermoPhase* _th(int n) {
    return Storage::__storage->__thtable[n];
}

static bool checkSpecies(int i, int k) {
    try {
        if (k < 0 || k >= _mix(i)->nSpecies()) 
            throw CanteraError("checkSpecies",
                "illegal species index ("+int2str(k)+") ");
        return true;
    }
    catch (CanteraError) {
        return false;
    }
}

static bool checkElement(int i, int m) {
    try {
        if (m < 0 || m >= _mix(i)->nElements()) 
            throw CanteraError("checkElement",
                "illegal element index ("+int2str(m)+") ");
        return true;
    }
    catch (CanteraError) {
        return false;
    }
}

static bool checkPhase(int i, int n) {
    try {
        if (n < 0 || n >= int(_mix(i)->nPhases())) 
            throw CanteraError("checkPhase",
                "illegal phase index ("+int2str(n)+") ");
        return true;
    }
    catch (CanteraError) {
        return false;
    }
}

namespace Cantera {
    int _equilflag(const char* xy);
}

extern "C" {  

    int DLL_EXPORT mix_new() {
        mix_t* m = new MultiPhase;
        return Cabinet<mix_t>::cabinet()->add(m);
    }

    int DLL_EXPORT mix_del(int i) {
        Cabinet<mix_t>::cabinet()->del(i);
        return 0;
    }

    int DLL_EXPORT mix_copy(int i) {
        return Cabinet<mix_t>::cabinet()->newCopy(i);
    }

    int DLL_EXPORT mix_assign(int i, int j) {
        return Cabinet<mix_t>::cabinet()->assign(i,j);
    }

    int DLL_EXPORT mix_addPhase(int i, int j, double moles) {
        _mix(i)->addPhase(_th(j), moles);
        return 0;
    }

    int DLL_EXPORT mix_init(int i) {
        _mix(i)->init();
        return 0;
    }

    int DLL_EXPORT mix_nElements(int i) {
        return _mix(i)->nElements();
     }

    int DLL_EXPORT mix_elementIndex(int i, char* name) {
        return _mix(i)->elementIndex(string(name));
    }

    int DLL_EXPORT mix_nSpecies(int i) {
        return _mix(i)->nSpecies();
    }

    int DLL_EXPORT mix_speciesIndex(int i, int k, int p) {
        return _mix(i)->speciesIndex(k, p);
    }

    doublereal DLL_EXPORT mix_nAtoms(int i, int k, int m) {
        bool ok = (checkSpecies(i,k) && checkElement(i,m));
        if (ok) 
            return _mix(i)->nAtoms(k,m);
        else
            return DERR;
    }

    double DLL_EXPORT mix_nPhases(int i) {
        return _mix(i)->nPhases();
    }

    doublereal DLL_EXPORT mix_phaseMoles(int i, int n) {
        if (!checkPhase(i, n)) return DERR;
        return _mix(i)->phaseMoles(n);
    }

    int DLL_EXPORT mix_setPhaseMoles(int i, int n, double v) {
        if (!checkPhase(i, n)) return ERR;
        if (v < 0.0) return -1;
        _mix(i)->setPhaseMoles(n, v);
        return 0;
    }

    int DLL_EXPORT mix_setMoles(int i, int nlen, double* n) {
        try {
            if (nlen < _mix(i)->nSpecies()) 
                throw CanteraError("setMoles","array size too small.");
            _mix(i)->setMoles(n);
            return 0;
        }
        catch (CanteraError) {
            return ERR;
        }
    }


    int DLL_EXPORT mix_setMolesByName(int i, char* n) {
        try {
            _mix(i)->setMolesByName(string(n));
            return 0;
        }
        catch (CanteraError) { return -1; }
    }

    int DLL_EXPORT mix_setTemperature(int i, double t) {
        if (t < 0.0) return -1;
        _mix(i)->setTemperature(t);
        return 0;
    }

    doublereal DLL_EXPORT mix_temperature(int i) {
        return _mix(i)->temperature();
    }

    doublereal DLL_EXPORT mix_minTemp(int i) {
        return _mix(i)->minTemp();
    }

    doublereal DLL_EXPORT mix_maxTemp(int i) {
        return _mix(i)->maxTemp();
    }

    doublereal DLL_EXPORT mix_charge(int i) {
        return _mix(i)->charge();
    }

    doublereal DLL_EXPORT mix_phaseCharge(int i, int p) {
        if (!checkPhase(i,p)) return DERR;
        return _mix(i)->phaseCharge(p);
    }

    int DLL_EXPORT mix_setPressure(int i, double p) {
        if (p < 0.0) return -1;
        _mix(i)->setPressure(p);
        return 0;
    }

    doublereal DLL_EXPORT mix_pressure(int i) {
        return _mix(i)->pressure();
    }

    doublereal DLL_EXPORT mix_speciesMoles(int i, int k) {
        if (!checkSpecies(i,k)) return DERR;
        return _mix(i)->speciesMoles(k);
    }

    doublereal DLL_EXPORT mix_elementMoles(int i, int m) {
        if (!checkElement(i,m)) return DERR;
        return _mix(i)->elementMoles(m);
    }

    
  doublereal DLL_EXPORT mix_equilibrate(int i, char* XY,
					doublereal rtol, int maxsteps, 
					int maxiter, int loglevel) { 
    try {
      return equilibrate(*_mix(i), XY, 
			 rtol, maxsteps, maxiter, loglevel);
    }
    catch (CanteraError) {
      return DERR;
    }
  }


  doublereal DLL_EXPORT mix_vcs_equilibrate(int i, char* XY, int estimateEquil,
					    int printLvl, int solver,
					    doublereal rtol, int maxsteps,
					    int maxiter, int loglevel) { 
    try {
#ifdef WITH_VCSNONIDEAL
      int retn = vcs_equilibrate(*_mix(i), XY, estimateEquil, printLvl, solver,
				 rtol, maxsteps, maxiter, loglevel);
#else
      int retn = -1;
      throw CanteraError("mix_vcs_equilibrate", 
      "The VCS NonIdeal equilibrium solver isn't compiled in\n"
      " To use this feature add export WITH_VCS_NONIDEAL='y' to the preconfig file");
#endif
      return (double) retn;
    }
    catch (CanteraError) {
      return DERR;
    }
  }

    int DLL_EXPORT mix_getChemPotentials(int i, int lenmu, double* mu) {
        try {
            if (lenmu < _mix(i)->nSpecies()) 
                throw CanteraError("getChemPotentials","array too small");
            _mix(i)->getChemPotentials(mu);
            return 0;
        }
        catch (CanteraError) {
            return -1;
        }
    }

    int DLL_EXPORT mix_getValidChemPotentials(int i, double bad_mu, 
        int standard, int lenmu, double* mu) {
        bool st = (standard == 1);
        try {
            if (lenmu < _mix(i)->nSpecies()) 
                throw CanteraError("getChemPotentials","array too small");
            _mix(i)->getValidChemPotentials(bad_mu, mu, st);
            return 0;
        }
        catch (CanteraError) {
            return -1;
        }
    }

    double DLL_EXPORT mix_enthalpy(int i)  {
        return _mix(i)->enthalpy();
    }

    double DLL_EXPORT mix_entropy(int i)  {
        return _mix(i)->entropy();
    }
        

    double DLL_EXPORT mix_gibbs(int i) {
        return _mix(i)->gibbs();
    }

    double DLL_EXPORT mix_cp(int i) {
        return _mix(i)->cp();
    }

    double DLL_EXPORT mix_volume(int i) {
        return _mix(i)->volume();
    }

    int DLL_EXPORT mix_speciesPhaseIndex(int i, int k) {
        return _mix(i)->speciesPhaseIndex(k);
    }

    double DLL_EXPORT mix_moleFraction(int i, int k) {
        return _mix(i)->moleFraction(k);
    }

}
