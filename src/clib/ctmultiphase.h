/**
 * @file ctmultiphase.h
 */
#ifndef CTC_MULTIPHASE_H
#define CTC_MULTIPHASE_H

#include "clib_defs.h"

extern "C" {
    CANTERA_CAPI int mix_new();
    CANTERA_CAPI int mix_del(int i);
    CANTERA_CAPI int mix_copy(int i);
    CANTERA_CAPI int mix_assign(int i, int j);
    CANTERA_CAPI int mix_addPhase(int i, int j, double moles);
    CANTERA_CAPI int mix_init(int i);
    CANTERA_CAPI size_t mix_nElements(int i);
    CANTERA_CAPI size_t mix_elementIndex(int i, char* name);
    CANTERA_CAPI size_t mix_speciesIndex(int i, int k, int p);
    CANTERA_CAPI size_t mix_nSpecies(int i);
    CANTERA_CAPI int mix_setTemperature(int i, double t);
    CANTERA_CAPI double mix_temperature(int i);
    CANTERA_CAPI double mix_minTemp(int i);
    CANTERA_CAPI double mix_maxTemp(int i);
    CANTERA_CAPI double mix_charge(int i);
    CANTERA_CAPI double mix_phaseCharge(int i, int p);
    CANTERA_CAPI int mix_setPressure(int i, double p);
    CANTERA_CAPI double mix_pressure(int i);
    CANTERA_CAPI double mix_nAtoms(int i, int k, int m);
    CANTERA_CAPI size_t mix_nPhases(int i);
    CANTERA_CAPI double mix_phaseMoles(int i, int n);
    CANTERA_CAPI int mix_setPhaseMoles(int i, int n, double v);
    CANTERA_CAPI int mix_setMoles(int i, size_t nlen, double* n);
    CANTERA_CAPI int mix_setMolesByName(int i, char* n);
    CANTERA_CAPI double mix_speciesMoles(int i, int k);
    CANTERA_CAPI double mix_elementMoles(int i, int m);
    CANTERA_CAPI double mix_equilibrate(int i, char* XY, double err,
                                        int maxsteps, int maxiter, int loglevel);
    CANTERA_CAPI double mix_vcs_equilibrate(int i, char* XY, int estimateEquil,
                                            int printLvl, int solver,
                                            double rtol, int maxsteps,
                                            int maxiter, int loglevel);
    CANTERA_CAPI int mix_getChemPotentials(int i, size_t lenmu, double* mu);
    CANTERA_CAPI int mix_getValidChemPotentials(int i, double bad_mu,
            int standard, size_t lenmu,
            double* mu);
    CANTERA_CAPI double mix_enthalpy(int i);
    CANTERA_CAPI double mix_entropy(int i);
    CANTERA_CAPI double mix_gibbs(int i);
    CANTERA_CAPI double mix_cp(int i);
    CANTERA_CAPI double mix_volume(int i);

    CANTERA_CAPI size_t mix_speciesPhaseIndex(int i, int k);
    CANTERA_CAPI double mix_moleFraction(int i, int k);
}
#endif
