#ifndef CTC_MULTIPHASE_H
#define CTC_MULTIPHASE_H

#include "clib_defs.h"

extern "C" {

    int DLL_IMPORT mix_new();
    int DLL_IMPORT mix_del(int i);
    int DLL_IMPORT mix_copy(int i);
    int DLL_IMPORT mix_assign(int i, int j);
    int DLL_IMPORT mix_addPhase(int i, int j, double moles);
    int DLL_IMPORT mix_nElements(int i);
    int DLL_IMPORT mix_elementIndex(int i, char* name);
    int DLL_IMPORT mix_speciesIndex(int i, int k, int p);
    int DLL_IMPORT mix_nSpecies(int i);
    int DLL_IMPORT mix_setTemperature(int i, double t);
    double DLL_IMPORT mix_temperature(int i);
    int DLL_IMPORT mix_setPressure(int i, double p);
    double DLL_IMPORT mix_pressure(int i);
    double DLL_IMPORT mix_nAtoms(int i, int k, int m);
    double DLL_IMPORT mix_phaseMoles(int i, int n);
    int DLL_IMPORT mix_setPhaseMoles(int i, int n, double v);
    int DLL_IMPORT mix_setMoles(int i, int nlen, double* n);
        int DLL_IMPORT mix_setMolesByName(int i, char* n);
    double DLL_IMPORT mix_speciesMoles(int i, int k);
    double DLL_IMPORT mix_elementMoles(int i, int m);
    double DLL_IMPORT mix_equilibrate(int i, char* XY, 
        double err, int maxiter);
    int DLL_IMPORT mix_getChemPotentials(int i, int lenmu, double* mu);
}
#endif
