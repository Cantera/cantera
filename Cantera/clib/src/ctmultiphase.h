/**
 * @file ctmultiphase.h
 */
/*
 *      $Id: ctmultiphase.h,v 1.7 2009/07/11 17:16:09 hkmoffa Exp $
 */

#ifndef CTC_MULTIPHASE_H
#define CTC_MULTIPHASE_H

#include "clib_defs.h"

extern "C" {

  EEXXTT int DLL_CPREFIX mix_new();
  EEXXTT int DLL_CPREFIX mix_del(int i);
  EEXXTT int DLL_CPREFIX mix_copy(int i);
  EEXXTT int DLL_CPREFIX mix_assign(int i, int j);
  EEXXTT int DLL_CPREFIX mix_addPhase(int i, int j, double moles);
  EEXXTT int DLL_CPREFIX mix_init(int i);
  EEXXTT int DLL_CPREFIX mix_nElements(int i);
  EEXXTT int DLL_CPREFIX mix_elementIndex(int i, char* name);
  EEXXTT int DLL_CPREFIX mix_speciesIndex(int i, int k, int p);
  EEXXTT int DLL_CPREFIX mix_nSpecies(int i);
  EEXXTT int DLL_CPREFIX mix_setTemperature(int i, double t);
  EEXXTT double DLL_CPREFIX mix_temperature(int i);
  EEXXTT double DLL_CPREFIX mix_minTemp(int i);
  EEXXTT double DLL_CPREFIX mix_maxTemp(int i);
  EEXXTT double DLL_CPREFIX mix_charge(int i);
  EEXXTT double DLL_CPREFIX mix_phaseCharge(int i, int p);
  EEXXTT int DLL_CPREFIX mix_setPressure(int i, double p);
  EEXXTT double DLL_CPREFIX mix_pressure(int i);
  EEXXTT double DLL_CPREFIX mix_nAtoms(int i, int k, int m);
  EEXXTT double DLL_CPREFIX mix_nPhases(int i);
  EEXXTT double DLL_CPREFIX mix_phaseMoles(int i, int n);
  EEXXTT int DLL_CPREFIX mix_setPhaseMoles(int i, int n, double v);
  EEXXTT int DLL_CPREFIX mix_setMoles(int i, int nlen, double* n);
  EEXXTT int DLL_CPREFIX mix_setMolesByName(int i, char* n);
  EEXXTT double DLL_CPREFIX mix_speciesMoles(int i, int k);
  EEXXTT double DLL_CPREFIX mix_elementMoles(int i, int m);
  EEXXTT double DLL_CPREFIX mix_equilibrate(int i, char* XY, 
				    double err, int maxsteps, int maxiter, int loglevel);
  EEXXTT double DLL_EXPORT mix_vcs_equilibrate(int i, char* XY, int estimateEquil,
					int printLvl, int solver,
					double rtol, int maxsteps,
					int maxiter, int loglevel);
  EEXXTT int DLL_CPREFIX mix_getChemPotentials(int i, int lenmu, double* mu);
  EEXXTT int DLL_CPREFIX mix_getValidChemPotentials(int i, double bad_mu, 
					    int standard, int lenmu, double* mu);

  EEXXTT double DLL_CPREFIX mix_enthalpy(int i);
  EEXXTT double DLL_CPREFIX mix_entropy(int i);
  EEXXTT double DLL_CPREFIX mix_gibbs(int i);
  EEXXTT double DLL_CPREFIX mix_cp(int i);
  EEXXTT double DLL_CPREFIX mix_volume(int i);

  EEXXTT int DLL_CPREFIX mix_speciesPhaseIndex(int i, int k);
  EEXXTT double DLL_CPREFIX mix_moleFraction(int i, int k);

}
#endif
