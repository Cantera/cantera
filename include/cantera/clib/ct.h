/**
 * @file ct.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CTC_CT_H
#define CTC_CT_H

#include "clib_defs.h"

#ifdef __cplusplus
extern "C" {
#endif

    CANTERA_CAPI int ct_appdelete();

    CANTERA_CAPI int thermo_newFromFile(const char* filename, const char* phasename);
    CANTERA_CAPI int thermo_newFromXML(int mxml);
    CANTERA_CAPI int thermo_del(int n);
    CANTERA_CAPI size_t thermo_nElements(int n);
    CANTERA_CAPI size_t thermo_nSpecies(int n);
    CANTERA_CAPI double thermo_temperature(int n);
    CANTERA_CAPI int thermo_setTemperature(int n, double t);
    CANTERA_CAPI double thermo_density(int n);
    CANTERA_CAPI int thermo_setDensity(int n, double rho);
    CANTERA_CAPI double thermo_molarDensity(int n);
    CANTERA_CAPI int thermo_setMolarDensity(int n, double ndens);
    CANTERA_CAPI double thermo_meanMolecularWeight(int n);
    CANTERA_CAPI double thermo_moleFraction(int n, size_t k);
    CANTERA_CAPI double thermo_massFraction(int n, size_t k);
    CANTERA_CAPI int thermo_getMoleFractions(int n, size_t lenx, double* x);
    CANTERA_CAPI int thermo_getMassFractions(int n, size_t leny, double* y);
    CANTERA_CAPI int thermo_setMoleFractions(int n, size_t lenx,
                                            double* x, int norm);
    CANTERA_CAPI int thermo_setMassFractions(int n, size_t leny,
                                            double* y, int norm);
    CANTERA_CAPI int thermo_setMoleFractionsByName(int n, const char* x);
    CANTERA_CAPI int thermo_setMassFractionsByName(int n, const char* y);
    CANTERA_CAPI int thermo_getAtomicWeights(int n, size_t lenm, double* atw);
    CANTERA_CAPI int thermo_getMolecularWeights(int n, size_t lenm, double* mw);
    CANTERA_CAPI int thermo_getCharges(int n, size_t lenm, double* sc);
    CANTERA_CAPI int thermo_getElementName(int n, size_t k, size_t lennm, char* nm);
    CANTERA_CAPI int thermo_getSpeciesName(int n, size_t m, size_t lennm, char* nm);
    CANTERA_CAPI int thermo_getName(int n, size_t lennm, char* nm);
    CANTERA_CAPI int thermo_setName(int n, const char* nm);
    CANTERA_CAPI size_t thermo_elementIndex(int n, const char* nm);
    CANTERA_CAPI size_t thermo_speciesIndex(int n, const char* nm);
    CANTERA_CAPI int thermo_report(int nth,
                                  int ibuf, char* buf, int show_thermo);
    CANTERA_CAPI int thermo_print(int nth, int show_thermo, double threshold);
    CANTERA_CAPI double thermo_nAtoms(int n, size_t k, size_t m);
    CANTERA_CAPI int thermo_addElement(int n, const char* name, double weight);
    CANTERA_CAPI int thermo_getEosType(int n, size_t leneos, char* eos);
    CANTERA_CAPI double thermo_refPressure(int n);
    CANTERA_CAPI double thermo_minTemp(int n, int k);
    CANTERA_CAPI double thermo_maxTemp(int n, int k);
    CANTERA_CAPI double thermo_enthalpy_mole(int n);
    CANTERA_CAPI double thermo_intEnergy_mole(int n);
    CANTERA_CAPI double thermo_entropy_mole(int n);
    CANTERA_CAPI double thermo_gibbs_mole(int n);
    CANTERA_CAPI double thermo_cp_mole(int n);
    CANTERA_CAPI double thermo_cv_mole(int n);
    CANTERA_CAPI double thermo_pressure(int n);
    CANTERA_CAPI int thermo_setPressure(int n, double p);
    CANTERA_CAPI double thermo_enthalpy_mass(int n);
    CANTERA_CAPI double thermo_intEnergy_mass(int n);
    CANTERA_CAPI double thermo_entropy_mass(int n);
    CANTERA_CAPI double thermo_gibbs_mass(int n);
    CANTERA_CAPI double thermo_cp_mass(int n);
    CANTERA_CAPI double thermo_cv_mass(int n);
    CANTERA_CAPI double thermo_electricPotential(int n);
    CANTERA_CAPI double thermo_thermalExpansionCoeff(int n);
    CANTERA_CAPI double thermo_isothermalCompressibility(int n);
    CANTERA_CAPI int thermo_chemPotentials(int n, size_t lenm, double* murt);
    CANTERA_CAPI int thermo_getEnthalpies_RT(int n, size_t lenm, double* h_rt);
    CANTERA_CAPI int thermo_getEntropies_R(int n, size_t lenm, double* s_r);
    CANTERA_CAPI int thermo_getCp_R(int n, size_t lenm, double* cp_r);
    CANTERA_CAPI int thermo_setElectricPotential(int n, double v);
    CANTERA_CAPI int thermo_set_RP(int n, double* vals);
    CANTERA_CAPI int thermo_set_HP(int n, double* vals);
    CANTERA_CAPI int thermo_set_UV(int n, double* vals);
    CANTERA_CAPI int thermo_set_SV(int n, double* vals);
    CANTERA_CAPI int thermo_set_SP(int n, double* vals);
    CANTERA_CAPI int thermo_set_ST(int n, double* vals);
    CANTERA_CAPI int thermo_set_TV(int n, double* vals);
    CANTERA_CAPI int thermo_set_PV(int n, double* vals);
    CANTERA_CAPI int thermo_set_UP(int n, double* vals);
    CANTERA_CAPI int thermo_set_VH(int n, double* vals);
    CANTERA_CAPI int thermo_set_TH(int n, double* vals);
    CANTERA_CAPI int thermo_set_SH(int n, double* vals);
    CANTERA_CAPI int thermo_equilibrate(int n, const char* XY, int solver,
                                        double rtol, int maxsteps, int maxiter,
                                        int loglevel);

    CANTERA_CAPI double thermo_critTemperature(int n);
    CANTERA_CAPI double thermo_critPressure(int n);
    CANTERA_CAPI double thermo_critDensity(int n);
    CANTERA_CAPI double thermo_vaporFraction(int n);
    CANTERA_CAPI double thermo_satTemperature(int n, double p);
    CANTERA_CAPI double thermo_satPressure(int n, double t);
    CANTERA_CAPI int thermo_setState_Psat(int n, double p, double x);
    CANTERA_CAPI int thermo_setState_Tsat(int n, double t, double x);

    CANTERA_CAPI int kin_newFromFile(const char* filename, const char* phasename,
                                     int reactingPhase, int neighbor1, int neighbor2,
                                     int neighbor3, int neighbor4);
    CANTERA_CAPI int kin_newFromXML(int mxml, int iphase,
                                    int neighbor1, int neighbor2, int neighbor3,
                                    int neighbor4);
    CANTERA_CAPI int kin_del(int n);
    CANTERA_CAPI size_t kin_nSpecies(int n);
    CANTERA_CAPI size_t kin_nReactions(int n);
    CANTERA_CAPI size_t kin_nPhases(int n);
    CANTERA_CAPI size_t kin_phaseIndex(int n, const char* ph);
    CANTERA_CAPI size_t kin_reactionPhaseIndex(int n);
    CANTERA_CAPI double kin_reactantStoichCoeff(int n, int i, int k);
    CANTERA_CAPI double kin_productStoichCoeff(int n, int i, int k);
    CANTERA_CAPI int kin_reactionType(int n, int i);
    CANTERA_CAPI int kin_getFwdRatesOfProgress(int n, size_t len, double* fwdROP);
    CANTERA_CAPI int kin_getRevRatesOfProgress(int n, size_t len, double* revROP);
    CANTERA_CAPI int kin_getNetRatesOfProgress(int n, size_t len, double* netROP);
    CANTERA_CAPI int kin_getEquilibriumConstants(int n, size_t len, double* kc);

    CANTERA_CAPI int kin_getFwdRateConstants(int n, size_t len, double* kfwd);
    CANTERA_CAPI int kin_getRevRateConstants(int n, int doIrreversible, size_t len, double* krev);
    CANTERA_CAPI int kin_getDelta(int n, int job, size_t len, double* delta);
    CANTERA_CAPI int kin_getCreationRates(int n, size_t len, double* cdot);
    CANTERA_CAPI int kin_getDestructionRates(int n, size_t len, double* ddot);
    CANTERA_CAPI int kin_getNetProductionRates(int n, size_t len, double* wdot);
    CANTERA_CAPI int kin_getSourceTerms(int n, size_t len, double* ydot);
    CANTERA_CAPI double kin_multiplier(int n, int i);
    CANTERA_CAPI int kin_getReactionString(int n, int i, int len, char* buf);
    CANTERA_CAPI int kin_setMultiplier(int n, int i, double v);

    CANTERA_CAPI int kin_isReversible(int n, int i);
    CANTERA_CAPI int kin_getType(int n, size_t len, char* name);
    CANTERA_CAPI size_t kin_start(int n, int p);
    CANTERA_CAPI size_t kin_speciesIndex(int n, const char* nm, const char* ph);
    CANTERA_CAPI int kin_advanceCoverages(int n, double tstep);
    CANTERA_CAPI size_t kin_phase(int n, size_t i);

    CANTERA_CAPI int trans_newDefault(int th, int loglevel);
    CANTERA_CAPI int trans_new(const char* model, int th, int loglevel);
    CANTERA_CAPI int trans_del(int n);
    CANTERA_CAPI double trans_viscosity(int n);
    CANTERA_CAPI double trans_electricalConductivity(int n);
    CANTERA_CAPI double trans_thermalConductivity(int n);
    CANTERA_CAPI int trans_getThermalDiffCoeffs(int n, int ldt, double* dt);
    CANTERA_CAPI int trans_getMixDiffCoeffs(int n, int ld, double* d);
    CANTERA_CAPI int trans_getBinDiffCoeffs(int n, int ld, double* d);
    CANTERA_CAPI int trans_getMultiDiffCoeffs(int n, int ld, double* d);
    CANTERA_CAPI int trans_setParameters(int n, int type, int k, double* d);
    CANTERA_CAPI int trans_getMolarFluxes(int n, const double* state1,
                                          const double* state2, double delta, double* fluxes);
    CANTERA_CAPI int trans_getMassFluxes(int n, const double* state1,
                                         const double* state2, double delta, double* fluxes);

    CANTERA_CAPI int ct_getCanteraError(int buflen, char* buf);
    CANTERA_CAPI int ct_setLogWriter(void* logger);
    CANTERA_CAPI int ct_addCanteraDirectory(size_t buflen, const char* buf);
    CANTERA_CAPI int ct_getDataDirectories(int buflen, char* buf, const char* sep);
    CANTERA_CAPI int ct_getCanteraVersion(int buflen, char* buf);
    CANTERA_CAPI int ct_getGitCommit(int buflen, char* buf);
    CANTERA_CAPI int ct_suppress_thermo_warnings(int suppress);
    CANTERA_CAPI int ct_clearStorage();

    CANTERA_CAPI int ct_ck2cti(const char* in_file, const char* db_file,
                               const char* tr_file, const char* id_tag,
                               int debug, int validate);

#ifdef __cplusplus
}
#endif

#endif
