/**
 * @file ct.h
 */
#ifndef CTC_CT_H
#define CTC_CT_H

#include "clib_defs.h"
#include "cantera/base/config.h"

extern "C" {
    CANTERA_CAPI int ct_appdelete();
    CANTERA_CAPI size_t phase_nElements(int n);
    CANTERA_CAPI size_t phase_nSpecies(int n);
    CANTERA_CAPI double phase_temperature(int n);
    CANTERA_CAPI int phase_setTemperature(int n, double t);
    CANTERA_CAPI double phase_density(int n);
    CANTERA_CAPI int phase_setDensity(int n, double rho);
    CANTERA_CAPI double phase_molarDensity(int n);
    CANTERA_CAPI int phase_setMolarDensity(int n, double ndens);
    CANTERA_CAPI double phase_meanMolecularWeight(int n);
    CANTERA_CAPI double phase_moleFraction(int n, size_t k);
    CANTERA_CAPI double phase_massFraction(int n, size_t k);
    CANTERA_CAPI int phase_getMoleFractions(int n, size_t lenx, double* x);
    CANTERA_CAPI int phase_getMassFractions(int n, size_t leny, double* y);
    CANTERA_CAPI int phase_setMoleFractions(int n, size_t lenx,
                                            double* x, int norm);
    CANTERA_CAPI int phase_setMassFractions(int n, size_t leny,
                                            double* y, int norm);
    CANTERA_CAPI int phase_setMoleFractionsByName(int n, char* x);
    CANTERA_CAPI int phase_setMassFractionsByName(int n, char* y);
    CANTERA_CAPI int phase_getAtomicWeights(int n, size_t lenm, double* atw);
    CANTERA_CAPI int phase_getMolecularWeights(int n, size_t lenm, double* mw);
    CANTERA_CAPI int phase_getElementName(int n, size_t k, size_t lennm, char* nm);
    CANTERA_CAPI int phase_getSpeciesName(int n, size_t m, size_t lennm, char* nm);
    CANTERA_CAPI int phase_getName(int n, size_t lennm, char* nm);
    CANTERA_CAPI int phase_setName(int n, const char* nm);
    CANTERA_CAPI size_t phase_elementIndex(int n, char* nm);
    CANTERA_CAPI size_t phase_speciesIndex(int n, char* nm);
    CANTERA_CAPI int phase_report(int nth,
                                  int ibuf, char* buf, int show_thermo);
    CANTERA_CAPI int write_phase(int nth, int show_thermo);

    CANTERA_CAPI double phase_nAtoms(int n, size_t k, size_t m);

    CANTERA_CAPI int phase_addElement(int n, char* name, double weight);

    CANTERA_CAPI int newThermoFromXML(int mxml);
    CANTERA_CAPI size_t th_nSpecies(size_t n);
    CANTERA_CAPI int th_eosType(int n);
    CANTERA_CAPI double th_refPressure(int n);
    CANTERA_CAPI double th_minTemp(int n, int k=-1);
    CANTERA_CAPI double th_maxTemp(int n, int k=-1);
    CANTERA_CAPI double th_enthalpy_mole(int n);
    CANTERA_CAPI double th_intEnergy_mole(int n);
    CANTERA_CAPI double th_entropy_mole(int n);
    CANTERA_CAPI double th_gibbs_mole(int n);
    CANTERA_CAPI double th_cp_mole(int n);
    CANTERA_CAPI double th_cv_mole(int n);
    CANTERA_CAPI double th_pressure(int n);
    CANTERA_CAPI int th_setPressure(int n, double p);
    CANTERA_CAPI double th_enthalpy_mass(int n);
    CANTERA_CAPI double th_intEnergy_mass(int n);
    CANTERA_CAPI double th_entropy_mass(int n);
    CANTERA_CAPI double th_gibbs_mass(int n);
    CANTERA_CAPI double th_cp_mass(int n);
    CANTERA_CAPI double th_cv_mass(int n);
    CANTERA_CAPI double th_electricPotential(int n);
    CANTERA_CAPI int th_chemPotentials(int n, size_t lenm, double* murt);
    CANTERA_CAPI int th_elementPotentials(int n, size_t lenm, double* lambda);
    CANTERA_CAPI int th_getEnthalpies_RT(int n, size_t lenm, double* h_rt);
    CANTERA_CAPI int th_getEntropies_R(int n, size_t lenm, double* s_r);
    CANTERA_CAPI int th_getCp_R(int n, size_t lenm, double* cp_r);
    CANTERA_CAPI int th_setElectricPotential(int n, double v);
    CANTERA_CAPI int th_set_HP(int n, double* vals);
    CANTERA_CAPI int th_set_UV(int n, double* vals);
    CANTERA_CAPI int th_set_SV(int n, double* vals);
    CANTERA_CAPI int th_set_SP(int n, double* vals);
    CANTERA_CAPI int th_equil(int n, char* XY, int solver,
                              double rtol, int maxsteps, int maxiter, int loglevel);

    CANTERA_CAPI double th_critTemperature(int n);
    CANTERA_CAPI double th_critPressure(int n);
    CANTERA_CAPI double th_critDensity(int n);
    CANTERA_CAPI double th_vaporFraction(int n);
    CANTERA_CAPI double th_satTemperature(int n, double p);
    CANTERA_CAPI double th_satPressure(int n, double t);
    CANTERA_CAPI int th_setState_Psat(int n, double p, double x);
    CANTERA_CAPI int th_setState_Tsat(int n, double t, double x);

    CANTERA_CAPI size_t newKineticsFromXML(int mxml, int iphase,
                                           int neighbor1=-1, int neighbor2=-1, int neighbor3=-1,
                                           int neighbor4=-1);
    CANTERA_CAPI int installRxnArrays(int pxml, int ikin,
                                      char* default_phase);
    CANTERA_CAPI size_t kin_nSpecies(int n);
    CANTERA_CAPI size_t kin_nReactions(int n);
    CANTERA_CAPI size_t kin_nPhases(int n);
    CANTERA_CAPI size_t kin_phaseIndex(int n, char* ph);
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
    CANTERA_CAPI int kin_getActivationEnergies(int n, size_t len, double* E);
    CANTERA_CAPI int kin_getDelta(int n, int job, size_t len, double* delta);
    CANTERA_CAPI int kin_getCreationRates(int n, size_t len, double* cdot);
    CANTERA_CAPI int kin_getDestructionRates(int n, size_t len, double* ddot);
    CANTERA_CAPI int kin_getNetProductionRates(int n, size_t len, double* wdot);
    CANTERA_CAPI int kin_getSourceTerms(int n, size_t len, double* ydot);
    CANTERA_CAPI double kin_multiplier(int n, int i);
    CANTERA_CAPI int kin_getReactionString(int n, int i, int len, char* buf);
    CANTERA_CAPI int kin_setMultiplier(int n, int i, double v);

    CANTERA_CAPI int kin_isReversible(int n, int i);
    CANTERA_CAPI int kin_type(int n);
    CANTERA_CAPI size_t kin_start(int n, int p);
    CANTERA_CAPI size_t kin_speciesIndex(int n, const char* nm, const char* ph);
    CANTERA_CAPI int kin_advanceCoverages(int n, double tstep);
    CANTERA_CAPI size_t kin_phase(int n, size_t i);

    CANTERA_CAPI size_t newTransport(char* model,
                                     int th, int loglevel);
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

    CANTERA_CAPI int import_phase(int nth, int nxml, char* id);
    CANTERA_CAPI int import_kinetics(int nxml, char* id,
                                     int nphases, int* ith, int nkin);
    CANTERA_CAPI int getCanteraError(int buflen, char* buf);
    CANTERA_CAPI int showCanteraErrors();
    CANTERA_CAPI int write_HTML_log(const char* file);
    CANTERA_CAPI int setLogWriter(void* logger);
    CANTERA_CAPI int addCanteraDirectory(size_t buflen, char* buf);
    CANTERA_CAPI int clearStorage();
    CANTERA_CAPI int delThermo(int n);
    CANTERA_CAPI int delKinetics(int n);
    CANTERA_CAPI int delTransport(int n);
    CANTERA_CAPI int readlog(int n, char* buf);
    CANTERA_CAPI int buildSolutionFromXML(char* src, int ixml, char* id,
                                          int ith, int ikin);

    CANTERA_CAPI int ck_to_cti(char* in_file, char* db_file,
                               char* tr_file, char* id_tag, int debug, int validate);
    CANTERA_CAPI int writelogfile(char* logfile);
}

#endif
