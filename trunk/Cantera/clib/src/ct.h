#ifndef CTC_CT_H
#define CTC_CT_H

#include "clib_defs.h"
#include "../../src/config.h"

extern "C" {

    int DLL_IMPORT phase_nElements(int n);
    int DLL_IMPORT phase_nSpecies(int n);
    double DLL_IMPORT phase_temperature(int n);
    int DLL_IMPORT phase_setTemperature(int n, double t);
    double DLL_IMPORT phase_density(int n);
    int DLL_IMPORT phase_setDensity(int n, double rho);
    double DLL_IMPORT phase_molarDensity(int n);
    int DLL_IMPORT phase_setMolarDensity(int n, double ndens);
    double DLL_IMPORT phase_meanMolecularWeight(int n);
    double DLL_IMPORT phase_moleFraction(int n, int k);
    double DLL_IMPORT phase_massFraction(int n, int k);
    int DLL_IMPORT phase_getMoleFractions(int n, int lenx, double* x);
    int DLL_IMPORT phase_getMassFractions(int n, int leny, double* y);
    int DLL_IMPORT phase_setMoleFractions(int n, int lenx, 
        double* x, int norm);
    int DLL_IMPORT phase_setMassFractions(int n, int leny, 
        double* y, int norm);
    int DLL_IMPORT phase_setMoleFractionsByName(int n, char* x);
    int DLL_IMPORT phase_setMassFractionsByName(int n, char* y);
    int DLL_IMPORT phase_getAtomicWeights(int n, int lenm, double* atw);
    int DLL_IMPORT phase_getMolecularWeights(int n, int lenm, double* mw);
    int DLL_IMPORT phase_getElementName(int n, int k, int lennm, char* nm);
    int DLL_IMPORT phase_getSpeciesName(int n, int m, int lennm, char* nm);
    int DLL_IMPORT phase_getName(int n, int lennm, char* nm);
    int DLL_IMPORT phase_setName(int n, const char* nm);
    int DLL_IMPORT phase_elementIndex(int n, char* nm);
    int DLL_IMPORT phase_speciesIndex(int n, char* nm);
    int DLL_IMPORT phase_report(int nth, 
        int ibuf, char* buf, int show_thermo);
    int DLL_EXPORT write_phase(int nth, int show_thermo);

    double DLL_IMPORT phase_nAtoms(int n, int k, int m);

    int DLL_IMPORT phase_addElement(int n, char* name, double weight);
    int DLL_IMPORT phase_addSpecies(int n, char* name, int phase,
        int ncomp, double* comp, int thermoType, int ncoeffs, 
        double* coeffs, double minTemp, double maxTemp, double refPressure,
        double charge, double weight);

  //int DLL_IMPORT newThermo(char* model);
    int DLL_IMPORT newThermoFromXML(int mxml);
    int DLL_IMPORT th_thermoIndex(char* id);
    int DLL_IMPORT th_phase(int n);
    int DLL_IMPORT th_nSpecies(int n);
    int DLL_IMPORT th_eosType(int n);
    double DLL_IMPORT th_refPressure(int n);
    double DLL_IMPORT th_minTemp(int n, int k=-1);
    double DLL_IMPORT th_maxTemp(int n, int k=-1);
    double DLL_IMPORT th_enthalpy_mole(int n);
    double DLL_IMPORT th_intEnergy_mole(int n);
    double DLL_IMPORT th_entropy_mole(int n);
    double DLL_IMPORT th_gibbs_mole(int n);
    double DLL_IMPORT th_cp_mole(int n);
    double DLL_IMPORT th_cv_mole(int n);
    double DLL_IMPORT th_pressure(int n);
    int DLL_IMPORT th_setPressure(int n, double p);
    double DLL_IMPORT th_enthalpy_mass(int n);
    double DLL_IMPORT th_intEnergy_mass(int n);
    double DLL_IMPORT th_entropy_mass(int n);
    double DLL_IMPORT th_gibbs_mass(int n);
    double DLL_IMPORT th_cp_mass(int n);
    double DLL_IMPORT th_cv_mass(int n);
    double DLL_IMPORT th_electricPotential(int n);
    int DLL_IMPORT th_chemPotentials(int n, int lenm, double* murt);
    int DLL_IMPORT th_elementPotentials(int n, int lenm, double* lambda);
    int DLL_IMPORT th_getEnthalpies_RT(int n, int lenm, double* h_rt);
    int DLL_IMPORT th_getEntropies_R(int n, int lenm, double* s_r);
    int DLL_IMPORT th_getCp_R(int n, int lenm, double* cp_r);
    int DLL_IMPORT th_setElectricPotential(int n, double v);
    int DLL_IMPORT get_eos(char* fname, char* phase_id);

    int DLL_IMPORT th_set_HP(int n, double* vals);
    int DLL_IMPORT th_set_UV(int n, double* vals);
    int DLL_IMPORT th_set_SV(int n, double* vals);
    int DLL_IMPORT th_set_SP(int n, double* vals);
    int DLL_IMPORT th_equil(int n, char* XY, int solver,
        double rtol, int maxsteps, int maxiter, int loglevel);
    
    double DLL_IMPORT th_critTemperature(int n);
    double DLL_IMPORT th_critPressure(int n);
    double DLL_IMPORT th_critDensity(int n);
    double DLL_IMPORT th_vaporFraction(int n);
    double DLL_IMPORT th_satTemperature(int n, double p);
    double DLL_IMPORT th_satPressure(int n, double t);
    int DLL_IMPORT th_setState_Psat(int n, double p, double x);
    int DLL_IMPORT th_setState_Tsat(int n, double t, double x);

    int DLL_IMPORT newKineticsFromXML(int mxml, int iphase, 
        int neighbor1=-1, int neighbor2=-1, int neighbor3=-1,
        int neighbor4=-1);
    int DLL_IMPORT installRxnArrays(int pxml, int ikin, 
        char* default_phase);
    int DLL_IMPORT kin_nSpecies(int n);
    int DLL_IMPORT kin_nReactions(int n);
    double DLL_IMPORT kin_reactantStoichCoeff(int n, int i, int k);
    double DLL_IMPORT kin_productStoichCoeff(int n, int i, int k);
    int DLL_IMPORT kin_reactionType(int n, int i);
    int DLL_IMPORT kin_getFwdRatesOfProgress(int n, int len, double* fwdROP);
    int DLL_IMPORT kin_getRevRatesOfProgress(int n, int len, double* revROP);
    int DLL_IMPORT kin_getNetRatesOfProgress(int n, int len, double* netROP);
    int DLL_IMPORT kin_getEquilibriumConstants(int n, int len, double* kc);

    int DLL_IMPORT kin_getFwdRateConstants(int n, int len, double* kfwd);
    int DLL_IMPORT kin_getRevRateConstants(int n, int doIrreversible, int len, double* krev);
    int DLL_IMPORT kin_getActivationEnergies(int n, int len, double* E);
    int DLL_IMPORT kin_getDelta(int n, int job, int len, double* delta);
    int DLL_IMPORT kin_getCreationRates(int n, int len, double* cdot);
    int DLL_IMPORT kin_getDestructionRates(int n, int len, double* ddot);
    int DLL_IMPORT kin_getNetProductionRates(int n, int len, double* wdot);
    int DLL_IMPORT kin_getSourceTerms(int n, int len, double* ydot);
    double DLL_IMPORT kin_multiplier(int n, int i);
    int DLL_IMPORT kin_getReactionString(int n, int i, int len, char* buf);
    int DLL_IMPORT kin_setMultiplier(int n, int i, double v);

    int DLL_IMPORT kin_isReversible(int n, int i);
    int DLL_IMPORT kin_type(int n);
    int DLL_IMPORT kin_start(int n, int p);
    int DLL_IMPORT kin_speciesIndex(int n, const char* nm, const char* ph);
    int DLL_IMPORT kin_advanceCoverages(int n, double tstep);
    int DLL_IMPORT kin_phase(int n, int i);

    int DLL_IMPORT newTransport(char* model,  
        int th, int loglevel);
    double DLL_IMPORT trans_viscosity(int n);
    double DLL_IMPORT trans_thermalConductivity(int n);
    int DLL_IMPORT trans_getThermalDiffCoeffs(int n, int ldt, double* dt);
    int DLL_IMPORT trans_getMixDiffCoeffs(int n, int ld, double* d);
    int DLL_IMPORT trans_getBinDiffCoeffs(int n, int ld, double* d);
    int DLL_IMPORT trans_getMultiDiffCoeffs(int n, int ld, double* d);
    int DLL_IMPORT trans_setParameters(int n, int type, int k, double* d);
    int DLL_IMPORT trans_getMolarFluxes(int n, const double* state1,
        const double* state2, double delta, double* fluxes);

    int DLL_IMPORT import_phase(int nth, int nxml, char* id);
    int DLL_IMPORT import_kinetics(int nxml, char* id, 
        int nphases, int* ith, int nkin);
    int DLL_IMPORT import_from_file(int nth, int nkin, char* fname, char* db, 
        char* id, int validate, double threshold);
    int DLL_IMPORT getCanteraError(int buflen, char* buf);
    int DLL_IMPORT showCanteraErrors();
    int DLL_IMPORT write_HTML_log(char* file);
    int DLL_IMPORT setLogWriter(void* logger);
    int DLL_IMPORT addCanteraDirectory(int buflen, char* buf);
    int DLL_IMPORT clearStorage();
    int DLL_IMPORT delPhase(int n);
    int DLL_IMPORT delThermo(int n);
    int DLL_IMPORT delKinetics(int n);
    int DLL_IMPORT delTransport(int n);
    int DLL_IMPORT readlog(int n, char* buf);
    int DLL_IMPORT buildSolutionFromXML(char* src, int ixml, char* id, 
        int ith, int ikin);

    int DLL_IMPORT ck_to_cti(char* in_file, char* db_file,
        char* tr_file, char* id_tag, int debug, int validate);
    int DLL_IMPORT writelogfile(char* logfile);
}

#endif
