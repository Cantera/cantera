/**
 * @file ct.h
 */
/*
 *      $Id: ct.h,v 1.24 2009/07/11 17:16:09 hkmoffa Exp $
 */

#ifndef CTC_CT_H
#define CTC_CT_H

#include "clib_defs.h"

#ifdef CANTERA_USE_INTERNAL
#include "config.h"
#else
#include "cantera/config.h"
#endif


extern "C" {
    EEXXTT int DLL_CPREFIX ct_appdelete();

    EEXXTT int DLL_CPREFIX phase_nElements(int n);
    EEXXTT int DLL_CPREFIX phase_nSpecies(int n);
    EEXXTT double DLL_CPREFIX phase_temperature(int n);
    EEXXTT int DLL_CPREFIX phase_setTemperature(int n, double t);
    EEXXTT double DLL_CPREFIX phase_density(int n);
    EEXXTT int DLL_CPREFIX phase_setDensity(int n, double rho);
    EEXXTT double DLL_CPREFIX phase_molarDensity(int n);
    EEXXTT int DLL_CPREFIX phase_setMolarDensity(int n, double ndens);
    EEXXTT double DLL_CPREFIX phase_meanMolecularWeight(int n);
    EEXXTT double DLL_CPREFIX phase_moleFraction(int n, int k);
    EEXXTT double DLL_CPREFIX phase_massFraction(int n, int k);
    EEXXTT int DLL_CPREFIX phase_getMoleFractions(int n, int lenx, double* x);
    EEXXTT int DLL_CPREFIX phase_getMassFractions(int n, int leny, double* y);
    EEXXTT int DLL_CPREFIX phase_setMoleFractions(int n, int lenx, 
        double* x, int norm);
    EEXXTT int DLL_CPREFIX phase_setMassFractions(int n, int leny, 
        double* y, int norm);
    EEXXTT int DLL_CPREFIX phase_setMoleFractionsByName(int n, char* x);
    EEXXTT int DLL_CPREFIX phase_setMassFractionsByName(int n, char* y);
    EEXXTT int DLL_CPREFIX phase_getAtomicWeights(int n, int lenm, double* atw);
    EEXXTT int DLL_CPREFIX phase_getMolecularWeights(int n, int lenm, double* mw);
    EEXXTT int DLL_CPREFIX phase_getElementName(int n, int k, int lennm, char* nm);
    EEXXTT int DLL_CPREFIX phase_getSpeciesName(int n, int m, int lennm, char* nm);
    EEXXTT int DLL_CPREFIX phase_getName(int n, int lennm, char* nm);
    EEXXTT int DLL_CPREFIX phase_setName(int n, const char* nm);
    EEXXTT int DLL_CPREFIX phase_elementIndex(int n, char* nm);
    EEXXTT int DLL_CPREFIX phase_speciesIndex(int n, char* nm);
    EEXXTT int DLL_CPREFIX phase_report(int nth, 
        int ibuf, char* buf, int show_thermo);
    EEXXTT int DLL_EXPORT write_phase(int nth, int show_thermo);

    EEXXTT double DLL_CPREFIX phase_nAtoms(int n, int k, int m);

    EEXXTT int DLL_CPREFIX phase_addElement(int n, char* name, double weight);
    EEXXTT int DLL_CPREFIX phase_addSpecies(int n, char* name, int phase,
        int ncomp, double* comp, int thermoType, int ncoeffs, 
        double* coeffs, double minTemp, double maxTemp, double refPressure,
        double charge, double weight);

  //int DLL_CPREFIX newThermo(char* model);
    EEXXTT int DLL_CPREFIX newThermoFromXML(int mxml);
    EEXXTT int DLL_CPREFIX th_thermoIndex(char* id);
    EEXXTT int DLL_CPREFIX th_phase(int n);
    EEXXTT int DLL_CPREFIX th_nSpecies(int n);
    EEXXTT int DLL_CPREFIX th_eosType(int n);
    EEXXTT double DLL_CPREFIX th_refPressure(int n);
    EEXXTT double DLL_CPREFIX th_minTemp(int n, int k=-1);
    EEXXTT double DLL_CPREFIX th_maxTemp(int n, int k=-1);
    EEXXTT double DLL_CPREFIX th_enthalpy_mole(int n);
    EEXXTT double DLL_CPREFIX th_intEnergy_mole(int n);
    EEXXTT double DLL_CPREFIX th_entropy_mole(int n);
    EEXXTT double DLL_CPREFIX th_gibbs_mole(int n);
    EEXXTT double DLL_CPREFIX th_cp_mole(int n);
    EEXXTT double DLL_CPREFIX th_cv_mole(int n);
    EEXXTT double DLL_CPREFIX th_pressure(int n);
    EEXXTT int DLL_CPREFIX th_setPressure(int n, double p);
    EEXXTT double DLL_CPREFIX th_enthalpy_mass(int n);
    EEXXTT double DLL_CPREFIX th_intEnergy_mass(int n);
    EEXXTT double DLL_CPREFIX th_entropy_mass(int n);
    EEXXTT double DLL_CPREFIX th_gibbs_mass(int n);
    EEXXTT double DLL_CPREFIX th_cp_mass(int n);
    EEXXTT double DLL_CPREFIX th_cv_mass(int n);
    EEXXTT double DLL_CPREFIX th_electricPotential(int n);
    EEXXTT int DLL_CPREFIX th_chemPotentials(int n, int lenm, double* murt);
    EEXXTT int DLL_CPREFIX th_elementPotentials(int n, int lenm, double* lambda);
    EEXXTT int DLL_CPREFIX th_getEnthalpies_RT(int n, int lenm, double* h_rt);
    EEXXTT int DLL_CPREFIX th_getEntropies_R(int n, int lenm, double* s_r);
    EEXXTT int DLL_CPREFIX th_getCp_R(int n, int lenm, double* cp_r);
    EEXXTT int DLL_CPREFIX th_setElectricPotential(int n, double v);
    EEXXTT int DLL_CPREFIX get_eos(char* fname, char* phase_id);

    EEXXTT int DLL_CPREFIX th_set_HP(int n, double* vals);
    EEXXTT int DLL_CPREFIX th_set_UV(int n, double* vals);
    EEXXTT int DLL_CPREFIX th_set_SV(int n, double* vals);
    EEXXTT int DLL_CPREFIX th_set_SP(int n, double* vals);
    EEXXTT int DLL_CPREFIX th_equil(int n, char* XY, int solver,
        double rtol, int maxsteps, int maxiter, int loglevel);
    
    EEXXTT double DLL_CPREFIX th_critTemperature(int n);
    EEXXTT double DLL_CPREFIX th_critPressure(int n);
    EEXXTT double DLL_CPREFIX th_critDensity(int n);
    EEXXTT double DLL_CPREFIX th_vaporFraction(int n);
    EEXXTT double DLL_CPREFIX th_satTemperature(int n, double p);
    EEXXTT double DLL_CPREFIX th_satPressure(int n, double t);
    EEXXTT int DLL_CPREFIX th_setState_Psat(int n, double p, double x);
    EEXXTT int DLL_CPREFIX th_setState_Tsat(int n, double t, double x);

    EEXXTT int DLL_CPREFIX newKineticsFromXML(int mxml, int iphase, 
        int neighbor1=-1, int neighbor2=-1, int neighbor3=-1,
        int neighbor4=-1);
    EEXXTT int DLL_CPREFIX installRxnArrays(int pxml, int ikin, 
        char* default_phase);
    EEXXTT int DLL_CPREFIX kin_nSpecies(int n);
    EEXXTT int DLL_CPREFIX kin_nReactions(int n);
    EEXXTT int DLL_CPREFIX kin_nPhases(int n);
    EEXXTT int DLL_CPREFIX kin_phaseIndex(int n, char* ph);
    EEXXTT int DLL_CPREFIX kin_reactionPhaseIndex(int n);
    EEXXTT double DLL_CPREFIX kin_reactantStoichCoeff(int n, int i, int k);
    EEXXTT double DLL_CPREFIX kin_productStoichCoeff(int n, int i, int k);
    EEXXTT int DLL_CPREFIX kin_reactionType(int n, int i);
    EEXXTT int DLL_CPREFIX kin_getFwdRatesOfProgress(int n, int len, double* fwdROP);
    EEXXTT int DLL_CPREFIX kin_getRevRatesOfProgress(int n, int len, double* revROP);
    EEXXTT int DLL_CPREFIX kin_getNetRatesOfProgress(int n, int len, double* netROP);
    EEXXTT int DLL_CPREFIX kin_getEquilibriumConstants(int n, int len, double* kc);

    EEXXTT int DLL_CPREFIX kin_getFwdRateConstants(int n, int len, double* kfwd);
    EEXXTT int DLL_CPREFIX kin_getRevRateConstants(int n, int doIrreversible, int len, double* krev);
    EEXXTT int DLL_CPREFIX kin_getActivationEnergies(int n, int len, double* E);
    EEXXTT int DLL_CPREFIX kin_getDelta(int n, int job, int len, double* delta);
    EEXXTT int DLL_CPREFIX kin_getCreationRates(int n, int len, double* cdot);
    EEXXTT int DLL_CPREFIX kin_getDestructionRates(int n, int len, double* ddot);
    EEXXTT int DLL_CPREFIX kin_getNetProductionRates(int n, int len, double* wdot);
    EEXXTT int DLL_CPREFIX kin_getSourceTerms(int n, int len, double* ydot);
    EEXXTT double DLL_CPREFIX kin_multiplier(int n, int i);
    EEXXTT int DLL_CPREFIX kin_getReactionString(int n, int i, int len, char* buf);
    EEXXTT int DLL_CPREFIX kin_setMultiplier(int n, int i, double v);

    EEXXTT int DLL_CPREFIX kin_isReversible(int n, int i);
    EEXXTT int DLL_CPREFIX kin_type(int n);
    EEXXTT int DLL_CPREFIX kin_start(int n, int p);
    EEXXTT int DLL_CPREFIX kin_speciesIndex(int n, const char* nm, const char* ph);
    EEXXTT int DLL_CPREFIX kin_advanceCoverages(int n, double tstep);
    EEXXTT int DLL_CPREFIX kin_phase(int n, int i);

    EEXXTT int DLL_CPREFIX newTransport(char* model,  
        int th, int loglevel);
    EEXXTT double DLL_CPREFIX trans_viscosity(int n);
    EEXXTT double DLL_CPREFIX trans_thermalConductivity(int n);
    EEXXTT int DLL_CPREFIX trans_getThermalDiffCoeffs(int n, int ldt, double* dt);
    EEXXTT int DLL_CPREFIX trans_getMixDiffCoeffs(int n, int ld, double* d);
    EEXXTT int DLL_CPREFIX trans_getBinDiffCoeffs(int n, int ld, double* d);
    EEXXTT int DLL_CPREFIX trans_getMultiDiffCoeffs(int n, int ld, double* d);
    EEXXTT int DLL_CPREFIX trans_setParameters(int n, int type, int k, double* d);
    EEXXTT int DLL_CPREFIX trans_getMolarFluxes(int n, const double* state1,
        const double* state2, double delta, double* fluxes);

    EEXXTT int DLL_CPREFIX import_phase(int nth, int nxml, char* id);
    EEXXTT int DLL_CPREFIX import_kinetics(int nxml, char* id, 
        int nphases, int* ith, int nkin);
    EEXXTT int DLL_CPREFIX import_from_file(int nth, int nkin, char* fname, char* db, 
        char* id, int validate, double threshold);
    EEXXTT int DLL_CPREFIX getCanteraError(int buflen, char* buf);
    EEXXTT int DLL_CPREFIX showCanteraErrors();
    EEXXTT int DLL_CPREFIX write_HTML_log(char* file);
    EEXXTT int DLL_CPREFIX setLogWriter(void* logger);
    EEXXTT int DLL_CPREFIX addCanteraDirectory(int buflen, char* buf);
    EEXXTT int DLL_CPREFIX clearStorage();
    EEXXTT int DLL_CPREFIX delPhase(int n);
    EEXXTT int DLL_CPREFIX delThermo(int n);
    EEXXTT int DLL_CPREFIX delKinetics(int n);
    EEXXTT int DLL_CPREFIX delTransport(int n);
    EEXXTT int DLL_CPREFIX readlog(int n, char* buf);
    EEXXTT int DLL_CPREFIX buildSolutionFromXML(char* src, int ixml, char* id, 
        int ith, int ikin);

    EEXXTT int DLL_CPREFIX ck_to_cti(char* in_file, char* db_file,
        char* tr_file, char* id_tag, int debug, int validate);
    EEXXTT int DLL_CPREFIX writelogfile(char* logfile);
}

#endif
