#ifndef CTC_ONEDIM_H
#define CTC_ONEDIM_H

#include "clib_defs.h"
#include "cantera/config.h"

extern "C" {

    int DLL_IMPORT domain_clear();
    int DLL_IMPORT domain_del(int i);
    int DLL_IMPORT domain_type(int i);
    int DLL_IMPORT domain_index(int i);
    int DLL_IMPORT domain_nComponents(int i);
    int DLL_IMPORT domain_nPoints(int i);
    int DLL_IMPORT domain_componentName(int i, int n, int sz, char* nameout);
    int DLL_IMPORT domain_componentIndex(int i, char* name);
    int DLL_IMPORT domain_setBounds(int i, int n, double lower, 
        double upper);
    double DLL_EXPORT domain_lowerBound(int i, int n);
    double DLL_EXPORT domain_upperBound(int i, int n);
    int DLL_IMPORT domain_setTolerances(int i, int n, double rtol, 
        double atol, int itime);
    double DLL_IMPORT domain_rtol(int i, int n);
    double DLL_IMPORT domain_atol(int i, int n);
    int DLL_IMPORT domain_setupGrid(int i, int npts, double* grid);
    int DLL_IMPORT domain_setID(int i, char* id);
    int DLL_IMPORT domain_setDesc(int i, char* desc);
    double DLL_IMPORT domain_grid(int i, int n);

    int DLL_IMPORT bdry_setMdot(int i, double mdot);
    int DLL_IMPORT bdry_setTemperature(int i, double t);
    int DLL_IMPORT bdry_setMoleFractions(int i, char* x);
    double DLL_IMPORT bdry_temperature(int i);
    double DLL_IMPORT bdry_massFraction(int i, int k);
    double DLL_IMPORT bdry_mdot(int i);

    int DLL_IMPORT reactingsurf_setkineticsmgr(int i, int j);
    int DLL_IMPORT reactingsurf_enableCoverageEqs(int i, int onoff);

    int DLL_IMPORT inlet_new();
    int DLL_IMPORT outlet_new();
    int DLL_IMPORT outletres_new();
    int DLL_IMPORT symm_new();
    int DLL_IMPORT surf_new();
    int DLL_IMPORT reactingsurf_new();

    int DLL_IMPORT inlet_setSpreadRate(int i, double v);

    int DLL_IMPORT stflow_new(int iph, int ikin, int itr, int itype=1);
    int DLL_IMPORT stflow_setTransport(int i, int itr, int iSoret);
    int DLL_IMPORT stflow_enableSoret(int i, int iSoret);
    int DLL_IMPORT stflow_setPressure(int i, double p);
    int DLL_IMPORT stflow_setFixedTempProfile(int i, int n, double* pos, 
        int m, double* temp);
    int DLL_IMPORT stflow_solveSpeciesEqs(int i, int flag);
    int DLL_IMPORT stflow_solveEnergyEqn(int i, int flag);

    int DLL_IMPORT sim1D_clear();
    int DLL_IMPORT sim1D_new(int nd, int* domains);
    int DLL_IMPORT sim1D_del(int i);
    int DLL_IMPORT sim1D_setValue(int i, int dom, int comp, int localPoint, double value);
    int DLL_IMPORT sim1D_setProfile(int i, int dom, int comp, 
        int np, double* pos, int nv, double* v);
    int DLL_IMPORT sim1D_setFlatProfile(int i, int dom, int comp, double v);
    int DLL_IMPORT sim1D_showSolution(int i, char* fname);
    int DLL_IMPORT sim1D_setTimeStep(int i, double stepsize, int ns, integer* nsteps);
    int DLL_IMPORT sim1D_getInitialSoln(int i);
    int DLL_IMPORT sim1D_solve(int i, int loglevel, int refine_grid);
    int DLL_IMPORT sim1D_refine(int i, int loglevel);
    int DLL_IMPORT sim1D_setRefineCriteria(int i, int dom, double ratio,
        double slope, double curve, double prune);
    int DLL_IMPORT sim1D_save(int i, char* fname, char* id, 
        char* desc);
    int DLL_IMPORT sim1D_restore(int i, char* fname, char* id);
    int DLL_IMPORT sim1D_writeStats(int i);
    int DLL_IMPORT sim1D_domainIndex(int i, char* name);
    double DLL_IMPORT sim1D_value(int i, int idom, int icomp, int localPoint);
    double DLL_IMPORT sim1D_workValue(int i, int idom, 
        int icomp, int localPoint);
    int DLL_IMPORT sim1D_eval(int i, double rdt, int count);
    int DLL_IMPORT sim1D_setMaxJacAge(int i, int ss_age, int ts_age);
    int DLL_IMPORT sim1D_timeStepFactor(int i, double tfactor);
    int DLL_IMPORT sim1D_setTimeStepLimits(int i, double tsmin, double tsmax);
    int DLL_IMPORT sim1D_setFixedTemperature(int i, double temp);
    int DLL_IMPORT sim1D_evalSSJacobian(int i);
    double DLL_IMPORT sim1D_jacobian(int i, int m, int n);
    int DLL_IMPORT sim1D_size(int i);
}


#endif
