/**
 * @file ctonedim.h
 */
#ifndef CTC_ONEDIM_H
#define CTC_ONEDIM_H

#include "clib_defs.h"
#include "cantera/base/config.h"

extern "C" {
    CANTERA_CAPI int domain_clear();
    CANTERA_CAPI int domain_del(int i);
    CANTERA_CAPI int domain_type(int i);
    CANTERA_CAPI size_t domain_index(int i);
    CANTERA_CAPI size_t domain_nComponents(int i);
    CANTERA_CAPI size_t domain_nPoints(int i);
    CANTERA_CAPI int domain_componentName(int i, int n, int sz, char* nameout);
    CANTERA_CAPI size_t domain_componentIndex(int i, char* name);
    CANTERA_CAPI int domain_setBounds(int i, int n, double lower,
                                      double upper);
    CANTERA_CAPI double domain_lowerBound(int i, int n);
    CANTERA_CAPI double domain_upperBound(int i, int n);
    CANTERA_CAPI int domain_setTolerances(int i, int n, double rtol,
                                          double atol, int itime);
    CANTERA_CAPI double domain_rtol(int i, int n);
    CANTERA_CAPI double domain_atol(int i, int n);
    CANTERA_CAPI int domain_setupGrid(int i, size_t npts, double* grid);
    CANTERA_CAPI int domain_setID(int i, char* id);
    CANTERA_CAPI int domain_setDesc(int i, char* desc);
    CANTERA_CAPI double domain_grid(int i, int n);

    CANTERA_CAPI int bdry_setMdot(int i, double mdot);
    CANTERA_CAPI int bdry_setTemperature(int i, double t);
    CANTERA_CAPI int bdry_setMoleFractions(int i, char* x);
    CANTERA_CAPI double bdry_temperature(int i);
    CANTERA_CAPI double bdry_massFraction(int i, int k);
    CANTERA_CAPI double bdry_mdot(int i);

    CANTERA_CAPI int reactingsurf_setkineticsmgr(int i, int j);
    CANTERA_CAPI int reactingsurf_enableCoverageEqs(int i, int onoff);

    CANTERA_CAPI int inlet_new();
    CANTERA_CAPI int outlet_new();
    CANTERA_CAPI int outletres_new();
    CANTERA_CAPI int symm_new();
    CANTERA_CAPI int surf_new();
    CANTERA_CAPI int reactingsurf_new();

    CANTERA_CAPI int inlet_setSpreadRate(int i, double v);

    CANTERA_CAPI int stflow_new(int iph, int ikin, int itr, int itype=1);
    CANTERA_CAPI int stflow_setTransport(int i, int itr, int iSoret);
    CANTERA_CAPI int stflow_enableSoret(int i, int iSoret);
    CANTERA_CAPI int stflow_setPressure(int i, double p);
    CANTERA_CAPI double stflow_pressure(int i);
    CANTERA_CAPI int stflow_setFixedTempProfile(int i, size_t n, double* pos,
            size_t m, double* temp);
    CANTERA_CAPI  int stflow_solveSpeciesEqs(int i, int flag);
    CANTERA_CAPI int stflow_solveEnergyEqn(int i, int flag);

    CANTERA_CAPI int sim1D_clear();
    CANTERA_CAPI int sim1D_new(size_t nd, int* domains);
    CANTERA_CAPI int sim1D_del(int i);
    CANTERA_CAPI int sim1D_setValue(int i, int dom, int comp, int localPoint, double value);
    CANTERA_CAPI int sim1D_setProfile(int i, int dom, int comp,
                                      size_t np, double* pos, size_t nv, double* v);
    CANTERA_CAPI int sim1D_setFlatProfile(int i, int dom, int comp, double v);
    CANTERA_CAPI int sim1D_showSolution(int i, char* fname);
    CANTERA_CAPI int sim1D_setTimeStep(int i, double stepsize, size_t ns, integer* nsteps);
    CANTERA_CAPI int sim1D_getInitialSoln(int i);
    CANTERA_CAPI int sim1D_solve(int i, int loglevel, int refine_grid);
    CANTERA_CAPI int sim1D_refine(int i, int loglevel);
    CANTERA_CAPI int sim1D_setRefineCriteria(int i, int dom, double ratio,
            double slope, double curve, double prune);
    CANTERA_CAPI int sim1D_setGridMin(int i, int dom, double gridmin);
    CANTERA_CAPI int sim1D_save(int i, char* fname, char* id,
                                char* desc);
    CANTERA_CAPI int sim1D_restore(int i, char* fname, char* id);
    CANTERA_CAPI int sim1D_writeStats(int i, int printTime = 1);
    CANTERA_CAPI int sim1D_domainIndex(int i, char* name);
    CANTERA_CAPI double sim1D_value(int i, int idom, int icomp, int localPoint);
    CANTERA_CAPI double sim1D_workValue(int i, int idom,
                                        int icomp, int localPoint);
    CANTERA_CAPI int sim1D_eval(int i, double rdt, int count);
    CANTERA_CAPI int sim1D_setMaxJacAge(int i, int ss_age, int ts_age);
    CANTERA_CAPI int sim1D_timeStepFactor(int i, double tfactor);
    CANTERA_CAPI int sim1D_setTimeStepLimits(int i, double tsmin, double tsmax);
    CANTERA_CAPI int sim1D_setFixedTemperature(int i, double temp);
    CANTERA_CAPI int sim1D_evalSSJacobian(int i);
    CANTERA_CAPI double sim1D_jacobian(int i, int m, int n);
    CANTERA_CAPI size_t sim1D_size(int i);
}

#endif
