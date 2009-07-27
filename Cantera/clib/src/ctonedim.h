/**
 * @file ctonedim.h
 */
/*
 *      $Id: ctonedim.h,v 1.16 2009/07/11 17:16:09 hkmoffa Exp $
 */

#ifndef CTC_ONEDIM_H
#define CTC_ONEDIM_H

#include "clib_defs.h"

#ifdef CANTERA_USE_INTERNAL
#include "config.h"
#else
#include "cantera/config.h"
#endif

extern "C" {

    EEXXTT int DLL_CPREFIX domain_clear();
    EEXXTT int DLL_CPREFIX domain_del(int i);
    EEXXTT int DLL_CPREFIX domain_type(int i);
    EEXXTT int DLL_CPREFIX domain_index(int i);
    EEXXTT int DLL_CPREFIX domain_nComponents(int i);
    EEXXTT int DLL_CPREFIX domain_nPoints(int i);
    EEXXTT int DLL_CPREFIX domain_componentName(int i, int n, int sz, char* nameout);
    EEXXTT int DLL_CPREFIX domain_componentIndex(int i, char* name);
    EEXXTT int DLL_CPREFIX domain_setBounds(int i, int n, double lower, 
        double upper);
    EEXXTT double DLL_EXPORT domain_lowerBound(int i, int n);
    EEXXTT double DLL_EXPORT domain_upperBound(int i, int n);
    EEXXTT int DLL_CPREFIX domain_setTolerances(int i, int n, double rtol, 
        double atol, int itime);
    EEXXTT double DLL_CPREFIX domain_rtol(int i, int n);
    EEXXTT double DLL_CPREFIX domain_atol(int i, int n);
    EEXXTT int DLL_CPREFIX domain_setupGrid(int i, int npts, double* grid);
    EEXXTT int DLL_CPREFIX domain_setID(int i, char* id);
    EEXXTT int DLL_CPREFIX domain_setDesc(int i, char* desc);
    EEXXTT double DLL_CPREFIX domain_grid(int i, int n);

    EEXXTT int DLL_CPREFIX bdry_setMdot(int i, double mdot);
    EEXXTT int DLL_CPREFIX bdry_setTemperature(int i, double t);
    EEXXTT int DLL_CPREFIX bdry_setMoleFractions(int i, char* x);
    EEXXTT double DLL_CPREFIX bdry_temperature(int i);
    EEXXTT double DLL_CPREFIX bdry_massFraction(int i, int k);
    EEXXTT double DLL_CPREFIX bdry_mdot(int i);

    EEXXTT int DLL_CPREFIX reactingsurf_setkineticsmgr(int i, int j);
    EEXXTT int DLL_CPREFIX reactingsurf_enableCoverageEqs(int i, int onoff);

    EEXXTT int DLL_CPREFIX inlet_new();
    EEXXTT int DLL_CPREFIX outlet_new();
    EEXXTT int DLL_CPREFIX outletres_new();
    EEXXTT int DLL_CPREFIX symm_new();
    EEXXTT int DLL_CPREFIX surf_new();
    EEXXTT int DLL_CPREFIX reactingsurf_new();

    EEXXTT int DLL_CPREFIX inlet_setSpreadRate(int i, double v);

    EEXXTT int DLL_CPREFIX stflow_new(int iph, int ikin, int itr, int itype=1);
    EEXXTT int DLL_CPREFIX stflow_setTransport(int i, int itr, int iSoret);
    EEXXTT int DLL_CPREFIX stflow_enableSoret(int i, int iSoret);
    EEXXTT int DLL_CPREFIX stflow_setPressure(int i, double p);
    EEXXTT int DLL_CPREFIX stflow_setFixedTempProfile(int i, int n, double* pos, 
        int m, double* temp);
    EEXXTT  int DLL_CPREFIX stflow_solveSpeciesEqs(int i, int flag);
    EEXXTT int DLL_CPREFIX stflow_solveEnergyEqn(int i, int flag);

    EEXXTT int DLL_CPREFIX sim1D_clear();
    EEXXTT int DLL_CPREFIX sim1D_new(int nd, int* domains);
    EEXXTT int DLL_CPREFIX sim1D_del(int i);
    EEXXTT int DLL_CPREFIX sim1D_setValue(int i, int dom, int comp, int localPoint, double value);
    EEXXTT int DLL_CPREFIX sim1D_setProfile(int i, int dom, int comp, 
        int np, double* pos, int nv, double* v);
    EEXXTT int DLL_CPREFIX sim1D_setFlatProfile(int i, int dom, int comp, double v);
    EEXXTT int DLL_CPREFIX sim1D_showSolution(int i, char* fname);
    EEXXTT int DLL_CPREFIX sim1D_setTimeStep(int i, double stepsize, int ns, integer* nsteps);
    EEXXTT int DLL_CPREFIX sim1D_getInitialSoln(int i);
    EEXXTT int DLL_CPREFIX sim1D_solve(int i, int loglevel, int refine_grid);
    EEXXTT int DLL_CPREFIX sim1D_refine(int i, int loglevel);
    EEXXTT int DLL_CPREFIX sim1D_setRefineCriteria(int i, int dom, double ratio,
        double slope, double curve, double prune);
    EEXXTT int DLL_CPREFIX sim1D_save(int i, char* fname, char* id, 
        char* desc);
    EEXXTT int DLL_CPREFIX sim1D_restore(int i, char* fname, char* id);
    EEXXTT int DLL_CPREFIX sim1D_writeStats(int i);
    EEXXTT int DLL_CPREFIX sim1D_domainIndex(int i, char* name);
    EEXXTT double DLL_CPREFIX sim1D_value(int i, int idom, int icomp, int localPoint);
    EEXXTT double DLL_CPREFIX sim1D_workValue(int i, int idom, 
        int icomp, int localPoint);
    EEXXTT int DLL_CPREFIX sim1D_eval(int i, double rdt, int count);
    EEXXTT int DLL_CPREFIX sim1D_setMaxJacAge(int i, int ss_age, int ts_age);
    EEXXTT int DLL_CPREFIX sim1D_timeStepFactor(int i, double tfactor);
    EEXXTT int DLL_CPREFIX sim1D_setTimeStepLimits(int i, double tsmin, double tsmax);
    EEXXTT int DLL_CPREFIX sim1D_setFixedTemperature(int i, double temp);
    EEXXTT int DLL_CPREFIX sim1D_evalSSJacobian(int i);
    EEXXTT double DLL_CPREFIX sim1D_jacobian(int i, int m, int n);
    EEXXTT int DLL_CPREFIX sim1D_size(int i);
}


#endif
