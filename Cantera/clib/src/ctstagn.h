/**
 * @file ctstagn.h
 */
/*
 *      $Id: ctstagn.h,v 1.2 2009/07/11 17:16:09 hkmoffa Exp $
 */

#ifndef CTC_STAGN_H
#define CTC_STAGN_H

// Cantera includes
//#include "stagn.h"

//#include "Cabinet.h"
//#include "Storage.h"
#include "clib_defs.h"

//inline StFlow* _flow(int i) {
//    return Cabinet<StFlow>::cabinet()->item(i);
//}

extern "C" {

    int DLL_IMPORT flow_new(int type, int iph, int np);
    int DLL_IMPORT flow_del(int i);
    int DLL_IMPORT flow_copy(int i);
    int DLL_IMPORT flow_assign(int i, int j);
    int DLL_IMPORT flow_setupgrid(int i, int npts, double* grid);
    int DLL_EXPORT flow_setthermo(int i, int k);
    int DLL_IMPORT flow_setkinetics(int i, int k);
    int DLL_IMPORT flow_settransport(int i, int k, int soret);
    int DLL_IMPORT flow_solveenergyeqn(int i, int j);
    int DLL_IMPORT flow_fixtemperature(int i, int j);
    int DLL_IMPORT flow_setenergyfactor(int i, double e);
    int DLL_IMPORT flow_fixspecies(int i, int j);
    int DLL_IMPORT flow_solvespecies(int i, int j);
    //    int DLL_IMPORT flow_integratechem(int i, double* x, double dt);
    int DLL_IMPORT flow_settemperature(int i, int j, double t);
    int DLL_IMPORT flow_setpressure(int i, double p);
    int DLL_IMPORT flow_setmassfraction(int i, int j, int k, double t);
    int DLL_IMPORT flow_outputtec(int i, double* x, char* fname, 
        char* title, int zone);
    int DLL_IMPORT flow_showsolution(int i, char* fname, double* x);
    int DLL_IMPORT flow_settolerances(int i, int nr, 
        double* rtol, int na, double* atol);
    int DLL_IMPORT flow_resize(int i, int points);
    int DLL_IMPORT flow_setsteadymode(int i);
    int DLL_IMPORT flow_settransientmode(int i, double dt, double* x);

    int DLL_IMPORT flow_restore(int i, int job, char* fname, char* id, 
        int& size_z, double* z, int& size_soln, double* soln);
    int DLL_IMPORT flow_setfixedpoint(int i, int j0, double t0);
    int DLL_IMPORT flow_setboundaries(int i, int nleft, int nright);
    int DLL_IMPORT bdry_new(int type, int iph, int kin);
    int DLL_IMPORT bdry_del(int i);
    int DLL_IMPORT bdry_copy(int i);
    int DLL_IMPORT bdry_assign(int i, int j);
    int DLL_IMPORT bdry_set(int i, int n, double* v);

    int DLL_IMPORT onedim_new(int nd, int* domains, int* types);
    int DLL_IMPORT onedim_del(int i);
    int DLL_IMPORT onedim_addFlow(int i, int n);
    //int DLL_IMPORT onedim_addSurf(int i, int n);
    int DLL_EXPORT onedim_eval(int i, double* x0, double* r);
    int DLL_IMPORT onedim_solve(int i, double* x0, double* x1, int loglevel);
    double DLL_IMPORT onedim_ssnorm(int i, double* x0, double* x1);
    int DLL_IMPORT onedim_setsteadymode(int i);
    int DLL_IMPORT onedim_settransientmode(int i, double dt, double* x);
    int DLL_IMPORT onedim_setnewtonoptions(int i, int maxage);
    int DLL_IMPORT onedim_resize(int i);
    int DLL_IMPORT onedim_writeStats(int i);
    double DLL_IMPORT onedim_timestep(int i, int nsteps, double dt,
        double* x, double* xnew, int loglevel);
    int DLL_IMPORT onedim_save(int i, char* fname, char* id, char* desc, double* soln);

}

#endif
