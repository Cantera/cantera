/**
 * @file ctreactor.h
 */
/*
 *      $Id: ctreactor.h,v 1.7 2009/07/11 17:16:09 hkmoffa Exp $
 */

#ifndef CTC_REACTOR_H
#define CTC_REACTOR_H

#include "clib_defs.h"

extern "C" {  

    EEXXTT int DLL_CPREFIX reactor_new(int type);
    EEXXTT int DLL_CPREFIX reactor_del(int i);
    EEXXTT int DLL_CPREFIX reactor_copy(int i);
    EEXXTT int DLL_CPREFIX reactor_assign(int i, int j);
    EEXXTT int DLL_CPREFIX reactor_setInitialVolume(int i, double v);
    EEXXTT int DLL_CPREFIX reactor_setInitialTime(int i, double t);
    EEXXTT int DLL_CPREFIX reactor_setEnergy(int i, int eflag);
    EEXXTT int DLL_CPREFIX reactor_setThermoMgr(int i, int n);
    EEXXTT int DLL_CPREFIX reactor_setKineticsMgr(int i, int n);
    EEXXTT int DLL_CPREFIX reactor_advance(int i, double t);
    EEXXTT double DLL_CPREFIX reactor_step(int i, double t);
    EEXXTT double DLL_CPREFIX reactor_time(int i);
    EEXXTT double DLL_CPREFIX reactor_mass(int i);
    EEXXTT double DLL_CPREFIX reactor_volume(int i);
    EEXXTT double DLL_CPREFIX reactor_density(int i);
    EEXXTT double DLL_CPREFIX reactor_temperature(int i);
    EEXXTT double DLL_CPREFIX reactor_enthalpy_mass(int i);
    EEXXTT double DLL_CPREFIX reactor_intEnergy_mass(int i);
    EEXXTT double DLL_CPREFIX reactor_pressure(int i);
    EEXXTT double DLL_CPREFIX reactor_massFraction(int i, int k);
    EEXXTT int DLL_CPREFIX reactor_nSensParams(int i);
    EEXXTT int DLL_CPREFIX reactor_addSensitivityReaction(int i, int rxn);
    EEXXTT int DLL_CPREFIX flowReactor_setMassFlowRate(int i, double mdot);

    EEXXTT int DLL_CPREFIX reactornet_new();
    EEXXTT int DLL_CPREFIX reactornet_del(int i);
    EEXXTT int DLL_CPREFIX reactornet_copy(int i);
    EEXXTT int DLL_CPREFIX reactornet_assign(int i, int j);
    EEXXTT int DLL_CPREFIX reactornet_setInitialTime(int i, double t);
    EEXXTT int DLL_CPREFIX reactornet_setMaxTimeStep(int i, double maxstep);
    EEXXTT int DLL_CPREFIX reactornet_setTolerances(int i, double rtol, double atol);
    EEXXTT int DLL_CPREFIX reactornet_setSensitivityTolerances(int i, double rtol, double atol);
    EEXXTT int DLL_CPREFIX reactornet_addreactor(int i, int n);
    EEXXTT int DLL_CPREFIX reactornet_advance(int i, double t);
    EEXXTT double DLL_CPREFIX reactornet_step(int i, double t);
    EEXXTT double DLL_CPREFIX reactornet_time(int i);
    EEXXTT double DLL_CPREFIX reactornet_rtol(int i);
    EEXXTT double DLL_CPREFIX reactornet_atol(int i);
    EEXXTT double DLL_CPREFIX reactornet_sensitivity(int i, char* v, int p, int r);

    EEXXTT int DLL_CPREFIX flowdev_new(int type);
    EEXXTT int DLL_CPREFIX flowdev_del(int i);
    EEXXTT int DLL_CPREFIX flowdev_install(int i, int n, int m);
    EEXXTT int DLL_CPREFIX flowdev_setMaster(int i, int n);
    EEXXTT double DLL_CPREFIX flowdev_massFlowRate(int i, double time);
    EEXXTT int DLL_CPREFIX flowdev_setMassFlowRate(int i, double mdot);
    EEXXTT int DLL_CPREFIX flowdev_setParameters(int i, int n, double* v);
    EEXXTT int DLL_CPREFIX flowdev_setFunction(int i, int n);
    EEXXTT int DLL_CPREFIX flowdev_ready(int i);

     EEXXTT int DLL_CPREFIX wall_new(int type);
     EEXXTT int DLL_CPREFIX wall_del(int i);
     EEXXTT int DLL_CPREFIX wall_copy(int i);
     EEXXTT int DLL_CPREFIX wall_assign(int i, int j);
     EEXXTT int DLL_CPREFIX wall_install(int i, int n, int m);
     EEXXTT int DLL_CPREFIX wall_setkinetics(int i, int n, int m);
     EEXXTT double DLL_CPREFIX wall_vdot(int i, double t);
     EEXXTT double DLL_CPREFIX wall_Q(int i, double t);
     EEXXTT double DLL_CPREFIX wall_area(int i);
     EEXXTT int DLL_CPREFIX wall_setArea(int i, double v);
     EEXXTT int DLL_CPREFIX wall_setThermalResistance(int i, double rth);
     EEXXTT int DLL_CPREFIX wall_setHeatTransferCoeff(int i, double u);
     EEXXTT int DLL_CPREFIX wall_setHeatFlux(int i, int n);
     EEXXTT int DLL_CPREFIX wall_setExpansionRateCoeff(int i, double k);
     EEXXTT int DLL_CPREFIX wall_setVelocity(int i, int n);
     EEXXTT int DLL_CPREFIX wall_setEmissivity(int i, double epsilon);
     EEXXTT int DLL_CPREFIX wall_ready(int i);
     EEXXTT int DLL_CPREFIX wall_addSensitivityReaction(int i, int lr, int rxn);

}

#endif
