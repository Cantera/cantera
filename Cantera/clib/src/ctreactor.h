#ifndef CTC_REACTOR_H
#define CTC_REACTOR_H

#include "clib_defs.h"

extern "C" {  

    int DLL_IMPORT reactor_new(int type);
    int DLL_IMPORT reactor_del(int i);
    int DLL_IMPORT reactor_copy(int i);
    int DLL_IMPORT reactor_assign(int i, int j);
    int DLL_IMPORT reactor_setInitialVolume(int i, double v);
    int DLL_IMPORT reactor_setInitialTime(int i, double t);
    int DLL_IMPORT reactor_setEnergy(int i, int eflag);
    int DLL_IMPORT reactor_setThermoMgr(int i, int n);
    int DLL_IMPORT reactor_setKineticsMgr(int i, int n);
    int DLL_IMPORT reactor_advance(int i, double t);
    double DLL_IMPORT reactor_step(int i, double t);
    double DLL_IMPORT reactor_time(int i);
    double DLL_IMPORT reactor_mass(int i);
    double DLL_IMPORT reactor_volume(int i);
    double DLL_IMPORT reactor_density(int i);
    double DLL_IMPORT reactor_temperature(int i);
    double DLL_IMPORT reactor_enthalpy_mass(int i);
    double DLL_IMPORT reactor_intEnergy_mass(int i);
    double DLL_IMPORT reactor_pressure(int i);
    double DLL_IMPORT reactor_massFraction(int i, int k);
    int DLL_IMPORT reactor_nSensParams(int i);
    int DLL_IMPORT reactor_addSensitivityReaction(int i, int rxn);
    int DLL_IMPORT flowReactor_setMassFlowRate(int i, double mdot);

    int DLL_IMPORT reactornet_new();
    int DLL_IMPORT reactornet_del(int i);
    int DLL_IMPORT reactornet_copy(int i);
    int DLL_IMPORT reactornet_assign(int i, int j);
    int DLL_IMPORT reactornet_setInitialTime(int i, double t);
    int DLL_IMPORT reactornet_setMaxTimeStep(int i, double maxstep);
    int DLL_IMPORT reactornet_setTolerances(int i, double rtol, double atol);
    int DLL_IMPORT reactornet_setSensitivityTolerances(int i, double rtol, double atol);
    int DLL_IMPORT reactornet_addreactor(int i, int n);
    int DLL_IMPORT reactornet_advance(int i, double t);
    double DLL_IMPORT reactornet_step(int i, double t);
    double DLL_IMPORT reactornet_time(int i);
    double DLL_IMPORT reactornet_rtol(int i);
    double DLL_IMPORT reactornet_atol(int i);
    double DLL_IMPORT reactornet_sensitivity(int i, char* v, int p, int r);

    int DLL_IMPORT flowdev_new(int type);
    int DLL_IMPORT flowdev_del(int i);
    int DLL_IMPORT flowdev_install(int i, int n, int m);
    int DLL_IMPORT flowdev_setMaster(int i, int n);
    double DLL_IMPORT flowdev_massFlowRate(int i, double time);
    int DLL_IMPORT flowdev_setMassFlowRate(int i, double mdot);
    int DLL_IMPORT flowdev_setParameters(int i, int n, double* v);
    int DLL_IMPORT flowdev_setFunction(int i, int n);
    int DLL_IMPORT flowdev_ready(int i);

     int DLL_IMPORT wall_new(int type);
     int DLL_IMPORT wall_del(int i);
     int DLL_IMPORT wall_copy(int i);
     int DLL_IMPORT wall_assign(int i, int j);
     int DLL_IMPORT wall_install(int i, int n, int m);
     int DLL_IMPORT wall_setkinetics(int i, int n, int m);
     double DLL_IMPORT wall_vdot(int i, double t);
     double DLL_IMPORT wall_Q(int i, double t);
     double DLL_IMPORT wall_area(int i);
     int DLL_IMPORT wall_setArea(int i, double v);
     int DLL_IMPORT wall_setThermalResistance(int i, double rth);
     int DLL_IMPORT wall_setHeatTransferCoeff(int i, double u);
     int DLL_IMPORT wall_setHeatFlux(int i, int n);
     int DLL_IMPORT wall_setExpansionRateCoeff(int i, double k);
     int DLL_IMPORT wall_setVelocity(int i, int n);
    int DLL_IMPORT wall_setEmissivity(int i, double epsilon);
     int DLL_IMPORT wall_ready(int i);
    int DLL_IMPORT wall_addSensitivityReaction(int i, int lr, int rxn);

}

#endif
