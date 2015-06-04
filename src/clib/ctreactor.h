/**
 * @file ctreactor.h
 */
#ifndef CTC_REACTOR_H
#define CTC_REACTOR_H

#include "clib_defs.h"

extern "C" {
    CANTERA_CAPI int reactor_new(int type);
    CANTERA_CAPI int reactor_del(int i);
    CANTERA_CAPI int reactor_copy(int i);
    CANTERA_CAPI int reactor_assign(int i, int j);
    CANTERA_CAPI int reactor_setInitialVolume(int i, double v);
    CANTERA_CAPI int reactor_setEnergy(int i, int eflag);
    CANTERA_CAPI int reactor_setThermoMgr(int i, int n);
    CANTERA_CAPI int reactor_setKineticsMgr(int i, int n);
    CANTERA_CAPI double reactor_mass(int i);
    CANTERA_CAPI double reactor_volume(int i);
    CANTERA_CAPI double reactor_density(int i);
    CANTERA_CAPI double reactor_temperature(int i);
    CANTERA_CAPI double reactor_enthalpy_mass(int i);
    CANTERA_CAPI double reactor_intEnergy_mass(int i);
    CANTERA_CAPI double reactor_pressure(int i);
    CANTERA_CAPI double reactor_massFraction(int i, int k);
    CANTERA_CAPI size_t reactor_nSensParams(int i);
    CANTERA_CAPI int reactor_addSensitivityReaction(int i, int rxn);
    CANTERA_CAPI int flowReactor_setMassFlowRate(int i, double mdot);

    CANTERA_CAPI int reactornet_new();
    CANTERA_CAPI int reactornet_del(int i);
    CANTERA_CAPI int reactornet_copy(int i);
    CANTERA_CAPI int reactornet_assign(int i, int j);
    CANTERA_CAPI int reactornet_setInitialTime(int i, double t);
    CANTERA_CAPI int reactornet_setMaxTimeStep(int i, double maxstep);
    CANTERA_CAPI int reactornet_setTolerances(int i, double rtol, double atol);
    CANTERA_CAPI int reactornet_setSensitivityTolerances(int i, double rtol, double atol);
    CANTERA_CAPI int reactornet_addreactor(int i, int n);
    CANTERA_CAPI int reactornet_advance(int i, double t);
    CANTERA_CAPI double reactornet_step(int i, double t);
    CANTERA_CAPI double reactornet_time(int i);
    CANTERA_CAPI double reactornet_rtol(int i);
    CANTERA_CAPI double reactornet_atol(int i);
    CANTERA_CAPI double reactornet_sensitivity(int i, char* v, int p, int r);

    CANTERA_CAPI int flowdev_new(int type);
    CANTERA_CAPI int flowdev_del(int i);
    CANTERA_CAPI int flowdev_install(int i, int n, int m);
    CANTERA_CAPI int flowdev_setMaster(int i, int n);
    CANTERA_CAPI double flowdev_massFlowRate(int i, double time);
    CANTERA_CAPI int flowdev_setMassFlowRate(int i, double mdot);
    CANTERA_CAPI int flowdev_setParameters(int i, int n, double* v);
    CANTERA_CAPI int flowdev_setFunction(int i, int n);
    CANTERA_CAPI int flowdev_ready(int i);

    CANTERA_CAPI int wall_new(int type);
    CANTERA_CAPI int wall_del(int i);
    CANTERA_CAPI int wall_copy(int i);
    CANTERA_CAPI int wall_assign(int i, int j);
    CANTERA_CAPI int wall_install(int i, int n, int m);
    CANTERA_CAPI int wall_setkinetics(int i, int n, int m);
    CANTERA_CAPI double wall_vdot(int i, double t);
    CANTERA_CAPI double wall_Q(int i, double t);
    CANTERA_CAPI double wall_area(int i);
    CANTERA_CAPI int wall_setArea(int i, double v);
    CANTERA_CAPI int wall_setThermalResistance(int i, double rth);
    CANTERA_CAPI int wall_setHeatTransferCoeff(int i, double u);
    CANTERA_CAPI int wall_setHeatFlux(int i, int n);
    CANTERA_CAPI int wall_setExpansionRateCoeff(int i, double k);
    CANTERA_CAPI int wall_setVelocity(int i, int n);
    CANTERA_CAPI int wall_setEmissivity(int i, double epsilon);
    CANTERA_CAPI int wall_ready(int i);
    CANTERA_CAPI int wall_addSensitivityReaction(int i, int lr, int rxn);

    CANTERA_CAPI int clear_reactors();
}

#endif
