/**
 * @file ctreactor.h
 *
 * @warning  This module is an experimental part of the %Cantera API and
 *      may be changed or removed without notice.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CTC_REACTOR_H
#define CTC_REACTOR_H

#include "clib_defs.h"

#ifdef __cplusplus
extern "C" {
#endif

    CANTERA_CAPI int reactor_new(const char* type, int n, const char* name);
    CANTERA_CAPI int reactor_del(int i);
    CANTERA_CAPI int reactor_name(int i, int len, char* nbuf);
    CANTERA_CAPI int reactor_setName(int i, const char* name);
    CANTERA_CAPI int reactor_setInitialVolume(int i, double v);
    CANTERA_CAPI int reactor_setChemistry(int i, int cflag);
    CANTERA_CAPI int reactor_setEnergy(int i, int eflag);
    CANTERA_CAPI int reactor_insert(int i, int n);
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
    CANTERA_CAPI int reactor_setMassFlowRate(int i, double mdot);

    CANTERA_CAPI int reactornet_new();
    CANTERA_CAPI int reactornet_del(int i);
    CANTERA_CAPI int reactornet_setInitialTime(int i, double t);
    CANTERA_CAPI int reactornet_setMaxTimeStep(int i, double maxstep);
    CANTERA_CAPI int reactornet_setTolerances(int i, double rtol, double atol);
    CANTERA_CAPI int reactornet_setSensitivityTolerances(int i, double rtol, double atol);
    CANTERA_CAPI int reactornet_addreactor(int i, int n);
    CANTERA_CAPI int reactornet_advance(int i, double t);
    CANTERA_CAPI double reactornet_step(int i);
    CANTERA_CAPI double reactornet_time(int i);
    CANTERA_CAPI double reactornet_rtol(int i);
    CANTERA_CAPI double reactornet_atol(int i);
    CANTERA_CAPI double reactornet_sensitivity(int i, const char* v, int p, int r);

    CANTERA_CAPI int flowdev_new(const char* type, const char* name);
    CANTERA_CAPI int flowdev_del(int i);
    CANTERA_CAPI int flowdev_name(int i, int len, char* nbuf);
    CANTERA_CAPI int flowdev_setName(int i, const char* name);
    CANTERA_CAPI int flowdev_install(int i, int n, int m);
    CANTERA_CAPI int flowdev_setPrimary(int i, int n);
    CANTERA_CAPI double flowdev_massFlowRate(int i);
    CANTERA_CAPI int flowdev_setMassFlowCoeff(int i, double v);
    CANTERA_CAPI int flowdev_setValveCoeff(int i, double v);
    CANTERA_CAPI int flowdev_setPressureCoeff(int i, double v);
    CANTERA_CAPI int flowdev_setPressureFunction(int i, int n);
    CANTERA_CAPI int flowdev_setTimeFunction(int i, int n);

    CANTERA_CAPI int wall_new(const char* type, const char* name);
    CANTERA_CAPI int wall_del(int i);
    CANTERA_CAPI int wall_name(int i, int len, char* nbuf);
    CANTERA_CAPI int wall_setName(int i, const char* name);
    CANTERA_CAPI int wall_install(int i, int n, int m);
    CANTERA_CAPI double wall_expansionRate(int i);
    CANTERA_CAPI double wall_heatRate(int i);
    CANTERA_CAPI double wall_area(int i);
    CANTERA_CAPI int wall_setArea(int i, double v);
    CANTERA_CAPI int wall_setThermalResistance(int i, double rth);
    CANTERA_CAPI int wall_setHeatTransferCoeff(int i, double u);
    CANTERA_CAPI int wall_setHeatFlux(int i, int n);
    CANTERA_CAPI int wall_setExpansionRateCoeff(int i, double k);
    CANTERA_CAPI int wall_setVelocity(int i, int n);
    CANTERA_CAPI int wall_setEmissivity(int i, double epsilon);
    CANTERA_CAPI int wall_ready(int i);

    CANTERA_CAPI int reactorsurface_new(const char* name);
    CANTERA_CAPI int reactorsurface_del(int i);
    CANTERA_CAPI int reactorsurface_name(int i, int len, char* nbuf);
    CANTERA_CAPI int reactorsurface_setName(int i, const char* name);
    CANTERA_CAPI int reactorsurface_install(int i, int n);
    CANTERA_CAPI int reactorsurface_setkinetics(int i, int n);
    CANTERA_CAPI double reactorsurface_area(int i);
    CANTERA_CAPI int reactorsurface_setArea(int i, double v);
    CANTERA_CAPI int reactorsurface_addSensitivityReaction(int i, int rxn);

    CANTERA_CAPI int ct_clearReactors();

#ifdef __cplusplus
}
#endif

#endif
