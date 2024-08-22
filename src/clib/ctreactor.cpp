/**
 * @file ctreactor.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/clib/ctreactor.h"

// Cantera includes
#include "cantera/numerics/Func1.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/zerodim.h"
#include "cantera/base/stringUtils.h"
#include "clib_utils.h"

using namespace Cantera;

typedef Cabinet<ReactorBase> ReactorCabinet;
typedef Cabinet<ReactorNet> NetworkCabinet;
typedef Cabinet<Connector> ConnectorCabinet;
typedef Cabinet<Func1> FuncCabinet;
typedef Cabinet<ThermoPhase> ThermoCabinet;
typedef Cabinet<Kinetics> KineticsCabinet;
typedef Cabinet<Solution> SolutionCabinet;
typedef Cabinet<ReactorSurface> ReactorSurfaceCabinet;

template<> ReactorCabinet* ReactorCabinet::s_storage = 0;
template<> NetworkCabinet* NetworkCabinet::s_storage = 0;
template<> ConnectorCabinet* ConnectorCabinet::s_storage = 0;
template<> ReactorSurfaceCabinet* ReactorSurfaceCabinet::s_storage = 0;
template<> FuncCabinet* FuncCabinet::s_storage; // defined in ctfunc.cpp
template<> ThermoCabinet* ThermoCabinet::s_storage; // defined in ct.cpp
template<> KineticsCabinet* KineticsCabinet::s_storage; // defined in ct.cpp
template<> SolutionCabinet* SolutionCabinet::s_storage; // defined in ct.cpp

extern "C" {

    // reactor

    int reactor_new(const char* type, int n, const char* name)
    {
        try {
            return ReactorCabinet::add(newReactor(type, SolutionCabinet::at(n), name));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int reactor_del(int i)
    {
        try {
            ReactorCabinet::del(i);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int reactor_name(int i, int len, char* nbuf)
    {
        try {
            return static_cast<int>(
                copyString(ReactorCabinet::at(i)->name(), nbuf, len));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int reactor_setName(int i, const char* name)
    {
        try {
            ReactorCabinet::at(i)->setName(name);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int reactor_setInitialVolume(int i, double v)
    {
        try {
            ReactorCabinet::at(i)->setInitialVolume(v);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int reactor_setThermoMgr(int i, int n)
    {
        try {
            ReactorCabinet::at(i)->setThermoMgr(*ThermoCabinet::at(n));
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int reactor_setKineticsMgr(int i, int n)
    {
        try {
            ReactorCabinet::at(i)->setKineticsMgr(*KineticsCabinet::at(n));
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int reactor_insert(int i, int n)
    {
        try {
            ReactorCabinet::as<Reactor>(i)->insert(SolutionCabinet::at(n));
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    double reactor_mass(int i)
    {
        try {
            return ReactorCabinet::at(i)->mass();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double reactor_volume(int i)
    {
        try {
            return ReactorCabinet::at(i)->volume();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double reactor_density(int i)
    {
        try {
            return ReactorCabinet::at(i)->density();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double reactor_temperature(int i)
    {
        try {
            return ReactorCabinet::at(i)->temperature();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double reactor_enthalpy_mass(int i)
    {
        try {
            return ReactorCabinet::at(i)->enthalpy_mass();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double reactor_intEnergy_mass(int i)
    {
        try {
            return ReactorCabinet::at(i)->intEnergy_mass();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double reactor_pressure(int i)
    {
        try {
            return ReactorCabinet::at(i)->pressure();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double reactor_massFraction(int i, int k)
    {
        try {
            return ReactorCabinet::at(i)->massFraction(k);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    int reactor_setChemistry(int i, int cflag)
    {
        try {
            ReactorCabinet::as<Reactor>(i)->setChemistry(cflag != 0);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int reactor_setEnergy(int i, int eflag)
    {
        try {
            ReactorCabinet::as<Reactor>(i)->setEnergy(eflag);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int reactor_setMassFlowRate(int i, double mdot)
    {
        try {
            ReactorCabinet::as<FlowReactor>(i)->setMassFlowRate(mdot);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    size_t reactor_nSensParams(int i)
    {
        try {
            return ReactorCabinet::as<Reactor>(i)->nSensParams();
        } catch (...) {
            return handleAllExceptions(npos, npos);
        }
    }

    int reactor_addSensitivityReaction(int i, int rxn)
    {
        try {
            ReactorCabinet::as<Reactor>(i)->addSensitivityReaction(rxn);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    // reactor networks

    int reactornet_new()
    {
        try {
            return NetworkCabinet::add(make_shared<ReactorNet>());
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int reactornet_del(int i)
    {
        try {
            NetworkCabinet::del(i);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int reactornet_setInitialTime(int i, double t)
    {
        try {
            NetworkCabinet::at(i)->setInitialTime(t);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int reactornet_setMaxTimeStep(int i, double maxstep)
    {
        try {
            NetworkCabinet::at(i)->setMaxTimeStep(maxstep);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int reactornet_setTolerances(int i, double rtol, double atol)
    {
        try {
            NetworkCabinet::at(i)->setTolerances(rtol, atol);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int reactornet_setSensitivityTolerances(int i, double rtol, double atol)
    {
        try {
            NetworkCabinet::at(i)->setSensitivityTolerances(rtol, atol);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int reactornet_addreactor(int i, int n)
    {
        try {
            NetworkCabinet::at(i)->addReactor(
                dynamic_cast<Reactor&>(*ReactorCabinet::at(n)));
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int reactornet_advance(int i, double t)
    {
        try {
            NetworkCabinet::at(i)->advance(t);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    double reactornet_step(int i)
    {
        try {
            return NetworkCabinet::at(i)->step();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double reactornet_time(int i)
    {
        try {
            return NetworkCabinet::at(i)->time();
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    double reactornet_rtol(int i)
    {
        try {
            return NetworkCabinet::at(i)->rtol();
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    double reactornet_atol(int i)
    {
        try {
            return NetworkCabinet::at(i)->atol();
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    double reactornet_sensitivity(int i, const char* v, int p, int r)
    {
        try {
            return NetworkCabinet::at(i)->sensitivity(v, p, r);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    // flow devices

    int flowdev_new(const char* type, const char* name)
    {
        try {
            return ConnectorCabinet::add(newFlowDevice(type, name));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int flowdev_del(int i)
    {
        try {
            ConnectorCabinet::del(i);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int flowdev_name(int i, int len, char* nbuf)
    {
        try {
            return static_cast<int>(
                copyString(ConnectorCabinet::at(i)->name(), nbuf, len));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int flowdev_setName(int i, const char* name)
    {
        try {
            ConnectorCabinet::at(i)->setName(name);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int flowdev_install(int i, int n, int m)
    {
        try {
            bool ok = ConnectorCabinet::as<FlowDevice>(i)->install(
                *ReactorCabinet::at(n), *ReactorCabinet::at(m));
            if (!ok) {
                throw CanteraError("flowdev_install",
                                   "Could not install flow device.");
            }
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int flowdev_setPrimary(int i, int n)
    {
        try {
            ConnectorCabinet::as<PressureController>(i)->setPrimary(
                ConnectorCabinet::as<FlowDevice>(n).get());
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    double flowdev_massFlowRate(int i)
    {
        try {
            return ConnectorCabinet::as<FlowDevice>(i)->massFlowRate();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    int flowdev_setMassFlowCoeff(int i, double v)
    {
        try {
            ConnectorCabinet::as<MassFlowController>(i)->setMassFlowCoeff(v);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int flowdev_setValveCoeff(int i, double v)
    {
        try {
            ConnectorCabinet::as<Valve>(i)->setValveCoeff(v);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int flowdev_setPressureCoeff(int i, double v)
    {
        try {
            ConnectorCabinet::as<PressureController>(i)->setPressureCoeff(v);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int flowdev_setPressureFunction(int i, int n)
    {
        try {
            ConnectorCabinet::as<FlowDevice>(i)->setPressureFunction(
                FuncCabinet::at(n).get());
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int flowdev_setTimeFunction(int i, int n)
    {
        try {
            ConnectorCabinet::as<FlowDevice>(i)->setTimeFunction(
                FuncCabinet::at(n).get());
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    /////////////    Walls   ///////////////////////

    int wall_new(const char* type, const char* name)
    {
        try {
            return ConnectorCabinet::add(newWall(type, name));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int wall_del(int i)
    {
        try {
            ConnectorCabinet::del(i);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int wall_name(int i, int len, char* nbuf)
    {
        try {
            return static_cast<int>(
                copyString(ConnectorCabinet::at(i)->name(), nbuf, len));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int wall_setName(int i, const char* name)
    {
        try {
            ConnectorCabinet::at(i)->setName(name);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int wall_install(int i, int n, int m)
    {
        try {
            ConnectorCabinet::as<Wall>(i)->install(
                *ReactorCabinet::at(n), *ReactorCabinet::at(m));
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    double wall_expansionRate(int i)
    {
        try {
            return ConnectorCabinet::as<Wall>(i)->expansionRate();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double wall_heatRate(int i)
    {
        try {
            return ConnectorCabinet::as<Wall>(i)->heatRate();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double wall_area(int i)
    {
        try {
            return ConnectorCabinet::as<Wall>(i)->area();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    int wall_setArea(int i, double v)
    {
        try {
            ConnectorCabinet::as<Wall>(i)->setArea(v);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int wall_setThermalResistance(int i, double rth)
    {
        try {
            ConnectorCabinet::as<Wall>(i)->setThermalResistance(rth);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int wall_setHeatTransferCoeff(int i, double u)
    {
        try {
            ConnectorCabinet::as<Wall>(i)->setHeatTransferCoeff(u);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int wall_setHeatFlux(int i, int n)
    {
        try {
            ConnectorCabinet::as<Wall>(i)->setHeatFlux(FuncCabinet::at(n).get());
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int wall_setExpansionRateCoeff(int i, double k)
    {
        try {
            ConnectorCabinet::as<Wall>(i)->setExpansionRateCoeff(k);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int wall_setVelocity(int i, int n)
    {
        try {
            ConnectorCabinet::as<Wall>(i)->setVelocity(FuncCabinet::at(n).get());
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int wall_setEmissivity(int i, double epsilon)
    {
        try {
            ConnectorCabinet::as<Wall>(i)->setEmissivity(epsilon);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int wall_ready(int i)
    {
        try {
            return int(ConnectorCabinet::as<Wall>(i)->ready());
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    // ReactorSurface

    int reactorsurface_new(const char* name)
    {
        try {
            return ReactorSurfaceCabinet::add(make_shared<ReactorSurface>(name));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int reactorsurface_del(int i)
    {
        try {
            ReactorSurfaceCabinet::del(i);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int reactorsurface_name(int i, int len, char* nbuf)
    {
        try {
            return static_cast<int>(
                copyString(ReactorSurfaceCabinet::at(i)->name(), nbuf, len));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int reactorsurface_setName(int i, const char* name)
    {
        try {
            ReactorSurfaceCabinet::at(i)->setName(name);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int reactorsurface_install(int i, int n)
    {
        try {
            ReactorCabinet::at(n)->addSurface(ReactorSurfaceCabinet::at(i).get());
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int reactorsurface_setkinetics(int i, int n)
    {
        try {
            ReactorSurfaceCabinet::at(i)->setKinetics(KineticsCabinet::at(n).get());
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    double reactorsurface_area(int i)
    {
        try {
            return ReactorSurfaceCabinet::at(i)->area();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    int reactorsurface_setArea(int i, double v)
    {
        try {
            ReactorSurfaceCabinet::at(i)->setArea(v);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int reactorsurface_addSensitivityReaction(int i, int rxn)
    {
        try {
            ReactorSurfaceCabinet::at(i)->addSensitivityReaction(rxn);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int ct_clearReactors()
    {
        try {
            ReactorCabinet::clear();
            NetworkCabinet::clear();
            ConnectorCabinet::clear();
            ReactorSurfaceCabinet::clear();
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }
}
