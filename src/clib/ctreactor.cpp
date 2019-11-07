/**
 * @file ctreactor.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#define CANTERA_USE_INTERNAL
#include "cantera/clib/ctreactor.h"

// Cantera includes
#include "cantera/numerics/Func1.h"
#include "cantera/zerodim.h"
#include "Cabinet.h"

using namespace Cantera;
using namespace std;

typedef Cabinet<ReactorBase> ReactorCabinet;
typedef Cabinet<ReactorNet> NetworkCabinet;
typedef Cabinet<FlowDevice> FlowDeviceCabinet;
typedef Cabinet<WallBase> WallCabinet;
typedef Cabinet<Func1> FuncCabinet;
typedef Cabinet<ThermoPhase> ThermoCabinet;
typedef Cabinet<Kinetics> KineticsCabinet;
typedef Cabinet<ReactorSurface> ReactorSurfaceCabinet;

template<> ReactorCabinet* ReactorCabinet::s_storage = 0;
template<> NetworkCabinet* NetworkCabinet::s_storage = 0;
template<> FlowDeviceCabinet* FlowDeviceCabinet::s_storage = 0;
template<> WallCabinet* WallCabinet::s_storage = 0;
template<> ReactorSurfaceCabinet* ReactorSurfaceCabinet::s_storage = 0;

extern "C" {

    // reactor

    //! @deprecated To be changed after Cantera 2.5.
    int reactor_new(int type)
    {
        warn_deprecated("reactor_new(int)",
                        "To be changed after Cantera 2.5. "
                        "Argument changed to string instead of int; use"
                        "reactor_new2(char*) during transition.");
        try {
            ReactorBase* r = ReactorFactory::factory()->newReactor(type);
            return ReactorCabinet::add(r);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int reactor_new2(const char* type)
    {
        try {
            ReactorBase* r = ReactorFactory::factory()->newReactor(type);
            return ReactorCabinet::add(r);
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

    int reactor_setInitialVolume(int i, double v)
    {
        try {
            ReactorCabinet::item(i).setInitialVolume(v);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int reactor_setThermoMgr(int i, int n)
    {
        try {
            ReactorCabinet::item(i).setThermoMgr(ThermoCabinet::item(n));
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int reactor_setKineticsMgr(int i, int n)
    {
        try {
            ReactorCabinet::item(i).setKineticsMgr(KineticsCabinet::item(n));
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    double reactor_mass(int i)
    {
        try {
            return ReactorCabinet::item(i).mass();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double reactor_volume(int i)
    {
        try {
            return ReactorCabinet::item(i).volume();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double reactor_density(int i)
    {
        try {
            return ReactorCabinet::item(i).density();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double reactor_temperature(int i)
    {
        try {
            return ReactorCabinet::item(i).temperature();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double reactor_enthalpy_mass(int i)
    {
        try {
            return ReactorCabinet::item(i).enthalpy_mass();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double reactor_intEnergy_mass(int i)
    {
        try {
            return ReactorCabinet::item(i).intEnergy_mass();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double reactor_pressure(int i)
    {
        try {
            return ReactorCabinet::item(i).pressure();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double reactor_massFraction(int i, int k)
    {
        try {
            return ReactorCabinet::item(i).massFraction(k);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    int reactor_setChemistry(int i, int cflag)
    {
        try {
            ReactorCabinet::get<Reactor>(i).setChemistry(cflag != 0);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int reactor_setEnergy(int i, int eflag)
    {
        try {
            ReactorCabinet::get<Reactor>(i).setEnergy(eflag);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int flowReactor_setMassFlowRate(int i, double mdot)
    {
        try {
            ReactorCabinet::get<FlowReactor>(i).setMassFlowRate(mdot);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    size_t reactor_nSensParams(int i)
    {
        try {
            return ReactorCabinet::get<Reactor>(i).nSensParams();
        } catch (...) {
            return handleAllExceptions(npos, npos);
        }
    }

    int reactor_addSensitivityReaction(int i, int rxn)
    {
        try {
            ReactorCabinet::get<Reactor>(i).addSensitivityReaction(rxn);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    // reactor networks

    int reactornet_new()
    {
        try {
            return NetworkCabinet::add(new ReactorNet());
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
            NetworkCabinet::item(i).setInitialTime(t);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int reactornet_setMaxTimeStep(int i, double maxstep)
    {
        try {
            NetworkCabinet::item(i).setMaxTimeStep(maxstep);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int reactornet_setTolerances(int i, double rtol, double atol)
    {
        try {
            NetworkCabinet::item(i).setTolerances(rtol, atol);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int reactornet_setSensitivityTolerances(int i, double rtol, double atol)
    {
        try {
            NetworkCabinet::item(i).setSensitivityTolerances(rtol, atol);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int reactornet_addreactor(int i, int n)
    {
        try {
            NetworkCabinet::item(i).addReactor(
                dynamic_cast<Reactor&>(ReactorCabinet::item(n)));
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int reactornet_advance(int i, double t)
    {
        try {
            NetworkCabinet::item(i).advance(t);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    double reactornet_step(int i)
    {
        try {
            return NetworkCabinet::item(i).step();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double reactornet_time(int i)
    {
        try {
            return NetworkCabinet::item(i).time();
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    double reactornet_rtol(int i)
    {
        try {
            return NetworkCabinet::item(i).rtol();
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    double reactornet_atol(int i)
    {
        try {
            return NetworkCabinet::item(i).atol();
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    double reactornet_sensitivity(int i, const char* v, int p, int r)
    {
        try {
            return NetworkCabinet::item(i).sensitivity(v, p, r);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    // flow devices

    int flowdev_new(int type)
    {
        warn_deprecated("flowdev_new(int)",
                        "To be changed after Cantera 2.5. "
                        "Argument changed to string instead of int; use"
                        "flowdev_new2(char*) during transition.");
        try {
            FlowDevice* f = FlowDeviceFactory::factory()->newFlowDevice(type);
            return FlowDeviceCabinet::add(f);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int flowdev_new2(const char* type)
    {
        try {
            FlowDevice* f = FlowDeviceFactory::factory()->newFlowDevice(type);
            return FlowDeviceCabinet::add(f);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int flowdev_del(int i)
    {
        try {
            FlowDeviceCabinet::del(i);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int flowdev_install(int i, int n, int m)
    {
        try {
            bool ok = FlowDeviceCabinet::item(i).install(ReactorCabinet::item(n),
                      ReactorCabinet::item(m));
            if (!ok) {
                throw CanteraError("flowdev_install",
                                   "Could not install flow device.");
            }
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int flowdev_setMaster(int i, int n)
    {
        try {
            FlowDeviceCabinet::get<PressureController>(i).setMaster(
                &FlowDeviceCabinet::item(n));
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    double flowdev_massFlowRate(int i, double time)
    {
        try {
            return FlowDeviceCabinet::item(i).massFlowRate(time);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    int flowdev_setMassFlowRate(int i, double mdot)
    {
        /* @deprecated To be removed after Cantera 2.5. */
        try {
            FlowDeviceCabinet::item(i).setMassFlowRate(mdot);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int flowdev_setParameters(int i, int n, const double* v)
    {
        /* @deprecated To be removed after Cantera 2.5. */
        try {
            FlowDeviceCabinet::item(i).setParameters(n, v);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int flowdev_setMassFlowCoeff(int i, double v)
    {
        try {
            FlowDeviceCabinet::get<MassFlowController>(i).setMassFlowCoeff(v);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int flowdev_setValveCoeff(int i, double v)
    {
        try {
            FlowDeviceCabinet::get<Valve>(i).setValveCoeff(v);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int flowdev_setPressureCoeff(int i, double v)
    {
        try {
            FlowDeviceCabinet::get<PressureController>(i).setPressureCoeff(v);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int flowdev_setFunction(int i, int n)
    {
        /* @deprecated To be removed after Cantera 2.5. */
        try {
            FlowDeviceCabinet::item(i).setFunction(&FuncCabinet::item(n));
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int flowdev_setPressureFunction(int i, int n)
    {
        try {
            FlowDeviceCabinet::item(i).setPressureFunction(&FuncCabinet::item(n));
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int flowdev_setTimeFunction(int i, int n)
    {
        try {
            FlowDeviceCabinet::item(i).setTimeFunction(&FuncCabinet::item(n));
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    /////////////    Walls   ///////////////////////

    int wall_new(int type)
    {
        warn_deprecated("wall_new(int)",
                        "To be changed after Cantera 2.5. "
                        "Argument changed to string instead of int; use"
                        "wall_new2(char*) during transition.");
        try {
            WallBase* w = WallFactory::factory()->newWall(type);
            return WallCabinet::add(w);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int wall_new2(const char* type)
    {
        try {
            WallBase* w = WallFactory::factory()->newWall(type);
            return WallCabinet::add(w);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int wall_del(int i)
    {
        try {
            WallCabinet::del(i);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int wall_install(int i, int n, int m)
    {
        try {
            WallCabinet::item(i).install(ReactorCabinet::item(n),
                                         ReactorCabinet::item(m));
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    double wall_vdot(int i, double t)
    {
        try {
            return WallCabinet::item(i).vdot(t);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double wall_Q(int i, double t)
    {
        try {
            return WallCabinet::item(i).Q(t);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double wall_area(int i)
    {
        try {
            return WallCabinet::item(i).area();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    int wall_setArea(int i, double v)
    {
        try {
            WallCabinet::item(i).setArea(v);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int wall_setThermalResistance(int i, double rth)
    {
        try {
            WallCabinet::get<Wall>(i).setThermalResistance(rth);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int wall_setHeatTransferCoeff(int i, double u)
    {
        try {
            WallCabinet::get<Wall>(i).setHeatTransferCoeff(u);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int wall_setHeatFlux(int i, int n)
    {
        try {
            WallCabinet::get<Wall>(i).setHeatFlux(&FuncCabinet::item(n));
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int wall_setExpansionRateCoeff(int i, double k)
    {
        try {
            WallCabinet::get<Wall>(i).setExpansionRateCoeff(k);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int wall_setVelocity(int i, int n)
    {
        try {
            WallCabinet::get<Wall>(i).setVelocity(&FuncCabinet::item(n));
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int wall_setEmissivity(int i, double epsilon)
    {
        try {
            WallCabinet::get<Wall>(i).setEmissivity(epsilon);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int wall_ready(int i)
    {
        try {
            return int(WallCabinet::item(i).ready());
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    // ReactorSurface

    int reactorsurface_new(int type)
    {
        try {
            return ReactorSurfaceCabinet::add(new ReactorSurface());
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

    int reactorsurface_install(int i, int n)
    {
        try {
            ReactorCabinet::item(n).addSurface(&ReactorSurfaceCabinet::item(i));
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int reactorsurface_setkinetics(int i, int n)
    {
        try {
            ReactorSurfaceCabinet::item(i).setKinetics(&KineticsCabinet::item(n));
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    double reactorsurface_area(int i)
    {
        try {
            return ReactorSurfaceCabinet::item(i).area();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    int reactorsurface_setArea(int i, double v)
    {
        try {
            ReactorSurfaceCabinet::item(i).setArea(v);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int reactorsurface_addSensitivityReaction(int i, int rxn)
    {
        try {
            ReactorSurfaceCabinet::item(i).addSensitivityReaction(rxn);
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
            FlowDeviceCabinet::clear();
            WallCabinet::clear();
            ReactorSurfaceCabinet::clear();
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }
}
