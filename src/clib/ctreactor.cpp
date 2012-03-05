/**
 * @file ctreactor.cpp
 */
#define CANTERA_USE_INTERNAL
#include "ctreactor.h"

// Cantera includes
#include "cantera/zeroD/Reactor.h"
#include "cantera/zeroD/FlowReactor.h"
#include "cantera/zeroD/ConstPressureReactor.h"
#include "cantera/zeroD/ReactorNet.h"
#include "cantera/zeroD/Reservoir.h"
#include "cantera/zeroD/Wall.h"
#include "cantera/zeroD/flowControllers.h"
#include "Cabinet.h"

using namespace Cantera;
using namespace std;

typedef Cabinet<ReactorBase> ReactorCabinet;
typedef Cabinet<ReactorNet> NetworkCabinet;
typedef Cabinet<FlowDevice> FlowDeviceCabinet;
typedef Cabinet<Wall> WallCabinet;
typedef Cabinet<Func1> FuncCabinet;
typedef Cabinet<ThermoPhase> ThermoCabinet;
typedef Cabinet<Kinetics> KineticsCabinet;

template<> ReactorCabinet* ReactorCabinet::__storage = 0;
template<> NetworkCabinet* NetworkCabinet::__storage = 0;
template<> FlowDeviceCabinet* FlowDeviceCabinet::__storage = 0;
template<> WallCabinet* WallCabinet::__storage = 0;

extern "C" {

    // reactor

    int reactor_new(int type)
    {
        ReactorBase* r=0;
        if (type == ReactorType) {
            r = new Reactor();
        } else if (type == FlowReactorType) {
            r = new FlowReactor();
        } else if (type == ConstPressureReactorType) {
            r = new ConstPressureReactor();
        } else if (type == ReservoirType) {
            r = new Reservoir();
        } else {
            r = new ReactorBase();
        }
        return ReactorCabinet::add(r);
    }

    int reactor_del(int i)
    {
        ReactorCabinet::del(i);
        return 0;
    }

    int reactor_copy(int i)
    {
        return ReactorCabinet::newCopy(i);
    }

    int reactor_assign(int i, int j)
    {
        return ReactorCabinet::assign(i,j);
    }

    int reactor_setInitialVolume(int i, double v)
    {
        ReactorCabinet::item(i).setInitialVolume(v);
        return 0;
    }

    int reactor_setInitialTime(int i, double t)
    {
        ReactorCabinet::item(i).setInitialTime(t);
        return 0;
    }

    int reactor_setThermoMgr(int i, int n)
    {
        ReactorCabinet::item(i).setThermoMgr(ThermoCabinet::item(n));
        return 0;
    }

    int reactor_setKineticsMgr(int i, int n)
    {
        ReactorBase* r = &ReactorCabinet::item(i);
        if (r->type() >= ReactorType) {
            ((Reactor*)r)->setKineticsMgr(KineticsCabinet::item(n));
        }
        return 0;
    }

    int reactor_advance(int i, double t)
    {
        try {
            ReactorCabinet::item(i).advance(t);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    double reactor_step(int i, double t)
    {
        return ReactorCabinet::item(i).step(t);
    }

    double reactor_time(int i)
    {
        return ReactorCabinet::item(i).time();
    }

    double reactor_mass(int i)
    {
        return ReactorCabinet::item(i).mass();
    }

    double reactor_volume(int i)
    {
        return ReactorCabinet::item(i).volume();
    }

    double reactor_density(int i)
    {
        return ReactorCabinet::item(i).density();
    }

    double reactor_temperature(int i)
    {
        return ReactorCabinet::item(i).temperature();
    }

    double reactor_enthalpy_mass(int i)
    {
        return ReactorCabinet::item(i).enthalpy_mass();
    }

    double reactor_intEnergy_mass(int i)
    {
        return ReactorCabinet::item(i).intEnergy_mass();
    }

    double reactor_pressure(int i)
    {
        return ReactorCabinet::item(i).pressure();
    }

    double reactor_massFraction(int i, int k)
    {
        return ReactorCabinet::item(i).massFraction(k);
    }

    int reactor_setEnergy(int i, int eflag)
    {
        ReactorBase* r = &ReactorCabinet::item(i);
        if (r->type() >= ReactorType) {
            ((Reactor*)r)->setEnergy(eflag);
        }
        return 0;
    }

    int flowReactor_setMassFlowRate(int i, double mdot)
    {
        ReactorBase* r = &ReactorCabinet::item(i);
        if (r->type() >= ReactorType) {
            ((FlowReactor*)r)->setMassFlowRate(mdot);
        }
        return 0;
    }

    size_t reactor_nSensParams(int i)
    {
        ReactorBase* r = &ReactorCabinet::item(i);
        if (r->type() >= ReactorType) {
            return ((Reactor*)r)->nSensParams();
        } else {
            std::cout << "type problem..." << r->type() << std::endl;
            return 0;
        }
    }

    int reactor_addSensitivityReaction(int i, int rxn)
    {
        ReactorBase* r = &ReactorCabinet::item(i);
        ((Reactor*)r)->addSensitivityReaction(rxn);
        return 0;
    }


    // reactor networks

    int reactornet_new()
    {
        ReactorNet* r = new ReactorNet();
        return NetworkCabinet::add(r);
    }

    int reactornet_del(int i)
    {
        try {
            NetworkCabinet::del(i);
            return 0;
        } catch (...) {
            return -1;
        }
    }

    int reactornet_copy(int i)
    {
        return NetworkCabinet::newCopy(i);
    }

    int reactornet_assign(int i, int j)
    {
        return NetworkCabinet::assign(i,j);
    }

    int reactornet_setInitialTime(int i, double t)
    {
        NetworkCabinet::item(i).setInitialTime(t);
        return 0;
    }

    int reactornet_setMaxTimeStep(int i, double maxstep)
    {
        NetworkCabinet::item(i).setMaxTimeStep(maxstep);
        return 0;
    }

    int reactornet_setTolerances(int i, double rtol, double atol)
    {
        NetworkCabinet::item(i).setTolerances(rtol, atol);
        return 0;
    }

    int reactornet_setSensitivityTolerances(int i, double rtol, double atol)
    {
        NetworkCabinet::item(i).setSensitivityTolerances(rtol, atol);
        return 0;
    }

    int reactornet_addreactor(int i, int n)
    {
        try {
            NetworkCabinet::item(i).addReactor(&ReactorCabinet::item(n));
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

    double reactornet_step(int i, double t)
    {
        try {
            return NetworkCabinet::item(i).step(t);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    double reactornet_time(int i)
    {
        return NetworkCabinet::item(i).time();
    }

    double reactornet_rtol(int i)
    {
        return NetworkCabinet::item(i).rtol();
    }

    double reactornet_atol(int i)
    {
        return NetworkCabinet::item(i).atol();
    }

    double reactornet_sensitivity(int i, char* v, int p, int r)
    {
        return NetworkCabinet::item(i).sensitivity(v, p, r);
    }


    // flow devices

    int flowdev_new(int type)
    {
        FlowDevice* r;
        switch (type) {
        case MFC_Type:
            r = new MassFlowController();
            break;
        case PressureController_Type:
            r = new PressureController();
            break;
        case Valve_Type:
            r = new Valve();
            break;
        default:
            r = new FlowDevice();
        }
        return FlowDeviceCabinet::add(r);
    }

    int flowdev_del(int i)
    {
        FlowDeviceCabinet::del(i);
        return 0;
    }

    int flowdev_install(int i, int n, int m)
    {
        try {
            bool ok = FlowDeviceCabinet::item(i).install(ReactorCabinet::item(n),
                                                         ReactorCabinet::item(m));
            if (!ok) {
                throw CanteraError("install","Could not install flow device.");
            }
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int flowdev_setMaster(int i, int n)
    {
        if (FlowDeviceCabinet::item(i).type() == PressureController_Type) {
            dynamic_cast<PressureController&>(FlowDeviceCabinet::item(i)).setMaster(
                &FlowDeviceCabinet::item(n));
        }
        return 0;
    }

    double flowdev_massFlowRate(int i, double time)
    {
        return FlowDeviceCabinet::item(i).massFlowRate(time);
    }

    int flowdev_setMassFlowRate(int i, double mdot)
    {
        FlowDeviceCabinet::item(i).setMassFlowRate(mdot);
        return 0;
    }

    int flowdev_setParameters(int i, int n, double* v)
    {
        FlowDeviceCabinet::item(i).setParameters(n, v);
        return 0;
    }

    int flowdev_setFunction(int i, int n)
    {
        FlowDeviceCabinet::item(i).setFunction(&FuncCabinet::item(n));
        return 0;
    }

    int flowdev_ready(int i)
    {
        bool ok = FlowDeviceCabinet::item(i).ready();
        if (ok) {
            return 1;
        }
        return 0;
    }


    /////////////    Walls   ///////////////////////


    int wall_new(int type)
    {
        Wall* r;
        r = new Wall();
        return WallCabinet::add(r);
    }

    int wall_del(int i)
    {
        WallCabinet::del(i);
        return 0;
    }

    int wall_copy(int i)
    {
        return WallCabinet::newCopy(i);
    }

    int wall_assign(int i, int j)
    {
        return WallCabinet::assign(i,j);
    }

    int wall_install(int i, int n, int m)
    {
        WallCabinet::item(i).install(ReactorCabinet::item(n),
                                     ReactorCabinet::item(m));
        return 0;
    }

    int wall_setkinetics(int i, int n, int m)
    {
        Kinetics* left=0, *right=0;
        if (n > 0)
            if (KineticsCabinet::item(n).type() == cInterfaceKinetics) {
                left = &KineticsCabinet::item(n);
            }
        if (m > 0)
            if (KineticsCabinet::item(m).type() == cInterfaceKinetics) {
                right = &KineticsCabinet::item(m);
            }
        WallCabinet::item(i).setKinetics(left, right);
        return 0;
    }

    double wall_vdot(int i, double t)
    {
        return WallCabinet::item(i).vdot(t);
    }

    double wall_Q(int i, double t)
    {
        return WallCabinet::item(i).Q(t);
    }

    double wall_area(int i)
    {
        return WallCabinet::item(i).area();
    }

    int wall_setArea(int i, double v)
    {
        WallCabinet::item(i).setArea(v);
        return 0;
    }

    int wall_setThermalResistance(int i, double rth)
    {
        WallCabinet::item(i).setThermalResistance(rth);
        return 0;
    }

    int wall_setHeatTransferCoeff(int i, double u)
    {
        WallCabinet::item(i).setHeatTransferCoeff(u);
        return 0;
    }

    int wall_setHeatFlux(int i, int n)
    {
        WallCabinet::item(i).setHeatFlux(&FuncCabinet::item(n));
        return 0;
    }

    int wall_setExpansionRateCoeff(int i, double k)
    {
        WallCabinet::item(i).setExpansionRateCoeff(k);
        return 0;
    }

    int wall_setVelocity(int i, int n)
    {
        WallCabinet::item(i).setVelocity(&FuncCabinet::item(n));
        return 0;
    }

    int wall_setEmissivity(int i, double epsilon)
    {
        WallCabinet::item(i).setEmissivity(epsilon);
        return 0;
    }

    int wall_ready(int i)
    {
        if (WallCabinet::item(i).ready()) {
            return 1;
        } else {
            return 0;
        }
    }

    int wall_addSensitivityReaction(int i, int lr, int rxn)
    {
        WallCabinet::item(i).addSensitivityReaction(lr, rxn);
        return 0;
    }
}
