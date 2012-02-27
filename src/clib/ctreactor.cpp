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
#include "Storage.h"

using namespace Cantera;
using namespace std;

typedef Cabinet<ReactorBase> ReactorCabinet;
typedef Cabinet<ReactorNet> NetworkCabinet;
typedef Cabinet<FlowDevice> FlowDeviceCabinet;
typedef Cabinet<Wall> WallCabinet;
typedef Cabinet<Func1> FuncCabinet;

inline Kinetics* _kin(int n)
{
    return Storage::__storage->__ktable[n];
}

inline ThermoPhase* _th(int n)
{
    return Storage::__storage->__thtable[n];
}

template<> ReactorCabinet* ReactorCabinet::__storage = 0;
template<> NetworkCabinet* NetworkCabinet::__storage = 0;
template<> FlowDeviceCabinet* FlowDeviceCabinet::__storage = 0;
template<> WallCabinet* WallCabinet::__storage = 0;

extern "C" {

    // reactor

    int DLL_EXPORT reactor_new(int type)
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

    int DLL_EXPORT reactor_del(int i)
    {
        ReactorCabinet::del(i);
        return 0;
    }

    int DLL_EXPORT reactor_copy(int i)
    {
        return ReactorCabinet::newCopy(i);
    }

    int DLL_EXPORT reactor_assign(int i, int j)
    {
        return ReactorCabinet::assign(i,j);
    }

    int DLL_EXPORT reactor_setInitialVolume(int i, double v)
    {
        ReactorCabinet::item(i).setInitialVolume(v);
        return 0;
    }

    int DLL_EXPORT reactor_setInitialTime(int i, double t)
    {
        ReactorCabinet::item(i).setInitialTime(t);
        return 0;
    }

    int DLL_EXPORT reactor_setThermoMgr(int i, int n)
    {
        ReactorCabinet::item(i).setThermoMgr(*_th(n));
        return 0;
    }

    int DLL_EXPORT reactor_setKineticsMgr(int i, int n)
    {
        ReactorBase* r = &ReactorCabinet::item(i);
        if (r->type() >= ReactorType) {
            ((Reactor*)r)->setKineticsMgr(*_kin(n));
        }
        return 0;
    }

    int DLL_EXPORT reactor_advance(int i, double t)
    {
        try {
            ReactorCabinet::item(i).advance(t);
            return 0;
        } catch (CanteraError) {
            return -1;
        }
    }

    double DLL_EXPORT reactor_step(int i, double t)
    {
        return ReactorCabinet::item(i).step(t);
    }

    double DLL_EXPORT reactor_time(int i)
    {
        return ReactorCabinet::item(i).time();
    }

    double DLL_EXPORT reactor_mass(int i)
    {
        return ReactorCabinet::item(i).mass();
    }

    double DLL_EXPORT reactor_volume(int i)
    {
        return ReactorCabinet::item(i).volume();
    }

    double DLL_EXPORT reactor_density(int i)
    {
        return ReactorCabinet::item(i).density();
    }

    double DLL_EXPORT reactor_temperature(int i)
    {
        return ReactorCabinet::item(i).temperature();
    }

    double DLL_EXPORT reactor_enthalpy_mass(int i)
    {
        return ReactorCabinet::item(i).enthalpy_mass();
    }

    double DLL_EXPORT reactor_intEnergy_mass(int i)
    {
        return ReactorCabinet::item(i).intEnergy_mass();
    }

    double DLL_EXPORT reactor_pressure(int i)
    {
        return ReactorCabinet::item(i).pressure();
    }

    double DLL_EXPORT reactor_massFraction(int i, int k)
    {
        return ReactorCabinet::item(i).massFraction(k);
    }

    int DLL_EXPORT reactor_setEnergy(int i, int eflag)
    {
        ReactorBase* r = &ReactorCabinet::item(i);
        if (r->type() >= ReactorType) {
            ((Reactor*)r)->setEnergy(eflag);
        }
        return 0;
    }

    int DLL_EXPORT flowReactor_setMassFlowRate(int i, double mdot)
    {
        ReactorBase* r = &ReactorCabinet::item(i);
        if (r->type() >= ReactorType) {
            ((FlowReactor*)r)->setMassFlowRate(mdot);
        }
        return 0;
    }

    size_t DLL_EXPORT reactor_nSensParams(int i)
    {
        ReactorBase* r = &ReactorCabinet::item(i);
        if (r->type() >= ReactorType) {
            return ((Reactor*)r)->nSensParams();
        } else {
            std::cout << "type problem..." << r->type() << std::endl;
            return 0;
        }
    }

    int DLL_EXPORT reactor_addSensitivityReaction(int i, int rxn)
    {
        ReactorBase* r = &ReactorCabinet::item(i);
        ((Reactor*)r)->addSensitivityReaction(rxn);
        return 0;
    }


    // reactor networks

    int DLL_EXPORT reactornet_new()
    {
        ReactorNet* r = new ReactorNet();
        return NetworkCabinet::add(r);
    }

    int DLL_EXPORT reactornet_del(int i)
    {
        try {
            NetworkCabinet::del(i);
            return 0;
        } catch (...) {
            return -1;
        }
    }

    int DLL_EXPORT reactornet_copy(int i)
    {
        return NetworkCabinet::newCopy(i);
    }

    int DLL_EXPORT reactornet_assign(int i, int j)
    {
        return NetworkCabinet::assign(i,j);
    }

    int DLL_EXPORT reactornet_setInitialTime(int i, double t)
    {
        NetworkCabinet::item(i).setInitialTime(t);
        return 0;
    }

    int DLL_EXPORT reactornet_setMaxTimeStep(int i, double maxstep)
    {
        NetworkCabinet::item(i).setMaxTimeStep(maxstep);
        return 0;
    }

    int DLL_EXPORT reactornet_setTolerances(int i, double rtol, double atol)
    {
        NetworkCabinet::item(i).setTolerances(rtol, atol);
        return 0;
    }

    int DLL_EXPORT reactornet_setSensitivityTolerances(int i, double rtol, double atol)
    {
        NetworkCabinet::item(i).setSensitivityTolerances(rtol, atol);
        return 0;
    }

    int DLL_EXPORT reactornet_addreactor(int i, int n)
    {
        try {
            NetworkCabinet::item(i).addReactor(&ReactorCabinet::item(n));
            return 0;
        } catch (CanteraError) {
            return -1;
        }
    }

    int DLL_EXPORT reactornet_advance(int i, double t)
    {
        try {
            NetworkCabinet::item(i).advance(t);
            return 0;
        } catch (...) {
            return -1;
        }
    }

    double DLL_EXPORT reactornet_step(int i, double t)
    {
        try {
            return NetworkCabinet::item(i).step(t);
        } catch (...) {
            return DERR;
        }
    }

    double DLL_EXPORT reactornet_time(int i)
    {
        return NetworkCabinet::item(i).time();
    }

    double DLL_EXPORT reactornet_rtol(int i)
    {
        return NetworkCabinet::item(i).rtol();
    }

    double DLL_EXPORT reactornet_atol(int i)
    {
        return NetworkCabinet::item(i).atol();
    }

    double DLL_EXPORT reactornet_sensitivity(int i, char* v, int p, int r)
    {
        return NetworkCabinet::item(i).sensitivity(v, p, r);
    }


    // flow devices

    int DLL_EXPORT flowdev_new(int type)
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

    int DLL_EXPORT flowdev_del(int i)
    {
        FlowDeviceCabinet::del(i);
        return 0;
    }

    int DLL_EXPORT flowdev_install(int i, int n, int m)
    {
        try {
            bool ok = FlowDeviceCabinet::item(i).install(ReactorCabinet::item(n),
                                                         ReactorCabinet::item(m));
            if (!ok) {
                throw CanteraError("install","Could not install flow device.");
            }
            return 0;
        } catch (CanteraError) {
            return -1;
        }
    }

    int DLL_EXPORT flowdev_setMaster(int i, int n)
    {
        if (FlowDeviceCabinet::item(i).type() == PressureController_Type) {
            dynamic_cast<PressureController&>(FlowDeviceCabinet::item(i)).setMaster(
                &FlowDeviceCabinet::item(n));
        }
        return 0;
    }

    double DLL_EXPORT flowdev_massFlowRate(int i, double time)
    {
        return FlowDeviceCabinet::item(i).massFlowRate(time);
    }

    int DLL_EXPORT flowdev_setMassFlowRate(int i, double mdot)
    {
        FlowDeviceCabinet::item(i).setMassFlowRate(mdot);
        return 0;
    }

    int DLL_EXPORT flowdev_setParameters(int i, int n, double* v)
    {
        FlowDeviceCabinet::item(i).setParameters(n, v);
        return 0;
    }

    int DLL_EXPORT flowdev_setFunction(int i, int n)
    {
        FlowDeviceCabinet::item(i).setFunction(&FuncCabinet::item(n));
        return 0;
    }

    int DLL_EXPORT flowdev_ready(int i)
    {
        bool ok = FlowDeviceCabinet::item(i).ready();
        if (ok) {
            return 1;
        }
        return 0;
    }


    /////////////    Walls   ///////////////////////


    int DLL_EXPORT wall_new(int type)
    {
        Wall* r;
        r = new Wall();
        return WallCabinet::add(r);
    }

    int DLL_EXPORT wall_del(int i)
    {
        WallCabinet::del(i);
        return 0;
    }

    int DLL_EXPORT wall_copy(int i)
    {
        return WallCabinet::newCopy(i);
    }

    int DLL_EXPORT wall_assign(int i, int j)
    {
        return WallCabinet::assign(i,j);
    }

    int DLL_EXPORT wall_install(int i, int n, int m)
    {
        WallCabinet::item(i).install(ReactorCabinet::item(n),
                                     ReactorCabinet::item(m));
        return 0;
    }

    int DLL_EXPORT wall_setkinetics(int i, int n, int m)
    {
        Kinetics* left=0, *right=0;
        if (n > 0)
            if (_kin(n)->type() == cInterfaceKinetics) {
                left = _kin(n);
            }
        if (m > 0)
            if (_kin(m)->type() == cInterfaceKinetics) {
                right = _kin(m);
            }
        WallCabinet::item(i).setKinetics(left, right);
        return 0;
    }

    double DLL_EXPORT wall_vdot(int i, double t)
    {
        return WallCabinet::item(i).vdot(t);
    }

    double DLL_EXPORT wall_Q(int i, double t)
    {
        return WallCabinet::item(i).Q(t);
    }

    double DLL_EXPORT wall_area(int i)
    {
        return WallCabinet::item(i).area();
    }

    int DLL_EXPORT wall_setArea(int i, double v)
    {
        WallCabinet::item(i).setArea(v);
        return 0;
    }

    int DLL_EXPORT wall_setThermalResistance(int i, double rth)
    {
        WallCabinet::item(i).setThermalResistance(rth);
        return 0;
    }

    int DLL_EXPORT wall_setHeatTransferCoeff(int i, double u)
    {
        WallCabinet::item(i).setHeatTransferCoeff(u);
        return 0;
    }

    int DLL_EXPORT wall_setHeatFlux(int i, int n)
    {
        WallCabinet::item(i).setHeatFlux(&FuncCabinet::item(n));
        return 0;
    }

    int DLL_EXPORT wall_setExpansionRateCoeff(int i, double k)
    {
        WallCabinet::item(i).setExpansionRateCoeff(k);
        return 0;
    }

    int DLL_EXPORT wall_setVelocity(int i, int n)
    {
        WallCabinet::item(i).setVelocity(&FuncCabinet::item(n));
        return 0;
    }

    int DLL_EXPORT wall_setEmissivity(int i, double epsilon)
    {
        WallCabinet::item(i).setEmissivity(epsilon);
        return 0;
    }

    int DLL_EXPORT wall_ready(int i)
    {
        if (WallCabinet::item(i).ready()) {
            return 1;
        } else {
            return 0;
        }
    }

    int DLL_EXPORT wall_addSensitivityReaction(int i, int lr, int rxn)
    {
        WallCabinet::item(i).addSensitivityReaction(lr, rxn);
        return 0;
    }
}
