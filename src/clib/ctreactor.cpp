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
#include "cantera/zeroD/ReactorFactory.h"
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

    int reactor_copy(int i)
    {
        try {
            return ReactorCabinet::newCopy(i);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int reactor_assign(int i, int j)
    {
        try {
            return ReactorCabinet::assign(i,j);
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
            // @todo This should not fail silently
            if (ReactorCabinet::item(i).type() >= ReactorType) {
                ReactorCabinet::get<Reactor>(i).setKineticsMgr(KineticsCabinet::item(n));
            }
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

    int reactor_setEnergy(int i, int eflag)
    {
        try {
            // @todo This should not fail silently
            if (ReactorCabinet::item(i).type() >= ReactorType) {
                ReactorCabinet::get<Reactor>(i).setEnergy(eflag);
            }
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

    int reactornet_copy(int i)
    {
        try {
            return NetworkCabinet::newCopy(i);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int reactornet_assign(int i, int j)
    {
        try {
            return NetworkCabinet::assign(i,j);
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

    double reactornet_sensitivity(int i, char* v, int p, int r)
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
        try {
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
                throw CanteraError("install","Could not install flow device.");
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
        try {
            FlowDeviceCabinet::item(i).setMassFlowRate(mdot);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int flowdev_setParameters(int i, int n, double* v)
    {
        try {
            FlowDeviceCabinet::item(i).setParameters(n, v);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int flowdev_setFunction(int i, int n)
    {
        try {
            FlowDeviceCabinet::item(i).setFunction(&FuncCabinet::item(n));
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int flowdev_ready(int i)
    {
        try {
            return int(FlowDeviceCabinet::item(i).ready());
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }


    /////////////    Walls   ///////////////////////


    int wall_new(int type)
    {
        try {
            return WallCabinet::add(new Wall());
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

    int wall_copy(int i)
    {
        try {
            return WallCabinet::newCopy(i);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int wall_assign(int i, int j)
    {
        try {
            return WallCabinet::assign(i,j);
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

    int wall_setkinetics(int i, int n, int m)
    {
        try {
            Kinetics* left=0, *right=0;
            if (n > 0 && KineticsCabinet::item(n).type() == cInterfaceKinetics) {
                left = &KineticsCabinet::item(n);
            }
            if (m > 0 && KineticsCabinet::item(m).type() == cInterfaceKinetics) {
                right = &KineticsCabinet::item(m);
            }
            WallCabinet::item(i).setKinetics(left, right);
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
            WallCabinet::item(i).setThermalResistance(rth);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int wall_setHeatTransferCoeff(int i, double u)
    {
        try {
            WallCabinet::item(i).setHeatTransferCoeff(u);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int wall_setHeatFlux(int i, int n)
    {
        try {
            WallCabinet::item(i).setHeatFlux(&FuncCabinet::item(n));
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int wall_setExpansionRateCoeff(int i, double k)
    {
        try {
            WallCabinet::item(i).setExpansionRateCoeff(k);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int wall_setVelocity(int i, int n)
    {
        try {
            WallCabinet::item(i).setVelocity(&FuncCabinet::item(n));
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int wall_setEmissivity(int i, double epsilon)
    {
        try {
            WallCabinet::item(i).setEmissivity(epsilon);
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

    int wall_addSensitivityReaction(int i, int lr, int rxn)
    {
        try {
            WallCabinet::item(i).addSensitivityReaction(lr, rxn);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }
}
