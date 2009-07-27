/**
 * @file ctreactor.cpp
 */
/*
 *      $Id: ctreactor.cpp,v 1.19 2009/07/11 17:16:09 hkmoffa Exp $
 */


#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#define CANTERA_USE_INTERNAL
#include "ctreactor.h"

// Cantera includes
#include "Reactor.h"
#include "FlowReactor.h"
#include "ConstPressureReactor.h"
#include "ReactorNet.h"
#include "Reservoir.h"
#include "Wall.h"
#include "flowControllers.h"

#include "Cabinet.h"
#include "Storage.h"

using namespace CanteraZeroD;

typedef ReactorBase reactor_t;
typedef ReactorNet  reactornet_t;
typedef FlowDevice  flowdev_t;
typedef Wall        wall_t;

template<> Cabinet<reactor_t>*    Cabinet<reactor_t>::__storage = 0;
template<> Cabinet<reactornet_t>*    Cabinet<reactornet_t>::__storage = 0;
template<> Cabinet<flowdev_t>*    Cabinet<flowdev_t>::__storage = 0;
template<> Cabinet<wall_t>*       Cabinet<wall_t>::__storage = 0;

inline reactor_t* _reactor(int i) {
    return Cabinet<reactor_t>::cabinet()->item(i);
}

inline reactornet_t* _reactornet(int i) {
    return Cabinet<reactornet_t>::cabinet()->item(i);
}

inline flowdev_t* _flowdev(int i) {
    return Cabinet<flowdev_t>::cabinet()->item(i);
}

inline wall_t* _wall(int i) {
    return Cabinet<wall_t>::cabinet()->item(i);
}

inline Kinetics* _kin(int n) {
    return Storage::__storage->__ktable[n];
}

inline ThermoPhase* _th(int n) {
    return Storage::__storage->__thtable[n];
}


inline Func1* _func(int i) {
    return Cabinet<Func1>::cabinet()->item(i);
} 

extern "C" {  


    // reactor

    int DLL_EXPORT reactor_new(int type) {
        reactor_t* r=0;
        if (type == ReactorType)
            r = new Reactor();
        else if (type == FlowReactorType)
            r = new FlowReactor();
        else if (type == ConstPressureReactorType)
            r = new ConstPressureReactor();
        else if (type == ReservoirType)
            r = new Reservoir();
        else 
            r = new ReactorBase();
        return Cabinet<reactor_t>::cabinet()->add(r);
    }

    int DLL_EXPORT reactor_del(int i) {
        Cabinet<reactor_t>::cabinet()->del(i);
        return 0;
    }

    int DLL_EXPORT reactor_copy(int i) {
        return Cabinet<reactor_t>::cabinet()->newCopy(i);
    }

    int DLL_EXPORT reactor_assign(int i, int j) {
        return Cabinet<reactor_t>::cabinet()->assign(i,j);
    }

    int DLL_EXPORT reactor_setInitialVolume(int i, double v) {
        _reactor(i)->setInitialVolume(v);
        return 0;
    }

    int DLL_EXPORT reactor_setInitialTime(int i, double t) {
       _reactor(i)->setInitialTime(t);
       return 0;
    }

    int DLL_EXPORT reactor_setThermoMgr(int i, int n) {
        _reactor(i)->setThermoMgr(*_th(n));
        return 0;
    }

    int DLL_EXPORT reactor_setKineticsMgr(int i, int n) {
        reactor_t* r = _reactor(i);
        if (r->type() >= ReactorType) 
            ((Reactor*)r)->setKineticsMgr(*_kin(n));
        return 0;
    }

    int DLL_EXPORT reactor_advance(int i, double t) {
       try {
           _reactor(i)->advance(t);
           return 0;
       }
       catch (CanteraError) {return -1;}
    }

       double DLL_EXPORT reactor_step(int i, double t) {
       return _reactor(i)->step(t);
    }

    double DLL_EXPORT reactor_time(int i) {
        return _reactor(i)->time();
    }

    double DLL_EXPORT reactor_mass(int i) {
        return _reactor(i)->mass();
    }

    double DLL_EXPORT reactor_volume(int i) {
        return _reactor(i)->volume();
    }

    double DLL_EXPORT reactor_density(int i) {
        return _reactor(i)->density();
    }

    double DLL_EXPORT reactor_temperature(int i) {
        return _reactor(i)->temperature();
    }

    double DLL_EXPORT reactor_enthalpy_mass(int i) {
        return _reactor(i)->enthalpy_mass();
    }

    double DLL_EXPORT reactor_intEnergy_mass(int i) {
        return _reactor(i)->intEnergy_mass();
    }

    double DLL_EXPORT reactor_pressure(int i) {
        return _reactor(i)->pressure();
    }

    double DLL_EXPORT reactor_massFraction(int i, int k) {
        return _reactor(i)->massFraction(k);
    }

    int DLL_EXPORT reactor_setEnergy(int i, int eflag) {
        reactor_t* r = _reactor(i);
        if (r->type() >= ReactorType) ((Reactor*)r)->setEnergy(eflag);
        return 0;
    }

    int DLL_EXPORT flowReactor_setMassFlowRate(int i, double mdot) {
        reactor_t* r = _reactor(i);
        if (r->type() >= ReactorType) ((FlowReactor*)r)->setMassFlowRate(mdot);
        return 0;
    }

    int DLL_EXPORT reactor_nSensParams(int i) {
        reactor_t* r = _reactor(i);
        if (r->type() >= ReactorType) 
            return ((Reactor*)r)->nSensParams();
        else {
            cout << "type problem..." << r->type() << endl;
            return 0;
        }
    }
 
    int DLL_EXPORT reactor_addSensitivityReaction(int i, int rxn) {
        reactor_t* r = _reactor(i);
        ((Reactor*)r)->addSensitivityReaction(rxn);
        return 0;
    }


    // reactor networks

    int DLL_EXPORT reactornet_new() {
        ReactorNet* r = new ReactorNet();
        return Cabinet<reactornet_t>::cabinet()->add(r);
    }

    int DLL_EXPORT reactornet_del(int i) {
        try {
            Cabinet<reactornet_t>::cabinet()->del(i);
            return 0;
        }
        catch (...) {
            return -1;
        }
    }

    int DLL_EXPORT reactornet_copy(int i) {
        return Cabinet<reactornet_t>::cabinet()->newCopy(i);
    }

    int DLL_EXPORT reactornet_assign(int i, int j) {
        return Cabinet<reactornet_t>::cabinet()->assign(i,j);
    }

    int DLL_EXPORT reactornet_setInitialTime(int i, double t) {
        _reactornet(i)->setInitialTime(t);
        return 0;
    }

    int DLL_EXPORT reactornet_setMaxTimeStep(int i, double maxstep) {
        _reactornet(i)->setMaxTimeStep(maxstep);
        return 0;
    }

    int DLL_EXPORT reactornet_setTolerances(int i, double rtol, double atol) {
        _reactornet(i)->setTolerances(rtol, atol);
        return 0;
    }

    int DLL_EXPORT reactornet_setSensitivityTolerances(int i, double rtol, double atol) {
        _reactornet(i)->setSensitivityTolerances(rtol, atol);
        return 0;
    }

    int DLL_EXPORT reactornet_addreactor(int i, int n) {
        try {
            _reactornet(i)->addReactor(_reactor(n));
            return 0;
        }
        catch (CanteraError) {
            return -1;
        }
    }

    int DLL_EXPORT reactornet_advance(int i, double t) {
        try {
            _reactornet(i)->advance(t);
            return 0;
        }
        catch (...) {return -1;}
    }

    double DLL_EXPORT reactornet_step(int i, double t) {
        try {
            return _reactornet(i)->step(t);
        }
        catch (...) {
            return DERR;
        }
    }

    double DLL_EXPORT reactornet_time(int i) {
        return _reactornet(i)->time();
    }

    double DLL_EXPORT reactornet_rtol(int i) {
        return _reactornet(i)->rtol();
    }

    double DLL_EXPORT reactornet_atol(int i) {
        return _reactornet(i)->atol();
    }

    double DLL_EXPORT reactornet_sensitivity(int i, char* v, int p, int r) {
        return _reactornet(i)->sensitivity(v, p, r);
    }


    // flow devices

    int DLL_EXPORT flowdev_new(int type) {
        flowdev_t* r;
        switch (type) {
        case MFC_Type:
            r = new MassFlowController(); break;
        case PressureController_Type:
            r = new PressureController(); break;
        case Valve_Type:
            r = new Valve(); break;
        default:
            r = new FlowDevice();
        }
        return Cabinet<flowdev_t>::cabinet()->add(r);
    }

    int DLL_EXPORT flowdev_del(int i) {
        Cabinet<flowdev_t>::cabinet()->del(i);
        return 0;
    }

    int DLL_EXPORT flowdev_install(int i, int n, int m) {
        try {
            bool ok = _flowdev(i)->install(*_reactor(n), *_reactor(m) );
            if (!ok) throw CanteraError("install","Could not install flow device.");
            return 0;
        }
        catch (CanteraError) {
            return -1;
        }
    }

    int DLL_EXPORT flowdev_setMaster(int i, int n) {
        if (_flowdev(i)->type() == PressureController_Type) {
            ((PressureController*)_flowdev(i))->setMaster(_flowdev(n));
        }
        return 0;
    }

    double DLL_EXPORT flowdev_massFlowRate(int i, double time) {
        return _flowdev(i)->massFlowRate(time);
    }

    int DLL_EXPORT flowdev_setMassFlowRate(int i, double mdot) {
        _flowdev(i)->setMassFlowRate(mdot);
        return 0;
    }

    int DLL_EXPORT flowdev_setParameters(int i, int n, double* v) {
        _flowdev(i)->setParameters(n, v);
        return 0;
    }

    int DLL_EXPORT flowdev_setFunction(int i, int n) {
        _flowdev(i)->setFunction(_func(n));
        return 0;
    }

    int DLL_EXPORT flowdev_ready(int i) {
        bool ok = _flowdev(i)->ready();
        if (ok) return 1;
        return 0;
    }


    /////////////    Walls   ///////////////////////


    int DLL_EXPORT wall_new(int type) {
        wall_t* r;
        r = new Wall();
        return Cabinet<wall_t>::cabinet()->add(r);
    }

    int DLL_EXPORT wall_del(int i) {
        Cabinet<wall_t>::cabinet()->del(i);
        return 0;
    }

    int DLL_EXPORT wall_copy(int i) {
        return Cabinet<wall_t>::cabinet()->newCopy(i);
    }

    int DLL_EXPORT wall_assign(int i, int j) {
        return Cabinet<wall_t>::cabinet()->assign(i,j);
    }

    int DLL_EXPORT wall_install(int i, int n, int m) {
        _wall(i)->install(*_reactor(n), *_reactor(m) );
        return 0;
    }

    int DLL_EXPORT wall_setkinetics(int i, int n, int m) {
        Kinetics *left=0, *right=0;
        if (n > 0) 
            if (_kin(n)->type() == cInterfaceKinetics)
                left = _kin(n);
        if (m > 0) 
            if (_kin(m)->type() == cInterfaceKinetics)
                right = _kin(m);
        _wall(i)->setKinetics(left, right);
        return 0;
    }

    double DLL_EXPORT wall_vdot(int i, double t) {
        return _wall(i)->vdot(t);
    }

    double DLL_EXPORT wall_Q(int i, double t) {
        return _wall(i)->Q(t);
    }

    double DLL_EXPORT wall_area(int i) {
        return _wall(i)->area();
    }

    int DLL_EXPORT wall_setArea(int i, double v) {
        _wall(i)->setArea(v);
        return 0;
    }

    int DLL_EXPORT wall_setThermalResistance(int i, double rth) {
        _wall(i)->setThermalResistance(rth);
        return 0;
    }

    int DLL_EXPORT wall_setHeatTransferCoeff(int i, double u) {
        _wall(i)->setHeatTransferCoeff(u);
        return 0;
    }

    int DLL_EXPORT wall_setHeatFlux(int i, int n) {
        _wall(i)->setHeatFlux(_func(n));
        return 0;
    }

    int DLL_EXPORT wall_setExpansionRateCoeff(int i, double k) {
        _wall(i)->setExpansionRateCoeff(k);
        return 0;
    }

    int DLL_EXPORT wall_setVelocity(int i, int n) {
        _wall(i)->setVelocity(_func(n));
        return 0;
    }

    int DLL_EXPORT wall_setEmissivity(int i, double epsilon) {
        _wall(i)->setEmissivity(epsilon);
        return 0;
    }

    int DLL_EXPORT wall_ready(int i) {
        if (_wall(i)->ready()) return 1;
        else return 0;
    }

    int DLL_EXPORT wall_addSensitivityReaction(int i, int lr, int rxn) {
        _wall(i)->addSensitivityReaction(lr, rxn);
        return 0;
    }


}
