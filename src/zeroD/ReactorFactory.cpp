//! @file ReactorFactory.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/ReactorFactory.h"
#include "cantera/zeroD/Reservoir.h"
#include "cantera/zeroD/Reactor.h"
#include "cantera/zeroD/MoleReactor.h"
#include "cantera/zeroD/FlowReactor.h"
#include "cantera/zeroD/ConstPressureReactor.h"
#include "cantera/zeroD/ConstPressureMoleReactor.h"
#include "cantera/zeroD/IdealGasReactor.h"
#include "cantera/zeroD/IdealGasMoleReactor.h"
#include "cantera/zeroD/IdealGasConstPressureReactor.h"
#include "cantera/zeroD/ReactorDelegator.h"
#include "cantera/zeroD/IdealGasConstPressureMoleReactor.h"

namespace Cantera
{

ReactorFactory* ReactorFactory::s_factory = 0;
std::mutex ReactorFactory::reactor_mutex;

ReactorFactory::ReactorFactory()
{
    reg("Reservoir",
        [](shared_ptr<Solution> sol, bool clone, const string& name)
        { return new Reservoir(sol, clone, name); });
    reg("Reactor",
        [](shared_ptr<Solution> sol, bool clone, const string& name)
        { return new Reactor(sol, clone, name); });
    reg("ConstPressureReactor",
        [](shared_ptr<Solution> sol, bool clone, const string& name)
        { return new ConstPressureReactor(sol, clone, name); });
    reg("FlowReactor",
        [](shared_ptr<Solution> sol, bool clone, const string& name)
        { return new FlowReactor(sol, clone, name); });
    reg("IdealGasReactor",
        [](shared_ptr<Solution> sol, bool clone, const string& name)
        { return new IdealGasReactor(sol, clone, name); });
    reg("IdealGasConstPressureReactor",
        [](shared_ptr<Solution> sol, bool clone, const string& name)
        { return new IdealGasConstPressureReactor(sol, clone, name); });
    reg("ExtensibleReactor",
        [](shared_ptr<Solution> sol, bool clone, const string& name)
        { return new ReactorDelegator<Reactor>(sol, clone, name); });
    reg("ExtensibleIdealGasReactor",
        [](shared_ptr<Solution> sol, bool clone, const string& name)
        { return new ReactorDelegator<IdealGasReactor>(sol, clone, name); });
    reg("ExtensibleConstPressureReactor",
        [](shared_ptr<Solution> sol, bool clone, const string& name)
        { return new ReactorDelegator<ConstPressureReactor>(sol, clone, name); });
    reg("ExtensibleIdealGasConstPressureReactor",
        [](shared_ptr<Solution> sol, bool clone, const string& name)
        { return new ReactorDelegator<IdealGasConstPressureReactor>(sol, clone, name); });
    reg("ExtensibleMoleReactor",
        [](shared_ptr<Solution> sol, bool clone, const string& name)
        { return new ReactorDelegator<MoleReactor>(sol, clone, name); });
    reg("ExtensibleConstPressureMoleReactor",
        [](shared_ptr<Solution> sol, bool clone, const string& name)
        { return new ReactorDelegator<ConstPressureMoleReactor>(sol, clone, name); });
    reg("ExtensibleIdealGasMoleReactor",
        [](shared_ptr<Solution> sol, bool clone, const string& name)
        { return new ReactorDelegator<IdealGasMoleReactor>(sol, clone, name); });
    reg("ExtensibleIdealGasConstPressureMoleReactor",
        [](shared_ptr<Solution> sol, bool clone, const string& name)
        { return new ReactorDelegator<IdealGasConstPressureMoleReactor>(sol, clone, name); });
    reg("IdealGasConstPressureMoleReactor",
        [](shared_ptr<Solution> sol, bool clone, const string& name)
        { return new IdealGasConstPressureMoleReactor(sol, clone, name); });
    reg("IdealGasMoleReactor",
        [](shared_ptr<Solution> sol, bool clone, const string& name)
        { return new IdealGasMoleReactor(sol, clone, name); });
    reg("ConstPressureMoleReactor",
        [](shared_ptr<Solution> sol, bool clone, const string& name)
        { return new ConstPressureMoleReactor(sol, clone, name); });
    reg("MoleReactor",
        [](shared_ptr<Solution> sol, bool clone, const string& name)
        { return new MoleReactor(sol, clone, name); });
}

ReactorFactory* ReactorFactory::factory() {
    std::unique_lock<std::mutex> lock(reactor_mutex);
    if (!s_factory) {
        s_factory = new ReactorFactory;
    }
    return s_factory;
}

void ReactorFactory::deleteFactory() {
    std::unique_lock<std::mutex> lock(reactor_mutex);
    delete s_factory;
    s_factory = 0;
}

// ---------- ReactorSurfaceFactory methods ----------

ReactorSurfaceFactory* ReactorSurfaceFactory::s_factory = 0;
std::mutex ReactorSurfaceFactory::s_mutex;

ReactorSurfaceFactory* ReactorSurfaceFactory::factory() {
    std::unique_lock<std::mutex> lock(s_mutex);
    if (!s_factory) {
        s_factory = new ReactorSurfaceFactory;
    }
    return s_factory;
}

void ReactorSurfaceFactory::deleteFactory() {
    std::unique_lock<std::mutex> lock(s_mutex);
    delete s_factory;
    s_factory = 0;
}

ReactorSurfaceFactory::ReactorSurfaceFactory()
{
    reg("ReactorSurface",
        [](shared_ptr<Solution> surf, const vector<shared_ptr<ReactorBase>>& reactors,
           bool clone, const string& name)
        { return new ReactorSurface(surf, reactors, clone, name); });
    reg("MoleReactorSurface",
        [](shared_ptr<Solution> surf, const vector<shared_ptr<ReactorBase>>& reactors,
           bool clone, const string& name)
        { return new MoleReactorSurface(surf, reactors, clone, name); });
    reg("FlowReactorSurface",
        [](shared_ptr<Solution> surf, const vector<shared_ptr<ReactorBase>>& reactors,
           bool clone, const string& name)
        { return new FlowReactorSurface(surf, reactors, clone, name); });
    reg("ExtensibleReactorSurface",
        [](shared_ptr<Solution> surf, const vector<shared_ptr<ReactorBase>>& reactors,
           bool clone, const string& name)
        { return new ReactorDelegator<ReactorSurface>(surf, reactors, clone, name); });
    reg("ExtensibleMoleReactorSurface",
        [](shared_ptr<Solution> surf, const vector<shared_ptr<ReactorBase>>& reactors,
           bool clone, const string& name)
        { return new ReactorDelegator<MoleReactorSurface>(surf, reactors, clone, name); });
    reg("ExtensibleFlowReactorSurface",
        [](shared_ptr<Solution> surf, const vector<shared_ptr<ReactorBase>>& reactors,
           bool clone, const string& name)
        { return new ReactorDelegator<FlowReactorSurface>(surf, reactors, clone, name); });
}

// ---------- free functions ----------

shared_ptr<ReactorBase> newReactorBase(
    const string& model, shared_ptr<Solution> phase, bool clone, const string& name)
{
    return shared_ptr<ReactorBase>(
        ReactorFactory::factory()->create(model, phase, clone, name));
}

shared_ptr<Reactor> newReactor(
    const string& model, shared_ptr<Solution> phase, bool clone, const string& name)
{
    auto reactor = std::dynamic_pointer_cast<Reactor>(
        newReactorBase(model, phase, clone, name));
    if (!reactor) {
        throw CanteraError("newReactor4",
            "Model type '{}' does not specify a bulk reactor.", model);
    }
    return reactor;
}

shared_ptr<Reservoir> newReservoir(
    shared_ptr<Solution> phase, bool clone, const string& name)
{
    auto reservoir = std::dynamic_pointer_cast<Reservoir>(
        newReactorBase("Reservoir", phase, clone, name));
    if (!reservoir) {
        throw CanteraError("newReservoir",
            "Caught unexpected inconsistency in factory method.");
    }
    return reservoir;
}

shared_ptr<ReactorSurface> newReactorSurface(shared_ptr<Solution> phase,
    const vector<shared_ptr<ReactorBase>>& reactors, bool clone, const string& name)
{
    if (reactors.empty()) {
        throw CanteraError("newReactorSurface",
            "At least one adjacent reactor must be specified.");
    }
    string model = "ReactorSurface";
    if (std::dynamic_pointer_cast<FlowReactor>(reactors[0])) {
        model = "FlowReactorSurface";
    }
    for (const auto& r : reactors) {
        if (std::dynamic_pointer_cast<MoleReactor>(r)) {
            model = "MoleReactorSurface";
            break;
        }
    }
    return shared_ptr<ReactorSurface>(ReactorSurfaceFactory::factory()->create(
        model, phase, reactors, clone, name));
}

shared_ptr<ReactorSurface> newReactorSurface(const string& model,
    shared_ptr<Solution> phase, const vector<shared_ptr<ReactorBase>>& reactors,
    bool clone, const string& name)
{
    return shared_ptr<ReactorSurface>(ReactorSurfaceFactory::factory()->create(
        model, phase, reactors, clone, name));
}

}
