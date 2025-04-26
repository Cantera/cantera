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
        [](shared_ptr<Solution> sol, const string& name)
        { return new Reservoir(sol, name); });
    reg("Reactor",
        [](shared_ptr<Solution> sol, const string& name)
        { return new Reactor(sol, name); });
    reg("ConstPressureReactor",
        [](shared_ptr<Solution> sol, const string& name)
        { return new ConstPressureReactor(sol, name); });
    reg("FlowReactor",
        [](shared_ptr<Solution> sol, const string& name)
        { return new FlowReactor(sol, name); });
    reg("IdealGasReactor",
        [](shared_ptr<Solution> sol, const string& name)
        { return new IdealGasReactor(sol, name); });
    reg("IdealGasConstPressureReactor",
        [](shared_ptr<Solution> sol, const string& name)
        { return new IdealGasConstPressureReactor(sol, name); });
    reg("ExtensibleReactor",
        [](shared_ptr<Solution> sol, const string& name)
        { return new ReactorDelegator<Reactor>(sol, name); });
    reg("ExtensibleIdealGasReactor",
        [](shared_ptr<Solution> sol, const string& name)
        { return new ReactorDelegator<IdealGasReactor>(sol, name); });
    reg("ExtensibleConstPressureReactor",
        [](shared_ptr<Solution> sol, const string& name)
        { return new ReactorDelegator<ConstPressureReactor>(sol, name); });
    reg("ExtensibleIdealGasConstPressureReactor",
        [](shared_ptr<Solution> sol, const string& name)
        { return new ReactorDelegator<IdealGasConstPressureReactor>(sol, name); });
    reg("ExtensibleMoleReactor",
        [](shared_ptr<Solution> sol, const string& name)
        { return new ReactorDelegator<MoleReactor>(sol, name); });
    reg("ExtensibleConstPressureMoleReactor",
        [](shared_ptr<Solution> sol, const string& name)
        { return new ReactorDelegator<ConstPressureMoleReactor>(sol, name); });
    reg("ExtensibleIdealGasMoleReactor",
        [](shared_ptr<Solution> sol, const string& name)
        { return new ReactorDelegator<IdealGasMoleReactor>(sol, name); });
    reg("ExtensibleIdealGasConstPressureMoleReactor",
        [](shared_ptr<Solution> sol, const string& name)
        { return new ReactorDelegator<IdealGasConstPressureMoleReactor>(sol, name); });
    reg("IdealGasConstPressureMoleReactor",
        [](shared_ptr<Solution> sol, const string& name)
        { return new IdealGasConstPressureMoleReactor(sol, name); });
    reg("IdealGasMoleReactor",
        [](shared_ptr<Solution> sol, const string& name)
        { return new IdealGasMoleReactor(sol, name); });
    reg("ConstPressureMoleReactor",
        [](shared_ptr<Solution> sol, const string& name)
        { return new ConstPressureMoleReactor(sol, name); });
    reg("MoleReactor",
        [](shared_ptr<Solution> sol, const string& name)
        { return new MoleReactor(sol, name); });
    reg("ReactorSurface",
        [](shared_ptr<Solution> sol, const string& name)
        { return new ReactorSurface(sol, name); });
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

shared_ptr<ReactorBase> newReactorBase(
    const string& model, shared_ptr<Solution> contents, const string& name)
{
    return shared_ptr<ReactorBase>(
        ReactorFactory::factory()->create(model, contents, name));
}

shared_ptr<ReactorBase> newReactor(
    const string& model, shared_ptr<Solution> contents, const string& name)
{
    warn_deprecated("newReactor", "Behavior changes after Cantera 3.2, when a "
        "'shared_ptr<Reactor>' is returned.\nFor new behavior, use 'newReactor4'.");
    return newReactorBase(model, contents, name);
}

shared_ptr<Reactor> newReactor4(
    const string& model, shared_ptr<Solution> contents, const string& name)
{
    auto reactor = std::dynamic_pointer_cast<Reactor>(
        newReactorBase(model, contents, name));
    if (!reactor) {
        throw CanteraError("newReactor4",
            "Model type '{}' does not specify a bulk reactor.", model);
    }
    return reactor;
}

shared_ptr<Reservoir> newReservoir(shared_ptr<Solution> contents, const string& name)
{
    auto reservoir = std::dynamic_pointer_cast<Reservoir>(
        newReactorBase("Reservoir", contents, name));
    if (!reservoir) {
        throw CanteraError("newReservoir",
            "Caught unexpected inconsistency in factory method.");
    }
    return reservoir;
}

}
