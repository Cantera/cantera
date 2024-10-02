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
        { return Reservoir::create(sol, name); });
    reg("Reactor",
        [](shared_ptr<Solution> sol, const string& name)
        { return Reactor::create(sol, name); });
    reg("ConstPressureReactor",
        [](shared_ptr<Solution> sol, const string& name)
        { return ConstPressureReactor::create(sol, name); });
    reg("FlowReactor",
        [](shared_ptr<Solution> sol, const string& name)
        { return FlowReactor::create(sol, name); });
    reg("IdealGasReactor",
        [](shared_ptr<Solution> sol, const string& name)
        { return IdealGasReactor::create(sol, name); });
    reg("IdealGasConstPressureReactor",
        [](shared_ptr<Solution> sol, const string& name)
        { return IdealGasConstPressureReactor::create(sol, name); });
    reg("ExtensibleReactor",
        [](shared_ptr<Solution> sol, const string& name)
        { return ReactorDelegator<Reactor>::create(sol, name); });
    reg("ExtensibleIdealGasReactor",
        [](shared_ptr<Solution> sol, const string& name)
        { return ReactorDelegator<IdealGasReactor>::create(sol, name); });
    reg("ExtensibleConstPressureReactor",
        [](shared_ptr<Solution> sol, const string& name)
        { return ReactorDelegator<ConstPressureReactor>::create(sol, name); });
    reg("ExtensibleIdealGasConstPressureReactor",
        [](shared_ptr<Solution> sol, const string& name)
        { return ReactorDelegator<IdealGasConstPressureReactor>::create(sol, name); });
    reg("ExtensibleMoleReactor",
        [](shared_ptr<Solution> sol, const string& name)
        { return ReactorDelegator<MoleReactor>::create(sol, name); });
    reg("ExtensibleConstPressureMoleReactor",
        [](shared_ptr<Solution> sol, const string& name)
        { return ReactorDelegator<ConstPressureMoleReactor>::create(sol, name); });
    reg("ExtensibleIdealGasMoleReactor",
        [](shared_ptr<Solution> sol, const string& name)
        { return ReactorDelegator<IdealGasMoleReactor>::create(sol, name); });
    reg("ExtensibleIdealGasConstPressureMoleReactor",
        [](shared_ptr<Solution> sol, const string& name)
        { return
            ReactorDelegator<IdealGasConstPressureMoleReactor>::create(sol, name); });
    reg("IdealGasConstPressureMoleReactor",
        [](shared_ptr<Solution> sol, const string& name)
        { return IdealGasConstPressureMoleReactor::create(sol, name); });
    reg("IdealGasMoleReactor",
        [](shared_ptr<Solution> sol, const string& name)
        { return IdealGasMoleReactor::create(sol, name); });
    reg("ConstPressureMoleReactor",
        [](shared_ptr<Solution> sol, const string& name)
        { return ConstPressureMoleReactor::create(sol, name); });
    reg("MoleReactor",
        [](shared_ptr<Solution> sol, const string& name)
        { return MoleReactor::create(sol, name); });
    reg("ReactorSurface",
        [](shared_ptr<Solution> sol, const string& name)
        { return ReactorSurface::create(sol, name); });

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

shared_ptr<ReactorNode> newReactorNode(
    const string& model, shared_ptr<Solution> contents, const string& name)
{
    return ReactorFactory::factory()->create(model, contents, name);
}

shared_ptr<ReactorNode> newReactorNode(const string& model)
{
    warn_deprecated("newReactorNode",
        "Transitional method to be removed after Cantera 3.1. Use newReactorNode with "
        "contents instead.");
    return ReactorFactory::factory()->create(model, nullptr, "");
}

shared_ptr<ReactorBase> newReactor3(const string& model)
{
    warn_deprecated("newReactor3",
        "Superseded by newReactorNode with contents; to be removed after Cantera 3.1.");
    auto reactor = std::dynamic_pointer_cast<ReactorBase>(
        newReactorNode(model, nullptr, ""));
    if (!reactor) {
        throw CanteraError("newReactor",
            "Detected incompatible ReactorBase type '{}'", model);
    }
    return reactor;
}

}
