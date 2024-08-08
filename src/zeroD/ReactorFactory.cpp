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
    reg("Reservoir", [](const string& name) { return new Reservoir(name); });
    reg("Reactor", [](const string& name) { return new Reactor(name); });
    reg("ConstPressureReactor", [](const string& name) { return new ConstPressureReactor(name); });
    reg("FlowReactor", [](const string& name) { return new FlowReactor(name); });
    reg("IdealGasReactor", [](const string& name) { return new IdealGasReactor(name); });
    reg("IdealGasConstPressureReactor", [](const string& name) { return new IdealGasConstPressureReactor(name); });
    reg("ExtensibleReactor", [](const string& name) { return new ReactorDelegator<Reactor>(name); });
    reg("ExtensibleIdealGasReactor",
        [](const string& name) { return new ReactorDelegator<IdealGasReactor>(name); });
    reg("ExtensibleConstPressureReactor",
        [](const string& name) { return new ReactorDelegator<ConstPressureReactor>(name); });
    reg("ExtensibleIdealGasConstPressureReactor",
        [](const string& name) { return new ReactorDelegator<IdealGasConstPressureReactor>(name); });
    reg("ExtensibleMoleReactor",
        [](const string& name) { return new ReactorDelegator<MoleReactor>(name); });
    reg("ExtensibleConstPressureMoleReactor",
        [](const string& name) { return new ReactorDelegator<ConstPressureMoleReactor>(name); });
    reg("ExtensibleIdealGasMoleReactor",
        [](const string& name) { return new ReactorDelegator<IdealGasMoleReactor>(name); });
    reg("ExtensibleIdealGasConstPressureMoleReactor",
        [](const string& name) { return new ReactorDelegator<IdealGasConstPressureMoleReactor>(name); });
    reg("IdealGasConstPressureMoleReactor",
        [](const string& name) { return new IdealGasConstPressureMoleReactor(name); });
    reg("IdealGasMoleReactor", [](const string& name) { return new IdealGasMoleReactor(name); });
    reg("ConstPressureMoleReactor", [](const string& name) { return new ConstPressureMoleReactor(name); });
    reg("MoleReactor", [](const string& name) { return new MoleReactor(name); });
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

shared_ptr<ReactorBase> newReactor(const string& model)
{
    warn_deprecated("newReactor",
        "Creation of empty reactor objects is deprecated in Cantera 3.1 and will be \n"
        "removed thereafter; reactor contents should be provided in the constructor.");
    return shared_ptr<ReactorBase>(ReactorFactory::factory()->create(model, ""));
}

shared_ptr<ReactorBase> newReactor(
    const string& model, shared_ptr<Solution> contents, const string& name)
{
    // once empty reactors are no longer supported, the create factory method should
    // support passing a Solution object
    auto ret = shared_ptr<ReactorBase>(ReactorFactory::factory()->create(model, name));
    ret->setSolution(contents);
    return ret;
}

shared_ptr<ReactorBase> newReactor3(const string& model)
{
    warn_deprecated("newReactor3", "To be removed after Cantera 3.1.");
    return newReactor(model);
}

}
