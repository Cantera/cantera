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
    reg("Reservoir", []() { return new Reservoir(); });
    reg("Reactor", []() { return new Reactor(); });
    reg("ConstPressureReactor", []() { return new ConstPressureReactor(); });
    reg("FlowReactor", []() { return new FlowReactor(); });
    reg("IdealGasReactor", []() { return new IdealGasReactor(); });
    reg("IdealGasConstPressureReactor", []() { return new IdealGasConstPressureReactor(); });
    reg("ExtensibleReactor", []() { return new ReactorDelegator<Reactor>(); });
    reg("ExtensibleIdealGasReactor",
        []() { return new ReactorDelegator<IdealGasReactor>(); });
    reg("ExtensibleConstPressureReactor",
        []() { return new ReactorDelegator<ConstPressureReactor>(); });
    reg("ExtensibleIdealGasConstPressureReactor",
        []() { return new ReactorDelegator<IdealGasConstPressureReactor>(); });
    reg("ExtensibleMoleReactor",
        []() { return new ReactorDelegator<MoleReactor>(); });
    reg("ExtensibleConstPressureMoleReactor",
        []() { return new ReactorDelegator<ConstPressureMoleReactor>(); });
    reg("ExtensibleIdealGasMoleReactor",
        []() { return new ReactorDelegator<IdealGasMoleReactor>(); });
    reg("ExtensibleIdealGasConstPressureMoleReactor",
        []() { return new ReactorDelegator<IdealGasConstPressureMoleReactor>(); });
    reg("IdealGasConstPressureMoleReactor", []() { return new
        IdealGasConstPressureMoleReactor(); });
    reg("IdealGasMoleReactor", []() { return new IdealGasMoleReactor(); });
    reg("ConstPressureMoleReactor", []() { return new ConstPressureMoleReactor(); });
    reg("MoleReactor", []() { return new MoleReactor(); });
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
    return shared_ptr<ReactorBase>(ReactorFactory::factory()->create(model));
}

shared_ptr<ReactorBase> newReactor(
    const string& model, shared_ptr<Solution> contents, const string& name)
{
    // once empty reactors are no longer supported, the create factory method should
    // support passing a Solution object and a name
    auto ret = shared_ptr<ReactorBase>(ReactorFactory::factory()->create(model));
    ret->setSolution(contents);
    ret->setName(name);
    return ret;
}

shared_ptr<ReactorBase> newReactor3(const string& model)
{
    warn_deprecated("newReactor3", "To be removed after Cantera 3.1.");
    return newReactor(model);
}

}
