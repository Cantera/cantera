//! @file ReactorFactory.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/ReactorFactory.h"
#include "cantera/zeroD/Reservoir.h"
#include "cantera/zeroD/Reactor.h"
#include "cantera/zeroD/FlowReactor.h"
#include "cantera/zeroD/ConstPressureReactor.h"
#include "cantera/zeroD/IdealGasReactor.h"
#include "cantera/zeroD/IdealGasMoleReactor.h"
#include "cantera/zeroD/IdealGasConstPressureReactor.h"
#include "cantera/zeroD/ReactorDelegator.h"
#include "cantera/zeroD/IdealGasConstPressureMoleReactor.h"

using namespace std;
namespace Cantera
{

class Reservoir;

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
    reg("IdealGasConstPressureMoleReactor", []() { return new IdealGasConstPressureMoleReactor(); });
    reg("IdealGasMoleReactor", []() { return new IdealGasMoleReactor(); });
}

ReactorBase* ReactorFactory::newReactor(const std::string& reactorType)
{
    return create(reactorType);
}

}
