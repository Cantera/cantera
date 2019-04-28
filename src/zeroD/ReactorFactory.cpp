//! @file ReactorFactory.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/ReactorFactory.h"
#include "cantera/zeroD/Reservoir.h"
#include "cantera/zeroD/Reactor.h"
#include "cantera/zeroD/FlowReactor.h"
#include "cantera/zeroD/ConstPressureReactor.h"
#include "cantera/zeroD/IdealGasReactor.h"
#include "cantera/zeroD/IdealGasConstPressureReactor.h"

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

    // only used by clib
    reg_type("Reservoir", ReservoirType);
    reg_type("Reactor", ReactorType);
    reg_type("ConstPressureReactor", ConstPressureReactorType);
    reg_type("FlowReactor", FlowReactorType);
    reg_type("IdealGasReactor", IdealGasReactorType);
    reg_type("IdealGasConstPressureReactor", IdealGasConstPressureReactorType);
}

ReactorBase* ReactorFactory::newReactor(const std::string& reactorType)
{
    return create(reactorType);
}

ReactorBase* ReactorFactory::newReactor(int ir)
{
    try {
        return create(m_types.at(ir));
    } catch (out_of_range&) {
        throw CanteraError("ReactorFactory::newReactor",
                           "unknown reactor type!");
    }
}

}
