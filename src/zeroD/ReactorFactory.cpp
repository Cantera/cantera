//! @file ReactorFactory.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

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
}

ReactorBase* ReactorFactory::newReactor(const std::string& reactorType)
{
    return create(reactorType);
}


ReactorBase* ReactorFactory::newReactor(int ir)
{
    static const unordered_map<int, string> types {
        {ReservoirType, "Reservoir"},
        {ReactorType, "Reactor"},
        {ConstPressureReactorType, "ConstPressureReactor"},
        {FlowReactorType, "FlowReactor"},
        {IdealGasReactorType, "IdealGasReactor"},
        {IdealGasConstPressureReactorType, "IdealGasConstPressureReactor"}
    };

    try {
        return create(types.at(ir));
    } catch (out_of_range&) {
        throw CanteraError("ReactorFactory::newReactor",
                           "unknown reactor type!");
    }
}

}
