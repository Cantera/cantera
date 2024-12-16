//! @file FlowDeviceFactory.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/FlowDeviceFactory.h"
#include "cantera/zeroD/flowControllers.h"

namespace Cantera
{

FlowDeviceFactory* FlowDeviceFactory::s_factory = 0;
std::mutex FlowDeviceFactory::flowDevice_mutex;

FlowDeviceFactory::FlowDeviceFactory()
{
    reg("MassFlowController", [](const string& name) {
        return new MassFlowController(name);
    });
    reg("PressureController", [](const string& name) {
        return new PressureController(name);
    });
    reg("Valve", [](const string& name) {
        return new Valve(name);
    });
}

FlowDeviceFactory* FlowDeviceFactory::factory() {
    std::unique_lock<std::mutex> lock(flowDevice_mutex);
    if (!s_factory) {
        s_factory = new FlowDeviceFactory;
    }
    return s_factory;
}

void FlowDeviceFactory::deleteFactory() {
    std::unique_lock<std::mutex> lock(flowDevice_mutex);
    delete s_factory;
    s_factory = 0;
}

shared_ptr<FlowDevice> newFlowDevice(const string& model, const string& name)
{
    return shared_ptr<FlowDevice>(FlowDeviceFactory::factory()->create(model, name));
}

}
