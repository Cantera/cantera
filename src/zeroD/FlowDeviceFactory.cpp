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
    reg("MassFlowController", []() { return new MassFlowController(); });
    reg("PressureController", []() { return new PressureController(); });
    reg("Valve", []() { return new Valve(); });
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


FlowDevice* FlowDeviceFactory::newFlowDevice(const std::string& flowDeviceType)
{
    return create(flowDeviceType);
}

FlowDevice* newFlowDevice(const string& model)
{
    return FlowDeviceFactory::factory()->newFlowDevice(model);
}

}
