//! @file FlowDeviceFactory.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/FlowDeviceFactory.h"
#include "cantera/zeroD/flowControllers.h"

using namespace std;
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

FlowDevice* FlowDeviceFactory::newFlowDevice(const std::string& flowDeviceType)
{
    return create(flowDeviceType);
}

}
