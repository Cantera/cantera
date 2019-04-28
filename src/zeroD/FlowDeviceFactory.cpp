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

    // only used by clib
    reg_type("MassFlowController", MFC_Type);
    reg_type("PressureController", PressureController_Type);
    reg_type("Valve", Valve_Type);
}

FlowDevice* FlowDeviceFactory::newFlowDevice(const std::string& flowDeviceType)
{
    return create(flowDeviceType);
}

FlowDevice* FlowDeviceFactory::newFlowDevice(int ir)
{
    try {
        return create(m_types.at(ir));
    } catch (out_of_range&) {
        throw CanteraError("FlowDeviceFactory::newFlowDevice",
                           "unknown flowDevice type!");
    }
}

}
