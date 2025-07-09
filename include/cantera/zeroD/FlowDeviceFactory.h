//! @file FlowDeviceFactory.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef FLOWDEVICE_FACTORY_H
#define FLOWDEVICE_FACTORY_H

#include "cantera/base/FactoryBase.h"
#include "cantera/zeroD/FlowDevice.h"

namespace Cantera
{

//! Factory class to create FlowDevice objects.
//!
//! This class is mainly used via the newFlowDevice() function, for example:
//!
//! ```cpp
//!     shared_ptr<FlowDevice> mfc = newFlowDevice("MassFlowController");
//! ```
class FlowDeviceFactory : public Factory<FlowDevice, const string&>
{
public:
    static FlowDeviceFactory* factory();

    void deleteFactory() override;

private:
    static FlowDeviceFactory* s_factory;
    static std::mutex flowDevice_mutex;
    FlowDeviceFactory();
};

//! @defgroup flowDeviceGroup Flow Devices
//! Flow device objects connect zero-dimensional reactors.
//! FlowDevice objects should be instantiated via the newFlowDevice function, for
//! example:
//!
//! ```cpp
//!     shared_ptr<FlowDevice> mfc = newFlowDevice("MassFlowController", "my_mfc");
//! ```
//! @ingroup zerodGroup
//! @{

//! Create a FlowDevice object of the specified type
//! @since Starting in %Cantera 3.1, this method returns a `shared_ptr<FlowDevice>`
shared_ptr<FlowDevice> newFlowDevice(const string& model, const string& name="(none)");

//! Create a FlowDevice object of the specified type
//! @since New in %Cantera 3.0.
//! @deprecated Replaced by newFlowDevice. To be removed after %Cantera 3.1.
shared_ptr<FlowDevice> newFlowDevice3(const string& model);

//! @}
}

#endif
