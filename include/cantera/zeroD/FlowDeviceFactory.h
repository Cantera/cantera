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
//! This class is mainly used via the newFlowDevice3() function, for example:
//!
//! ```cpp
//!     shared_ptr<FlowDevice> mfc = newFlowDevice3("MassFlowController");
//! ```
class FlowDeviceFactory : public Factory<FlowDevice>
{
public:
    static FlowDeviceFactory* factory();

    void deleteFactory() override;

    //! Create a new flow device by type name.
    /*!
     * @param flowDeviceType the type to be created.
     * @deprecated To be removed after %Cantera 3.0; replaceable by newFlowDevice3.
     */
    FlowDevice* newFlowDevice(const string& flowDeviceType);

private:
    static FlowDeviceFactory* s_factory;
    static std::mutex flowDevice_mutex;
    FlowDeviceFactory();
};

//! @defgroup flowDeviceGroup Flow Devices
//! Flow device objects connect zero-dimensional reactors.
//! FlowDevice objects should be instantiated via the newFlowDevice3() function, for
//! example:
//!
//! ```cpp
//!     shared_ptr<FlowDevice> mfc = newFlowDevice3("MassFlowController");
//! ```
//! @ingroup zerodGroup
//! @{

//! Create a FlowDevice object of the specified type
//! @deprecated To be changed after %Cantera 3.0; for new behavior, see
//! newFlowDevice3().
FlowDevice* newFlowDevice(const string& model);

//! Create a FlowDevice object of the specified type
//! @since New in %Cantera 3.0.
//! @todo Transition back to newFlowDevice() after %Cantera 3.0
shared_ptr<FlowDevice> newFlowDevice3(const string& model);

//! @}
}

#endif
