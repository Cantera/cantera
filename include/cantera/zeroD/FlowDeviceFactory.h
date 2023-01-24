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
//!     unique_ptr<FlowDevice> mfc(newFlowDevice("MassFlowController"));
//! ```
//!
//! @ingroup ZeroD
class FlowDeviceFactory : public Factory<FlowDevice>
{
public:
    static FlowDeviceFactory* factory();

    virtual void deleteFactory();

    //! Create a new flow device by type name.
    /*!
     * @param flowDeviceType the type to be created.
     */
    virtual FlowDevice* newFlowDevice(const std::string& flowDeviceType);

private:
    static FlowDeviceFactory* s_factory;
    static std::mutex flowDevice_mutex;
    FlowDeviceFactory();
};

//! Create a FlowDevice object of the specified type
//! @ingroup ZeroD
FlowDevice* newFlowDevice(const string& model);

}

#endif
