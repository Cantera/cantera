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
     * @deprecated  To be removed after Cantera 3.0; replaceable by newFlowDevice3.
     */
    virtual FlowDevice* newFlowDevice(const std::string& flowDeviceType);

private:
    static FlowDeviceFactory* s_factory;
    static std::mutex flowDevice_mutex;
    FlowDeviceFactory();
};

//! Create a FlowDevice object of the specified type
//! @deprecated  To be changed after Cantera 3.0; for new behavior, see newFlowDevice3.
//! @ingroup ZeroD
FlowDevice* newFlowDevice(const string& model);

//! Create a FlowDevice object of the specified type
//! @ingroup ZeroD
//! @since New in Cantera 3.0.
//! @todo Transition back to newFlowDevice after Cantera 3.0
shared_ptr<FlowDevice> newFlowDevice3(const string& model);

}

#endif
