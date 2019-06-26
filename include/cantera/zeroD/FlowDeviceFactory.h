//! @file FlowDeviceFactory.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef FLOWDEVICE_FACTORY_H
#define FLOWDEVICE_FACTORY_H

#include "cantera/base/FactoryBase.h"
#include "cantera/zeroD/FlowDevice.h"

namespace Cantera
{

class FlowDeviceFactory : public Factory<FlowDevice>
{
public:
    static FlowDeviceFactory* factory() {
        std::unique_lock<std::mutex> lock(flowDevice_mutex);
        if (!s_factory) {
            s_factory = new FlowDeviceFactory;
        }
        return s_factory;
    }

    virtual void deleteFactory() {
        std::unique_lock<std::mutex> lock(flowDevice_mutex);
        delete s_factory;
        s_factory = 0;
    }

    //! Create a new flow device by type identifier.
    /*!
     * @param n the type to be created.
     */
    virtual FlowDevice* newFlowDevice(int n);

    //! Create a new flow device by type name.
    /*!
     * @param flowDeviceType the type to be created.
     */
    virtual FlowDevice* newFlowDevice(const std::string& flowDeviceType);

    //! Register a new flow device type identifier.
    /*!
     * @param name the name of the flow device type.
     * @param type the type identifier of the flow device.
     * Integer type identifiers are used by clib and matlab interfaces.
     *
     * @deprecated To be removed after Cantera 2.5.
     */
    void reg_type(const std::string& name, const int type) {
        m_types[type] = name;
    }

protected:
    //! Map containing flow device type identifier / type name pairs.
    //! @deprecated To be removed after Cantera 2.5.
    std::unordered_map<int, std::string> m_types;

private:
    static FlowDeviceFactory* s_factory;
    static std::mutex flowDevice_mutex;
    FlowDeviceFactory();
};

//! Create a FlowDevice object of the specified type
inline FlowDevice* newFlowDevice(const std::string& model)
{
    return FlowDeviceFactory::factory()->newFlowDevice(model);
}

}

#endif
