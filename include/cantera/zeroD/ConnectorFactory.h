//! @file ConnectorFactory.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CONNECTOR_FACTORY_H
#define CONNECTOR_FACTORY_H

#include "cantera/base/FactoryBase.h"
#include "cantera/zeroD/Connector.h"

namespace Cantera
{

class FlowDevice;
class WallBase;

//! Factory class to create Connector objects.
//!
//! This class is mainly used via the newConnector() function, for example:
//!
//! ```cpp
//!     shared_ptr<Connector> valve = newConnector("Valve", r0, r1, "my_valve");
//! ```
//!
//! where `r0` and `r1` are reactor objects.
class ConnectorFactory :
    public Factory<Connector, shared_ptr<ReactorBase>, shared_ptr<ReactorBase>, const string&>
{
public:
    static ConnectorFactory* factory();

    void deleteFactory() override;

private:
    static ConnectorFactory* s_factory;
    static std::mutex connector_mutex;
    ConnectorFactory();
};

//! @defgroup connectorGroup Connectors
//! %Connector objects connect zero-dimensional reactors.
//! Connector objects should be instantiated via the newConnector() function, for
//! example:
//!
//! ```cpp
//!     shared_ptr<Connector> valve = newConnector("Valve", r0, r1, "my_valve");
//! ```
//!
//! where `r0` and `r1` are reactor objects.
//!
//! @since New in %Cantera 3.1.
//!
//! @ingroup zerodGroup
//! @{

//! Create a Connector object of the specified type
//! @param model  String specifying reactor type.
//! @param r0  First reactor.
//! @param r1  Second reactor.
//! @param name  Name of the connector.
//! @since New in %Cantera 3.1.
shared_ptr<Connector> newConnector(const string& model,
                                   shared_ptr<ReactorBase> r0,
                                   shared_ptr<ReactorBase> r1,
                                   const string& name="(none)");

//! Create a FlowDevice object of the specified type
//! @since Starting in %Cantera 3.1, this method returns a `shared_ptr<FlowDevice>`
shared_ptr<FlowDevice> newFlowDevice(const string& model, const string& name="(none)");

//! Create a FlowDevice object of the specified type
//! @since New in %Cantera 3.0.
//! @deprecated Replaced by newFlowDevice. To be removed after %Cantera 3.1.
shared_ptr<FlowDevice> newFlowDevice3(const string& model);

//! Create a WallBase object of the specified type
//! @since Starting in %Cantera 3.1, this method returns a `shared_ptr<WallBase>`
shared_ptr<WallBase> newWall(const string& model, const string& name="(none)");

//! Create a WallBase object of the specified type
//! @since New in %Cantera 3.0.
//! @deprecated Replaced by newWall. To be removed after %Cantera 3.1.
shared_ptr<WallBase> newWall3(const string& model);

//! @}
}

#endif
