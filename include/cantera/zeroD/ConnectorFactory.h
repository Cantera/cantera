//! @file ConnectorFactory.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CONNECTOR_FACTORY_H
#define CONNECTOR_FACTORY_H

#include "cantera/base/FactoryBase.h"
#include "cantera/zeroD/ConnectorNode.h"

namespace Cantera
{

class FlowDevice;
class WallBase;

//! Factory class to create ConnectorNode objects.
//!
//! This class is mainly used via the newConnectorNode() function, for example:
//!
//! ```cpp
//!     shared_ptr<ConnectorNode> valve = newConnectorNode("Valve", r0, r1, "my_valve");
//! ```
//!
//! where `r0` and `r1` are reactor objects.
class ConnectorFactory :
    public Factory<ConnectorNode,
                   shared_ptr<ReactorBase>, shared_ptr<ReactorBase>, const string&>
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
//! %ConnectorNode objects connect zero-dimensional reactors.
//! ConnectorNode objects should be instantiated via the newConnectorNode() function,
//! for example:
//!
//! ```cpp
//!     shared_ptr<ConnectorNode> valve = newConnectorNode("Valve", r0, r1, "my_valve");
//! ```
//!
//! where `r0` and `r1` are reactor objects.
//!
//! @since New in %Cantera 3.2.
//!
//! @ingroup zerodGroup
//! @{

//! Create a ConnectorNode object of the specified type
//! @param model  String specifying reactor type.
//! @param r0  First reactor.
//! @param r1  Second reactor.
//! @param name  Name of the connector.
//! @since New in %Cantera 3.2.
shared_ptr<ConnectorNode> newConnectorNode(const string& model,
                                           shared_ptr<ReactorBase> r0,
                                           shared_ptr<ReactorBase> r1,
                                           const string& name="(none)");

//! Create a FlowDevice object of the specified type
//! @since Starting in %Cantera 3.1, this method returns a `shared_ptr<FlowDevice>`
//! @deprecated To be removed after %Cantera 3.2. Use version that provides reactors
//!     as parameter instead.
shared_ptr<FlowDevice> newFlowDevice(const string& model, const string& name="(none)");

//! Create a FlowDevice object of the specified type.
//! @copydetails newConnectorNode
shared_ptr<FlowDevice> newFlowDevice(const string& model,
                                     shared_ptr<ReactorBase> r0,
                                     shared_ptr<ReactorBase> r1,
                                     const string& name="(none)");

//! Create a WallBase object of the specified type
//! @since Starting in %Cantera 3.1, this method returns a `shared_ptr<WallBase>`
//! @deprecated To be removed after %Cantera 3.2. Use version that provides reactors
//!     as parameter instead.
shared_ptr<WallBase> newWall(const string& model, const string& name="(none)");

//! Create a WallBase object of the specified type.
//! @copydetails newConnectorNode
shared_ptr<WallBase> newWall(const string& model,
                             shared_ptr<ReactorBase> r0,
                             shared_ptr<ReactorBase> r1,
                             const string& name="(none)");

//! @}
}

#endif
