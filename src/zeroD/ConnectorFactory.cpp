//! @file ConnectorFactory.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/ConnectorFactory.h"
#include "cantera/zeroD/flowControllers.h"
#include "cantera/zeroD/Wall.h"

namespace Cantera
{

ConnectorFactory* ConnectorFactory::s_factory = 0;
std::mutex ConnectorFactory::connector_mutex;

ConnectorFactory::ConnectorFactory()
{
    reg("MassFlowController",
        [](shared_ptr<ReactorBase> r0, shared_ptr<ReactorBase> r1, const string& name)
        { return new MassFlowController(r0, r1, name); });
    reg("PressureController",
        [](shared_ptr<ReactorBase> r0, shared_ptr<ReactorBase> r1, const string& name)
        { return new PressureController(r0, r1, name); });
    reg("Valve",
        [](shared_ptr<ReactorBase> r0, shared_ptr<ReactorBase> r1, const string& name)
        { return new Valve(r0, r1, name); });
    reg("Wall",
        [](shared_ptr<ReactorBase> r0, shared_ptr<ReactorBase> r1, const string& name)
        { return new Wall(r0, r1, name); });
}

ConnectorFactory* ConnectorFactory::factory() {
    std::unique_lock<std::mutex> lock(connector_mutex);
    if (!s_factory) {
        s_factory = new ConnectorFactory;
    }
    return s_factory;
}

void ConnectorFactory::deleteFactory() {
    std::unique_lock<std::mutex> lock(connector_mutex);
    delete s_factory;
    s_factory = 0;
}

shared_ptr<ConnectorNode> newConnectorNode(
    const string& model,
    shared_ptr<ReactorBase> r0, shared_ptr<ReactorBase> r1, const string& name)
{
    return shared_ptr<ConnectorNode>(
        ConnectorFactory::factory()->create(model, r0, r1, name));
}

}
