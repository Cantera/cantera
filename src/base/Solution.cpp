/**
 *  @file Solution.cpp
 *   Definition file for class Solution.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/Solution.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/kinetics/KineticsFactory.h"
#include "cantera/transport/TransportBase.h"
#include "cantera/transport/TransportFactory.h"

namespace Cantera
{

Solution::Solution() {}

  Solution::Solution(const std::string& infile, std::string name, std::string transport,
                     std::vector<shared_ptr<Solution> > adjacent) {

    // thermo phase
    m_thermo = shared_ptr<ThermoPhase>(newPhase(infile, name));

    // kinetics
    std::vector<ThermoPhase*> phases;
    for (auto & adj : adjacent) {
        phases.push_back(adj->thermo().get());
    }
    phases.push_back(m_thermo.get());
    m_kinetics = std::move(newKinetics(phases, infile, name));

    // transport
    if (transport=="") {
        m_transport = shared_ptr<Transport>(newDefaultTransportMgr(m_thermo.get()));
    } else if (transport!="None") {
        m_transport = shared_ptr<Transport>(newTransportMgr(transport, m_thermo.get()));
    }
}

std::string Solution::name() const {
    if (m_thermo) {
        return m_thermo->name();
    } else {
        throw CanteraError("Solution::name()",
                           "Requires associated 'ThermoPhase'");
    }
}

void Solution::setName(const std::string& name) {
    if (m_thermo) {
        m_thermo->setName(name);
    } else {
        throw CanteraError("Solution::setName()",
                           "Requires associated 'ThermoPhase'");
    }
}

void Solution::setThermoPhase(shared_ptr<ThermoPhase> thermo) {
    m_thermo = thermo;
    if (m_thermo) {
        m_thermo->setRoot(shared_from_this());
    }
}

void Solution::setKinetics(shared_ptr<Kinetics> kinetics) {
    m_kinetics = kinetics;
    if (m_kinetics) {
        m_kinetics->setRoot(shared_from_this());
    }
}

void Solution::setTransport(shared_ptr<Transport> transport) {
    m_transport = transport;
    if (m_transport) {
        m_transport->setRoot(shared_from_this());
    }
}

} // namespace Cantera
