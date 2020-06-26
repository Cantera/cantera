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
#include "cantera/base/stringUtils.h"

namespace Cantera
{

Solution::Solution() :
    m_description("") {
}

std::string Solution::name() const {
    if (m_thermo) {
        return m_thermo->name();
    } else {
        throw CanteraError("Solution::name",
                           "Requires associated 'ThermoPhase'");
    }
}

void Solution::setName(const std::string& name) {
    if (m_thermo) {
        m_thermo->setName(name);
    } else {
        throw CanteraError("Solution::setName",
                           "Requires associated 'ThermoPhase'");
    }
}

void Solution::setThermo(shared_ptr<ThermoPhase> thermo) {
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

AnyMap Solution::parameters(bool withInput) const
{
    AnyMap out = m_thermo->parameters(false);
    if (m_kinetics) {
        out.update(m_kinetics->parameters());
    }
    if (m_transport) {
        out.update(m_transport->parameters());
    }
    if (withInput) {
        out.update(m_thermo->input());
    }
    return out;
}

shared_ptr<Solution> newSolution(const std::string& infile,
                                 const std::string& name,
                                 const std::string& transport,
                                 const std::vector<shared_ptr<Solution>>& adjacent) {

    // instantiate Solution object
    auto sol = Solution::create();

    // description
    size_t dot = infile.find_last_of(".");
    std::string extension;
    if (dot != npos) {
        extension = toLowerCopy(infile.substr(dot+1));
    }
    if (extension == "yml" || extension == "yaml") {
        AnyMap root = AnyMap::fromYamlFile(infile);
        sol->setDescription(root.getString("description", ""));
    }

    // thermo phase
    sol->setThermo(shared_ptr<ThermoPhase>(newPhase(infile, name)));

    // kinetics
    std::vector<ThermoPhase*> phases;
    phases.push_back(sol->thermo().get());
    for (auto& adj : adjacent) {
        phases.push_back(adj->thermo().get());
    }
    sol->setKinetics(newKinetics(phases, infile, name));

    // transport
    if (transport == "") {
        sol->setTransport(shared_ptr<Transport>(newDefaultTransportMgr(sol->thermo().get())));
    } else if (transport != "None") {
        sol->setTransport(shared_ptr<Transport>(newTransportMgr(transport, sol->thermo().get())));
    }

    return sol;
}

} // namespace Cantera
