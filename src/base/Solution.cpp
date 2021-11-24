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

Solution::Solution() {}

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
    AnyValue empty("<NULL>");
    if (m_kinetics) {
        out.update(m_kinetics->parameters());
    }
    if (!m_transport) {
        out["transport"] = empty;
    } else if (m_transport->transportType() == "None") {
        out["transport"] = empty;
    } else {
        out.update(m_transport->parameters());
    }
    if (withInput) {
        auto transport = out["transport"];
        AnyMap input = m_thermo->input();
        out.update(input);
        if (input.hasKey("transport")) {
            // revert changes / ensure that correct model is referenced
            out["transport"] = transport;
        }
    }
    if (out["transport"] == empty) {
        out.erase("transport");
    }
    return out;
}

const AnyMap& Solution::header() const
{
    return m_header;
}

AnyMap& Solution::header()
{
    return m_header;
}

const std::string Solution::source() const {
    AnyValue source = m_header.getMetadata("filename");
    return source.empty() ? "<unknown>" : source.asString();
}

void Solution::setSource(const std::string& source) {
    AnyValue filename(source);
    m_header.setMetadata("filename", filename);
}

shared_ptr<Solution> newSolution(const std::string& infile,
                                 const std::string& name,
                                 const std::string& transport,
                                 const std::vector<shared_ptr<Solution>>& adjacent)
{
    // get file extension
    size_t dot = infile.find_last_of(".");
    std::string extension;
    if (dot != npos) {
        extension = toLowerCopy(infile.substr(dot+1));
    }

    if (extension == "yml" || extension == "yaml") {
        // load YAML file
        auto rootNode = AnyMap::fromYamlFile(infile);
        AnyMap& phaseNode = rootNode["phases"].getMapWhere("name", name);
        auto sol = newSolution(phaseNode, rootNode, transport, adjacent);
        sol->setSource(infile);
        return sol;
    }

    // instantiate Solution object
    auto sol = Solution::create();
    sol->setSource(infile);

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
        sol->setTransport(shared_ptr<Transport>(
            newDefaultTransportMgr(sol->thermo().get())));
    } else if (transport == "None") {
        sol->setTransport(shared_ptr<Transport>(newTransportMgr("None")));
    } else {
        sol->setTransport(shared_ptr<Transport>(
            newTransportMgr(transport, sol->thermo().get())));
    }

    return sol;
}

std::vector<std::string> getExcludes(const AnyValue& node)
{
    std::vector<std::string> exclude{"phases", "species", "reactions"};

    // retrieve alternate names of reactions sections from phases entries
    if (node.is<std::vector<AnyMap>>()) {

        for (const auto& phaseNode : node.as<std::vector<AnyMap>>()) {

            if (phaseNode.hasKey("reactions")) {
                const auto& reactionsNode = phaseNode.at("reactions");
                if (reactionsNode.is<std::string>()) {
                    // this may include 'none', 'all', or similar
                    exclude.push_back(reactionsNode.asString());
                } else if (reactionsNode.is<std::vector<std::string>>()) {
                    // List of sections from which all species should be added
                    for (const auto& item : reactionsNode.as<std::vector<std::string>>()) {
                        exclude.push_back(item);
                    }
                } else if (reactionsNode.is<std::vector<AnyMap>>()) {
                    // Mapping of rules to apply for each specified section containing
                    // reactions
                    for (const auto& item : reactionsNode.as<std::vector<AnyMap>>()) {
                        exclude.push_back(item.begin()->first);
                    }
                }
            }
        }
    }

    return exclude;
}

shared_ptr<Solution> newSolution(AnyMap& phaseNode,
                                 const AnyMap& rootNode,
                                 const std::string& transport,
                                 const std::vector<shared_ptr<Solution>>& adjacent)
{
    // instantiate Solution object
    auto sol = Solution::create();
    sol->setSource("custom YAML");

    // thermo phase
    sol->setThermo(shared_ptr<ThermoPhase>(newPhase(phaseNode, rootNode)));

    // kinetics
    std::vector<ThermoPhase*> phases;
    phases.push_back(sol->thermo().get());
    for (auto& adj : adjacent) {
        phases.push_back(adj->thermo().get());
    }
    sol->setKinetics(newKinetics(phases, phaseNode, rootNode));

    // transport
    if (transport == "") {
        sol->setTransport(shared_ptr<Transport>(
            newDefaultTransportMgr(sol->thermo().get())));
    } else if (transport == "None") {
        sol->setTransport(shared_ptr<Transport>(newTransportMgr("None")));
    } else {
        sol->setTransport(shared_ptr<Transport>(
            newTransportMgr(transport, sol->thermo().get())));
    }

    // find root-level fields that are not related to phases, species or reactions
    const auto& phasesNode = rootNode.at("phases");
    auto exclude = getExcludes(phasesNode);

    // save root-level information (YAML header)
    AnyMap header;
    for (const auto& item : rootNode) {
        std::string key = item.first;
        if (find(exclude.begin(), exclude.end(), key) == exclude.end()) {
            header[key] = item.second;
        }
    }
    header.setUnits(rootNode.units());
    sol->header() = header;

    return sol;
}

} // namespace Cantera
