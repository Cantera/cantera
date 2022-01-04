/**
 *  @file Solution.cpp
 *   Definition file for class Solution.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/Solution.h"
#include "cantera/base/Interface.h"
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
}

void Solution::setKinetics(shared_ptr<Kinetics> kinetics) {
    m_kinetics = kinetics;
}

void Solution::setTransport(shared_ptr<Transport> transport) {
    m_transport = transport;
}

void Solution::addAdjacent(shared_ptr<Solution> adjacent) {
    if (m_adjacentByName.count(adjacent->name())) {
        throw CanteraError("Solution::addAdjacent",
            "Solution '{}' already contains an adjacent phase named '{}'.",
            name(), adjacent->name());
    }
    if (m_thermo && adjacent->thermo()
        && adjacent->thermo()->nDim() <= m_thermo->nDim())
    {
        throw CanteraError("Solution::addAdjacent",
            "Adjacent phases should have higher dimensionality than the reacting ",
            "phase.\n'{}' is {}-dimensional while '{}' is {}-dimensional",
            adjacent->thermo()->name(), adjacent->thermo()->nDim(),
            m_thermo->name(), m_thermo->nDim());
    }
    m_adjacent.push_back(adjacent);
    m_adjacentByName[adjacent->name()] = adjacent;
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

shared_ptr<Solution> newSolution(const std::string& infile, const std::string& name,
    const std::string& transport, const std::vector<std::string>& adjacent)
{
    // @todo Remove file extension check after Cantera 2.6
    // get file extension
    size_t dot = infile.find_last_of(".");
    std::string extension;
    if (dot != npos) {
        extension = toLowerCopy(infile.substr(dot+1));
    }

    if (extension == "xml" || extension == "cti") {
        throw CanteraError("newSolution(string infile, string name, string transport, "
            "vector<string> adjacent)",
            "This constructor is only compatible with YAML input files");
    }

    auto rootNode = AnyMap::fromYamlFile(infile);
    AnyMap& phaseNode = rootNode["phases"].getMapWhere("name", name);

    std::vector<shared_ptr<Solution>> adjPhases;
    // Create explicitly-specified adjacent bulk phases
    for (auto& name : adjacent) {
        auto& adjNode = rootNode["phases"].getMapWhere("name", name);
        adjPhases.push_back(newSolution(adjNode, rootNode));
    }
    return newSolution(phaseNode, rootNode, transport, adjPhases);
}

shared_ptr<Solution> newSolution(const AnyMap& phaseNode,
                                 const AnyMap& rootNode,
                                 const std::string& transport,
                                 const std::vector<shared_ptr<Solution>>& adjacent,
                                 const std::map<std::string, shared_ptr<Solution>>& related)
{
    // thermo phase
    auto thermo = shared_ptr<ThermoPhase>(newPhase(phaseNode, rootNode));

    // instantiate Solution object of the correct derived type
    shared_ptr<Solution> sol;
    switch (thermo->nDim()) {
    case 2:
        sol = Interface::create();
        break;
    default:
        sol = Solution::create();
    }
    sol->setSource("custom YAML");
    sol->setThermo(thermo);

    // Add explicitly-specified adjacent phases
    for (auto& adj : adjacent) {
        sol->addAdjacent(adj);
    }

    // If no adjacent phases were explicitly specified, look for them in the interface
    // phase definition
    if (adjacent.empty() && phaseNode.hasKey("adjacent-phases")) {
        auto all_related = related;
        for (auto& phase : adjacent) {
            all_related[phase->name()] = phase;
        }

        // Helper function for adding individual phases
        auto addPhase = [&](const AnyValue& phases, const AnyMap& root,
                            const std::string& name)
        {
            if (!all_related.count(name)) {
                // Create a new phase only if there isn't already one with the same name
                auto adj = newSolution(phases.getMapWhere("name", name), root,
                                    "", {}, all_related);
                all_related[name] = adj;
                for (size_t i = 0; i < adj->nAdjacent(); i++) {
                    all_related[adj->adjacent(i)->name()] = adj->adjacent(i);
                }
            }
            sol->addAdjacent(all_related[name]);
        };

        auto& adjPhases = phaseNode["adjacent-phases"];
        if (adjPhases.is<std::vector<std::string>>()) {
            // 'adjacent' is a list of bulk phases from the current input file
            for (auto& phase : adjPhases.as<std::vector<std::string>>()) {
                addPhase(rootNode["phases"], rootNode, phase);
            }
        } else if (adjPhases.is<std::vector<AnyMap>>()) {
            // Each element of 'adjacent' is a map with one item, where the key is
            // a section in this file or another YAML file, and the value is a list of
            // phase names to read from that section
            for (auto& item : adjPhases.asVector<AnyMap>()) {
                const std::string& source = item.begin()->first;
                const auto& names = item.begin()->second.asVector<std::string>();
                const auto& slash = boost::ifind_last(source, "/");
                if (slash) {
                    // source is a different input file
                    std::string fileName(source.begin(), slash.begin());
                    std::string node(slash.end(), source.end());
                    AnyMap phaseSource = AnyMap::fromYamlFile(fileName,
                        rootNode.getString("__file__", ""));
                    for (auto& phase : names) {
                        addPhase(phaseSource[node], phaseSource, phase);
                    }
                } else if (rootNode.hasKey(source)) {
                    // source is in the current file
                    for (auto& phase : names) {
                        addPhase(rootNode[source], rootNode, phase);
                    }
                } else {
                    throw InputFileError("newSolution", adjPhases,
                        "Could not find a phases section named '{}'.", source);
                }
            }
        } else {
            throw InputFileError("addAdjacentPhases", adjPhases,
                "Could not parse adjacent phase declaration of type '{}'",
                adjPhases.type_str());
        }
    }

    // kinetics
    std::vector<ThermoPhase*> phases;
    phases.push_back(sol->thermo().get());
    for (size_t i = 0; i < sol->nAdjacent(); i++) {
        phases.push_back(sol->adjacent(i)->thermo().get());
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

    // save root-level information (YAML header)
    AnyMap header;
    for (const auto& item : rootNode.ordered()) {
        std::string key = item.first;
        if (key == "phases") {
            // header ends with "phases" field
            break;
        } else if (key != "units") {
            header[key] = item.second;
        }
    }
    sol->header() = header;

    return sol;
}

} // namespace Cantera
