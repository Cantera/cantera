/**
 *  @file Solution.cpp
 *   Definition file for class Solution.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/Solution.h"
#include "cantera/base/Interface.h"
#include "cantera/base/ExtensionManager.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/kinetics/KineticsFactory.h"
#include "cantera/transport/Transport.h"
#include "cantera/transport/TransportFactory.h"
#include "cantera/base/stringUtils.h"

#include <boost/algorithm/string.hpp>

namespace Cantera
{

shared_ptr<Solution> Solution::clone(const vector<shared_ptr<Solution>>& adjacent,
    bool withKinetics, bool withTransport) const
{
    shared_ptr<Solution> out = create();
    out->setThermo(m_thermo->clone());
    vector<shared_ptr<ThermoPhase>> kinPhases;
    if (withKinetics) {
        kinPhases.push_back(out->thermo());
        if (!adjacent.empty()) {
            // Use the provided adjacent phases
            for (auto& soln : adjacent) {
                kinPhases.push_back(soln->thermo());
                out->addAdjacent(soln);
            }
        } else {
            // Clone new adjacent phases
            for (size_t i = 1; i < m_kinetics->nPhases(); i++) {
                auto soln = m_kinetics->phase(i)->root()->clone();
                kinPhases.push_back(soln->thermo());
                out->addAdjacent(soln);
            }
        }
        out->setKinetics(m_kinetics->clone(kinPhases));
    } else {
        out->setKinetics(newKinetics("none"));
    }
    if (withTransport) {
        out->setTransport(m_transport->clone(out->thermo()));
    } else {
        out->setTransport(newTransport(m_thermo, "none"));
    }
    return out;
}

string Solution::name() const {
    if (m_thermo) {
        return m_thermo->name();
    } else {
        throw CanteraError("Solution::name",
                           "Requires associated 'ThermoPhase'");
    }
}

void Solution::setName(const string& name) {
    if (m_thermo) {
        m_thermo->setName(name);
    } else {
        throw CanteraError("Solution::setName",
                           "Requires associated 'ThermoPhase'");
    }
}

void Solution::setThermo(shared_ptr<ThermoPhase> thermo) {
    m_thermo = thermo;
    m_thermo->setSolution(weak_from_this());
    for (const auto& [id, callback] : m_changeCallbacks) {
        callback();
    }
}

void Solution::setKinetics(shared_ptr<Kinetics> kinetics) {
    if (kinetics == m_kinetics) {
        return;
    }
    m_kinetics = kinetics;
    if (m_kinetics) {
        m_kinetics->setRoot(shared_from_this());
    }
    for (const auto& [id, callback] : m_changeCallbacks) {
        callback();
    }
}

string Solution::transportModel()
{
    if (!m_transport) {
        throw CanteraError("Solution::transportModel",
            "The Transport object is not initialized.");
    }
    return m_transport->transportModel();
}

void Solution::setTransport(shared_ptr<Transport> transport) {
    if (transport == m_transport) {
        return;
    }
    m_transport = transport;
    for (const auto& [id, callback] : m_changeCallbacks) {
        callback();
    }
}

void Solution::setTransportModel(const string& model) {
    if (!m_thermo) {
        throw CanteraError("Solution::setTransportModel",
            "Unable to set Transport model without valid ThermoPhase object.");
    }
    if (m_transport && transportModel() == model) {
        return;
    }
    setTransport(newTransport(m_thermo, model));
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
    } else if (m_transport->transportModel() == "none") {
        out["transport"] = empty;
    } else {
        out.update(m_transport->parameters());
    }
    if (withInput) {
        auto transport = out["transport"];
        AnyMap input = m_thermo->input();
        if (input.hasKey("reactions")) {
            // all reactions are listed in the standard 'reactions' section
            input.erase("reactions");
        }
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

const string Solution::source() const {
    AnyValue source = m_header.getMetadata("filename");
    return source.empty() ? "<unknown>" : source.asString();
}

void Solution::setSource(const string& source) {
    AnyValue filename(source);
    m_header.setMetadata("filename", filename);
}

void Solution::holdExternalHandle(const string& name,
                                  shared_ptr<ExternalHandle> handle)
{
    m_externalHandles[name] = handle;
}

shared_ptr<ExternalHandle> Solution::getExternalHandle(const string& name) const
{
    if (m_externalHandles.count(name)) {
        return m_externalHandles.at(name);
    } else {
        return shared_ptr<ExternalHandle>();
    }
}

void Solution::registerChangedCallback(void *id, const function<void()>& callback)
{
    m_changeCallbacks[id] = callback;
}

void Solution::removeChangedCallback(void* id)
{
    m_changeCallbacks.erase(id);
}

shared_ptr<Solution> newSolution(const string &infile,
                                 const string &name,
                                 const string &transport,
                                 const vector<shared_ptr<Solution>> &adjacent)
{
    // get file extension
    size_t dot = infile.find_last_of(".");
    string extension;
    if (dot != npos) {
        extension = toLowerCopy(infile.substr(dot+1));
    }

    if (extension == "cti" || extension == "xml") {
        throw CanteraError("newSolution",
                           "The CTI and XML formats are no longer supported.");
    }

    // load YAML file
    auto rootNode = AnyMap::fromYamlFile(infile);
    const AnyMap& phaseNode = rootNode.at("phases").getMapWhere("name", name);
    auto sol = newSolution(phaseNode, rootNode, transport, adjacent);
    sol->setSource(infile);
    return sol;
}

shared_ptr<Solution> newSolution(const string& infile, const string& name,
    const string& transport, const vector<string>& adjacent)
{
    auto rootNode = AnyMap::fromYamlFile(infile);
    const AnyMap& phaseNode = rootNode.at("phases").getMapWhere("name", name);

    vector<shared_ptr<Solution>> adjPhases;
    // Create explicitly-specified adjacent bulk phases
    for (auto& name : adjacent) {
        const auto& adjNode = rootNode.at("phases").getMapWhere("name", name);
        adjPhases.push_back(newSolution(adjNode, rootNode));
    }
    return newSolution(phaseNode, rootNode, transport, adjPhases);
}

shared_ptr<Solution> newSolution(const AnyMap& phaseNode,
                                 const AnyMap& rootNode,
                                 const string& transport,
                                 const vector<shared_ptr<Solution>>& adjacent,
                                 const map<string, shared_ptr<Solution>>& related)
{
    // thermo phase
    auto thermo = newThermo(phaseNode, rootNode);

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
                            const string& name)
        {
            if (!all_related.count(name)) {
                // Create a new phase only if there isn't already one with the same name
                auto adj = newSolution(phases.getMapWhere("name", name), root,
                                       "default", {}, all_related);
                all_related[name] = adj;
                for (size_t i = 0; i < adj->nAdjacent(); i++) {
                    all_related[adj->adjacent(i)->name()] = adj->adjacent(i);
                }
            }
            sol->addAdjacent(all_related[name]);
        };

        auto& adjPhases = phaseNode["adjacent-phases"];
        if (adjPhases.is<vector<string>>()) {
            // 'adjacent' is a list of bulk phases from the current input file
            for (auto& phase : adjPhases.as<vector<string>>()) {
                addPhase(rootNode["phases"], rootNode, phase);
            }
        } else if (adjPhases.is<vector<AnyMap>>()) {
            // Each element of 'adjacent' is a map with one item, where the key is
            // a section in this file or another YAML file, and the value is a list of
            // phase names to read from that section
            for (auto& item : adjPhases.asVector<AnyMap>()) {
                const string& source = item.begin()->first;
                const auto& names = item.begin()->second.asVector<string>();
                const auto& slash = boost::ifind_last(source, "/");
                if (slash) {
                    // source is a different input file
                    string fileName(source.begin(), slash.begin());
                    string node(slash.end(), source.end());
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
    vector<shared_ptr<ThermoPhase>> phases;
    phases.push_back(sol->thermo());
    for (size_t i = 0; i < sol->nAdjacent(); i++) {
        phases.push_back(sol->adjacent(i)->thermo());
    }
    sol->setKinetics(newKinetics(phases, phaseNode, rootNode, sol));

    // set transport model by name
    sol->setTransportModel(transport);

    // save root-level information (YAML header)
    AnyMap header;
    for (const auto& [key, value] : rootNode.ordered()) {
        if (key == "phases") {
            // header ends with "phases" field
            break;
        }
        header[key] = value;
    }
    sol->header() = header;

    return sol;
}

} // namespace Cantera
