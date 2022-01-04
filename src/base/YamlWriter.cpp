// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/YamlWriter.h"
#include "cantera/base/AnyMap.h"
#include "cantera/base/Solution.h"
#include "cantera/base/stringUtils.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/Species.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/kinetics/Reaction.h"
#include "cantera/transport/TransportBase.h"

#include <set>
#include <fstream>
#include <chrono>

namespace Cantera {

YamlWriter::YamlWriter()
    : m_float_precision(15)
    , m_skip_user_defined(false)
{
}

void YamlWriter::setHeader(const AnyMap& header) {
    m_header = header;
}

void YamlWriter::addPhase(shared_ptr<Solution> soln, bool includeAdjacent) {
    for (auto& phase : m_phases) {
        if (phase->name() == soln->name()) {
            if (phase.get() == soln.get()) {
                // This phase has already been added, so nothing needs to be done. This
                // is expected in cases such as bulk phases adjacent to multiple
                // surface phases.
                return;
            } else {
                throw CanteraError("YamlWriter::addPhase",
                    "Duplicate phase name '{}'", soln->name());
            }
        }
    }
    m_phases.push_back(soln);
    if (includeAdjacent) {
        for (size_t i = 0; i < soln->nAdjacent(); i++) {
            addPhase(soln->adjacent(i));
        }
    }
}

void YamlWriter::addPhase(shared_ptr<ThermoPhase> thermo,
                          shared_ptr<Kinetics> kin,
                          shared_ptr<Transport> tran) {
    auto soln = Solution::create();
    soln->setThermo(thermo);
    soln->setKinetics(kin);
    soln->setTransport(tran);
    addPhase(soln);
}

std::string YamlWriter::toYamlString() const
{
    AnyMap output;
    bool hasDescription = m_header.hasKey("description");
    if (hasDescription) {
        output["description"] = m_header["description"];
    }
    output["generator"] = "YamlWriter";
    output["cantera-version"] = CANTERA_VERSION;
    output["git-commit"] = gitCommit();
    time_t now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    output["date"] = trimCopy(std::ctime(&now));
    if (hasDescription) {
        output["description"].setLoc(-6, 0);
    }
    output["generator"].setLoc(-5, 0);
    output["cantera-version"].setLoc(-4, 0);
    output["git-commit"].setLoc(-3, 0);
    output["date"].setLoc(-2, 0);

    // Add remaining header information, ignoring obsolete information
    std::set<std::string> exclude = {
        "description", "generator", "cantera-version", "git-commit", "date"};
    for (const auto& item : m_header) {
        std::string key = item.first;
        if (!exclude.count(key)) {
            output[key] = item.second;
        }
    }

    // Build phase definitions
    std::vector<AnyMap> phaseDefs(m_phases.size());
    size_t nspecies_total = 0;
    for (size_t i = 0; i < m_phases.size(); i++) {
        phaseDefs[i] = m_phases[i]->parameters(!m_skip_user_defined);
        if (m_phases[i]->nAdjacent()) {
            std::vector<std::string> adj_names;
            for (size_t j = 0; j < m_phases[i]->nAdjacent(); j++) {
                adj_names.push_back(m_phases[i]->adjacent(j)->name());
            }
            phaseDefs[i]["adjacent-phases"] = adj_names;
        }
        nspecies_total += m_phases[i]->thermo()->nSpecies();
    }
    output["phases"] = phaseDefs;

    // Build species definitions for all phases
    std::vector<AnyMap> speciesDefs;
    speciesDefs.reserve(nspecies_total);
    std::unordered_map<std::string, size_t> speciesDefIndex;
    for (const auto& phase : m_phases) {
        const auto thermo = phase->thermo();
        for (const auto& name : thermo->speciesNames()) {
            const auto& species = thermo->species(name);
            AnyMap speciesDef = species->parameters(thermo.get(), !m_skip_user_defined);

            if (speciesDefIndex.count(name) == 0) {
                speciesDefs.emplace_back(speciesDef);
                speciesDefIndex[name] = speciesDefs.size() - 1;
            } else if (speciesDefs[speciesDefIndex[name]] != speciesDef) {
                throw CanteraError("YamlWriter::toYamlString",
                    "Multiple species with different definitions are not "
                    "supported:\n>>>>>>\n{}\n======\n{}\n<<<<<<\n",
                    speciesDef.toYamlString(),
                    speciesDefs[speciesDefIndex[name]].toYamlString());
            }
        }
    }
    output["species"] = speciesDefs;

    // build reaction definitions for all phases
    std::map<std::string, std::vector<AnyMap>> allReactions;
    for (const auto& phase : m_phases) {
        const auto kin = phase->kinetics();
        if (!kin || !kin->nReactions()) {
            continue;
        }
        std::vector<AnyMap> reactions;
        for (size_t i = 0; i < kin->nReactions(); i++) {
            reactions.push_back(kin->reaction(i)->parameters(!m_skip_user_defined));
        }
        allReactions[phase->name()] = std::move(reactions);
    }

    // Figure out which phase definitions have identical sets of reactions,
    // and can share a reaction definition section

    // key: canonical phase in allReactions
    // value: phases using this reaction set
    std::map<std::string, std::vector<std::string>> phaseGroups;

    for (const auto& phase : m_phases) {
        const auto kin = phase->kinetics();
        std::string name = phase->name();
        if (!kin || !kin->nReactions()) {
            continue;
        }
        bool match = false;
        for (auto& group : phaseGroups) {
            if (allReactions[group.first] == allReactions[name]) {
                group.second.push_back(name);
                allReactions.erase(name);
                match = true;
                break;
            }
        }
        if (!match) {
            phaseGroups[name].push_back(name);
        }
    }

    // Generate the reactions section(s) in the output file
    if (phaseGroups.size() == 1) {
        output["reactions"] = std::move(allReactions[phaseGroups.begin()->first]);
    } else {
        for (const auto& group : phaseGroups) {
            std::string groupName;
            for (auto& name : group.second) {
                groupName += name + "-";
            }
            groupName += "reactions";
            output[groupName] = std::move(allReactions[group.first]);

            for (auto& name : group.second) {
                AnyMap& phaseDef = output["phases"].getMapWhere("name", name);
                phaseDef["reactions"] = std::vector<std::string>{groupName};
            }
        }
    }

    output.setMetadata("precision", AnyValue(m_float_precision));
    output.setUnits(m_output_units);
    return output.toYamlString();
}

void YamlWriter::toYamlFile(const std::string& filename) const
{
    std::ofstream out(filename);
    out << toYamlString();
}

void YamlWriter::setUnits(const std::map<std::string, std::string>& units)
{
    m_output_units = UnitSystem();
    m_output_units.setDefaults(units);
}

void YamlWriter::setUnitSystem(const UnitSystem& units)
{
    m_output_units = units;
}

}
