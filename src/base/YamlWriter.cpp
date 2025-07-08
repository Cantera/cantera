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
#include "cantera/transport/Transport.h"

#include <fstream>
#include <chrono>

namespace Cantera {

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

string YamlWriter::toYamlString() const
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
    set<string> exclude = {
        "description", "generator", "cantera-version", "git-commit", "date"};
    for (const auto& [key, value] : m_header) {
        if (!exclude.count(key)) {
            output[key] = value;
        }
    }

    // Build phase definitions
    vector<AnyMap> phaseDefs(m_phases.size());
    size_t nspecies_total = 0;
    for (size_t i = 0; i < m_phases.size(); i++) {
        phaseDefs[i] = m_phases[i]->parameters(!m_skip_user_defined);
        if (m_phases[i]->nAdjacent()) {
            vector<string> adj_names;
            for (size_t j = 0; j < m_phases[i]->nAdjacent(); j++) {
                adj_names.push_back(m_phases[i]->adjacent(j)->name());
            }
            phaseDefs[i]["adjacent-phases"] = adj_names;
        }
        nspecies_total += m_phases[i]->thermo()->nSpecies();
    }
    output["phases"] = phaseDefs;

    // Build species definitions for all phases
    vector<AnyMap> speciesDefs;
    speciesDefs.reserve(nspecies_total);
    std::unordered_map<string, size_t> speciesDefIndex;
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
    map<string, vector<AnyMap>> allReactions;
    for (const auto& phase : m_phases) {
        const auto kin = phase->kinetics();
        if (!kin || !kin->nReactions()) {
            continue;
        }
        // Update duplicate reaction labels in case of a dynamically-built mechanism
        kin->checkDuplicates(false, true);
        vector<AnyMap> reactions;
        for (size_t i = 0; i < kin->nReactions(); i++) {
            reactions.push_back(kin->reaction(i)->parameters(!m_skip_user_defined));
        }
        allReactions[phase->name()] = std::move(reactions);
    }

    // Figure out which phase definitions have identical sets of reactions,
    // and can share a reaction definition section

    // key: canonical phase in allReactions
    // value: phases using this reaction set
    map<string, vector<string>> phaseGroups;

    for (const auto& phase : m_phases) {
        const auto kin = phase->kinetics();
        string name = phase->name();
        if (!kin || !kin->nReactions()) {
            continue;
        }
        bool match = false;
        for (auto& [canonicalPhase, dependentPhases] : phaseGroups) {
            if (allReactions[canonicalPhase] == allReactions[name]) {
                dependentPhases.push_back(name);
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
        for (const auto& [canonicalPhase, dependentPhases] : phaseGroups) {
            string groupName;
            for (auto& name : dependentPhases) {
                groupName += name + "-";
            }
            groupName += "reactions";
            output[groupName] = std::move(allReactions[canonicalPhase]);

            for (auto& name : dependentPhases) {
                AnyMap& phaseDef = output["phases"].getMapWhere("name", name);
                phaseDef["reactions"] = vector<string>{groupName};
            }
        }
    }

    output.setMetadata("precision", AnyValue(m_float_precision));
    output.setUnits(m_output_units);
    return output.toYamlString();
}

void YamlWriter::toYamlFile(const string& filename) const
{
    std::ofstream out(filename);
    out << toYamlString();
}

void YamlWriter::setUnits(const map<string, string>& units)
{
    m_output_units = UnitSystem();
    m_output_units.setDefaults(units);
}

void YamlWriter::setUnitSystem(const UnitSystem& units)
{
    m_output_units = units;
}

}
