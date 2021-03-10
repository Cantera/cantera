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

#include <fstream>
#include <chrono>

namespace Cantera {

YamlWriter::YamlWriter()
    : m_float_precision(15)
    , m_skip_user_defined(false)
{
}

void YamlWriter::addPhase(shared_ptr<Solution> soln) {
    for (auto& phase : m_phases) {
        if (phase->name() == soln->name()) {
            throw CanteraError("YamlWriter::addPhase",
                "Duplicate phase name '{}'", soln->name());
        }
    }
    m_phases.push_back(soln);
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
    output["generator"] = "YamlWriter";
    output["cantera-version"] = CANTERA_VERSION;
    time_t now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    output["date"] = trimCopy(std::ctime(&now));

    // Build phase definitions
    std::vector<AnyMap> phaseDefs(m_phases.size());
    size_t nspecies_total = 0;
    for (size_t i = 0; i < m_phases.size(); i++) {
        m_phases[i]->thermo()->getParameters(phaseDefs[i]);
        nspecies_total += m_phases[i]->thermo()->nSpecies();
        const auto& kin = m_phases[i]->kinetics();
        if (kin) {
            kin->getParameters(phaseDefs[i]);
            if (phaseDefs[i].hasKey("kinetics") && kin->nReactions() == 0) {
                phaseDefs[i]["reactions"] = "none";
            }
        }
        const auto& tran = m_phases[i]->transport();
        if (tran) {
            tran->getParameters(phaseDefs[i]);
        }

        if (!m_skip_user_defined) {
            for (const auto& item : m_phases[i]->thermo()->input()) {
                if (!phaseDefs[i].hasKey(item.first)) {
                    phaseDefs[i][item.first] = item.second;
                }
            }
        }
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
            AnyMap speciesDef;
            species->getParameters(speciesDef, !m_skip_user_defined);

            thermo->getSpeciesParameters(name, speciesDef);
            if (!m_skip_user_defined
                && species->input.hasKey("equation-of-state")) {
                auto& eosIn = species->input["equation-of-state"].asVector<AnyMap>();
                for (const auto& eos : eosIn) {
                    auto& out = speciesDef["equation-of-state"].getMapWhere(
                        "model", eos["model"].asString(), true);
                    for (const auto& item : eos) {
                        if (!out.hasKey(item.first)) {
                            out[item.first] = item.second;
                        }
                    }
                }
            }

            if (!m_skip_user_defined) {
                for (const auto& item : species->input) {
                    if (!speciesDef.hasKey(item.first)) {
                        speciesDef[item.first] = item.second;
                    }
                }
            }

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
            const auto reaction = kin->reaction(i);
            AnyMap reactionDef;
            reaction->getParameters(reactionDef);
            if (!m_skip_user_defined) {
                for (const auto& item : reaction->input) {
                    if (!reactionDef.hasKey(item.first)) {
                        reactionDef[item.first] = item.second;
                    }
                }
            }

            reactions.push_back(std::move(reactionDef));
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

}
