// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/YamlWriter.h"
#include "cantera/base/AnyMap.h"
#include "cantera/base/Solution.h"
#include "cantera/base/stringUtils.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/Species.h"

#include <fstream>
#include <chrono>

namespace Cantera {

void YamlWriter::addPhase(shared_ptr<Solution> soln) {
    for (auto& phase : m_phases) {
        if (phase->name() == soln->name()) {
            throw CanteraError("YamlWriter::addPhase",
                "Duplicate phase name '{}'", soln->name());
        }
    }
    m_phases.push_back(soln);
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
    }
    output["phases"] = phaseDefs;

    // Build species definitions for all phases
    std::vector<AnyMap> speciesDefs;
    speciesDefs.reserve(nspecies_total);
    std::unordered_map<std::string, size_t> speciesDefIndex;
    for (const auto& phase : m_phases) {
        for (const auto& name : phase->thermo()->speciesNames()) {
            const auto& species = phase->thermo()->species(name);
            AnyMap speciesDef;
            species->getParameters(speciesDef);
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

    return output.toYamlString();
}

void YamlWriter::toYamlFile(const std::string& filename) const
{
    std::ofstream out(filename);
    out << toYamlString();
}

}
