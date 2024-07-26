// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/Elements.h"
#include "cantera/thermo/Species.h"
#include "cantera/thermo/SpeciesThermoInterpType.h"
#include "cantera/thermo/SpeciesThermoFactory.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/transport/TransportData.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/global.h"
#include <iostream>
#include <limits>

using namespace std;

namespace Cantera {

Species::Species(const string& name_, const Composition& comp_,
                 double charge_, double size_)
    : name(name_)
    , composition(comp_)
    , charge(charge_)
    , size(size_)
{
}

double Species::molecularWeight() {
    if (m_molecularWeight == Undef) {
        double weight = 0.0;
        const auto& elements = elementWeights();
        for (const auto& [eName, stoich] : composition) {
            auto search = elements.find(eName);
            if (search != elements.end()) {
                if (search->second < 0) {
                    throw CanteraError("setMolecularWeight",
                        "element '{}' has no stable isotopes", eName);
                }
                weight += search->second * stoich;
            }
        }
        setMolecularWeight(weight);
    }
    return m_molecularWeight;
}

void Species::setMolecularWeight(double weight) {
    if (m_molecularWeight != Undef) {
        double maxWeight = max(weight, m_molecularWeight);
        double weight_cmp = fabs(weight - m_molecularWeight) / maxWeight;
        if (weight_cmp > 1.0e-9) {
            warn_user(
                "Species::setMolecularWeight",
                "Molecular weight of species '{}' is changing from {} to {}.",
                this->name,
                m_molecularWeight,
                weight
            );
        }
    }

    m_molecularWeight = weight;
}

AnyMap Species::parameters(const ThermoPhase* phase, bool withInput) const
{
    AnyMap speciesNode;
    speciesNode["name"] = name;
    speciesNode["composition"] = composition;
    speciesNode["composition"].setFlowStyle();

    if (charge != 0) {
        speciesNode["charge"] = charge;
    } else {
        speciesNode.exclude("charge");
    }

    if (size != 1) {
        speciesNode["size"] = size;
    } else {
        speciesNode.exclude("size");
    }

    if (thermo) {
        AnyMap thermoNode = thermo->parameters(withInput);
        if (thermoNode.size()) {
            speciesNode["thermo"] = std::move(thermoNode);
        }
    }
    if (transport) {
        speciesNode["transport"] = transport->parameters(withInput);
    }
    if (phase) {
        phase->getSpeciesParameters(name, speciesNode);
    }
    if (withInput && input.hasKey("equation-of-state")) {
        auto& eosIn = input["equation-of-state"].asVector<AnyMap>();
        for (const auto& eos : eosIn) {
            auto& out = speciesNode["equation-of-state"].getMapWhere(
                "model", eos["model"].asString(), true);
            out.update(eos);
        }
    }
    if (withInput) {
        speciesNode.update(input);
    }
    return speciesNode;
}

unique_ptr<Species> newSpecies(const AnyMap& node)
{
    auto s = make_unique<Species>(node["name"].asString(),
                                  node["composition"].asMap<double>());

    if (node.hasKey("thermo")) {
        s->thermo = newSpeciesThermo(node["thermo"].as<AnyMap>());
    } else {
        s->thermo = make_shared<SpeciesThermoInterpType>();
    }

    s->size = node.getDouble("sites", 1.0);
    if (s->composition.find("E") != s->composition.end()) {
        s->charge = -s->composition["E"];
    }

    if (node.hasKey("transport")) {
        s->transport = newTransportData(node["transport"].as<AnyMap>());
        s->transport->validate(*s);
    }

    // Store input parameters in the "input" map, unless they are stored in a
    // child object
    const static set<string> known_keys{
        "thermo", "transport"
    };
    s->input.setUnits(node.units());
    for (const auto& [key, child] : node) {
        if (known_keys.count(key) == 0) {
            s->input[key] = child;
        }
    }
    s->input.applyUnits();
    s->input.copyMetadata(node);

    return s;
}

vector<shared_ptr<Species>> getSpecies(const AnyValue& items)
{
    vector<shared_ptr<Species>> all_species;
    for (const auto& node : items.asVector<AnyMap>()) {
        all_species.emplace_back(newSpecies(node));
    }
    return all_species;
}

}
