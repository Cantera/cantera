// This file is part of Cantera. See License.txt in the top-level directory or
// at https://www.cantera.org/license.txt for license and copyright information.

#include "cantera/plasma/ElectronCrossSection.h"

namespace Cantera {

ElectronCrossSection::ElectronCrossSection()
    : threshold(0.0)
{
}

ElectronCrossSection::~ElectronCrossSection()
{
}

void ElectronCrossSection::validate()
{
    if (kind == "effective") {
        if (data.size() > 0 && data[0].size() > 0) {
            if (data[0][0] != 0.0) {
                throw CanteraError("ElectronCrossSection::validate",
                    "Invalid energy value of type '{}' for '{}'. "
                    "Energy must starts at zero.", kind, target);
            }
        }
    } else if (kind == "elastic") {
        if (data[0][0] != 0.0) {
            throw CanteraError("ElectronCrossSection::validate",
                "Invalid energy value of type '{}' for '{}'. "
                "Energy must starts at zero.", kind, target);
        }
    } else if (kind != "ionization" && kind != "attachment" && kind != "excitation"){
        throw CanteraError("ElectronCrossSection::validate",
            "'{}' is an unknown type of cross section data.", kind);
    }
}

unique_ptr<ElectronCrossSection> newElectronCrossSection(const AnyMap& node)
{
    unique_ptr<ElectronCrossSection> ecs(new ElectronCrossSection());

    ecs->kind = node["kind"].asString();
    ecs->target = node["target"].asString();
    ecs->data = node["data"].asVector<vector_fp>();
    if (ecs->kind != "effective" && ecs->kind != "elastic") {
        ecs->threshold = node["threshold"].asDouble();
        ecs->product = node["product"].asString();
    }

    // Store all unparsed keys in the "extra" map
    const static std::set<std::string> known_keys{
        "kind", "target", "product", "data", "threshold"
    };

    for (const auto& item : node) {
        if (known_keys.count(item.first) == 0) {
            ecs->extra[item.first] = item.second;
        }
    }

    return ecs;
}

std::vector<shared_ptr<ElectronCrossSection>> getElectronCrossSection(const AnyValue& items)
{
    std::vector<shared_ptr<ElectronCrossSection> > all_cross_sections;
    for (const auto& node : items.asVector<AnyMap>()) {
        all_cross_sections.emplace_back(newElectronCrossSection(node));
    }
    return all_cross_sections;
}

}
