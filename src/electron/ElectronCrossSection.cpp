// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/electron/ElectronCrossSection.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/ctml.h"
#include <iostream>
#include <limits>
#include <set>

namespace Cantera {

ElectronCrossSection::ElectronCrossSection()
{
}

ElectronCrossSection::ElectronCrossSection(const std::string& kind_,
                                           const std::string& target_,
                                           const std::vector<std::vector<double>> data_,
                                           double mass_ratio_,
                                           double threshold_)
    : kind(kind_)
    , target(target_)
    , data(data_)
    , mass_ratio(mass_ratio_)
    , threshold(threshold_)
{
}

ElectronCrossSection::~ElectronCrossSection()
{
}

void ElectronCrossSection::validate()
{
    if (kind == "EFFECTIVE") {
        if (data[0][0] != 0.0) {
            throw CanteraError("ElectronCrossSection::validate",
                "Invalid energy value of type '{}' for '{}'. "
                "Energy must starts at zero.", kind, target);
        }
        if (mass_ratio >= 1.0 || mass_ratio < 0.0) {
            throw CanteraError("ElectronCrossSection::validate",
                "Invalid mass ratio of type '{}' for '{}'. "
                "Mass ratio of electron to target must be in the range of 0 to 1.", kind, target);
        }
    } else if (kind == "ELASTIC") {
        if (data[0][0] != 0.0) {
            throw CanteraError("ElectronCrossSection::validate",
                "Invalid energy value of type '{}' for '{}'. "
                "Energy must starts at zero.", kind, target);
        }
        if (mass_ratio >= 1.0 || mass_ratio < 0.0) {
            throw CanteraError("ElectronCrossSection::validate",
                "Invalid mass ratio of type '{}' for '{}'. "
                "Mass ratio of electron to target must be in the range of 0 to 1.", kind, target);
        }
    } else if (kind == "IONIZATION") {
        if (data[0][1] != 0.0) {
            throw CanteraError("ElectronCrossSection::validate",
                "Invalid cross section value of type '{}' for '{}'. "
                "Cross section must starts at zero.", kind, target);
        }
    } else if (kind == "ATTACHMENT") {
        if (data[0][1] != 0.0) {
            throw CanteraError("ElectronCrossSection::validate",
                "Invalid cross section value of type '{}' for '{}'. "
                "Cross section must starts at zero.", kind, target);
        }
    } else if (kind == "EXCITATION") {
        if (data[0][1] != 0.0) {
            throw CanteraError("ElectronCrossSection::validate",
                "Invalid cross section value of type '{}' for '{}'. "
                "Cross section must starts at zero.", kind, target);
        }
    } else {
        throw CanteraError("ElectronCrossSection::validate",
            "'{}' is an unknown type of cross section data.", kind);
    }
}

unique_ptr<ElectronCrossSection> newElectronCrossSection(const AnyMap& node)
{
    unique_ptr<ElectronCrossSection> ecs(new ElectronCrossSection(node["kind"].asString(),
                                      node["target"].asString(),
                                      node["data"].asVector<std::vector<double>>()));

    if (ecs->kind == "EFFECTIVE") {
        ecs->mass_ratio = node["mass_ratio"].asDouble();
    } else if (ecs->kind != "ATTACHMENT") {
        ecs->threshold = node["threshold"].asDouble();
    }

    // Store all unparsed keys in the "extra" map
    const static std::set<std::string> known_keys{
        "kind", "target", "data", "mass_ratio", "threshold"
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
