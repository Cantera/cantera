/**
 *  @file PlasmaElectronFactory.cpp
 *     Definitions for the factory class that can create known PlasmaElectron objects
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://www.cantera.org/license.txt for license and copyright information.

#include "cantera/plasma/PlasmaElectronFactory.h"
#include "cantera/plasma/PlasmaElectron.h"
#include "cantera/plasma/WeakIonGasElectron.h"


using namespace std;

namespace Cantera
{

PlasmaElectronFactory* PlasmaElectronFactory::s_factory = 0;
std::mutex PlasmaElectronFactory::electron_mutex;

PlasmaElectronFactory::PlasmaElectronFactory()
{
    reg("WeaklyIonizedGas", []() { return new WeakIonGasElectron(); });
}

PlasmaElectron* PlasmaElectronFactory::newPlasmaElectron(const std::string& model)
{
    return create(model);
}

unique_ptr<PlasmaElectron> newPlasmaElectron(const AnyMap& phaseNode, const AnyMap& rootNode, thermo_t* phase)
{
    unique_ptr<PlasmaElectron> electron(newPlasmaElectron(phaseNode["plasma"].asString()));
    if (phaseNode.hasKey("cross-sections")) {
        if (phaseNode["cross-sections"].is<vector<AnyMap>>()) {
            // Each item in 'cross-sections' is a map with one item, where the key is
            // a section in another YAML file, and the value is a
            // list of target species names to read from that section 
            for (const auto& crossSectionsNode : phaseNode["cross-sections"].asVector<AnyMap>()) {
                const string& source = crossSectionsNode.begin()->first;
                const auto& names = crossSectionsNode.begin()->second;
                const auto& slash = boost::ifind_last(source, "/");
                if (slash) {
                    // source is a different input file
                    std::string fileName(source.begin(), slash.begin());
                    std::string node(slash.end(), source.end());
                    AnyMap crossSections = AnyMap::fromYamlFile(fileName);
                    addElectronCrossSections(*electron, crossSections[node], names);
                } else {
                    throw InputFileError("newPlasmaElectron", crossSectionsNode,
                        "Could not find species section named '{}'", source);
                }
            }
        } else {
            throw InputFileError("newPlasmaElectron", phaseNode["cross-sections"],
                "Could not parse cross-sections declaration of type '{}'",
                phaseNode["cross-sections"].type_str());
        }
    }
    electron->init(phase);
    return electron;
}

void addElectronCrossSections(PlasmaElectron& electron, const AnyValue& crossSections, const AnyValue& names)
{
    if (names.is<vector<string>>()) {
        // 'names' is a list of target species names which should be found in 'cross-sections'
        for (const auto& name : names.asVector<string>()) {
            for (const auto& item : crossSections.asVector<AnyMap>()) {
                if (item["target"].asString() == name) {
                    electron.addElectronCrossSection(newElectronCrossSection(item));
                }
            }
        }
    } else if (names.is<string>() && names.asString() == "all") {
        // The keyword 'all' means to add all cross-sections from this source
        for (const auto& item : crossSections.asVector<AnyMap>()) {
            electron.addElectronCrossSection(newElectronCrossSection(item));
        }
    } else {
        throw InputFileError("addElectronCrossSections", names,
            "Could not parse cross-sections declaration of type '{}'", names.type_str());
    }
}

}
