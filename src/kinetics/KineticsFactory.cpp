/**
 *  @file KineticsFactory.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/KineticsFactory.h"
#include "cantera/kinetics/GasKinetics.h"
#include "cantera/kinetics/InterfaceKinetics.h"
#include "cantera/kinetics/EdgeKinetics.h"
#include "cantera/kinetics/importKinetics.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/base/xml.h"
#include "cantera/base/stringUtils.h"

using namespace std;

namespace Cantera
{

KineticsFactory* KineticsFactory::s_factory = 0;
std::mutex KineticsFactory::kinetics_mutex;

Kinetics* KineticsFactory::newKinetics(XML_Node& phaseData,
                                       vector<ThermoPhase*> th)
{
    // Look for a child of the XML element phase called "kinetics". It has an
    // attribute name "model". Store the value of that attribute in the variable
    // kintype
    string kintype = phaseData.child("kinetics")["model"];

    // Create a kinetics object of the desired type
    Kinetics* k = newKinetics(kintype);
    // Now that we have the kinetics manager, we can import the reaction
    // mechanism into it.
    importKinetics(phaseData, th, k);

    // Return the pointer to the kinetics manager
    return k;
}

KineticsFactory::KineticsFactory() {
    reg("none", []() { return new Kinetics(); });
    addAlias("none", "Kinetics");
    addAlias("none", "None");
    reg("gas", []() { return new GasKinetics(); });
    addAlias("gas", "gaskinetics");
    addAlias("gas", "Gas");
    reg("surface", []() { return new InterfaceKinetics(); });
    addAlias("surface", "interface");
    addAlias("surface", "Surf");
    addAlias("surface", "surf");
    reg("edge", []() { return new EdgeKinetics(); });
    addAlias("edge", "Edge");
}

Kinetics* KineticsFactory::newKinetics(const string& model)
{
    return create(toLowerCopy(model));
}

unique_ptr<Kinetics> newKinetics(const vector<ThermoPhase*>& phases,
                                 const AnyMap& phaseNode,
                                 const AnyMap& rootNode)
{
    std::string kinType = phaseNode.getString("kinetics", "none");
    kinType = KineticsFactory::factory()->canonicalize(kinType);
    if (kinType == "none") {
        // determine phase with minimum number of dimensions
        size_t nDim = 3;
        for (auto& phase : phases) {
            nDim = std::min(phase->nDim(), nDim);
        }
        // change kinetics type as necessary
        if (nDim == 2) {
            kinType = "surface";
        } else if (nDim == 1) {
            kinType = "edge";
        }
    }

    unique_ptr<Kinetics> kin(KineticsFactory::factory()->newKinetics(kinType));
    for (auto& phase : phases) {
        kin->addPhase(*phase);
    }
    kin->init();
    addReactions(*kin, phaseNode, rootNode);
    return kin;
}

unique_ptr<Kinetics> newKinetics(const std::vector<ThermoPhase*>& phases,
                                 const std::string& filename,
                                 const std::string& phase_name)
{
    size_t dot = filename.find_last_of(".");
    string extension;
    if (dot != npos) {
        extension = toLowerCopy(filename.substr(dot+1));
    }

    if (extension == "yml" || extension == "yaml") {
        AnyMap root = AnyMap::fromYamlFile(filename);
        AnyMap& phaseNode = root["phases"].getMapWhere("name", phase_name);
        return newKinetics(phases, phaseNode, root);
    } else {
        XML_Node* root = get_XML_File(filename);
        XML_Node* xphase = get_XML_NameID("phase", "#"+phase_name, root);
        if (!xphase) {
            throw CanteraError("newKinetics",
                "Couldn't find phase named '{}' in file '{}'.",
                phase_name, filename);
        }
        return unique_ptr<Kinetics>(newKineticsMgr(*xphase, phases));
    }
}

void addReactions(Kinetics& kin, const AnyMap& phaseNode, const AnyMap& rootNode)
{
    kin.skipUndeclaredThirdBodies(
        phaseNode.getBool("skip-undeclared-third-bodies", false));

    // Find sections containing reactions to add
    vector<string> sections, rules;

    if (phaseNode.hasKey("reactions")) {
        const auto& reactionsNode = phaseNode.at("reactions");
        if (reactionsNode.is<string>()) {
            if (rootNode.hasKey("reactions")) {
                // Specification of the rule for adding species from the default
                // 'reactions' section, if it exists
                sections.push_back("reactions");
                rules.push_back(reactionsNode.asString());
            } else if (reactionsNode.asString() != "none") {
                throw InputFileError("addReactions", reactionsNode,
                    "Phase entry implies existence of 'reactions' section "
                    "which does not exist in the current input file.");
            }
        } else if (reactionsNode.is<vector<string>>()) {
            // List of sections from which all species should be added
            for (const auto& item : reactionsNode.as<vector<string>>()) {
                sections.push_back(item);
                rules.push_back("all");
            }
        } else if (reactionsNode.is<vector<AnyMap>>()) {
            // Mapping of rules to apply for each specified section containing
            // reactions
            for (const auto& item : reactionsNode.as<vector<AnyMap>>()) {
                sections.push_back(item.begin()->first);
                rules.push_back(item.begin()->second.asString());
            }
        }
    } else if (kin.kineticsType() != "None") {
        if (!phaseNode.hasKey("kinetics")) {
            // Do nothing - default surface or edge kinetics require separate detection
            // while not adding reactions
        } else if (rootNode.hasKey("reactions")) {
            // Default behavior is to add all reactions from the 'reactions'
            // section, if a 'kinetics' model has been specified
            sections.push_back("reactions");
            rules.push_back("all");
        } else {
            throw InputFileError("addReactions", phaseNode,
                "Phase entry implies existence of 'reactions' section which "
                "does not exist in the current input file. Add the field "
                "'reactions: none' to the phase entry to specify a kinetics "
                "model with no reactions.");
        }
    }

    // Add reactions from each section
    fmt::memory_buffer add_rxn_err;
    for (size_t i = 0; i < sections.size(); i++) {
        if (rules[i] == "all") {
            kin.skipUndeclaredSpecies(false);
        } else if (rules[i] == "declared-species") {
            kin.skipUndeclaredSpecies(true);
        } else if (rules[i] == "none") {
            continue;
        } else {
            throw InputFileError("addReactions", phaseNode.at("reactions"),
                "Unknown rule '{}' for adding species from the '{}' section.",
                rules[i], sections[i]);
        }
        const auto& slash = boost::ifind_last(sections[i], "/");
        if (slash) {
            // specified section is in a different file
            string fileName (sections[i].begin(), slash.begin());
            string node(slash.end(), sections[i].end());
            AnyMap reactions = AnyMap::fromYamlFile(fileName,
                rootNode.getString("__file__", ""));
            for (const auto& R : reactions[node].asVector<AnyMap>()) {
                try {
                    kin.addReaction(newReaction(R, kin), false);
                } catch (CanteraError& err) {
                    fmt_append(add_rxn_err, "{}", err.what());
                }
            }
        } else {
            // specified section is in the current file
            for (const auto& R : rootNode.at(sections[i]).asVector<AnyMap>()) {
                #ifdef NDEBUG
                    try {
                        kin.addReaction(newReaction(R, kin), false);
                    } catch (CanteraError& err) {
                        fmt_append(add_rxn_err, "{}", err.what());
                    }
                #else
                    kin.addReaction(newReaction(R, kin), false);
                #endif
            }
        }
    }

    kin.checkDuplicates();
    if (add_rxn_err.size()) {
        throw CanteraError("addReactions", to_string(add_rxn_err));
    }

    kin.resizeReactions();
}

}
