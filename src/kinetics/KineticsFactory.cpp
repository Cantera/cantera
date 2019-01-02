/**
 *  @file KineticsFactory.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/KineticsFactory.h"
#include "cantera/kinetics/GasKinetics.h"
#include "cantera/kinetics/InterfaceKinetics.h"
#include "cantera/kinetics/EdgeKinetics.h"
#include "cantera/kinetics/importKinetics.h"
#include "cantera/base/xml.h"

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
    reg("gaskinetics", []() { return new GasKinetics(); });
    m_synonyms["gas"] = "gaskinetics";
    reg("interface", []() { return new InterfaceKinetics(); });
    reg("edge", []() { return new EdgeKinetics(); });
}

Kinetics* KineticsFactory::newKinetics(const string& model)
{
    return create(toLowerCopy(model));
}

unique_ptr<Kinetics> newKinetics(vector<ThermoPhase*>& phases,
                                 const AnyMap& phaseNode,
                                 const AnyMap& rootNode)
{
    unique_ptr<Kinetics> kin(KineticsFactory::factory()->newKinetics(
        phaseNode.getString("kinetics", "none")));
    for (auto& phase : phases) {
        kin->addPhase(*phase);
    }
    kin->init();
    if (kin->kineticsType() != "Kinetics") {
        addReactions(*kin, phaseNode, rootNode);
    }
    return kin;
}

void addReactions(Kinetics& kin, const AnyMap& phaseNode, const AnyMap& rootNode)
{
    // Find sections containing reactions to add
    vector<string> sections, rules;

    if (phaseNode.hasKey("reactions")) {
        const auto& reactionsNode = phaseNode.at("reactions");
        if (reactionsNode.is<string>()) {
            // Specification of the rule for adding species from the default
            // 'reactions' section
            sections.push_back("reactions");
            rules.push_back(reactionsNode.asString());
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
    } else {
        // Default behavior is to add all reactions from the 'reactions' section
        sections.push_back("reactions");
        rules.push_back("all");
    }

    // Add reactions from each section
    for (size_t i = 0; i < sections.size(); i++) {
        if (rules[i] == "all") {
            kin.skipUndeclaredSpecies(false);
            kin.skipUndeclaredThirdBodies(false);
        } else if (rules[i] == "declared-species") {
            kin.skipUndeclaredSpecies(true);
            kin.skipUndeclaredThirdBodies(true);
        } else if (rules[i] != "none") {
            throw CanteraError("setupKinetics", "Unknown rule '{}' for adding "
                "species from the '{}' section.", rules[i], sections[i]);
        }
        const auto& slash = boost::ifind_first(sections[i], "/");
        if (slash) {
            // specified section is in a different file
            string fileName (sections[i].begin(), slash.begin());
            string node(slash.end(), sections[i].end());
            AnyMap reactions = AnyMap::fromYamlFile(fileName);
            for (const auto& R : reactions[node].asVector<AnyMap>()) {
                kin.addReaction(newReaction(R, kin));
            }
        } else {
            // specified section is in the current file
            for (const auto& R : rootNode.at(sections[i]).asVector<AnyMap>()) {
                kin.addReaction(newReaction(R, kin));
            }
        }
    }
}

}
