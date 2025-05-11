/**
 *  @file KineticsFactory.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/KineticsFactory.h"
#include "cantera/kinetics/BulkKinetics.h"
#include "cantera/kinetics/InterfaceKinetics.h"
#include "cantera/kinetics/EdgeKinetics.h"
#include "cantera/kinetics/Reaction.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/Solution.h"
#include <boost/algorithm/string.hpp>

namespace Cantera
{

KineticsFactory* KineticsFactory::s_factory = 0;
std::mutex KineticsFactory::kinetics_mutex;

KineticsFactory::KineticsFactory() {
    reg("none", []() { return new Kinetics(); });
    addDeprecatedAlias("none", "Kinetics");
    addDeprecatedAlias("none", "None");
    reg("bulk", []() { return new BulkKinetics(); });
    addAlias("bulk", "gas");
    addDeprecatedAlias("bulk", "gaskinetics");
    addDeprecatedAlias("bulk", "Gas");
    reg("surface", []() { return new InterfaceKinetics(); });
    addAlias("surface", "interface");
    addDeprecatedAlias("surface", "Surf");
    addDeprecatedAlias("surface", "surf");
    reg("edge", []() { return new EdgeKinetics(); });
    addDeprecatedAlias("edge", "Edge");
}

KineticsFactory* KineticsFactory::factory() {
    std::unique_lock<std::mutex> lock(kinetics_mutex);
    if (!s_factory) {
        s_factory = new KineticsFactory;
    }
    return s_factory;
}

void KineticsFactory::deleteFactory() {
    std::unique_lock<std::mutex> lock(kinetics_mutex);
    delete s_factory;
    s_factory = 0;
}

Kinetics* KineticsFactory::newKinetics(const string& model)
{
    return create(toLowerCopy(model));
}

shared_ptr<Kinetics> newKinetics(const string& model)
{
    shared_ptr<Kinetics> kin(KineticsFactory::factory()->newKinetics(model));
    return kin;
}

shared_ptr<Kinetics> newKinetics(const vector<shared_ptr<ThermoPhase>>& phases,
                                 const AnyMap& phaseNode,
                                 const AnyMap& rootNode,
                                 shared_ptr<Solution> soln)
{
    string kinType = phaseNode.getString("kinetics", "none");
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

    shared_ptr<Kinetics> kin;
    if (soln && soln->thermo() && soln->thermo()->kinetics()) {
        // If kinetics was initiated in thermo already, use it directly
        kin = soln->thermo()->kinetics();
    } else {
        // Otherwise, create a new kinetics
        kin = std::shared_ptr<Kinetics>(KineticsFactory::factory()->newKinetics(kinType));
    }

    if (soln) {
        soln->setKinetics(kin);
    }
    for (auto& phase : phases) {
        kin->addThermo(phase);
    }

    if (!kin->ready()) {
        kin->init();
        addReactions(*kin, phaseNode, rootNode);
    }
    return kin;
}

shared_ptr<Kinetics> newKinetics(const vector<shared_ptr<ThermoPhase>>& phases,
                                 const string& filename)
{
    string reaction_phase = phases.at(0)->name();
    AnyMap root = AnyMap::fromYamlFile(filename);
    AnyMap& phaseNode = root["phases"].getMapWhere("name", reaction_phase);
    return newKinetics(phases, phaseNode, root);
}

vector<AnyMap> reactionsAnyMapList(Kinetics& kin, const AnyMap& phaseNode,
                                   const AnyMap& rootNode)
{
    kin.skipUndeclaredThirdBodies(
        phaseNode.getBool("skip-undeclared-third-bodies", false));
    kin.setExplicitThirdBodyDuplicateHandling(
        phaseNode.getString("explicit-third-body-duplicates", "warn"));

    loadExtensions(rootNode);

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
    } else if (kin.kineticsType() != "none") {
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
    vector<AnyMap> reactionsList;
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
            loadExtensions(reactions);

            for (const auto& R : reactions[node].asVector<AnyMap>()) {
                reactionsList.push_back(R);
            }
        } else {
            // specified section is in the current file
            for (const auto& R : rootNode.at(sections[i]).asVector<AnyMap>()) {
                reactionsList.push_back(R);
            }
        }
    }
    return reactionsList;
}

void addReactions(Kinetics& kin, vector<shared_ptr<Reaction>> rxnList)
{
    fmt::memory_buffer add_rxn_err;
    for (shared_ptr<Reaction> rxn : rxnList) {
        #ifdef NDEBUG
            try {
                kin.addReaction(rxn, false);
            } catch (CanteraError& err) {
                fmt_append(add_rxn_err, "{}", err.what());
            }
        #else
            kin.addReaction(rxn, false);
        #endif
    }

    if (add_rxn_err.size()) {
        throw CanteraError("addReactions", to_string(add_rxn_err));
    }
    kin.checkDuplicates();
    kin.resizeReactions();
}

void addReactions(Kinetics& kin, const AnyMap& phaseNode, const AnyMap& rootNode)
{
    fmt::memory_buffer add_rxn_err;
    for (AnyMap R : reactionsAnyMapList(kin, phaseNode, rootNode)) {
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

    if (add_rxn_err.size()) {
        throw CanteraError("addReactions", to_string(add_rxn_err));
    }
    kin.checkDuplicates();
    kin.resizeReactions();
}

}
