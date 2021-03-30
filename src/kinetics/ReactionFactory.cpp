 /**
 *  @file ReactionFactory.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/ThermoPhase.h"
#include "cantera/kinetics/Reaction.h"
#include "cantera/kinetics/ReactionFactory.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/base/ctml.h"
#include "cantera/base/AnyMap.h"
#include "cantera/base/stringUtils.h"

namespace Cantera
{

ReactionFactory* ReactionFactory::s_factory = 0;
std::mutex ReactionFactory::reaction_mutex;

ReactionFactory::ReactionFactory()
{
    // register elementary reactions
    reg("elementary", []() { return new ElementaryReaction(); });
    addAlias("elementary", "arrhenius");
    addAlias("elementary", "");
    reg_XML("elementary",
            [](Reaction* R, const XML_Node& node) {
                setupElementaryReaction(*(ElementaryReaction*)R, node);
            });
    reg_AnyMap("elementary",
               [](Reaction* R, const AnyMap& node, const Kinetics& kin) {
                   setupElementaryReaction(*(ElementaryReaction*)R, node, kin);
               });

    // register three-body reactions
    reg("three-body", []() { return new ThreeBodyReaction(); });
    addAlias("three-body", "threebody");
    addAlias("three-body", "three_body");
    reg_XML("three-body",
            [](Reaction* R, const XML_Node& node) {
                setupThreeBodyReaction(*(ThreeBodyReaction*)R, node);
            });
    reg_AnyMap("three-body",
               [](Reaction* R, const AnyMap& node, const Kinetics& kin) {
                   setupThreeBodyReaction(*(ThreeBodyReaction*)R, node, kin);
               });

    // register falloff reactions
    reg("falloff", []() { return new FalloffReaction(); });
    reg_XML("falloff",
            [](Reaction* R, const XML_Node& node) {
                setupFalloffReaction(*(FalloffReaction*)R, node);
            });
    reg_AnyMap("falloff",
               [](Reaction* R, const AnyMap& node, const Kinetics& kin) {
                   setupFalloffReaction(*(FalloffReaction*)R, node, kin);
               });

    // register falloff reactions
    reg("chemically-activated", []() { return new ChemicallyActivatedReaction(); });
    addAlias("chemically-activated", "chemact");
    addAlias("chemically-activated", "chemically_activated");
    reg_XML("chemically-activated",
            [](Reaction* R, const XML_Node& node) {
                setupChemicallyActivatedReaction(*(ChemicallyActivatedReaction*)R, node);
            });
    reg_AnyMap("chemically-activated",
               [](Reaction* R, const AnyMap& node, const Kinetics& kin) {
                   setupFalloffReaction(*(FalloffReaction*)R, node, kin);
               });

    // register pressure-depdendent-Arrhenius reactions
    reg("pressure-dependent-Arrhenius", []() { return new PlogReaction(); });
    addAlias("pressure-dependent-Arrhenius", "plog");
    addAlias("pressure-dependent-Arrhenius", "pdep_arrhenius");
    reg_XML("pressure-dependent-Arrhenius",
            [](Reaction* R, const XML_Node& node) {
                setupPlogReaction(*(PlogReaction*)R, node);
            });
    reg_AnyMap("pressure-dependent-Arrhenius",
               [](Reaction* R, const AnyMap& node, const Kinetics& kin) {
                   setupPlogReaction(*(PlogReaction*)R, node, kin);
               });

    // register Chebyshev reactions
    reg("Chebyshev", []() { return new ChebyshevReaction(); });
    addAlias("Chebyshev", "chebyshev");
    reg_XML("Chebyshev",
            [](Reaction* R, const XML_Node& node) {
                setupChebyshevReaction(*(ChebyshevReaction*)R, node);
            });
    reg_AnyMap("Chebyshev",
               [](Reaction* R, const AnyMap& node, const Kinetics& kin) {
                   setupChebyshevReaction(*(ChebyshevReaction*)R, node, kin);
               });

    // register custom reactions specified by Func1 objects
    reg("custom-rate-function", []() { return new CustomFunc1Reaction(); });
    reg_XML("custom-rate-function",
            [](Reaction* R, const XML_Node& node) {
                throw CanteraError("ReactionFactory", "Custom reactions based "
                    "on 'Func1' objects cannot be created from XML nodes'");
            });
    reg_AnyMap("custom-rate-function",
               [](Reaction* R, const AnyMap& node, const Kinetics& kin) {
                   ((Reaction2*)R)->setParameters(node, kin);
               });

    // register custom Python reactions
    reg("elementary-new", []() { return new TestReaction(); });
    reg_XML("elementary-new",
            [](Reaction* R, const XML_Node& node) {
                throw CanteraError("ReactionFactory", "Test reactions "
                                   "cannot be created from XML nodes'");
            });
    reg_AnyMap("elementary-new",
               [](Reaction* R, const AnyMap& node, const Kinetics& kin) {
                   (dynamic_cast<Reaction2*>(R))->setParameters(node, kin);
               });

    // register interface reactions
    reg("interface", []() { return new InterfaceReaction(); });
    addAlias("interface", "surface");
    addAlias("interface", "edge");
    reg_XML("interface",
            [](Reaction* R, const XML_Node& node) {
                setupInterfaceReaction(*(InterfaceReaction*)R, node);
            });
    reg_AnyMap("interface",
               [](Reaction* R, const AnyMap& node, const Kinetics& kin) {
                   setupInterfaceReaction(*(InterfaceReaction*)R, node, kin);
               });

    // register electrochemical reactions
    reg("electrochemical", []() { return new ElectrochemicalReaction(); });
    reg_XML("electrochemical",
            [](Reaction* R, const XML_Node& node) {
                setupElectrochemicalReaction(*(ElectrochemicalReaction*)R, node);
            });
    reg_AnyMap("electrochemical",
               [](Reaction* R, const AnyMap& node, const Kinetics& kin) {
                   setupElectrochemicalReaction(*(ElectrochemicalReaction*)R, node, kin);
               });

    // register Blowers Masel reactions
    reg("Blowers-Masel", []() { return new BlowersMaselReaction(); });
    reg_AnyMap("Blowers-Masel",
               [](Reaction* R, const AnyMap& node, const Kinetics& kin) {
                   setupBlowersMaselReaction(*(BlowersMaselReaction*)R, node, kin);
               });

    // register surface Blowers Masel reactions
    reg("surface-Blowers-Masel", []() { return new BlowersMaselInterfaceReaction(); });
    reg_AnyMap("surface-Blowers-Masel",
               [](Reaction* R, const AnyMap& node, const Kinetics& kin) {
                   setupBlowersMaselInterfaceReaction(*(BlowersMaselInterfaceReaction*)R, node, kin);
               });
}

bool isElectrochemicalReaction(Reaction& R, const Kinetics& kin)
{
    vector_fp e_counter(kin.nPhases(), 0.0);

    // Find the number of electrons in the products for each phase
    for (const auto& sp : R.products) {
        size_t kkin = kin.kineticsSpeciesIndex(sp.first);
        size_t i = kin.speciesPhaseIndex(kkin);
        size_t kphase = kin.thermo(i).speciesIndex(sp.first);
        e_counter[i] += sp.second * kin.thermo(i).charge(kphase);
    }

    // Subtract the number of electrons in the reactants for each phase
    for (const auto& sp : R.reactants) {
        size_t kkin = kin.kineticsSpeciesIndex(sp.first);
        size_t i = kin.speciesPhaseIndex(kkin);
        size_t kphase = kin.thermo(i).speciesIndex(sp.first);
        e_counter[i] -= sp.second * kin.thermo(i).charge(kphase);
    }

    // If the electrons change phases then the reaction is electrochemical
    for (double delta_e : e_counter) {
        if (std::abs(delta_e) > 1e-4) {
            return true;
        }
    }
    return false;
}

unique_ptr<Reaction> newReaction(const std::string& type)
{
    unique_ptr<Reaction> R(ReactionFactory::factory()->create(type));
    return R;
}

unique_ptr<Reaction> newReaction(const XML_Node& rxn_node)
{
    std::string type = toLowerCopy(rxn_node["type"]);

    // Modify the reaction type for interface reactions which contain
    // electrochemical reaction data
    if (rxn_node.child("rateCoeff").hasChild("electrochem")
        && (type == "edge" || type == "surface")) {
        type = "electrochemical";
    }

    Reaction* R;
    try {
        R = ReactionFactory::factory()->create(type);
    } catch (CanteraError& err) {
        throw CanteraError("newReaction",
            "Unknown reaction type '" + rxn_node["type"] + "'");
    }
    if (type != "electrochemical") {
        type = R->type();
    }
    ReactionFactory::factory()->setup_XML(type, R, rxn_node);

    return unique_ptr<Reaction>(R);
}

unique_ptr<Reaction> newReaction(const AnyMap& rxn_node,
                                 const Kinetics& kin)
{
    std::string type = "elementary";
    if (rxn_node.hasKey("type")) {
        type = rxn_node["type"].asString();
    }

    if (kin.thermo().nDim() < 3 && type == "elementary") {
        // See if this is an electrochemical reaction: type of
        // receiving reaction object is unimportant in this case
        ElementaryReaction testReaction;
        parseReactionEquation(testReaction, rxn_node["equation"], kin);
        if (isElectrochemicalReaction(testReaction, kin)) {
            type = "electrochemical";
        } else {
            type = "interface";
        }
    }
    if (kin.thermo().nDim() < 3 && type == "Blowers-Masel") {
        // Allow yaml file to specify "Blowers-Masel" for surface reactions
        type = "surface-Blowers-Masel";
    }

    Reaction* R;
    try {
        R = ReactionFactory::factory()->create(type);
    } catch (CanteraError& err) {
        throw InputFileError("ReactionFactory::newReaction", rxn_node["type"],
            "Unknown reaction type '{}'", type);
    }
    ReactionFactory::factory()->setup_AnyMap(type, R, rxn_node, kin);

    return unique_ptr<Reaction>(R);
}

}
