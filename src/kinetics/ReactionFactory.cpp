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
    reg("elementary", [](const AnyMap& node, const Kinetics& kin) {
        return new ElementaryReaction3(node, kin);
    });
    addAlias("elementary", "arrhenius");
    addAlias("elementary", "");

    // register elementary reactions (old framework)
    reg("elementary-legacy", [](const AnyMap& node, const Kinetics& kin) {
        ElementaryReaction* R = new ElementaryReaction2();
        if (!node.empty()) {
            setupElementaryReaction(*R, node, kin);
        }
        return R;
    });

    // register three-body reactions
    reg("three-body", [](const AnyMap& node, const Kinetics& kin) {
        return new ThreeBodyReaction3(node, kin);
    });
    addAlias("three-body", "threebody");
    addAlias("three-body", "three_body");

    // register three-body reactions (old framework)
    reg("three-body-legacy", [](const AnyMap& node, const Kinetics& kin) {
        ThreeBodyReaction* R = new ThreeBodyReaction2();
        if (!node.empty()) {
            setupThreeBodyReaction(*R, node, kin);
        }
        return R;
    });

    // register falloff reactions
    reg("falloff", [](const AnyMap& node, const Kinetics& kin) {
        return new FalloffReaction3(node, kin);
    });
    addAlias("falloff", "chemically-activated");

    // register falloff reactions (old framework)
    reg("falloff-legacy", [](const AnyMap& node, const Kinetics& kin) {
        FalloffReaction2* R = new FalloffReaction2();
        if (!node.empty()) {
            setupFalloffReaction(*R, node, kin);
        }
        return R;
    });

    // register falloff reactions
    reg("chemically-activated-legacy", [](const AnyMap& node, const Kinetics& kin) {
        FalloffReaction2* R = new ChemicallyActivatedReaction();
        if (!node.empty()) {
            setupFalloffReaction(*R, node, kin);
        }
        return R;
    });
    addAlias("chemically-activated-legacy", "chemact");
    addAlias("chemically-activated-legacy", "chemically_activated");

    // register pressure-dependent-Arrhenius reactions
    reg("pressure-dependent-Arrhenius", [](const AnyMap& node, const Kinetics& kin) {
        return new PlogReaction3(node, kin);
    });
    addAlias("pressure-dependent-Arrhenius", "plog");
    addAlias("pressure-dependent-Arrhenius", "pdep_arrhenius");

    // register pressure-dependent-Arrhenius reactions (old framework)
    reg("pressure-dependent-Arrhenius-legacy", [](const AnyMap& node, const Kinetics& kin) {
        PlogReaction* R = new PlogReaction2();
        if (!node.empty()) {
            setupPlogReaction(*R, node, kin);
        }
        return R;
    });

    // register Chebyshev reactions
    reg("Chebyshev", [](const AnyMap& node, const Kinetics& kin) {
        return new ChebyshevReaction3(node, kin);
    });
    addAlias("Chebyshev", "chebyshev");
    reg("Chebyshev-legacy", [](const AnyMap& node, const Kinetics& kin) {
        ChebyshevReaction* R = new ChebyshevReaction2();
        if (!node.empty()) {
            setupChebyshevReaction(*R, node, kin);
        }
        return R;
    });

    // register custom reactions specified by Func1 objects
    reg("custom-rate-function", [](const AnyMap& node, const Kinetics& kin) {
        return new CustomFunc1Reaction(node, kin);
    });

    // register interface reactions
    reg("interface", [](const AnyMap& node, const Kinetics& kin) {
        InterfaceReaction* R = new InterfaceReaction();
        if (!node.empty()) {
            setupInterfaceReaction(*R, node, kin);
        }
        return R;
    });
    addAlias("interface", "surface");
    addAlias("interface", "edge");

    // register electrochemical reactions
    reg("electrochemical", [](const AnyMap& node, const Kinetics& kin) {
        ElectrochemicalReaction* R = new ElectrochemicalReaction();
        if (!node.empty()) {
            setupElectrochemicalReaction(*R, node, kin);
        }
        return R;
    });

    // register Blowers Masel reactions
    reg("Blowers-Masel", [](const AnyMap& node, const Kinetics& kin) {
        return new BlowersMaselReaction(node, kin);
    });

    // register surface Blowers Masel reactions
    reg("surface-Blowers-Masel", [](const AnyMap& node, const Kinetics& kin) {
        BlowersMaselInterfaceReaction* R = new BlowersMaselInterfaceReaction();
        if (!node.empty()) {
            setupBlowersMaselInterfaceReaction(*R, node, kin);
        }
        return R;
    });
}

ReactionFactoryXML* ReactionFactoryXML::s_factory = 0;
std::mutex ReactionFactoryXML::reaction_mutex;

ReactionFactoryXML::ReactionFactoryXML()
{
    // register elementary reactions
    reg("elementary-legacy", [](const XML_Node& node) {
        Reaction* R = new ElementaryReaction2();
        setupElementaryReaction(*(ElementaryReaction2*)R, node);
        return R;
    });
    addAlias("elementary-legacy", "elementary");
    addAlias("elementary-legacy", "arrhenius");
    addAlias("elementary-legacy", "");

    // register three-body reactions
    reg("three-body-legacy", [](const XML_Node& node) {
        Reaction* R = new ThreeBodyReaction2();
        setupThreeBodyReaction(*(ThreeBodyReaction2*)R, node);
        return R;
    });
    addAlias("three-body-legacy", "three-body");
    addAlias("three-body-legacy", "threebody");
    addAlias("three-body-legacy", "three_body");

    // register falloff reactions
    reg("falloff-legacy", [](const XML_Node& node) {
        Reaction* R = new FalloffReaction2();
        setupFalloffReaction(*(FalloffReaction2*)R, node);
        return R;
    });
    addAlias("falloff-legacy", "falloff");

    // register falloff reactions
    reg("chemically-activated-legacy", [](const XML_Node& node) {
        Reaction* R = new ChemicallyActivatedReaction();
        setupChemicallyActivatedReaction(*(ChemicallyActivatedReaction*)R, node);
        return R;
    });
    addAlias("chemically-activated-legacy", "chemically-activated");
    addAlias("chemically-activated-legacy", "chemact");
    addAlias("chemically-activated-legacy", "chemically_activated");

    // register pressure-depdendent-Arrhenius reactions
    reg("pressure-dependent-Arrhenius-legacy", [](const XML_Node& node) {
        Reaction* R = new PlogReaction2();
        setupPlogReaction(*(PlogReaction2*)R, node);
        return R;
    });
    addAlias("pressure-dependent-Arrhenius-legacy", "pressure-dependent-Arrhenius");
    addAlias("pressure-dependent-Arrhenius-legacy", "plog");
    addAlias("pressure-dependent-Arrhenius-legacy", "pdep_arrhenius");

    // register Chebyshev reactions
    reg("Chebyshev-legacy", [](const XML_Node& node) {
        Reaction* R = new ChebyshevReaction2();
        setupChebyshevReaction(*(ChebyshevReaction2*)R, node);
        return R;
    });
    addAlias("Chebyshev-legacy", "chebyshev");

    // register interface reactions
    reg("interface", [](const XML_Node& node) {
        Reaction* R = new InterfaceReaction();
        setupInterfaceReaction(*(InterfaceReaction*)R, node);
        return R;
    });
    addAlias("interface", "surface");
    addAlias("interface", "edge");

    // register electrochemical reactions
    reg("electrochemical", [](const XML_Node& node) {
        Reaction* R = new ElectrochemicalReaction();
        setupElectrochemicalReaction(*(ElectrochemicalReaction*)R, node);
        return R;
    });
}

bool isThreeBody(const Reaction& R)
{
    // detect explicitly specified collision partner
    size_t found = 0;
    for (const auto& reac : R.reactants) {
        auto prod = R.products.find(reac.first);
        if (prod != R.products.end() &&
            trunc(reac.second) == reac.second && trunc(prod->second) == prod->second) {
            // candidate species with integer stoichiometric coefficients on both sides
            found += 1;
        }
    }
    if (found != 1) {
        return false;
    }

    // ensure that all reactants have integer stoichiometric coefficients
    size_t nreac = 0;
    for (const auto& reac : R.reactants) {
       if (trunc(reac.second) != reac.second) {
           return false;
       }
       nreac += reac.second;
    }

    // ensure that all products have integer stoichiometric coefficients
    size_t nprod = 0;
    for (const auto& prod : R.products) {
       if (trunc(prod.second) != prod.second) {
           return false;
       }
       nprod += prod.second;
    }

    // either reactant or product side involves exactly three species
    return (nreac == 3) || (nprod == 3);
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
    AnyMap rxn_node;
    Kinetics kin;
    unique_ptr<Reaction> R(ReactionFactory::factory()->create(type, rxn_node, kin));
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

    if (!(ReactionFactoryXML::factory()->exists(type))) {
        throw CanteraError("newReaction",
            "Unknown reaction type '" + rxn_node["type"] + "'");
    }
    if (type.empty()) {
        // Reaction type is not specified
        // See if this is a three-body reaction with a specified collision partner
        ElementaryReaction2 testReaction;
        setupReaction(testReaction, rxn_node);
        if (isThreeBody(testReaction)) {
            type = "three-body";
        }
    }
    Reaction* R = ReactionFactoryXML::factory()->create(type, rxn_node);
    return unique_ptr<Reaction>(R);
}

unique_ptr<Reaction> newReaction(const AnyMap& rxn_node, const Kinetics& kin)
{
    std::string type = "elementary";
    size_t nDim = kin.thermo(kin.reactionPhaseIndex()).nDim();
    if (rxn_node.hasKey("type")) {
        type = rxn_node["type"].asString();
    } else if (nDim == 3) {
        // Reaction type is not specified
        // See if this is a three-body reaction with a specified collision partner
        ElementaryReaction2 testReaction;
        parseReactionEquation(testReaction, rxn_node["equation"], kin);
        if (isThreeBody(testReaction)) {
            type = "three-body";
        }
    }

    if (nDim < 3 && type == "elementary") {
        // See if this is an electrochemical reaction: type of
        // receiving reaction object is unimportant in this case
        ElementaryReaction2 testReaction;
        parseReactionEquation(testReaction, rxn_node["equation"], kin);
        if (isElectrochemicalReaction(testReaction, kin)) {
            type = "electrochemical";
        } else {
            type = "interface";
        }
    }
    if (nDim < 3 && type == "Blowers-Masel") {
        // Allow yaml file to specify "Blowers-Masel" for surface reactions
        type = "surface-Blowers-Masel";
    }

    if (!(ReactionFactory::factory()->exists(type))) {
        throw InputFileError("ReactionFactory::newReaction", rxn_node["type"],
            "Unknown reaction type '{}'", type);
    }
    Reaction* R = ReactionFactory::factory()->create(type, rxn_node, kin);
    return unique_ptr<Reaction>(R);
}

}
