 /**
 *  @file ReactionFactory.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/ThermoPhase.h"
#include "cantera/kinetics/Reaction.h"
#include "ReactionFactory.h"
#include "cantera/kinetics/ReactionRateFactory.h"
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
    reg("reaction", [](const AnyMap& node, const Kinetics& kin) {
        return new Reaction(node, kin);
    });
    addAlias("reaction", "elementary");
    addAlias("reaction", "arrhenius");
    addAlias("reaction", "Arrhenius");
    addAlias("reaction", "");

    // register elementary reactions (old framework)
    reg("elementary-legacy", [](const AnyMap& node, const Kinetics& kin) {
        ElementaryReaction2* R = new ElementaryReaction2();
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
        ThreeBodyReaction2* R = new ThreeBodyReaction2();
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
        FalloffReaction2* R = new ChemicallyActivatedReaction2();
        if (!node.empty()) {
            setupFalloffReaction(*R, node, kin);
        }
        return R;
    });
    addAlias("chemically-activated-legacy", "chemact");
    addAlias("chemically-activated-legacy", "chemically_activated");

    addAlias("reaction", "pressure-dependent-Arrhenius");
    addAlias("reaction", "plog");
    addAlias("reaction", "pdep_arrhenius");

    // register pressure-dependent-Arrhenius reactions (old framework)
    reg("pressure-dependent-Arrhenius-legacy", [](const AnyMap& node, const Kinetics& kin) {
        PlogReaction2* R = new PlogReaction2();
        if (!node.empty()) {
            setupPlogReaction(*R, node, kin);
        }
        return R;
    });

    // register Chebyshev reactions
    addAlias("reaction", "Chebyshev");
    addAlias("reaction", "chebyshev");
    reg("Chebyshev-legacy", [](const AnyMap& node, const Kinetics& kin) {
        ChebyshevReaction2* R = new ChebyshevReaction2();
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
    reg("interface-legacy", [](const AnyMap& node, const Kinetics& kin) {
        InterfaceReaction2* R = new InterfaceReaction2();
        if (!node.empty()) {
            setupInterfaceReaction(*R, node, kin);
        }
        return R;
    });
    addAlias("interface-legacy", "interface");
    addAlias("interface-legacy", "surface");
    addAlias("interface-legacy", "edge");

    addAlias("reaction", "interface-Arrhenius");
    addAlias("reaction", "sticking-Arrhenius");

    // register electrochemical reactions
    reg("electrochemical-legacy", [](const AnyMap& node, const Kinetics& kin) {
        ElectrochemicalReaction2* R = new ElectrochemicalReaction2();
        if (!node.empty()) {
            setupElectrochemicalReaction(*R, node, kin);
        }
        return R;
    });
    addAlias("electrochemical-legacy", "electrochemical");

    addAlias("reaction", "two-temperature-plasma");

    addAlias("reaction", "Blowers-Masel");
    addAlias("reaction", "interface-Blowers-Masel");
    addAlias("reaction", "sticking-Blowers-Masel");
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
        Reaction* R = new ChemicallyActivatedReaction2();
        setupChemicallyActivatedReaction(*(ChemicallyActivatedReaction2*)R, node);
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
    reg("interface-legacy", [](const XML_Node& node) {
        Reaction* R = new InterfaceReaction2();
        setupInterfaceReaction(*(InterfaceReaction2*)R, node);
        return R;
    });
    addAlias("interface-legacy", "interface");
    addAlias("interface-legacy", "surface");
    addAlias("interface-legacy", "edge");

    // register electrochemical reactions
    reg("electrochemical-legacy", [](const XML_Node& node) {
        Reaction* R = new ElectrochemicalReaction2();
        setupElectrochemicalReaction(*(ElectrochemicalReaction2*)R, node);
        return R;
    });
    addAlias("electrochemical-legacy", "electrochemical");
}

}
