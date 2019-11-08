/**
 *  @file ReactionFactory.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/Reaction.h"
#include "cantera/kinetics/ReactionFactory.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/base/ctml.h"
#include "cantera/base/AnyMap.h"

namespace Cantera
{

ReactionFactory* ReactionFactory::s_factory = 0;
std::mutex ReactionFactory::reaction_mutex;

ReactionFactory::ReactionFactory()
{
    // register elementary reactions
    reg("elementary", []() { return new ElementaryReaction(); });
    m_synonyms["arrhenius"] = "elementary";
    m_synonyms[""] = "elementary";
    void setupElementaryReaction(ElementaryReaction&, const XML_Node&);
    reg_XML("elementary",
            [](Reaction* R, const XML_Node& node) {
                setupElementaryReaction(*(ElementaryReaction*)R, node);
            });
    void setupElementaryReaction(ElementaryReaction&, const AnyMap&,
                                 const Kinetics&);
    reg_AnyMap("elementary",
               [](Reaction* R, const AnyMap& node, const Kinetics& kin) {
                   setupElementaryReaction(*(ElementaryReaction*)R, node, kin);
               });

    // register three-body reactions
    reg("three-body", []() { return new ThreeBodyReaction(); });
    m_synonyms["threebody"] = "three-body";
    m_synonyms["three_body"] = "three-body";
    void setupThreeBodyReaction(ThreeBodyReaction&, const XML_Node&);
    reg_XML("three-body",
            [](Reaction* R, const XML_Node& node) {
                setupThreeBodyReaction(*(ThreeBodyReaction*)R, node);
            });
    void setupThreeBodyReaction(ThreeBodyReaction&, const AnyMap&,
                                const Kinetics&);
    reg_AnyMap("three-body",
               [](Reaction* R, const AnyMap& node, const Kinetics& kin) {
                   setupThreeBodyReaction(*(ThreeBodyReaction*)R, node, kin);
               });

    // register falloff reactions
    reg("falloff", []() { return new FalloffReaction(); });
    void setupFalloffReaction(FalloffReaction&, const XML_Node&);
    reg_XML("falloff",
            [](Reaction* R, const XML_Node& node) {
                setupFalloffReaction(*(FalloffReaction*)R, node);
            });
    void setupFalloffReaction(FalloffReaction&, const AnyMap&,
                              const Kinetics&);
    reg_AnyMap("falloff",
               [](Reaction* R, const AnyMap& node, const Kinetics& kin) {
                   setupFalloffReaction(*(FalloffReaction*)R, node, kin);
               });

    // register falloff reactions
    reg("chemically-activated", []() { return new ChemicallyActivatedReaction(); });
    m_synonyms["chemact"] = "chemically-activated";
    m_synonyms["chemically_activated"] = "chemically-activated";
    void setupChemicallyActivatedReaction(ChemicallyActivatedReaction&,
                                          const XML_Node&);
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
    m_synonyms["pressure-dependent-arrhenius"] = "pressure-dependent-Arrhenius";
    m_synonyms["plog"] = "pressure-dependent-Arrhenius";
    m_synonyms["pdep_arrhenius"] = "pressure-dependent-Arrhenius";
    void setupPlogReaction(PlogReaction&, const XML_Node&);
    reg_XML("pressure-dependent-Arrhenius",
            [](Reaction* R, const XML_Node& node) {
                setupPlogReaction(*(PlogReaction*)R, node);
            });
    void setupPlogReaction(PlogReaction&, const AnyMap&, const Kinetics&);
    reg_AnyMap("pressure-dependent-Arrhenius",
               [](Reaction* R, const AnyMap& node, const Kinetics& kin) {
                   setupPlogReaction(*(PlogReaction*)R, node, kin);
               });

    // register Chebyshev reactions
    reg("Chebyshev", []() { return new ChebyshevReaction(); });
    m_synonyms["chebyshev"] = "Chebyshev";
    void setupChebyshevReaction(ChebyshevReaction&, const XML_Node&);
    reg_XML("Chebyshev",
            [](Reaction* R, const XML_Node& node) {
                setupChebyshevReaction(*(ChebyshevReaction*)R, node);
            });
    void setupChebyshevReaction(ChebyshevReaction&, const AnyMap&,
                                const Kinetics&);
    reg_AnyMap("Chebyshev",
               [](Reaction* R, const AnyMap& node, const Kinetics& kin) {
                   setupChebyshevReaction(*(ChebyshevReaction*)R, node, kin);
               });

    // register interface reactions
    reg("interface", []() { return new InterfaceReaction(); });
    m_synonyms["surface"] = "interface";
    m_synonyms["edge"] = "interface";
    m_synonyms["global"] = "interface";
    void setupInterfaceReaction(InterfaceReaction&, const XML_Node&);
    reg_XML("interface",
            [](Reaction* R, const XML_Node& node) {
                setupInterfaceReaction(*(InterfaceReaction*)R, node);
            });
    void setupInterfaceReaction(InterfaceReaction&, const AnyMap&,
                                const Kinetics&);
    reg_AnyMap("interface",
               [](Reaction* R, const AnyMap& node, const Kinetics& kin) {
                   setupInterfaceReaction(*(InterfaceReaction*)R, node, kin);
               });

    // register electrochemical reactions
    reg("electrochemical", []() { return new ElectrochemicalReaction(); });
    m_synonyms["butlervolmer_noactivitycoeffs"] = "electrochemical";
    m_synonyms["butlervolmer"] = "electrochemical";
    m_synonyms["surfaceaffinity"] = "electrochemical";
    void setupElectrochemicalReaction(ElectrochemicalReaction&,
                                      const XML_Node&);
    reg_XML("electrochemical",
            [](Reaction* R, const XML_Node& node) {
                setupElectrochemicalReaction(*(ElectrochemicalReaction*)R, node);
            });
    void setupElectrochemicalReaction(ElectrochemicalReaction&,
                                      const AnyMap&, const Kinetics&);
    reg_AnyMap("electrochemical",
               [](Reaction* R, const AnyMap& node, const Kinetics& kin) {
                   setupElectrochemicalReaction(*(ElectrochemicalReaction*)R, node, kin);
               });
}

Reaction* ReactionFactory::newReaction(const std::string& type)
{
    Reaction* R = create(toLowerCopy(type));
    return R;
}

Reaction* ReactionFactory::newReaction(const XML_Node& rxn_node)
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
        R = create(type);
    } catch (CanteraError& err) {
        throw CanteraError("newReaction",
            "Unknown reaction type '" + rxn_node["type"] + "'");
    }
    setup_XML(R->type(), R, rxn_node);
    return R;
}

Reaction* ReactionFactory::newReaction(const AnyMap& node,
                                       const Kinetics& kin)
{
    std::string type = "elementary";
    if (node.hasKey("type")) {
        type = node["type"].asString();
    }

    if (kin.thermo().nDim() < 3) {
        // See if this is an electrochemical reaction
        Reaction testReaction(0);
        parseReactionEquation(testReaction, node["equation"], kin);
        if (isElectrochemicalReaction(testReaction, kin)) {
            type = "electrochemical";
        } else {
            type = "interface";
        }
    }

    Reaction* R;
    try {
        R = create(type);
    } catch (CanteraError& err) {
        throw InputFileError("ReactionFactory::newReaction", node["type"],
            "Unknown reaction type '{}'", type);
    }
    setup_AnyMap(type, R, node, kin);
    return R;
}

unique_ptr<Reaction> newReaction(const std::string& type)
{
    unique_ptr<Reaction> R(ReactionFactory::factory()->newReaction(type));
    return R;
}

unique_ptr<Reaction> newReaction(const XML_Node& rxn_node)
{
    unique_ptr<Reaction> R(ReactionFactory::factory()->newReaction(rxn_node));
    return R;
}

unique_ptr<Reaction> newReaction(const AnyMap& rxn_node,
                                 const Kinetics& kin)
{
    unique_ptr<Reaction> R(ReactionFactory::factory()->newReaction(rxn_node, kin));
    return R;
}

}
