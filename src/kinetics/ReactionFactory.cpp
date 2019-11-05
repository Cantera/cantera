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
    reg("elementary", []() { return new ElementaryReaction(); });
    m_synonyms["arrhenius"] = "elementary";
    m_synonyms[""] = "elementary";
    reg("three-body", []() { return new ThreeBodyReaction(); });
    m_synonyms["threebody"] = "three-body";
    m_synonyms["three_body"] = "three-body";
    reg("falloff", []() { return new FalloffReaction(); });
    reg("chemically-activated", []() { return new ChemicallyActivatedReaction(); });
    m_synonyms["chemact"] = "chemically-activated";
    m_synonyms["chemically_activated"] = "chemically-activated";
    reg("pressure-dependent-arrhenius", []() { return new PlogReaction(); });
    m_synonyms["plog"] = "pressure-dependent-arrhenius";
    m_synonyms["pdep_arrhenius"] = "pressure-dependent-arrhenius";
    reg("chebyshev", []() { return new ChebyshevReaction(); });
    reg("interface", []() { return new InterfaceReaction(); });
    m_synonyms["surface"] = "interface";
    m_synonyms["edge"] = "interface";
    m_synonyms["global"] = "interface";
    reg("electrochemical", []() { return new ElectrochemicalReaction(); });
    m_synonyms["butlervolmer_noactivitycoeffs"] = "electrochemical";
    m_synonyms["butlervolmer"] = "electrochemical";
    m_synonyms["surfaceaffinity"] = "electrochemical";
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
    R->setup(rxn_node);
    return R;
}

Reaction* ReactionFactory::newReaction(const AnyMap& node,
                                       const Kinetics& kin)
{
    std::string type = "elementary";
    if (node.hasKey("type")) {
        type = toLowerCopy(node["type"].asString());
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
        throw InputFileError("newReaction", node["type"],
            "Unknown reaction type '{}'", type);
    }
    R->setup(node, kin);
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
