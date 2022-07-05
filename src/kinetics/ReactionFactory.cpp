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

    // register three-body reactions
    reg("three-body", [](const AnyMap& node, const Kinetics& kin) {
        return new ThreeBodyReaction(node, kin);
    });
    addAlias("three-body", "threebody");
    addAlias("three-body", "three_body");

    // register falloff reactions
    reg("falloff", [](const AnyMap& node, const Kinetics& kin) {
        return new FalloffReaction(node, kin);
    });
    addAlias("falloff", "chemically-activated");
    addAlias("falloff", "chemact");
    addAlias("falloff", "chemically_activated");

    addAlias("reaction", "pressure-dependent-Arrhenius");
    addAlias("reaction", "plog");
    addAlias("reaction", "pdep_arrhenius");

    // register Chebyshev reactions
    addAlias("reaction", "Chebyshev");
    addAlias("reaction", "chebyshev");

    // register custom reactions specified by Func1 objects
    addAlias("reaction", "custom-rate-function");

    addAlias("reaction", "interface-Arrhenius");
    addAlias("reaction", "sticking-Arrhenius");

    // register legacy interface reaction names
    addAlias("reaction", "interface");
    addAlias("reaction", "surface");
    addAlias("reaction", "edge");
    addAlias("reaction", "electrochemical");

    addAlias("reaction", "two-temperature-plasma");

    addAlias("reaction", "Blowers-Masel");
    addAlias("reaction", "interface-Blowers-Masel");
    addAlias("reaction", "sticking-Blowers-Masel");
}

}
