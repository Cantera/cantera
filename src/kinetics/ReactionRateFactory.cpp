 /**
 *  @file ReactionRateFactory.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/ReactionRateFactory.h"
#include "cantera/kinetics/MultiRate.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/base/AnyMap.h"

namespace Cantera
{

ReactionRateFactory* ReactionRateFactory::s_factory = 0;
std::mutex ReactionRateFactory::rate_mutex;

ReactionRateFactory::ReactionRateFactory()
{
    // ArrheniusRate evaluator
    reg("Arrhenius", [](const AnyMap& node, const Units& rate_units) {
        return new ArrheniusRate(node, rate_units);
    });
    addAlias("Arrhenius", "");
    addAlias("Arrhenius", "elementary");
    addAlias("Arrhenius", "three-body");

    // BlowersMaselRate evaluator
    reg("Blowers-Masel", [](const AnyMap& node, const Units& rate_units) {
        return new BlowersMaselRate(node, rate_units);
    });

    // PlogRate evaluator
    reg("pressure-dependent-Arrhenius", [](const AnyMap& node, const Units& rate_units) {
        return new PlogRate(node, rate_units);
    });

    // ChebyshevRate evaluator
    reg("Chebyshev", [](const AnyMap& node, const Units& rate_units) {
        return new ChebyshevRate3(node, rate_units);
    });

    // CustomFunc1Rate evaluator
    reg("custom-rate-function", [](const AnyMap& node, const Units& rate_units) {
        return new CustomFunc1Rate(node, rate_units);
    });
}

shared_ptr<ReactionRateBase> newReactionRate(const std::string& type)
{
    return shared_ptr<ReactionRateBase> (
        ReactionRateFactory::factory()->create(type, AnyMap(), Units(0.0)));
}

shared_ptr<ReactionRateBase> newReactionRate(
    const AnyMap& rate_node, const Units& rate_units)
{
    std::string type = ""; // default is to create Arrhenius from empty
    if (rate_node.hasKey("type")) {
        type = rate_node["type"].asString();
    }

    if (!(ReactionRateFactory::factory()->exists(type))) {
        throw InputFileError("ReactionRateFactory::newReactionRate", rate_node,
            "Unknown reaction rate type '{}'", type);
    }

    return shared_ptr<ReactionRateBase> (
        ReactionRateFactory::factory()->create(type, rate_node, rate_units));
}

shared_ptr<ReactionRateBase> newReactionRate(const AnyMap& rate_node)
{
    const UnitSystem& system = rate_node.units();
    if (system.convertTo(1., "m") != 1. || system.convertTo(1., "kmol") != 1.) {
        throw InputFileError("ReactionRateFactory::newReactionRate",
            rate_node.at("__units__"),
            "Alternative units for 'length' or 'quantity` are not supported "
            "when creating\na standalone 'ReactionRate' object.");
    }
    AnyMap node(rate_node);
    return newReactionRate(node, Units(0.));
}

}
