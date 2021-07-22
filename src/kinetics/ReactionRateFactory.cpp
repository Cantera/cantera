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
    reg("ArrheniusRate", [](const AnyMap& node, const Units& rate_units) {
        return new ArrheniusRate(node, rate_units);
    });
    addAlias("ArrheniusRate", "");
    addAlias("ArrheniusRate", "elementary");
    addAlias("ArrheniusRate", "three-body");

    // PlogRate evaluator
    reg("PlogRate", [](const AnyMap& node, const Units& rate_units) {
        return new PlogRate(node, rate_units);
    });
    addAlias("PlogRate", "pressure-dependent-Arrhenius");

    // ChebyshevRate evaluator
    reg("ChebyshevRate", [](const AnyMap& node, const Units& rate_units) {
        return new ChebyshevRate3(node, rate_units);
    });
    addAlias("ChebyshevRate", "Chebyshev");

    // CustomFunc1Rate evaluator
    reg("custom-function", [](const AnyMap& node, const Units& rate_units) {
        return new CustomFunc1Rate(node, rate_units);
    });
    addAlias("custom-function", "custom-rate-function");
}

shared_ptr<ReactionRateBase> newReactionRate(const std::string& type)
{
    return shared_ptr<ReactionRateBase> (
        ReactionRateFactory::factory()->create(type, AnyMap(), Units(0.0)));
}

shared_ptr<ReactionRateBase> newReactionRate(
    const AnyMap& rate_node, const Units& rate_units)
{
    std::string type = "";
    if (rate_node.empty()) {
        throw InputFileError("ReactionRateFactory::newReactionRate", rate_node,
            "Received invalid empty node.");
    } else if (rate_node.hasKey("type")) {
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
    UnitSystem system = rate_node.units();
    if (system.convertTo(1., "m") != 1. || system.convertTo(1., "kmol") != 1.) {
        throw InputFileError("ReactionRateFactory::newReactionRate",
            rate_node.at("__units__"),
            "Alternative units for 'length' or 'quantity` are not supported "
            "when creating\na standalone 'ReactionRate' object.");
    }
    AnyMap node(rate_node);
    node["__standalone__"] = true;
    return newReactionRate(node, Units(1.));
}

std::string canonicalRateName(const std::string& type)
{
    if (ReactionRateFactory::factory()->exists(type)) {
        return ReactionRateFactory::factory()->canonicalize(type);
    }

    throw CanteraError("ReactionRateFactory::canonicalRateName",
        "Unknown reaction rate type alias '{}'.", type);
}

}
