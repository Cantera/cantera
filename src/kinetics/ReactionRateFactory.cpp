 /**
 *  @file ReactionRateFactory.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/ReactionRateFactory.h"
#include "cantera/kinetics/MultiRate.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/kinetics/Falloff.h"
#include "cantera/base/AnyMap.h"

namespace Cantera
{

ReactionRateFactory* ReactionRateFactory::s_factory = 0;
std::mutex ReactionRateFactory::rate_mutex;

ReactionRateFactory::ReactionRateFactory()
{
    // ArrheniusRate evaluator
    reg("Arrhenius", [](const AnyMap& node, const UnitStack& rate_units) {
        return new ArrheniusRate(node, rate_units);
    });
    addAlias("Arrhenius", "");
    addAlias("Arrhenius", "elementary");
    addAlias("Arrhenius", "three-body");

    // BlowersMaselRate evaluator
    reg("Blowers-Masel", [](const AnyMap& node, const UnitStack& rate_units) {
        return new BlowersMaselRate(node, rate_units);
    });

    // Lindemann falloff evaluator
    reg("Lindemann", [](const AnyMap& node, const UnitStack& rate_units) {
        return new LindemannRate(node, rate_units);
    });
    addAlias("Lindemann", "falloff");

    // Troe falloff evaluator
    reg("Troe", [](const AnyMap& node, const UnitStack& rate_units) {
        return new TroeRate(node, rate_units);
    });

    // SRI falloff evaluator
    reg("SRI", [](const AnyMap& node, const UnitStack& rate_units) {
        return new SriRate(node, rate_units);
    });

    // Tsang falloff evaluator
    reg("Tsang", [](const AnyMap& node, const UnitStack& rate_units) {
        return new TsangRate(node, rate_units);
    });

    // PlogRate evaluator
    reg("pressure-dependent-Arrhenius", [](const AnyMap& node, const UnitStack& rate_units) {
        return new PlogRate(node, rate_units);
    });

    // ChebyshevRate evaluator
    reg("Chebyshev", [](const AnyMap& node, const UnitStack& rate_units) {
        return new ChebyshevRate3(node, rate_units);
    });

    // CustomFunc1Rate evaluator
    reg("custom-rate-function", [](const AnyMap& node, const UnitStack& rate_units) {
        return new CustomFunc1Rate(node, rate_units);
    });
}

shared_ptr<ReactionRate> newReactionRate(const std::string& type)
{
    return shared_ptr<ReactionRate> (
        ReactionRateFactory::factory()->create(type, AnyMap(), UnitStack({})));
}

shared_ptr<ReactionRate> newReactionRate(
    const AnyMap& rate_node, const UnitStack& rate_units)
{
    std::string type = ""; // default is to create Arrhenius from empty
    if (rate_node.hasKey("type")) {
        type = rate_node["type"].asString();
    }

    if (type == "falloff" || type == "chemically-activated") {
        if (rate_node.hasKey("Troe")) {
            type = "Troe";
        } else if (rate_node.hasKey("SRI")) {
            type = "SRI";
        } else if (rate_node.hasKey("Tsang")) {
            type = "Tsang";
        } else {
            type = "Lindemann";
        }
    }

    if (!(ReactionRateFactory::factory()->exists(type))) {
        throw InputFileError("ReactionRateFactory::newReactionRate", rate_node,
            "Unknown reaction rate type '{}'", type);
    }

    return shared_ptr<ReactionRate> (
        ReactionRateFactory::factory()->create(type, rate_node, rate_units));
}

shared_ptr<ReactionRate> newReactionRate(const AnyMap& rate_node)
{
    const UnitSystem& system = rate_node.units();
    if (system.convertTo(1., "m") != 1. || system.convertTo(1., "kmol") != 1.) {
        throw InputFileError("ReactionRateFactory::newReactionRate",
            rate_node.at("__units__"),
            "Alternative units for 'length' or 'quantity` are not supported "
            "when creating\na standalone 'ReactionRate' object.");
    }
    AnyMap node(rate_node);
        return newReactionRate(node, UnitStack({}));
}

}
