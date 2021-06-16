 /**
 *  @file RateFactory.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/RateFactory.h"
#include "cantera/kinetics/MultiRate.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/base/AnyMap.h"

namespace Cantera
{

RateFactory* RateFactory::s_factory = 0;
std::mutex RateFactory::rate_mutex;

RateFactory::RateFactory()
{
    // ArrheniusRate evaluator
    reg("ArrheniusRate", [](const AnyMap& node, const Units& rate_units) {
        return new ArrheniusRate(node, rate_units);
    });
    addAlias("ArrheniusRate", "");
    addAlias("ArrheniusRate", "elementary");
    addAlias("ArrheniusRate", "arrhenius");
    addAlias("ArrheniusRate", "three-body");
    addAlias("ArrheniusRate", "threebody");
    addAlias("ArrheniusRate", "three_body");

    // PlogRate evaluator
    reg("PlogRate", [](const AnyMap& node, const Units& rate_units) {
        return new PlogRate(node, rate_units);
    });
    addAlias("PlogRate", "pressure-dependent-Arrhenius");
    addAlias("PlogRate", "plog");
    addAlias("PlogRate", "pressure-dependent-Arrhenius");

    // ChebyshevRate evaluator
    reg("ChebyshevRate", [](const AnyMap& node, const Units& rate_units) {
        return new ChebyshevRate3(node, rate_units);
    });
    addAlias("ChebyshevRate", "Chebyshev");
    addAlias("ChebyshevRate", "chebyshev");

    // CustomFunc1Rate evaluator
    reg("custom-function", [](const AnyMap& node, const Units& rate_units) {
        return new CustomFunc1Rate(node, rate_units);
    });
    addAlias("custom-function", "custom-rate-function");
}

shared_ptr<ReactionRateBase> newRate(const std::string& type)
{
    if (RateFactory::factory()->exists(type)) {
        return shared_ptr<ReactionRateBase> (
            RateFactory::factory()->create(type, AnyMap(), Units(0.0)));
    }
    return shared_ptr<ReactionRateBase> ();
}

shared_ptr<ReactionRateBase> newRate(const AnyMap& rate_node, const Units& rate_units)
{
    std::string type = "";
    if (rate_node.empty()) {
        throw InputFileError("RateFactory::newRate", rate_node,
            "Received invalid empty node.");
    } else if (rate_node.hasKey("type")) {
        type = rate_node["type"].asString();
    }

    if (!(RateFactory::factory()->exists(type))) {
        throw InputFileError("RateFactory::newRate", rate_node,
            "Unknown reaction rate type '{}'", type);
    }

    return shared_ptr<ReactionRateBase> (
        RateFactory::factory()->create(type, rate_node, rate_units));
}

shared_ptr<ReactionRateBase> newRate(const AnyMap& rate_node, const Kinetics& kin)
{
    if (rate_node.empty()) {
        return newRate(AnyMap(), Units(0.0));
    }
    size_t idx = kin.reactionPhaseIndex();
    Units rate_units = kin.thermo(idx).standardConcentrationUnits();
    return newRate(rate_node, rate_units);
}

std::string canonicalRateName(const std::string& type)
{
    if (RateFactory::factory()->exists(type)) {
        return RateFactory::factory()->canonicalize(type);
    } else {
        return "undefined";
    }
}

}
