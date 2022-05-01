/**
 *  @file Reaction.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/Reaction.h"
#include "ReactionFactory.h"
#include "cantera/kinetics/ReactionRateFactory.h"
#include "cantera/kinetics/FalloffFactory.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/base/ctml.h"
#include "cantera/base/Array.h"
#include "cantera/base/AnyMap.h"
#include "cantera/base/utilities.h"
#include "cantera/base/stringUtils.h"
#include <boost/algorithm/string/predicate.hpp>
#include <sstream>
#include <set>

#include <boost/algorithm/string.hpp>

namespace ba = boost::algorithm;

namespace Cantera
{

Reaction::Reaction()
    : reaction_type(NONE)
    , reversible(true)
    , duplicate(false)
    , allow_nonreactant_orders(false)
    , allow_negative_orders(false)
    , rate_units(0.0)
    , m_valid(true)
{
}

Reaction::Reaction(const Composition& reactants_,
                   const Composition& products_,
                   shared_ptr<ReactionRate> rate_)
    : reaction_type(NONE)
    , reactants(reactants_)
    , products(products_)
    , reversible(true)
    , duplicate(false)
    , allow_nonreactant_orders(false)
    , allow_negative_orders(false)
    , rate_units(0.0)
    , m_valid(true)
    , m_rate(rate_)
{
}

Reaction::Reaction(const AnyMap& node, const Kinetics& kin)
    : Reaction()
{
    setParameters(node, kin);
    if (kin.nPhases()) {
        size_t nDim = kin.thermo(kin.reactionPhaseIndex()).nDim();
        if (nDim == 3) {
            setRate(newReactionRate(node, calculateRateCoeffUnits3(kin)));
        } else {
            AnyMap rateNode = node;
            if (!rateNode.hasKey("type")) {
                // Reaction type is not specified
                rateNode["type"] = "Arrhenius";
            }
            std::string type = rateNode["type"].asString();
            if (rateNode.hasKey("rate-constant")) {
                if (!boost::algorithm::starts_with(type, "interface-")) {
                    rateNode["type"] = "interface-" + type;
                }
            } else if (node.hasKey("sticking-coefficient")) {
                if (!boost::algorithm::starts_with(type, "sticking-")) {
                    rateNode["type"] = "sticking-" + type;
                }
            } else {
                throw InputFileError("Reaction::Reaction", input,
                    "Unable to infer interface reaction type.");
            }
            setRate(newReactionRate(rateNode, calculateRateCoeffUnits3(kin)));
        }
    } else {
        // @deprecated This route is only used for legacy reaction types.
        setRate(newReactionRate(node));
    }
}

Reaction::Reaction(int type)
    : reaction_type(type)
    , reversible(true)
    , duplicate(false)
    , allow_nonreactant_orders(false)
    , allow_negative_orders(false)
    , rate_units(0.0)
    , m_valid(true)
{
    warn_deprecated("Reaction::Reaction()",
        "To be removed after Cantera 2.6. Use constructor without parameter "
        "'type' instead.");
}

Reaction::Reaction(int type, const Composition& reactants_,
                   const Composition& products_)
    : reaction_type(type)
    , reactants(reactants_)
    , products(products_)
    , reversible(true)
    , duplicate(false)
    , allow_nonreactant_orders(false)
    , allow_negative_orders(false)
    , rate_units(0.0)
    , m_valid(true)
{
    warn_deprecated("Reaction::Reaction()",
        "To be removed after Cantera 2.6. Use constructor without parameter "
        "'type' instead.");
}

void Reaction::validate()
{
    if (!allow_nonreactant_orders) {
        for (const auto& order : orders) {
            if (reactants.find(order.first) == reactants.end()) {
                throw InputFileError("Reaction::validate", input,
                    "Reaction order specified for non-reactant species '{}'",
                    order.first);
           }
        }
    }

    if (!allow_negative_orders) {
        for (const auto& order : orders) {
            if (order.second < 0.0) {
                throw InputFileError("Reaction::validate", input,
                    "Negative reaction order specified for species '{}'",
                    order.first);
            }
        }
    }

    // If reaction orders are specified, then this reaction does not follow
    // mass-action kinetics, and is not an elementary reaction. So check that it
    // is not reversible, since computing the reverse rate from thermochemistry
    // only works for elementary reactions.
    if (reversible && !orders.empty()) {
        throw InputFileError("Reaction::validate", input,
            "Reaction orders may only be given for irreversible reactions");
    }

    // Call validation of reaction rate evaluator
    if (!usesLegacy()) {
        m_rate->check(equation(), input);
        m_rate->validate(equation());
    }
}

AnyMap Reaction::parameters(bool withInput) const
{
    AnyMap out;
    getParameters(out);
    if (withInput) {
        out.update(input);
    }

    static bool reg = AnyMap::addOrderingRules("Reaction",
        {{"head", "type"},
         {"head", "equation"},
         {"tail", "duplicate"},
         {"tail", "orders"},
         {"tail", "negative-orders"},
         {"tail", "nonreactant-orders"}
        });
    if (reg) {
        out["__type__"] = "Reaction";
    }
    return out;
}

void Reaction::getParameters(AnyMap& reactionNode) const
{
    reactionNode["equation"] = equation();

    if (duplicate) {
        reactionNode["duplicate"] = true;
    }
    if (orders.size()) {
        reactionNode["orders"] = orders;
    }
    if (allow_negative_orders) {
        reactionNode["negative-orders"] = true;
    }
    if (allow_nonreactant_orders) {
        reactionNode["nonreactant-orders"] = true;
    }

    if (m_rate) {
        reactionNode.update(m_rate->parameters());

        // strip information not needed for reconstruction
        if (reactionNode.hasKey("type")) {
            std::string type = reactionNode["type"].asString();
            if (boost::algorithm::starts_with(type, "Arrhenius")) {
                reactionNode.erase("type");
            } else if (boost::algorithm::starts_with(type, "Blowers-Masel")) {
                reactionNode["type"] = "Blowers-Masel";
            }
        }
    }
}

void Reaction::setParameters(const AnyMap& node, const Kinetics& kin)
{
    if (node.empty()) {
        // empty node: used by newReaction() factory loader
        return;
    }

    input = node;
    input.copyMetadata(node);
    setEquation(node["equation"].asString(), &kin);
    // Non-stoichiometric reaction orders
    if (node.hasKey("orders")) {
        for (const auto& order : node["orders"].asMap<double>()) {
            orders[order.first] = order.second;
            if (kin.kineticsSpeciesIndex(order.first) == npos) {
                setValid(false);
            }
        }
    }

    // Flags
    id = node.getString("id", "");
    duplicate = node.getBool("duplicate", false);
    allow_negative_orders = node.getBool("negative-orders", false);
    allow_nonreactant_orders = node.getBool("nonreactant-orders", false);
}

void Reaction::setRate(shared_ptr<ReactionRate> rate)
{
    if (!rate) {
        // null pointer
        m_rate.reset();
    } else {
        m_rate = rate;
    }

    if (reactants.count("(+M)") && std::dynamic_pointer_cast<ChebyshevRate>(m_rate)) {
        warn_deprecated("Chebyshev reaction equation", input, "Specifying '(+M)' "
            "in the reaction equation for Chebyshev reactions is deprecated.");
        // remove optional third body notation
        reactants.erase("(+M)");
        products.erase("(+M)");
    }

    if (reactants.count("M") && std::dynamic_pointer_cast<PlogRate>(m_rate)) {
        throw InputFileError("Reaction::setRate", input, "Found superfluous 'M' in "
            "pressure-dependent-Arrhenius reaction.");
    }
}

std::string Reaction::reactantString() const
{
    std::ostringstream result;
    for (auto iter = reactants.begin(); iter != reactants.end(); ++iter) {
        if (iter != reactants.begin()) {
            result << " + ";
        }
        if (iter->second != 1.0) {
            result << iter->second << " ";
        }
        result << iter->first;
    }
  return result.str();
}

std::string Reaction::productString() const
{
    std::ostringstream result;
    for (auto iter = products.begin(); iter != products.end(); ++iter) {
        if (iter != products.begin()) {
            result << " + ";
        }
        if (iter->second != 1.0) {
            result << iter->second << " ";
        }
        result << iter->first;
    }
  return result.str();
}

std::string Reaction::equation() const
{
    if (reversible) {
        return reactantString() + " <=> " + productString();
    } else {
        return reactantString() + " => " + productString();
    }
}

void Reaction::setEquation(const std::string& equation, const Kinetics* kin)
{
    parseReactionEquation(*this, equation, input, kin);
}

std::string Reaction::type() const
{
    return "reaction";
}

void Reaction::calculateRateCoeffUnits(const Kinetics& kin)
{
    if (!valid()) {
        // If a reaction is invalid because of missing species in the Kinetics
        // object, determining the units of the rate coefficient is impossible.
        return;
    }

    // Determine the units of the rate coefficient
    Units rxn_phase_units = kin.thermo(kin.reactionPhaseIndex()).standardConcentrationUnits();
    rate_units = rxn_phase_units;
    rate_units *= Units(1.0, 0, 0, -1);
    for (const auto& order : orders) {
        const auto& phase = kin.speciesPhase(order.first);
        rate_units *= phase.standardConcentrationUnits().pow(-order.second);
    }
    for (const auto& stoich : reactants) {
        // Order for each reactant is the reactant stoichiometric coefficient,
        // unless already overridden by user-specified orders
        if (stoich.first == "M" || ba::starts_with(stoich.first, "(+")) {
            // calculateRateCoeffUnits may be called before these pseudo-species
            // have been stripped from the reactants
            continue;
        } else if (orders.find(stoich.first) == orders.end()) {
            const auto& phase = kin.speciesPhase(stoich.first);
            rate_units *= phase.standardConcentrationUnits().pow(-stoich.second);
        }
    }
}

UnitStack Reaction::calculateRateCoeffUnits3(const Kinetics& kin)
{
    if (!valid()) {
        // If a reaction is invalid because of missing species in the Kinetics
        // object, determining the units of the rate coefficient is impossible.
        return UnitStack({});
    }

    // Determine the units of the rate coefficient
    const auto& rxn_phase = kin.thermo(kin.reactionPhaseIndex());
    UnitStack rate_units(rxn_phase.standardConcentrationUnits());

    // Set output units to standardConcentrationUnits per second
    rate_units.join(1.);
    rate_units.update(Units(1.0, 0, 0, -1), 1.);

    for (const auto& order : orders) {
        const auto& phase = kin.speciesPhase(order.first);
        // Account for specified reaction orders
        rate_units.update(phase.standardConcentrationUnits(), -order.second);
    }
    for (const auto& stoich : reactants) {
        // Order for each reactant is the reactant stoichiometric coefficient,
        // unless already overridden by user-specified orders
        if (stoich.first == "M" || ba::starts_with(stoich.first, "(+")) {
            // calculateRateCoeffUnits may be called before these pseudo-species
            // have been stripped from the reactants
            continue;
        } else if (orders.find(stoich.first) == orders.end()) {
            const auto& phase = kin.speciesPhase(stoich.first);
            // Account for each reactant species
            rate_units.update(phase.standardConcentrationUnits(), -stoich.second);
        }
    }

    if (m_third_body) {
        // Account for third-body collision partner as the last entry
        rate_units.join(-1);
    }

    return rate_units;
}

void updateUndeclared(std::vector<std::string>& undeclared,
                      const Composition& comp, const Kinetics& kin)
{
    for (const auto& sp: comp) {
        if (kin.kineticsSpeciesIndex(sp.first) == npos) {
            undeclared.emplace_back(sp.first);
        }
    }
}

std::pair<std::vector<std::string>, bool> Reaction::undeclaredThirdBodies(
        const Kinetics& kin) const
{
    std::vector<std::string> undeclared;
    if (m_third_body) {
        updateUndeclared(undeclared, m_third_body->efficiencies, kin);
        return std::make_pair(undeclared, m_third_body->specified_collision_partner);
    }
    return std::make_pair(undeclared, false);
}

void Reaction::checkBalance(const Kinetics& kin) const
{
    Composition balr, balp;

    // iterate over products and reactants
    for (const auto& sp : products) {
        const ThermoPhase& ph = kin.speciesPhase(sp.first);
        size_t k = ph.speciesIndex(sp.first);
        double stoich = sp.second;
        for (size_t m = 0; m < ph.nElements(); m++) {
            balr[ph.elementName(m)] = 0.0; // so that balr contains all species
            balp[ph.elementName(m)] += stoich * ph.nAtoms(k, m);
        }
    }
    for (const auto& sp : reactants) {
        const ThermoPhase& ph = kin.speciesPhase(sp.first);
        size_t k = ph.speciesIndex(sp.first);
        double stoich = sp.second;
        for (size_t m = 0; m < ph.nElements(); m++) {
            balr[ph.elementName(m)] += stoich * ph.nAtoms(k, m);
        }
    }

    std::string msg;
    bool ok = true;
    for (const auto& el : balr) {
        const std::string& elem = el.first;
        double elemsum = balr[elem] + balp[elem];
        double elemdiff = fabs(balp[elem] - balr[elem]);
        if (elemsum > 0.0 && elemdiff / elemsum > 1e-4) {
            ok = false;
            msg += fmt::format("  {}           {}           {}\n",
                               elem, balr[elem], balp[elem]);
        }
    }
    if (!ok) {
        throw InputFileError("Reaction::checkBalance", input,
            "The following reaction is unbalanced: {}\n"
            "  Element    Reactants    Products\n{}",
            equation(), msg);
    }

    if (kin.thermo(kin.reactionPhaseIndex()).nDim() == 3) {
        return;
    }

    // Check that the number of surface sites is balanced
    double reac_sites = 0.0;
    double prod_sites = 0.0;
    auto& surf = dynamic_cast<const SurfPhase&>(kin.thermo(kin.surfacePhaseIndex()));
    for (const auto& reactant : reactants) {
        size_t k = surf.speciesIndex(reactant.first);
        if (k != npos) {
            reac_sites += reactant.second * surf.size(k);
        }
    }
    for (const auto& product : products) {
        size_t k = surf.speciesIndex(product.first);
        if (k != npos) {
            prod_sites += product.second * surf.size(k);
        }
    }
    if (fabs(reac_sites - prod_sites) > 1e-5 * (reac_sites + prod_sites)) {
        throw InputFileError("Reaction::checkBalance", input,
            "Number of surface sites not balanced in reaction {}.\n"
            "Reactant sites: {}\nProduct sites: {}",
            equation(), reac_sites, prod_sites);
    }
}

bool Reaction::checkSpecies(const Kinetics& kin) const
{
    // Check for undeclared species
    std::vector<std::string> undeclared;
    updateUndeclared(undeclared, reactants, kin);
    updateUndeclared(undeclared, products, kin);
    if (!undeclared.empty()) {
        if (kin.skipUndeclaredSpecies()) {
            return false;
        } else {
            throw InputFileError("Reaction::checkSpecies", input, "Reaction '{}'\n"
                "contains undeclared species: '{}'",
                equation(), boost::algorithm::join(undeclared, "', '"));
        }
    }

    undeclared.clear();
    updateUndeclared(undeclared, orders, kin);
    if (!undeclared.empty()) {
        if (kin.skipUndeclaredSpecies()) {
            return false;
        } else {
            if (input.hasKey("orders")) {
                throw InputFileError("Reaction::checkSpecies", input["orders"],
                    "Reaction '{}'\n"
                    "defines reaction orders for undeclared species: '{}'",
                    equation(), boost::algorithm::join(undeclared, "', '"));
            }
            // Error for empty input AnyMap (that is, XML)
            throw InputFileError("Reaction::checkSpecies", input, "Reaction '{}'\n"
                "defines reaction orders for undeclared species: '{}'",
                equation(), boost::algorithm::join(undeclared, "', '"));
        }
    }

    // Use helper function while there is no uniform handling of third bodies
    auto third = undeclaredThirdBodies(kin);
    undeclared = third.first;
    bool specified_collision_partner_ = third.second;
    if (!undeclared.empty()) {
        if (!kin.skipUndeclaredThirdBodies()) {
            if (input.hasKey("efficiencies")) {
                throw InputFileError("Reaction::checkSpecies", input["efficiencies"],
                    "Reaction '{}'\n"
                    "defines third-body efficiencies for undeclared species: '{}'",
                    equation(), boost::algorithm::join(undeclared, "', '"));
            }
            // Error for specified ThirdBody or empty input AnyMap
            throw InputFileError("Reaction::checkSpecies", input, "Reaction '{}'\n"
                "is a three-body reaction with undeclared species: '{}'",
                equation(), boost::algorithm::join(undeclared, "', '"));
        } else if (kin.skipUndeclaredSpecies() && specified_collision_partner_) {
            return false;
        }
    }

    checkBalance(kin);

    return true;
}

bool Reaction::usesElectrochemistry(const Kinetics& kin) const
{
    // Check electrochemistry
    vector_fp e_counter(kin.nPhases(), 0.0);

    // Find the number of electrons in the products for each phase
    for (const auto& sp : products) {
        size_t kkin = kin.kineticsSpeciesIndex(sp.first);
        size_t i = kin.speciesPhaseIndex(kkin);
        size_t kphase = kin.thermo(i).speciesIndex(sp.first);
        e_counter[i] += sp.second * kin.thermo(i).charge(kphase);
    }

    // Subtract the number of electrons in the reactants for each phase
    for (const auto& sp : reactants) {
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

ElementaryReaction2::ElementaryReaction2(const Composition& reactants_,
                                         const Composition products_,
                                         const Arrhenius2& rate_)
    : Reaction(reactants_, products_)
    , rate(rate_)
    , allow_negative_pre_exponential_factor(false)
{
    reaction_type = ELEMENTARY_RXN;
}

ElementaryReaction2::ElementaryReaction2()
    : Reaction()
    , allow_negative_pre_exponential_factor(false)
{
    reaction_type = ELEMENTARY_RXN;
}

void ElementaryReaction2::validate()
{
    Reaction::validate();
    if (!allow_negative_pre_exponential_factor &&
        rate.preExponentialFactor() < 0) {
        throw InputFileError("ElementaryReaction2::validate", input,
            "Undeclared negative pre-exponential factor found in reaction '"
            + equation() + "'");
    }
}

void ElementaryReaction2::getParameters(AnyMap& reactionNode) const
{
    Reaction::getParameters(reactionNode);
    if (allow_negative_pre_exponential_factor) {
        reactionNode["negative-A"] = true;
    }
    AnyMap rateNode;
    rate.getParameters(rateNode, rate_units);
    reactionNode["rate-constant"] = std::move(rateNode);
}

ThirdBody::ThirdBody(double default_eff)
    : default_efficiency(default_eff)
    , specified_collision_partner(false)
    , mass_action(true)
{
}

ThirdBody::ThirdBody(const AnyMap& node)
    : specified_collision_partner(false)
    , mass_action(true)
{
    setEfficiencies(node);
}

void ThirdBody::setEfficiencies(const AnyMap& node)
{
    default_efficiency = node.getDouble("default-efficiency", 1.0);
    if (node.hasKey("efficiencies")) {
        efficiencies = node["efficiencies"].asMap<double>();
    }
}

double ThirdBody::efficiency(const std::string& k) const
{
    return getValue(efficiencies, k, default_efficiency);
}

ThreeBodyReaction2::ThreeBodyReaction2()
{
    reaction_type = THREE_BODY_RXN;
}

ThreeBodyReaction2::ThreeBodyReaction2(const Composition& reactants_,
                                       const Composition& products_,
                                       const Arrhenius2& rate_,
                                       const ThirdBody& tbody)
    : ElementaryReaction2(reactants_, products_, rate_)
    , third_body(tbody)
{
    reaction_type = THREE_BODY_RXN;
}

std::string ThreeBodyReaction2::reactantString() const
{
    if (third_body.specified_collision_partner) {
        return ElementaryReaction2::reactantString() + " + "
            + third_body.efficiencies.begin()->first;
    } else {
        return ElementaryReaction2::reactantString() + " + M";
    }
}

std::string ThreeBodyReaction2::productString() const
{
    if (third_body.specified_collision_partner) {
        return ElementaryReaction2::productString() + " + "
            + third_body.efficiencies.begin()->first;
    } else {
        return ElementaryReaction2::productString() + " + M";
    }
}

void ThreeBodyReaction2::calculateRateCoeffUnits(const Kinetics& kin)
{
    ElementaryReaction2::calculateRateCoeffUnits(kin);
    bool specified_collision_partner_ = false;
    for (const auto& reac : reactants) {
        // While this reaction was already identified as a three-body reaction in a
        // pre-processing step, this method is often called before a three-body
        // reaction is fully instantiated. For the determination of the correct units,
        // it is necessary to check whether the reaction uses a generic 'M' or an
        // explicitly specified collision partner that may not have been deleted yet.
        if (reac.first != "M" && products.count(reac.first)) {
            // detected specified third-body collision partner
            specified_collision_partner_ = true;
        }
    }
    if (!specified_collision_partner_) {
        const ThermoPhase& rxn_phase = kin.thermo(kin.reactionPhaseIndex());
        rate_units *= rxn_phase.standardConcentrationUnits().pow(-1);
    }
}

void ThreeBodyReaction2::getParameters(AnyMap& reactionNode) const
{
    ElementaryReaction2::getParameters(reactionNode);
    if (!third_body.specified_collision_partner) {
        reactionNode["type"] = "three-body";
        reactionNode["efficiencies"] = third_body.efficiencies;
        reactionNode["efficiencies"].setFlowStyle();
        if (third_body.default_efficiency != 1.0) {
            reactionNode["default-efficiency"] = third_body.default_efficiency;
        }
    }
}

std::pair<std::vector<std::string>, bool> ThreeBodyReaction2::undeclaredThirdBodies(
        const Kinetics& kin) const
{
    std::vector<std::string> undeclared;
    updateUndeclared(undeclared, third_body.efficiencies, kin);
    return std::make_pair(undeclared, third_body.specified_collision_partner);
}

FalloffReaction2::FalloffReaction2()
    : Reaction()
    , falloff(new Lindemann())
    , allow_negative_pre_exponential_factor(false)
    , low_rate_units(0.0)
{
    reaction_type = FALLOFF_RXN;
}

FalloffReaction2::FalloffReaction2(
        const Composition& reactants_, const Composition& products_,
        const Arrhenius2& low_rate_, const Arrhenius2& high_rate_,
        const ThirdBody& tbody)
    : Reaction(reactants_, products_)
    , low_rate(low_rate_)
    , high_rate(high_rate_)
    , third_body(tbody)
    , falloff(new Lindemann())
    , allow_negative_pre_exponential_factor(false)
    , low_rate_units(0.0)
{
    reaction_type = FALLOFF_RXN;
}

std::string FalloffReaction2::reactantString() const
{
    if (third_body.default_efficiency == 0 &&
        third_body.efficiencies.size() == 1) {
        return Reaction::reactantString() + " (+" +
            third_body.efficiencies.begin()->first + ")";
    } else {
        return Reaction::reactantString() + " (+M)";
    }
}

std::string FalloffReaction2::productString() const
{
    if (third_body.default_efficiency == 0 &&
        third_body.efficiencies.size() == 1) {
        return Reaction::productString() + " (+" +
            third_body.efficiencies.begin()->first + ")";
    } else {
        return Reaction::productString() + " (+M)";
    }
}

void FalloffReaction2::validate()
{
    Reaction::validate();
    if (!allow_negative_pre_exponential_factor &&
        (low_rate.preExponentialFactor() < 0 ||
         high_rate.preExponentialFactor() < 0)) {
        throw InputFileError("FalloffReaction2::validate", input, "Negative "
            "pre-exponential factor found for reaction '" + equation() + "'");
    }
    if (low_rate.preExponentialFactor() * high_rate.preExponentialFactor() < 0) {
        throw InputFileError("FalloffReaction2::validate", input, "High and "
            "low rate pre-exponential factors must have the same sign."
            "Reaction: '{}'", equation());
    }
}

void FalloffReaction2::calculateRateCoeffUnits(const Kinetics& kin)
{
    Reaction::calculateRateCoeffUnits(kin);
    const ThermoPhase& rxn_phase = kin.thermo(kin.reactionPhaseIndex());
    low_rate_units = rate_units;
    low_rate_units *= rxn_phase.standardConcentrationUnits().pow(-1);
}

void FalloffReaction2::getParameters(AnyMap& reactionNode) const
{
    Reaction::getParameters(reactionNode);
    reactionNode["type"] = "falloff-legacy";
    AnyMap lowRateNode;
    low_rate.getParameters(lowRateNode, low_rate_units);
    reactionNode["low-P-rate-constant"] = std::move(lowRateNode);
    AnyMap highRateNode;
    high_rate.getParameters(highRateNode, rate_units);
    reactionNode["high-P-rate-constant"] = std::move(highRateNode);
    falloff->getParameters(reactionNode);

    reactionNode["efficiencies"] = third_body.efficiencies;
    reactionNode["efficiencies"].setFlowStyle();
    if (third_body.default_efficiency != 1.0) {
        reactionNode["default-efficiency"] = third_body.default_efficiency;
    }
}

std::pair<std::vector<std::string>, bool> FalloffReaction2::undeclaredThirdBodies(
        const Kinetics& kin) const
{
    std::vector<std::string> undeclared;
    updateUndeclared(undeclared, third_body.efficiencies, kin);
    return std::make_pair(undeclared, false);
}

ChemicallyActivatedReaction2::ChemicallyActivatedReaction2()
{
    reaction_type = CHEMACT_RXN;
}

ChemicallyActivatedReaction2::ChemicallyActivatedReaction2(
        const Composition& reactants_, const Composition& products_,
        const Arrhenius2& low_rate_, const Arrhenius2& high_rate_,
        const ThirdBody& tbody)
    : FalloffReaction2(reactants_, products_, low_rate_, high_rate_, tbody)
{
    reaction_type = CHEMACT_RXN;
}

void ChemicallyActivatedReaction2::calculateRateCoeffUnits(const Kinetics& kin)
{
    Reaction::calculateRateCoeffUnits(kin); // Skip FalloffReaction2
    const ThermoPhase& rxn_phase = kin.thermo(kin.reactionPhaseIndex());
    low_rate_units = rate_units;
    rate_units *= rxn_phase.standardConcentrationUnits();
}

void ChemicallyActivatedReaction2::getParameters(AnyMap& reactionNode) const
{
    FalloffReaction2::getParameters(reactionNode);
    reactionNode["type"] = "chemically-activated";
}

PlogReaction2::PlogReaction2()
    : Reaction()
{
    reaction_type = PLOG_RXN;
}

PlogReaction2::PlogReaction2(const Composition& reactants_,
                             const Composition& products_, const Plog& rate_)
    : Reaction(reactants_, products_)
    , rate(rate_)
{
    reaction_type = PLOG_RXN;
}

void PlogReaction2::getParameters(AnyMap& reactionNode) const
{
    Reaction::getParameters(reactionNode);
    reactionNode["type"] = "pressure-dependent-Arrhenius";
    rate.getParameters(reactionNode, rate_units);
}

ChebyshevReaction2::ChebyshevReaction2()
    : Reaction()
{
    reaction_type = CHEBYSHEV_RXN;
}

ChebyshevReaction2::ChebyshevReaction2(const Composition& reactants_,
                                       const Composition& products_,
                                       const ChebyshevRate& rate_)
    : Reaction(reactants_, products_)
    , rate(rate_)
{
    reaction_type = CHEBYSHEV_RXN;
}

void ChebyshevReaction2::getParameters(AnyMap& reactionNode) const
{
    Reaction::getParameters(reactionNode);
    reactionNode["type"] = "Chebyshev";
    rate.getParameters(reactionNode, rate_units);
}

InterfaceReaction2::InterfaceReaction2()
    : is_sticking_coefficient(false)
    , use_motz_wise_correction(false)
{
    reaction_type = INTERFACE_RXN;
}

InterfaceReaction2::InterfaceReaction2(const Composition& reactants_,
                                       const Composition& products_,
                                       const Arrhenius2& rate_,
                                       bool isStick)
    : ElementaryReaction2(reactants_, products_, rate_)
    , is_sticking_coefficient(isStick)
    , use_motz_wise_correction(false)
{
    reaction_type = INTERFACE_RXN;
}

void InterfaceReaction2::calculateRateCoeffUnits(const Kinetics& kin)
{
    ElementaryReaction2::calculateRateCoeffUnits(kin);
    if (is_sticking_coefficient || input.hasKey("sticking-coefficient")) {
        rate_units = Units(1.0); // sticking coefficients are dimensionless
    }
}

void InterfaceReaction2::getParameters(AnyMap& reactionNode) const
{
    ElementaryReaction2::getParameters(reactionNode);
    if (is_sticking_coefficient) {
        reactionNode["sticking-coefficient"] = std::move(reactionNode["rate-constant"]);
        reactionNode.erase("rate-constant");
    }
    if (use_motz_wise_correction) {
        reactionNode["Motz-Wise"] = true;
    }
    if (!sticking_species.empty()) {
        reactionNode["sticking-species"] = sticking_species;
    }
    if (!coverage_deps.empty()) {
        AnyMap deps;
        for (const auto& d : coverage_deps) {
            AnyMap dep;
            dep["a"] = d.second.a;
            dep["m"] = d.second.m;
            dep["E"].setQuantity(d.second.E, "K", true);
            deps[d.first] = std::move(dep);
        }
        reactionNode["coverage-dependencies"] = std::move(deps);
    }
}

void InterfaceReaction2::validate(Kinetics& kin)
{
    if (is_sticking_coefficient) {
        fmt::memory_buffer err_reactions;
        double T[] = {200.0, 500.0, 1000.0, 2000.0, 5000.0, 10000.0};
        for (size_t i=0; i < 6; i++) {
            double k = rate.updateRC(log(T[i]), 1/T[i]);
            if (k > 1) {
                fmt_append(err_reactions,
                    "\n Sticking coefficient is greater than 1 for reaction '{}'\n"
                    " at T = {:.1f}\n",
                    equation(), T[i]);
            }
        }
        if (err_reactions.size()) {
            warn_user("InterfaceReaction2::validate", to_string(err_reactions));
        }
    }
}

ElectrochemicalReaction2::ElectrochemicalReaction2()
    : beta(0.5)
    , exchange_current_density_formulation(false)
{
}

ElectrochemicalReaction2::ElectrochemicalReaction2(const Composition& reactants_,
                                                   const Composition& products_,
                                                   const Arrhenius2& rate_)
    : InterfaceReaction2(reactants_, products_, rate_)
    , beta(0.5)
    , exchange_current_density_formulation(false)
{
}

void ElectrochemicalReaction2::getParameters(AnyMap& reactionNode) const
{
    InterfaceReaction2::getParameters(reactionNode);
    if (beta != 0.5) {
        reactionNode["beta"] = beta;
    }
    if (exchange_current_density_formulation) {
        reactionNode["exchange-current-density-formulation"] = true;
    }
}

ThreeBodyReaction3::ThreeBodyReaction3()
{
    m_third_body.reset(new ThirdBody);
    setRate(newReactionRate(type()));
}

ThreeBodyReaction3::ThreeBodyReaction3(const Composition& reactants,
                                       const Composition& products,
                                       const ArrheniusRate& rate,
                                       const ThirdBody& tbody)
    : Reaction(reactants, products, make_shared<ArrheniusRate>(rate))
{
    m_third_body = std::make_shared<ThirdBody>(tbody);
}

ThreeBodyReaction3::ThreeBodyReaction3(const AnyMap& node, const Kinetics& kin)
{
    m_third_body.reset(new ThirdBody);
    if (!node.empty()) {
        setParameters(node, kin);
        setRate(newReactionRate(node, calculateRateCoeffUnits3(kin)));
    } else {
        setRate(newReactionRate(type()));
    }
}

bool ThreeBodyReaction3::detectEfficiencies()
{
    for (const auto& reac : reactants) {
        // detect explicitly specified collision partner
        if (products.count(reac.first)) {
            m_third_body->efficiencies[reac.first] = 1.;
        }
    }

    if (m_third_body->efficiencies.size() == 0) {
        return false;
    } else if (m_third_body->efficiencies.size() > 1) {
        throw CanteraError("ThreeBodyReaction3::detectEfficiencies",
            "Found more than one explicitly specified collision partner\n"
            "in reaction '{}'.", equation());
    }

    m_third_body->default_efficiency = 0.;
    m_third_body->specified_collision_partner = true;
    auto sp = m_third_body->efficiencies.begin();

    // adjust reactant coefficients
    auto reac = reactants.find(sp->first);
    if (trunc(reac->second) != 1) {
        reac->second -= 1.;
    } else {
        reactants.erase(reac);
    }

    // adjust product coefficients
    auto prod = products.find(sp->first);
    if (trunc(prod->second) != 1) {
        prod->second -= 1.;
    } else {
        products.erase(prod);
    }

    return true;
}

void ThreeBodyReaction3::setEquation(const std::string& equation, const Kinetics* kin)
{
    Reaction::setEquation(equation, kin);
    if (reactants.count("M") != 1 || products.count("M") != 1) {
        if (!detectEfficiencies()) {
            throw InputFileError("ThreeBodyReaction3::setParameters", input,
                "Reaction equation '{}' does not contain third body 'M'",
                equation);
        }
        return;
    }

    reactants.erase("M");
    products.erase("M");
}

void ThreeBodyReaction3::setParameters(const AnyMap& node, const Kinetics& kin)
{
    if (node.empty()) {
        // empty node: used by newReaction() factory loader
        return;
    }
    Reaction::setParameters(node, kin);
    if (!m_third_body->specified_collision_partner) {
        m_third_body->setEfficiencies(node);
    }
}

void ThreeBodyReaction3::getParameters(AnyMap& reactionNode) const
{
    Reaction::getParameters(reactionNode);
    if (!m_third_body->specified_collision_partner) {
        reactionNode["type"] = "three-body";
        reactionNode["efficiencies"] = m_third_body->efficiencies;
        reactionNode["efficiencies"].setFlowStyle();
        if (m_third_body->default_efficiency != 1.0) {
            reactionNode["default-efficiency"] = m_third_body->default_efficiency;
        }
    }
}

std::string ThreeBodyReaction3::reactantString() const
{
    if (m_third_body->specified_collision_partner) {
        return Reaction::reactantString() + " + "
            + m_third_body->efficiencies.begin()->first;
    } else {
        return Reaction::reactantString() + " + M";
    }
}

std::string ThreeBodyReaction3::productString() const
{
    if (m_third_body->specified_collision_partner) {
        return Reaction::productString() + " + "
            + m_third_body->efficiencies.begin()->first;
    } else {
        return Reaction::productString() + " + M";
    }
}

FalloffReaction3::FalloffReaction3()
    : Reaction()
{
    m_third_body.reset(new ThirdBody);
    m_third_body->mass_action = false;
    setRate(newReactionRate(type()));
}

FalloffReaction3::FalloffReaction3(const Composition& reactants,
                                   const Composition& products,
                                   const ReactionRate& rate,
                                   const ThirdBody& tbody)
    : Reaction(reactants, products)
{
    m_third_body = std::make_shared<ThirdBody>(tbody);
    m_third_body->mass_action = false;
    AnyMap node = rate.parameters();
    node.applyUnits();
    std::string rate_type = node["type"].asString();
    if (rate_type != "falloff" && rate_type != "chemically-activated") {
        // use node information to determine whether rate is a falloff rate
        throw CanteraError("FalloffReaction3::FalloffReaction3",
            "Incompatible types: '{}' is not a falloff rate object.", rate.type());
    }
    setRate(newReactionRate(node));
}

FalloffReaction3::FalloffReaction3(const AnyMap& node, const Kinetics& kin)
{
    m_third_body.reset(new ThirdBody);
    m_third_body->mass_action = false;
    if (!node.empty()) {
        setParameters(node, kin);
        setRate(newReactionRate(node, calculateRateCoeffUnits3(kin)));
    } else {
        setRate(newReactionRate(type()));
    }
}

std::string FalloffReaction3::type() const
{
    if (m_rate &&
        std::dynamic_pointer_cast<FalloffRate>(m_rate)->chemicallyActivated())
    {
        return "chemically-activated";
    }
    return "falloff";
}

std::string FalloffReaction3::reactantString() const
{
    if (m_third_body->specified_collision_partner) {
        return Reaction::reactantString() + " (+" +
            m_third_body->efficiencies.begin()->first + ")";
    } else {
        return Reaction::reactantString() + " (+M)";
    }
}

std::string FalloffReaction3::productString() const
{
    if (m_third_body->specified_collision_partner) {
        return Reaction::productString() + " (+" +
            m_third_body->efficiencies.begin()->first + ")";
    } else {
        return Reaction::productString() + " (+M)";
    }
}

void FalloffReaction3::setParameters(const AnyMap& node, const Kinetics& kin)
{
    if (node.empty()) {
        // empty node: used by newReaction() factory loader
        return;
    }
    Reaction::setParameters(node, kin);

    if (!m_third_body->specified_collision_partner) {
        m_third_body->setEfficiencies(node);
    }
}

void FalloffReaction3::setEquation(const std::string& equation, const Kinetics* kin)
{
    Reaction::setEquation(equation, kin);

    // Detect falloff third body based on partial setup;
    // parseReactionEquation (called via Reaction::setEquation) sets the
    // stoichiometric coefficient of the falloff species to -1.
    std::string third_body_str;
    std::string third_body;
    for (auto& reactant : reactants) {
        if (reactant.second == -1 && ba::starts_with(reactant.first, "(+")) {
            third_body_str = reactant.first;
            third_body = third_body_str.substr(2, third_body_str.size() - 3);
            break;
        }
    }

    // Equation must contain a third body, and it must appear on both sides
    if (third_body_str == "") {
        throw InputFileError("FalloffReaction3::setParameters", input,
            "Reactants for reaction '{}' do not contain a pressure-dependent "
            "third body", equation);
    }
    if (products.count(third_body_str) == 0) {
        throw InputFileError("FalloffReaction3::setParameters", input,
            "Unable to match third body '{}' in reactants and products of "
            "reaction '{}'", third_body, equation);
    }

    // Remove the dummy species
    reactants.erase(third_body_str);
    products.erase(third_body_str);

    if (third_body == "M") {
        m_third_body->specified_collision_partner = false;
    } else {
        // Specific species is listed as the third body
        m_third_body->default_efficiency = 0;
        m_third_body->efficiencies.emplace(third_body, 1.0);
        m_third_body->specified_collision_partner = true;
    }
}

void FalloffReaction3::getParameters(AnyMap& reactionNode) const
{
    Reaction::getParameters(reactionNode);
    if (m_third_body->specified_collision_partner) {
        // pass
    } else if (m_third_body->efficiencies.size()) {
        reactionNode["efficiencies"] = m_third_body->efficiencies;
        reactionNode["efficiencies"].setFlowStyle();
        if (m_third_body->default_efficiency != 1.0) {
            reactionNode["default-efficiency"] = m_third_body->default_efficiency;
        }
    }
}

CustomFunc1Reaction::CustomFunc1Reaction()
{
    setRate(newReactionRate(type()));
}

CustomFunc1Reaction::CustomFunc1Reaction(const Composition& reactants,
                                         const Composition& products,
                                         const CustomFunc1Rate& rate)
    : Reaction(reactants, products)
{
    m_rate.reset(new CustomFunc1Rate(rate));
}

CustomFunc1Reaction::CustomFunc1Reaction(const AnyMap& node, const Kinetics& kin)
{
    if (!node.empty()) {
        setParameters(node, kin);
        setRate(newReactionRate(node, calculateRateCoeffUnits3(kin)));
    } else {
        setRate(newReactionRate(type()));
    }
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
       nreac += static_cast<size_t>(reac.second);
    }

    // ensure that all products have integer stoichiometric coefficients
    size_t nprod = 0;
    for (const auto& prod : R.products) {
       if (trunc(prod.second) != prod.second) {
           return false;
       }
       nprod += static_cast<size_t>(prod.second);
    }

    // either reactant or product side involves exactly three species
    return (nreac == 3) || (nprod == 3);
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
        parseReactionEquation(testReaction, rxn_node["equation"].asString(),
                              rxn_node, &kin);
        if (isThreeBody(testReaction)) {
            type = "three-body";
        }
    }

    if (!(ReactionFactory::factory()->exists(type))) {
        throw InputFileError("ReactionFactory::newReaction", rxn_node["type"],
            "Unknown reaction type '{}'", type);
    }
    Reaction* R = ReactionFactory::factory()->create(type, rxn_node, kin);
    return unique_ptr<Reaction>(R);
}

Arrhenius2 readArrhenius(const XML_Node& arrhenius_node)
{
    return Arrhenius2(getFloat(arrhenius_node, "A", "toSI"),
                      getFloat(arrhenius_node, "b"),
                      getFloat(arrhenius_node, "E", "actEnergy") / GasConstant);
}

Arrhenius2 readArrhenius(const Reaction& R, const AnyValue& rate,
                         const Kinetics& kin, const UnitSystem& units,
                         int pressure_dependence=0)
{
    double A, b, Ta;
    Units rc_units = R.rate_units;
    if (pressure_dependence) {
        Units rxn_phase_units = kin.thermo(kin.reactionPhaseIndex()).standardConcentrationUnits();
        rc_units *= rxn_phase_units.pow(-pressure_dependence);
    }
    if (rate.is<AnyMap>()) {
        auto& rate_map = rate.as<AnyMap>();
        A = units.convert(rate_map["A"], rc_units);
        b = rate_map["b"].asDouble();
        Ta = units.convertActivationEnergy(rate_map["Ea"], "K");
    } else {
        auto& rate_vec = rate.asVector<AnyValue>(3);
        A = units.convert(rate_vec[0], rc_units);
        b = rate_vec[1].asDouble();
        Ta = units.convertActivationEnergy(rate_vec[2], "K");
    }
    return Arrhenius2(A, b, Ta);
}

//! Parse falloff parameters, given a rateCoeff node
/*!
 * @verbatim
 <falloff type="Troe"> 0.5 73.2 5000. 9999. </falloff>
 @endverbatim
*/
void readFalloff(FalloffReaction2& R, const XML_Node& rc_node)
{
    XML_Node& falloff = rc_node.child("falloff");
    std::vector<std::string> p;
    vector_fp falloff_parameters;
    getStringArray(falloff, p);
    size_t np = p.size();
    for (size_t n = 0; n < np; n++) {
        falloff_parameters.push_back(fpValueCheck(p[n]));
    }

    if (caseInsensitiveEquals(falloff["type"], "lindemann")) {
        if (np != 0) {
            throw CanteraError("readFalloff", "Lindemann parameterization "
                "takes no parameters, but {} were given", np);
        }
        R.falloff = newFalloff("Lindemann", falloff_parameters);
    } else if (caseInsensitiveEquals(falloff["type"], "troe")) {
        if (np != 3 && np != 4) {
            throw CanteraError("readFalloff", "Troe parameterization takes "
                "3 or 4 parameters, but {} were given", np);
        }
        R.falloff = newFalloff("Troe", falloff_parameters);
    } else if (caseInsensitiveEquals(falloff["type"], "sri")) {
        if (np != 3 && np != 5) {
            throw CanteraError("readFalloff", "SRI parameterization takes "
                "3 or 5 parameters, but {} were given", np);
        }
        R.falloff = newFalloff("SRI", falloff_parameters);
    } else if (caseInsensitiveEquals(falloff["type"], "tsang")) {
        if (np != 2) {
            throw CanteraError("readFalloff", "Tsang parameterization takes "
                "2 parameters, but {} were given", np);
        }
        R.falloff = newFalloff("Tsang", falloff_parameters);
    } else {
        throw CanteraError("readFalloff", "Unrecognized falloff type: '{}'",
                           falloff["type"]);
    }
}

void readFalloff(FalloffReaction2& R, const AnyMap& node)
{
    if (node.hasKey("Troe")) {
        auto& f = node["Troe"].as<AnyMap>();
        vector_fp params{
            f["A"].asDouble(),
            f["T3"].asDouble(),
            f["T1"].asDouble()
        };
        if (f.hasKey("T2")) {
            params.push_back(f["T2"].asDouble());
        }
        R.falloff = newFalloff("Troe", params);
    } else if (node.hasKey("SRI")) {
        auto& f = node["SRI"].as<AnyMap>();
        vector_fp params{
            f["A"].asDouble(),
            f["B"].asDouble(),
            f["C"].asDouble()
        };
        if (f.hasKey("D")) {
            params.push_back(f["D"].asDouble());
        }
        if (f.hasKey("E")) {
            params.push_back(f["E"].asDouble());
        }
        R.falloff = newFalloff("SRI", params);
    } else if (node.hasKey("Tsang")) {
        auto& f = node["Tsang"].as<AnyMap>();
        vector_fp params{
            f["A"].asDouble(),
            f["B"].asDouble()
        };
        R.falloff = newFalloff("Tsang", params);
    } else {
        R.falloff = newFalloff("Lindemann", {});
    }
}

void readEfficiencies(ThirdBody& tbody, const XML_Node& rc_node)
{
    if (!rc_node.hasChild("efficiencies")) {
        tbody.default_efficiency = 1.0;
        return;
    }
    const XML_Node& eff_node = rc_node.child("efficiencies");
    tbody.default_efficiency = fpValue(eff_node["default"]);
    tbody.efficiencies = parseCompString(eff_node.value());
}

void readEfficiencies(ThirdBody& tbody, const AnyMap& node)
{
    tbody.setEfficiencies(node);
}

bool detectEfficiencies(ThreeBodyReaction2& R)
{
    for (const auto& reac : R.reactants) {
        // detect explicitly specified collision partner
        if (R.products.count(reac.first)) {
            R.third_body.efficiencies[reac.first] = 1.;
        }
    }

    if (R.third_body.efficiencies.size() == 0) {
        return false;
    } else if (R.third_body.efficiencies.size() > 1) {
        throw CanteraError("detectEfficiencies",
            "Found more than one explicitly specified collision partner\n"
            "in reaction '{}'.", R.equation());
    }

    R.third_body.default_efficiency = 0.;
    R.third_body.specified_collision_partner = true;
    auto sp = R.third_body.efficiencies.begin();

    // adjust reactant coefficients
    auto reac = R.reactants.find(sp->first);
    if (trunc(reac->second) != 1) {
        reac->second -= 1.;
    } else {
        R.reactants.erase(reac);
    }

    // adjust product coefficients
    auto prod = R.products.find(sp->first);
    if (trunc(prod->second) != 1) {
        prod->second -= 1.;
    } else {
        R.products.erase(prod);
    }

    return true;
}

void setupReaction(Reaction& R, const XML_Node& rxn_node)
{
    // Reactant and product stoichiometries
    R.reactants = parseCompString(rxn_node.child("reactants").value());
    R.products = parseCompString(rxn_node.child("products").value());

    // Non-stoichiometric reaction orders
    std::vector<XML_Node*> orders = rxn_node.getChildren("order");
    for (size_t i = 0; i < orders.size(); i++) {
        R.orders[orders[i]->attrib("species")] = orders[i]->fp_value();
    }

    // Flags
    R.id = rxn_node.attrib("id");
    R.duplicate = rxn_node.hasAttrib("duplicate");
    const std::string& rev = rxn_node["reversible"];
    R.reversible = (rev == "true" || rev == "yes");
}

void parseReactionEquation(Reaction& R, const std::string& equation,
                           const AnyBase& reactionNode, const Kinetics* kin)
{
    // Parse the reaction equation to determine participating species and
    // stoichiometric coefficients
    std::vector<std::string> tokens;
    tokenizeString(equation, tokens);
    tokens.push_back("+"); // makes parsing last species not a special case

    size_t last_used = npos; // index of last-used token
    bool reactants = true;
    for (size_t i = 1; i < tokens.size(); i++) {
        if (tokens[i] == "+" || ba::starts_with(tokens[i], "(+") ||
            tokens[i] == "<=>" || tokens[i] == "=" || tokens[i] == "=>") {
            std::string species = tokens[i-1];

            double stoich;
            if (last_used != npos && tokens[last_used] == "(+") {
                // Falloff third body with space, such as "(+ M)"
                species = "(+" + species;
                stoich = -1;
            } else if (last_used == i-1 && ba::starts_with(species, "(+")
                       && ba::ends_with(species, ")")) {
                // Falloff 3rd body written without space, such as "(+M)"
                stoich = -1;
            } else if (last_used == i-2) { // Species with no stoich. coefficient
                stoich = 1.0;
            } else if (last_used == i-3) { // Stoich. coefficient and species
                try {
                    stoich = fpValueCheck(tokens[i-2]);
                } catch (CanteraError& err) {
                    throw InputFileError("parseReactionEquation", reactionNode,
                        err.getMessage());
                }
            } else {
                throw InputFileError("parseReactionEquation", reactionNode,
                    "Error parsing reaction string '{}'.\n"
                    "Current token: '{}'\nlast_used: '{}'",
                    equation, tokens[i],
                    (last_used == npos) ? "n/a" : tokens[last_used]);
            }
            if (!kin || (kin->kineticsSpeciesIndex(species) == npos
                         && stoich != -1 && species != "M"))
            {
                R.setValid(false);
            }

            if (reactants) {
                R.reactants[species] += stoich;
            } else {
                R.products[species] += stoich;
            }

            last_used = i;
        }

        // Tokens after this point are part of the products string
        if (tokens[i] == "<=>" || tokens[i] == "=") {
            R.reversible = true;
            reactants = false;
        } else if (tokens[i] == "=>") {
            R.reversible = false;
            reactants = false;
        }
    }
}

void setupReaction(Reaction& R, const AnyMap& node, const Kinetics& kin)
{
    parseReactionEquation(R, node["equation"].asString(), node, &kin);
    // Non-stoichiometric reaction orders
    if (node.hasKey("orders")) {
        for (const auto& order : node["orders"].asMap<double>()) {
            R.orders[order.first] = order.second;
            if (kin.kineticsSpeciesIndex(order.first) == npos) {
                R.setValid(false);
            }
        }
    }

    //Flags
    R.id = node.getString("id", "");
    R.duplicate = node.getBool("duplicate", false);
    R.allow_negative_orders = node.getBool("negative-orders", false);
    R.allow_nonreactant_orders = node.getBool("nonreactant-orders", false);

    R.input = node;
    R.calculateRateCoeffUnits(kin);
}

void setupElementaryReaction(ElementaryReaction2& R, const XML_Node& rxn_node)
{
    const XML_Node& rc_node = rxn_node.child("rateCoeff");
    if (rc_node.hasChild("Arrhenius")) {
        R.rate = readArrhenius(rc_node.child("Arrhenius"));
    } else if (rc_node.hasChild("Arrhenius_ExchangeCurrentDensity")) {
        R.rate = readArrhenius(rc_node.child("Arrhenius_ExchangeCurrentDensity"));
    } else {
        throw CanteraError("setupElementaryReaction", "Couldn't find Arrhenius node");
    }
    if (rxn_node["negative_A"] == "yes") {
        R.allow_negative_pre_exponential_factor = true;
    }
    if (rxn_node["negative_orders"] == "yes") {
        R.allow_negative_orders = true;
    }
    if (rxn_node["nonreactant_orders"] == "yes") {
        R.allow_nonreactant_orders = true;
    }
    setupReaction(R, rxn_node);
}

void setupElementaryReaction(ElementaryReaction2& R, const AnyMap& node,
                             const Kinetics& kin)
{
    setupReaction(R, node, kin);
    R.allow_negative_pre_exponential_factor = node.getBool("negative-A", false);
    R.rate = readArrhenius(R, node["rate-constant"], kin, node.units());
}

void setupThreeBodyReaction(ThreeBodyReaction2& R, const XML_Node& rxn_node)
{
    readEfficiencies(R.third_body, rxn_node.child("rateCoeff"));
    setupElementaryReaction(R, rxn_node);
    if (R.third_body.efficiencies.size() == 0) {
        detectEfficiencies(R);
    }
}

void setupThreeBodyReaction(ThreeBodyReaction2& R, const AnyMap& node,
                            const Kinetics& kin)
{
    setupElementaryReaction(R, node, kin);
    if (R.reactants.count("M") != 1 || R.products.count("M") != 1) {
        if (!detectEfficiencies(R)) {
            throw InputFileError("setupThreeBodyReaction", node["equation"],
                "Reaction equation '{}' does not contain third body 'M'",
                node["equation"].asString());
        }
    } else {
        R.reactants.erase("M");
        R.products.erase("M");
        readEfficiencies(R.third_body, node);
    }
}

void setupFalloffReaction(FalloffReaction2& R, const XML_Node& rxn_node)
{
    XML_Node& rc_node = rxn_node.child("rateCoeff");
    std::vector<XML_Node*> rates = rc_node.getChildren("Arrhenius");
    int nLow = 0;
    int nHigh = 0;
    for (size_t i = 0; i < rates.size(); i++) {
        XML_Node& node = *rates[i];
        if (node["name"] == "") {
            R.high_rate = readArrhenius(node);
            nHigh++;
        } else if (node["name"] == "k0") {
            R.low_rate = readArrhenius(node);
            nLow++;
        } else {
            throw CanteraError("setupFalloffReaction", "Found an Arrhenius XML "
                "node with an unexpected type '" + node["name"] + "'");
        }
    }
    if (nLow != 1 || nHigh != 1) {
        throw CanteraError("setupFalloffReaction", "Did not find the correct "
            "number of Arrhenius rate expressions");
    }
    if (rxn_node["negative_A"] == "yes") {
        R.allow_negative_pre_exponential_factor = true;
    }
    readFalloff(R, rc_node);
    readEfficiencies(R.third_body, rc_node);
    setupReaction(R, rxn_node);
}

void setupFalloffReaction(FalloffReaction2& R, const AnyMap& node,
                          const Kinetics& kin)
{
    setupReaction(R, node, kin);
    // setupReaction sets the stoichiometric coefficient for the falloff third
    // body to -1.
    std::string third_body;
    for (auto& reactant : R.reactants) {
        if (reactant.second == -1 && ba::starts_with(reactant.first, "(+")) {
            third_body = reactant.first;
            break;
        }
    }

    // Equation must contain a third body, and it must appear on both sides
    if (third_body == "") {
        throw InputFileError("setupFalloffReaction", node["equation"],
            "Reactants for reaction '{}' do not contain a pressure-dependent "
            "third body", node["equation"].asString());
    } else if (R.products.count(third_body) == 0) {
        throw InputFileError("setupFalloffReaction", node["equation"],
            "Unable to match third body '{}' in reactants and products of "
            "reaction '{}'", third_body, node["equation"].asString());
    }

    // Remove the dummy species
    R.reactants.erase(third_body);
    R.products.erase(third_body);

    R.allow_negative_pre_exponential_factor = node.getBool("negative-A", false);
    if (third_body == "(+M)") {
        readEfficiencies(R.third_body, node);
    } else {
        // Specific species is listed as the third body
        R.third_body.default_efficiency = 0;
        R.third_body.efficiencies[third_body.substr(2, third_body.size() - 3)] = 1.0;
    }

    R.low_rate = readArrhenius(R, node["low-P-rate-constant"], kin,
                                node.units(), 1);
    R.high_rate = readArrhenius(R, node["high-P-rate-constant"], kin,
                                node.units());

    readFalloff(R, node);
}

void setupChemicallyActivatedReaction(ChemicallyActivatedReaction2& R,
                                      const XML_Node& rxn_node)
{
    XML_Node& rc_node = rxn_node.child("rateCoeff");
    std::vector<XML_Node*> rates = rc_node.getChildren("Arrhenius");
    int nLow = 0;
    int nHigh = 0;
    for (size_t i = 0; i < rates.size(); i++) {
        XML_Node& node = *rates[i];
        if (node["name"] == "kHigh") {
            R.high_rate = readArrhenius(node);
            nHigh++;
        } else if (node["name"] == "") {
            R.low_rate = readArrhenius(node);
            nLow++;
        } else {
            throw CanteraError("setupChemicallyActivatedReaction", "Found an "
                "Arrhenius XML node with an unexpected type '" + node["name"] + "'");
        }
    }
    if (nLow != 1 || nHigh != 1) {
        throw CanteraError("setupChemicallyActivatedReaction", "Did not find "
            "the correct number of Arrhenius rate expressions");
    }
    readFalloff(R, rc_node);
    readEfficiencies(R.third_body, rc_node);
    setupReaction(R, rxn_node);
}

void setupPlogReaction(PlogReaction2& R, const XML_Node& rxn_node)
{
    XML_Node& rc = rxn_node.child("rateCoeff");
    std::multimap<double, Arrhenius2> rates;
    for (size_t m = 0; m < rc.nChildren(); m++) {
        const XML_Node& node = rc.child(m);
        rates.insert({getFloat(node, "P", "toSI"), readArrhenius(node)});
    }
    R.rate = Plog(rates);
    setupReaction(R, rxn_node);
}

void setupPlogReaction(PlogReaction2& R, const AnyMap& node, const Kinetics& kin)
{
    setupReaction(R, node, kin);
    std::multimap<double, Arrhenius2> rates;
    for (const auto& rate : node.at("rate-constants").asVector<AnyMap>()) {
        rates.insert({rate.convert("P", "Pa"),
                      readArrhenius(R, AnyValue(rate), kin, node.units())});
    }
    R.rate = Plog(rates);
}

void PlogReaction2::validate()
{
    Reaction::validate();
    rate.validate(equation());
}

void setupChebyshevReaction(ChebyshevReaction2& R, const XML_Node& rxn_node)
{
    XML_Node& rc = rxn_node.child("rateCoeff");
    const XML_Node& coeff_node = rc.child("floatArray");
    size_t nP = atoi(coeff_node["degreeP"].c_str());
    size_t nT = atoi(coeff_node["degreeT"].c_str());

    vector_fp coeffs_flat;
    getFloatArray(rc, coeffs_flat, false);
    Array2D coeffs(nT, nP);
    for (size_t t = 0; t < nT; t++) {
        for (size_t p = 0; p < nP; p++) {
            coeffs(t,p) = coeffs_flat[nP*t + p];
        }
    }
    R.rate = ChebyshevRate(getFloat(rc, "Tmin", "toSI"),
                           getFloat(rc, "Tmax", "toSI"),
                           getFloat(rc, "Pmin", "toSI"),
                           getFloat(rc, "Pmax", "toSI"),
                           coeffs);
    setupReaction(R, rxn_node);
}

void setupChebyshevReaction(ChebyshevReaction2&R, const AnyMap& node,
                            const Kinetics& kin)
{
    setupReaction(R, node, kin);
    R.reactants.erase("(+M)"); // remove optional third body notation
    R.products.erase("(+M)");
    const auto& T_range = node["temperature-range"].asVector<AnyValue>(2);
    const auto& P_range = node["pressure-range"].asVector<AnyValue>(2);
    auto& vcoeffs = node["data"].asVector<vector_fp>();
    Array2D coeffs(vcoeffs.size(), vcoeffs[0].size());
    for (size_t i = 0; i < coeffs.nRows(); i++) {
        if (vcoeffs[i].size() != vcoeffs[0].size()) {
            throw InputFileError("setupChebyshevReaction", node["data"],
                "Inconsistent number of coefficients in row {} of matrix", i + 1);
        }
        for (size_t j = 0; j < coeffs.nColumns(); j++) {
            coeffs(i, j) = vcoeffs[i][j];
        }
    }
    const UnitSystem& units = node.units();
    coeffs(0, 0) += std::log10(units.convertTo(1.0, R.rate_units));
    R.rate = ChebyshevRate(units.convert(T_range[0], "K"),
                           units.convert(T_range[1], "K"),
                           units.convert(P_range[0], "Pa"),
                           units.convert(P_range[1], "Pa"),
                           coeffs);
}

void setupInterfaceReaction(InterfaceReaction2& R, const XML_Node& rxn_node)
{
    XML_Node& arr = rxn_node.child("rateCoeff").child("Arrhenius");
    if (caseInsensitiveEquals(arr["type"], "stick")) {
        R.is_sticking_coefficient = true;
        R.sticking_species = arr["species"];

        if (caseInsensitiveEquals(arr["motz_wise"], "true")) {
            R.use_motz_wise_correction = true;
        } else if (caseInsensitiveEquals(arr["motz_wise"], "false")) {
            R.use_motz_wise_correction = false;
        } else {
            // Default value for all reactions
            XML_Node* parent = rxn_node.parent();
            if (parent && parent->name() == "reactionData"
                && caseInsensitiveEquals((*parent)["motz_wise"], "true")) {
                R.use_motz_wise_correction = true;
            }
        }
    }
    std::vector<XML_Node*> cov = arr.getChildren("coverage");
    for (const auto& node : cov) {
        CoverageDependency& cdep = R.coverage_deps[node->attrib("species")];
        cdep.a = getFloat(*node, "a", "toSI");
        cdep.m = getFloat(*node, "m");
        cdep.E = getFloat(*node, "e", "actEnergy") / GasConstant;
    }
    setupElementaryReaction(R, rxn_node);
}

void setupInterfaceReaction(InterfaceReaction2& R, const AnyMap& node,
                            const Kinetics& kin)
{
    setupReaction(R, node, kin);
    R.allow_negative_pre_exponential_factor = node.getBool("negative-A", false);

    if (node.hasKey("rate-constant")) {
        R.rate = readArrhenius(R, node["rate-constant"], kin, node.units());
    } else if (node.hasKey("sticking-coefficient")) {
        R.is_sticking_coefficient = true;
        R.rate = readArrhenius(R, node["sticking-coefficient"], kin, node.units());
        R.use_motz_wise_correction = node.getBool("Motz-Wise",
            kin.thermo().input().getBool("Motz-Wise", false));
        R.sticking_species = node.getString("sticking-species", "");
    } else {
        throw InputFileError("setupInterfaceReaction", node,
            "Reaction must include either a 'rate-constant' or"
            " 'sticking-coefficient' node.");
    }

    if (node.hasKey("coverage-dependencies")) {
        for (const auto& item : node["coverage-dependencies"].as<AnyMap>()) {
            double a, E, m;
            if (item.second.is<AnyMap>()) {
                auto& cov_map = item.second.as<AnyMap>();
                a = cov_map["a"].asDouble();
                m = cov_map["m"].asDouble();
                E = node.units().convertActivationEnergy(cov_map["E"], "K");
            } else {
                auto& cov_vec = item.second.asVector<AnyValue>(3);
                a = cov_vec[0].asDouble();
                m = cov_vec[1].asDouble();
                E = node.units().convertActivationEnergy(cov_vec[2], "K");
            }
            R.coverage_deps[item.first] = CoverageDependency(a, E, m);
        }
    }
}

void setupElectrochemicalReaction(ElectrochemicalReaction2& R,
                                  const XML_Node& rxn_node)
{
    XML_Node& rc = rxn_node.child("rateCoeff");
    std::string rc_type = toLowerCopy(rc["type"]);
    if (rc_type == "exchangecurrentdensity") {
        R.exchange_current_density_formulation = true;
    } else if (rc_type != "" && rc_type != "arrhenius") {
        throw CanteraError("setupElectrochemicalReaction",
            "Unknown rate coefficient type: '" + rc_type + "'");
    }
    if (rc.hasChild("Arrhenius_ExchangeCurrentDensity")) {
        R.exchange_current_density_formulation = true;
    }

    if (rc.hasChild("electrochem") && rc.child("electrochem").hasAttrib("beta")) {
        R.beta = fpValueCheck(rc.child("electrochem")["beta"]);
    }

    setupInterfaceReaction(R, rxn_node);

    // Override orders based on the <orders> node
    if (rxn_node.hasChild("orders")) {
        Composition orders = parseCompString(rxn_node.child("orders").value());
        for (const auto& order : orders) {
            R.orders[order.first] = order.second;
        }
    }
}

void setupElectrochemicalReaction(ElectrochemicalReaction2& R,
                                  const AnyMap& node, const Kinetics& kin)
{
    setupInterfaceReaction(R, node, kin);
    R.beta = node.getDouble("beta", 0.5);
    R.exchange_current_density_formulation = node.getBool(
        "exchange-current-density-formulation", false);
}

std::vector<shared_ptr<Reaction> > getReactions(const XML_Node& node)
{
    std::vector<shared_ptr<Reaction> > all_reactions;
    for (const auto& rxnnode : node.child("reactionData").getChildren("reaction")) {
        all_reactions.push_back(newReaction(*rxnnode));
    }
    return all_reactions;
}

std::vector<shared_ptr<Reaction>> getReactions(const AnyValue& items,
                                               Kinetics& kinetics)
{
    std::vector<shared_ptr<Reaction>> all_reactions;
    for (const auto& node : items.asVector<AnyMap>()) {
        shared_ptr<Reaction> R(newReaction(node, kinetics));
        R->validate();
        R->validate(kinetics);
        if (R->valid() && R->checkSpecies(kinetics)) {
            all_reactions.emplace_back(R);
        }
    }
    return all_reactions;
}

}
