/**
 *  @file Reaction.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/Reaction.h"
#include "cantera/kinetics/FalloffFactory.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/base/ctml.h"
#include "cantera/base/Array.h"
#include "cantera/base/AnyMap.h"
#include <sstream>
#include <set>

#include <boost/algorithm/string.hpp>

namespace ba = boost::algorithm;

namespace Cantera
{

Reaction::Reaction(int type)
    : reaction_type(type)
    , reversible(true)
    , duplicate(false)
    , allow_nonreactant_orders(false)
    , allow_negative_orders(false)
{
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
{
}

void Reaction::validate()
{
    if (!allow_nonreactant_orders) {
        for (const auto& order : orders) {
            if (reactants.find(order.first) == reactants.end()) {
                throw CanteraError("Reaction::validate", "Reaction order "
                    "specified for non-reactant species '" + order.first + "'");
            }
        }
    }

    if (!allow_negative_orders) {
        for (const auto& order : orders) {
            if (order.second < 0.0) {
                throw CanteraError("Reaction::validate", "Negative reaction "
                    "order specified for species '" + order.first + "'");
            }
        }
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

ElementaryReaction::ElementaryReaction(const Composition& reactants_,
                                       const Composition products_,
                                       const Arrhenius& rate_)
    : Reaction(ELEMENTARY_RXN, reactants_, products_)
    , rate(rate_)
    , allow_negative_pre_exponential_factor(false)
{
}

ElementaryReaction::ElementaryReaction()
    : Reaction(ELEMENTARY_RXN)
    , allow_negative_pre_exponential_factor(false)
{
}

void ElementaryReaction::validate()
{
    Reaction::validate();
    if (!allow_negative_pre_exponential_factor &&
        rate.preExponentialFactor() < 0) {
        throw CanteraError("ElementaryReaction::validate",
            "Undeclared negative pre-exponential factor found in reaction '"
            + equation() + "'");
    }
}

ThirdBody::ThirdBody(double default_eff)
    : default_efficiency(default_eff)
{
}

ThreeBodyReaction::ThreeBodyReaction()
{
    reaction_type = THREE_BODY_RXN;
}

ThreeBodyReaction::ThreeBodyReaction(const Composition& reactants_,
                                     const Composition& products_,
                                     const Arrhenius& rate_,
                                     const ThirdBody& tbody)
    : ElementaryReaction(reactants_, products_, rate_)
    , third_body(tbody)
{
    reaction_type = THREE_BODY_RXN;
}

std::string ThreeBodyReaction::reactantString() const {
    return ElementaryReaction::reactantString() + " + M";
}

std::string ThreeBodyReaction::productString() const {
    return ElementaryReaction::productString() + " + M";
}

FalloffReaction::FalloffReaction()
    : Reaction(FALLOFF_RXN)
    , falloff(new Falloff())
    , allow_negative_pre_exponential_factor(false)
{
}

FalloffReaction::FalloffReaction(
        const Composition& reactants_, const Composition& products_,
        const Arrhenius& low_rate_, const Arrhenius& high_rate_,
        const ThirdBody& tbody)
    : Reaction(FALLOFF_RXN, reactants_, products_)
    , low_rate(low_rate_)
    , high_rate(high_rate_)
    , third_body(tbody)
    , falloff(new Falloff())
{
}

std::string FalloffReaction::reactantString() const {
    if (third_body.default_efficiency == 0 &&
        third_body.efficiencies.size() == 1) {
        return Reaction::reactantString() + " (+" +
            third_body.efficiencies.begin()->first + ")";
    } else {
        return Reaction::reactantString() + " (+M)";
    }
}

std::string FalloffReaction::productString() const {
    if (third_body.default_efficiency == 0 &&
        third_body.efficiencies.size() == 1) {
        return Reaction::productString() + " (+" +
            third_body.efficiencies.begin()->first + ")";
    } else {
        return Reaction::productString() + " (+M)";
    }
}

void FalloffReaction::validate() {
    Reaction::validate();
    if (!allow_negative_pre_exponential_factor &&
        (low_rate.preExponentialFactor() < 0 ||
         high_rate.preExponentialFactor() < 0)) {
        throw CanteraError("FalloffReaction::validate", "Negative "
            "pre-exponential factor found for reaction '" + equation() + "'");
    }
    if (low_rate.preExponentialFactor() * high_rate.preExponentialFactor() < 0) {
        throw CanteraError("FalloffReaction::validate", "High and "
            "low rate pre-exponential factors must have the same sign."
            "Reaction: '{}'", equation());
    }
}

ChemicallyActivatedReaction::ChemicallyActivatedReaction()
{
    reaction_type = CHEMACT_RXN;
}

ChemicallyActivatedReaction::ChemicallyActivatedReaction(
        const Composition& reactants_, const Composition& products_,
        const Arrhenius& low_rate_, const Arrhenius& high_rate_,
        const ThirdBody& tbody)
    : FalloffReaction(reactants_, products_, low_rate_, high_rate_, tbody)
{
    reaction_type = CHEMACT_RXN;
}

PlogReaction::PlogReaction()
    : Reaction(PLOG_RXN)
{
}

PlogReaction::PlogReaction(const Composition& reactants_,
                           const Composition& products_, const Plog& rate_)
    : Reaction(PLOG_RXN, reactants_, products_)
    , rate(rate_)
{
}

ChebyshevReaction::ChebyshevReaction()
    : Reaction(CHEBYSHEV_RXN)
{
}

ChebyshevReaction::ChebyshevReaction(const Composition& reactants_,
                                     const Composition& products_,
                                     const ChebyshevRate& rate_)
    : Reaction(CHEBYSHEV_RXN, reactants_, products_)
    , rate(rate_)
{
}

InterfaceReaction::InterfaceReaction()
    : is_sticking_coefficient(false)
    , use_motz_wise_correction(false)
{
    reaction_type = INTERFACE_RXN;
}

InterfaceReaction::InterfaceReaction(const Composition& reactants_,
                                     const Composition& products_,
                                     const Arrhenius& rate_,
                                     bool isStick)
    : ElementaryReaction(reactants_, products_, rate_)
    , is_sticking_coefficient(isStick)
    , use_motz_wise_correction(false)
{
    reaction_type = INTERFACE_RXN;
}

ElectrochemicalReaction::ElectrochemicalReaction()
    : film_resistivity(0.0)
    , beta(0.5)
    , exchange_current_density_formulation(false)
{
}

ElectrochemicalReaction::ElectrochemicalReaction(const Composition& reactants_,
                                                 const Composition& products_,
                                                 const Arrhenius& rate_)
    : InterfaceReaction(reactants_, products_, rate_)
    , film_resistivity(0.0)
    , beta(0.5)
    , exchange_current_density_formulation(false)
{
}

Arrhenius readArrhenius(const XML_Node& arrhenius_node)
{
    return Arrhenius(getFloat(arrhenius_node, "A", "toSI"),
                     getFloat(arrhenius_node, "b"),
                     getFloat(arrhenius_node, "E", "actEnergy") / GasConstant);
}

Units rateCoeffUnits(const Reaction& R, const Kinetics& kin,
                     int pressure_dependence=0)
{
    if (R.reaction_type == INVALID_RXN) {
        // If a reaction is invalid because of missing species in the Kinetics
        // object, determining the units of the rate coefficient is impossible.
        return Units();
    } else if (R.reaction_type == INTERFACE_RXN
               && dynamic_cast<const InterfaceReaction&>(R).is_sticking_coefficient) {
        // Sticking coefficients are dimensionless
        return Units();
    }

    // Determine the units of the rate coefficient
    Units rxn_phase_units = kin.thermo(kin.reactionPhaseIndex()).standardConcentrationUnits();
    Units rcUnits = rxn_phase_units;
    rcUnits *= Units(1.0, 0, 0, -1);
    for (const auto& order : R.orders) {
        const auto& phase = kin.speciesPhase(order.first);
        rcUnits *= phase.standardConcentrationUnits().pow(-order.second);
    }
    for (const auto& stoich : R.reactants) {
        // Order for each reactant is the reactant stoichiometric coefficient,
        // unless already overridden by user-specified orders
        if (stoich.first == "M") {
            rcUnits *= rxn_phase_units.pow(-1);
        } else if (R.orders.find(stoich.first) == R.orders.end()) {
            const auto& phase = kin.speciesPhase(stoich.first);
            rcUnits *= phase.standardConcentrationUnits().pow(-stoich.second);
        }
    }

    // Incorporate pressure dependence for low-pressure falloff and high-
    // pressure chemically-activated reaction limits
    rcUnits *= rxn_phase_units.pow(-pressure_dependence);
    return rcUnits;
}

Arrhenius readArrhenius(const Reaction& R, const AnyValue& rate,
                        const Kinetics& kin, const UnitSystem& units,
                        int pressure_dependence=0)
{
    double A, b, Ta;
    Units rc_units = rateCoeffUnits(R, kin, pressure_dependence);
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
    return Arrhenius(A, b, Ta);
}

//! Parse falloff parameters, given a rateCoeff node
/*!
 * @verbatim
 <falloff type="Troe"> 0.5 73.2 5000. 9999. </falloff>
 @endverbatim
*/
void readFalloff(FalloffReaction& R, const XML_Node& rc_node)
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
    } else {
        throw CanteraError("readFalloff", "Unrecognized falloff type: '{}'",
                           falloff["type"]);
    }
}

void readFalloff(FalloffReaction& R, const AnyMap& node)
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
    tbody.default_efficiency = node.getDouble("default-efficiency", 1.0);
    if (node.hasKey("efficiencies")) {
        tbody.efficiencies = node["efficiencies"].asMap<double>();
    }
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

void parseReactionEquation(Reaction& R, const AnyValue& equation,
                           const Kinetics& kin) {
    // Parse the reaction equation to determine participating species and
    // stoichiometric coefficients
    std::vector<std::string> tokens;
    tokenizeString(equation.asString(), tokens);
    tokens.push_back("+"); // makes parsing last species not a special case

    size_t last_used = npos; // index of last-used token
    bool reactants = true;
    for (size_t i = 1; i < tokens.size(); i++) {
        if (tokens[i] == "+" || ba::starts_with(tokens[i], "(+") ||
            tokens[i] == "<=>" || tokens[i] == "=" || tokens[i] == "=>") {
            std::string species = tokens[i-1];

            double stoich;
            if (last_used != npos && tokens[last_used] == "(+") {
                // Falloff third body with space, e.g. "(+ M)"
                species = "(+" + species;
                stoich = -1;
            } else if (last_used == i-1 && ba::starts_with(species, "(+")
                       && ba::ends_with(species, ")")) {
                // Falloff 3rd body written without space, e.g. "(+M)"
                stoich = -1;
            } else if (last_used == i-2) { // Species with no stoich. coefficient
                stoich = 1.0;
            } else if (last_used == i-3) { // Stoich. coefficient and species
                try {
                    stoich = fpValueCheck(tokens[i-2]);
                } catch (CanteraError& err) {
                    throw InputFileError("parseReactionEquation", equation,
                        err.getMessage());
                }
            } else {
                throw InputFileError("parseReactionEquation", equation,
                    "Error parsing reaction string '{}'.\n"
                    "Current token: '{}'\nlast_used: '{}'",
                    equation.asString(), tokens[i],
                    (last_used == npos) ? "n/a" : tokens[last_used]);
            }
            if (kin.kineticsSpeciesIndex(species) == npos
                && stoich != -1 && species != "M") {
                R.reaction_type = INVALID_RXN;
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
    parseReactionEquation(R, node["equation"], kin);
    // Non-stoichiometric reaction orders
    std::map<std::string, double> orders;
    if (node.hasKey("orders")) {
        for (const auto& order : node["orders"].asMap<double>()) {
            R.orders[order.first] = order.second;
            if (kin.kineticsSpeciesIndex(order.first) == npos) {
                R.reaction_type = INVALID_RXN;
            }
        }
    }

    //Flags
    R.id = node.getString("id", "");
    R.duplicate = node.getBool("duplicate", false);
    R.allow_negative_orders = node.getBool("negative-orders", false);
    R.allow_nonreactant_orders = node.getBool("nonreactant-orders", false);

    R.input = node;
}

void setupElementaryReaction(ElementaryReaction& R, const XML_Node& rxn_node)
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

void setupElementaryReaction(ElementaryReaction& R, const AnyMap& node,
                             const Kinetics& kin)
{
    setupReaction(R, node, kin);
    R.allow_negative_pre_exponential_factor = node.getBool("negative-A", false);
    R.rate = readArrhenius(R, node["rate-constant"], kin, node.units());
}

void setupThreeBodyReaction(ThreeBodyReaction& R, const XML_Node& rxn_node)
{
    readEfficiencies(R.third_body, rxn_node.child("rateCoeff"));
    setupElementaryReaction(R, rxn_node);
}

void setupThreeBodyReaction(ThreeBodyReaction& R, const AnyMap& node,
                            const Kinetics& kin)
{
    setupElementaryReaction(R, node, kin);
    if (R.reactants.count("M") != 1 || R.products.count("M") != 1) {
        throw InputFileError("setupThreeBodyReaction", node["equation"],
            "Reaction equation '{}' does not contain third body 'M'",
            node["equation"].asString());
    }
    R.reactants.erase("M");
    R.products.erase("M");
    readEfficiencies(R.third_body, node);
}

void setupFalloffReaction(FalloffReaction& R, const XML_Node& rxn_node)
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

void setupFalloffReaction(FalloffReaction& R, const AnyMap& node,
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

    if (node["type"] == "falloff") {
        R.low_rate = readArrhenius(R, node["low-P-rate-constant"], kin,
                                   node.units(), 1);
        R.high_rate = readArrhenius(R, node["high-P-rate-constant"], kin,
                                    node.units());
    } else { // type == "chemically-activated"
        R.low_rate = readArrhenius(R, node["low-P-rate-constant"], kin,
                                   node.units());
        R.high_rate = readArrhenius(R, node["high-P-rate-constant"], kin,
                                    node.units(), -1);
    }

    readFalloff(R, node);
}

void setupChemicallyActivatedReaction(ChemicallyActivatedReaction& R,
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

void setupPlogReaction(PlogReaction& R, const XML_Node& rxn_node)
{
    XML_Node& rc = rxn_node.child("rateCoeff");
    std::multimap<double, Arrhenius> rates;
    for (size_t m = 0; m < rc.nChildren(); m++) {
        const XML_Node& node = rc.child(m);
        rates.insert({getFloat(node, "P", "toSI"), readArrhenius(node)});
    }
    R.rate = Plog(rates);
    setupReaction(R, rxn_node);
}

void setupPlogReaction(PlogReaction& R, const AnyMap& node, const Kinetics& kin)
{
    setupReaction(R, node, kin);
    std::multimap<double, Arrhenius> rates;
    for (const auto& rate : node.at("rate-constants").asVector<AnyMap>()) {
        rates.insert({rate.convert("P", "Pa"),
                      readArrhenius(R, AnyValue(rate), kin, node.units())});
    }
    R.rate = Plog(rates);
}

void PlogReaction::validate()
{
    Reaction::validate();
    rate.validate(equation());
}

void setupChebyshevReaction(ChebyshevReaction& R, const XML_Node& rxn_node)
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

void setupChebyshevReaction(ChebyshevReaction&R, const AnyMap& node,
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
    Units rcUnits = rateCoeffUnits(R, kin);
    coeffs(0, 0) += std::log10(units.convert(1.0, rcUnits));
    R.rate = ChebyshevRate(units.convert(T_range[0], "K"),
                           units.convert(T_range[1], "K"),
                           units.convert(P_range[0], "Pa"),
                           units.convert(P_range[1], "Pa"),
                           coeffs);
}

void setupInterfaceReaction(InterfaceReaction& R, const XML_Node& rxn_node)
{
    if (caseInsensitiveEquals(rxn_node["type"], "global")) {
        R.reaction_type = GLOBAL_RXN;
    }
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

void setupInterfaceReaction(InterfaceReaction& R, const AnyMap& node,
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

void setupElectrochemicalReaction(ElectrochemicalReaction& R,
                                  const XML_Node& rxn_node)
{
    // Fix reaction_type for some specialized reaction types
    std::string type = toLowerCopy(rxn_node["type"]);
    if (type == "butlervolmer") {
        R.reaction_type = BUTLERVOLMER_RXN;
        warn_deprecated("reaction type 'ButlerVolmer'",
            "To be removed after Cantera 2.5.");
    } else if (type == "butlervolmer_noactivitycoeffs") {
        R.reaction_type = BUTLERVOLMER_NOACTIVITYCOEFFS_RXN;
        warn_deprecated("reaction type 'butlervolmer_noactivitycoeffs'",
            "To be removed after Cantera 2.5.");
    } else if (type == "surfaceaffinity") {
        R.reaction_type = SURFACEAFFINITY_RXN;
        warn_deprecated("reaction type 'surfaceaffinity'",
            "To be removed after Cantera 2.5.");
    } else if (type == "global") {
        R.reaction_type = GLOBAL_RXN;
        warn_deprecated("reaction type 'global'",
            "To be removed after Cantera 2.5.");
    }

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

    if (rxn_node.hasChild("filmResistivity")) {
        warn_deprecated("reaction filmResistivity",
            "Not implemented. To be removed after Cantera 2.5.");
    }
    getOptionalFloat(rxn_node, "filmResistivity", R.film_resistivity);
    setupInterfaceReaction(R, rxn_node);

    // For Butler Volmer reactions, install the orders for the exchange current
    if (R.reaction_type == BUTLERVOLMER_NOACTIVITYCOEFFS_RXN ||
        R.reaction_type == BUTLERVOLMER_RXN) {
        if (!R.reversible) {
            throw CanteraError("setupElectrochemicalReaction",
                "A Butler-Volmer reaction must be reversible");
        }

        R.orders.clear();
        // Reaction orders based on species stoichiometric coefficients
        R.allow_nonreactant_orders = true;
        for (const auto& sp : R.reactants) {
            R.orders[sp.first] += sp.second * (1.0 - R.beta);
        }
        for (const auto& sp : R.products) {
            R.orders[sp.first] += sp.second * R.beta;
        }
    }

    // For affinity reactions, fill in the global reaction formulation terms
    if (rxn_node.hasChild("reactionOrderFormulation")) {
        warn_deprecated("reactionOrderFormulation",
            "To be removed after Cantera 2.5.");
        Composition initial_orders = R.orders;
        R.orders.clear();
        R.allow_nonreactant_orders = true;
        const XML_Node& rof_node = rxn_node.child("reactionOrderFormulation");
        if (caseInsensitiveEquals(rof_node["model"], "reactantorders")) {
            R.orders = initial_orders;
        } else if (caseInsensitiveEquals(rof_node["model"], "zeroorders")) {
            for (const auto& sp : R.reactants) {
                R.orders[sp.first] = 0.0;
            }
        } else if (caseInsensitiveEquals(rof_node["model"], "butlervolmerorders")) {
            // Reaction orders based on provided reaction orders
            for (const auto& sp : R.reactants) {
                double c = getValue(initial_orders, sp.first, sp.second);
                R.orders[sp.first] += c * (1.0 - R.beta);
            }
            for (const auto& sp : R.products) {
                double c = getValue(initial_orders, sp.first, sp.second);
                R.orders[sp.first] += c * R.beta;
            }
        } else {
            throw CanteraError("setupElectrochemicalReaction", "unknown model "
                    "for reactionOrderFormulation XML_Node: '" +
                    rof_node["model"] + "'");
       }
    }

    // Override orders based on the <orders> node
    if (rxn_node.hasChild("orders")) {
        Composition orders = parseCompString(rxn_node.child("orders").value());
        for (const auto& order : orders) {
            R.orders[order.first] = order.second;
        }
    }
}

void setupElectrochemicalReaction(ElectrochemicalReaction& R,
                                  const AnyMap& node, const Kinetics& kin)
{
    setupInterfaceReaction(R, node, kin);
    R.beta = node.getDouble("beta", 0.5);
    R.exchange_current_density_formulation = node.getBool(
        "exchange-current-density-formulation", false);
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

shared_ptr<Reaction> newReaction(const XML_Node& rxn_node)
{
    std::string type = toLowerCopy(rxn_node["type"]);

    // Modify the reaction type for interface reactions which contain
    // electrochemical reaction data
    if (rxn_node.child("rateCoeff").hasChild("electrochem")
        && (type == "edge" || type == "surface")) {
        type = "electrochemical";
    }

    // Create a new Reaction object of the appropriate type
    if (type == "elementary" || type == "arrhenius" || type == "") {
        auto R = make_shared<ElementaryReaction>();
        setupElementaryReaction(*R, rxn_node);
        return R;
    } else if (type == "threebody" || type == "three_body") {
        auto R = make_shared<ThreeBodyReaction>();
        setupThreeBodyReaction(*R, rxn_node);
        return R;
    } else if (type == "falloff") {
        auto R = make_shared<FalloffReaction>();
        setupFalloffReaction(*R, rxn_node);
        return R;
    } else if (type == "chemact" || type == "chemically_activated") {
        auto R = make_shared<ChemicallyActivatedReaction>();
        setupChemicallyActivatedReaction(*R, rxn_node);
        return R;
    } else if (type == "plog" || type == "pdep_arrhenius") {
        auto R = make_shared<PlogReaction>();
        setupPlogReaction(*R, rxn_node);
        return R;
    } else if (type == "chebyshev") {
        auto R = make_shared<ChebyshevReaction>();
        setupChebyshevReaction(*R, rxn_node);
        return R;
    } else if (type == "interface" || type == "surface" || type == "edge" ||
               type == "global") {
        auto R = make_shared<InterfaceReaction>();
        setupInterfaceReaction(*R, rxn_node);
        return R;
    } else if (type == "electrochemical" ||
               type == "butlervolmer_noactivitycoeffs" ||
               type == "butlervolmer" ||
               type == "surfaceaffinity") {
        auto R = make_shared<ElectrochemicalReaction>();
        setupElectrochemicalReaction(*R, rxn_node);
        return R;
    } else {
        throw CanteraError("newReaction",
            "Unknown reaction type '" + rxn_node["type"] + "'");
    }
}

unique_ptr<Reaction> newReaction(const AnyMap& node, const Kinetics& kin)
{
    std::string type = "elementary";
    if (node.hasKey("type")) {
        type = node["type"].asString();
    }

    if (kin.thermo(kin.reactionPhaseIndex()).nDim() < 3) {
        // See if this is an electrochemical reaction
        Reaction testReaction(0);
        parseReactionEquation(testReaction, node["equation"], kin);
        if (isElectrochemicalReaction(testReaction, kin)) {
            unique_ptr<ElectrochemicalReaction> R(new ElectrochemicalReaction());
            setupElectrochemicalReaction(*R, node, kin);
            return unique_ptr<Reaction>(move(R));
        } else {
            unique_ptr<InterfaceReaction> R(new InterfaceReaction());
            setupInterfaceReaction(*R, node, kin);
            return unique_ptr<Reaction>(move(R));
        }
    }

    if (type == "elementary") {
        unique_ptr<ElementaryReaction> R(new ElementaryReaction());
        setupElementaryReaction(*R, node, kin);
        return unique_ptr<Reaction>(move(R));
    } else if (type == "three-body") {
        unique_ptr<ThreeBodyReaction> R(new ThreeBodyReaction());
        setupThreeBodyReaction(*R, node, kin);
        return unique_ptr<Reaction>(move(R));
    } else if (type == "falloff") {
        unique_ptr<FalloffReaction> R(new FalloffReaction());
        setupFalloffReaction(*R, node, kin);
        return unique_ptr<Reaction>(move(R));
    } else if (type == "chemically-activated") {
        unique_ptr<ChemicallyActivatedReaction> R(new ChemicallyActivatedReaction());
        setupFalloffReaction(*R, node, kin);
        return unique_ptr<Reaction>(move(R));
    } else if (type == "pressure-dependent-Arrhenius") {
        unique_ptr<PlogReaction> R(new PlogReaction());
        setupPlogReaction(*R, node, kin);
        return unique_ptr<Reaction>(move(R));
    } else if (type == "Chebyshev") {
        unique_ptr<ChebyshevReaction> R(new ChebyshevReaction());
        setupChebyshevReaction(*R, node, kin);
        return unique_ptr<Reaction>(move(R));
    } else {
        throw InputFileError("newReaction", node["type"],
            "Unknown reaction type '{}'", type);
    }
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
        if (R->reaction_type != INVALID_RXN) {
            all_reactions.emplace_back(R);
        } else if (!kinetics.skipUndeclaredSpecies()) {
            throw InputFileError("getReactions", node,
                "Reaction '{}' contains undeclared species.", R->equation());
        }
    };
    return all_reactions;
}

}
