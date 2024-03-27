/**
 *  @file Reaction.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/Reaction.h"
#include "cantera/kinetics/ReactionRateFactory.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/base/Array.h"
#include "cantera/base/AnyMap.h"
#include "cantera/base/utilities.h"
#include "cantera/base/stringUtils.h"
#include <boost/algorithm/string/predicate.hpp>
#include <sstream>

#include <boost/algorithm/string.hpp>

namespace ba = boost::algorithm;

namespace Cantera
{

Reaction::Reaction(const Composition& reactants_,
                   const Composition& products_,
                   shared_ptr<ReactionRate> rate_,
                   shared_ptr<ThirdBody> tbody_)
    : reactants(reactants_)
    , products(products_)
    , m_from_composition(true)
    , m_third_body(tbody_)
{
    if (reactants.count("M") || products.count("M")) {
        throw CanteraError("Reaction::Reaction",
            "Third body 'M' must not be included in either reactant or product maps.");
    }
    setRate(rate_);

    // set flags ensuring correct serialization output
    Composition third;
    for (const auto& [name, stoich] : reactants) {
        if (products.count(name)) {
            third[name] = products.at(name) - stoich;
        }
    }
    if (tbody_) {
        string name = tbody_->name();
        if (reactants.count(name) && products.count(name)) {
            throw CanteraError("Reaction::Reaction",
                "'{}' not acting as third body collider must not be included in both "
                "reactant and product maps.", name);
        }
        if (name != "M") {
            m_third_body->explicit_3rd = true;
        }
    } else if (!tbody_ && third.size() == 1) {
        // implicit third body
        string name = third.begin()->first;
        m_third_body = make_shared<ThirdBody>(name);
        if (name != "M") {
            m_third_body->explicit_3rd = true;
        }
    }
    check();
}

Reaction::Reaction(const string& equation,
                   shared_ptr<ReactionRate> rate_,
                   shared_ptr<ThirdBody> tbody_)
    : m_third_body(tbody_)
{
    setRate(rate_);
    setEquation(equation);
    if (m_third_body && m_third_body->name() != "M") {
        m_third_body->explicit_3rd = true;
    }
    check();
}

Reaction::Reaction(const AnyMap& node, const Kinetics& kin)
{
    string rate_type = node.getString("type", "Arrhenius");
    if (!kin.nPhases()) {
        throw InputFileError("Reaction", node,
            "Cannot instantiate Reaction with empty Kinetics object.");
    }

    setParameters(node, kin);
    size_t nDim = kin.thermo(0).nDim();
    if (!valid()) {
        // If the reaction isn't valid (for example, contains undefined species),
        // setting up the rate constant won't work
        return;
    }
    if (nDim == 3) {
        if (ba::starts_with(rate_type, "three-body-")) {
            AnyMap rateNode = node;
            rateNode["type"] = rate_type.substr(11, rate_type.size() - 11);
            setRate(newReactionRate(rateNode, calculateRateCoeffUnits(kin)));
        } else {
            setRate(newReactionRate(node, calculateRateCoeffUnits(kin)));
        }
    } else {
        AnyMap rateNode = node;
        if (rateNode.hasKey("rate-constant")) {
            if (!ba::starts_with(rate_type, "interface-")) {
                rateNode["type"] = "interface-" + rate_type;
            }
        } else if (node.hasKey("sticking-coefficient")) {
            if (!ba::starts_with(rate_type, "sticking-")) {
                rateNode["type"] = "sticking-" + rate_type;
            }
        } else {
            throw InputFileError("Reaction::Reaction", input,
                "Unable to infer interface reaction type.");
        }
        setRate(newReactionRate(rateNode, calculateRateCoeffUnits(kin)));
    }
    check();
}

void Reaction::check()
{
    if (!allow_nonreactant_orders) {
        for (const auto& [name, order] : orders) {
            if (reactants.find(name) == reactants.end()) {
                throw InputFileError("Reaction::validate", input,
                    "Reaction order specified for non-reactant species '{}'", name);
           }
        }
    }

    if (!allow_negative_orders) {
        for (const auto& [name, order] : orders) {
            if (order < 0.0) {
                throw InputFileError("Reaction::validate", input,
                    "Negative reaction order specified for species '{}'", name);
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

    if (!m_rate) {
        return;
    }

    // Check reaction rate evaluator to ensure changes introduced after object
    // instantiation are considered.
    m_rate->check(equation());

    string rate_type = m_rate->type();
    if (m_third_body) {
        if (rate_type == "falloff" || rate_type == "chemically-activated") {
            if (m_third_body->mass_action && !m_from_composition) {
                throw InputFileError("Reaction::setRate", input,
                    "Third-body collider does not use '(+{})' notation.",
                    m_third_body->name());
            }
            m_third_body->mass_action = false;
        } else if (rate_type == "Chebyshev") {
            if (m_third_body->name() == "M") {
                warn_deprecated("Chebyshev reaction equation", input, "Specifying 'M' "
                    "in the reaction equation for Chebyshev reactions is deprecated.");
                m_third_body.reset();
            }
        } else if (rate_type == "pressure-dependent-Arrhenius") {
            if (m_third_body->name() == "M") {
                throw InputFileError("Reaction::setRate", input,
                    "Found superfluous '{}' in pressure-dependent-Arrhenius reaction.",
                    m_third_body->name());
            }
        }
    } else {
        if (rate_type == "falloff" || rate_type == "chemically-activated") {
            if (!m_from_composition) {
                throw InputFileError("Reaction::setRate", input,
                    "Reaction equation for falloff reaction '{}'\n does not "
                    "contain valid pressure-dependent third body", equation());
            }
            m_third_body = make_shared<ThirdBody>("(+M)");
        }
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
    if (!m_rate) {
        throw CanteraError("Reaction::getParameters",
            "Serialization of empty Reaction object is not supported.");
    }

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

    reactionNode.update(m_rate->parameters());

    // strip information not needed for reconstruction
    string rtype = reactionNode["type"].asString();
    if (rtype == "pressure-dependent-Arrhenius") {
        // skip
    } else if (m_explicit_type && ba::ends_with(rtype, "Arrhenius")) {
        // retain type information
        if (m_third_body) {
            reactionNode["type"] = "three-body";
        } else {
            reactionNode["type"] = "elementary";
        }
    } else if (ba::ends_with(rtype, "Arrhenius")) {
        reactionNode.erase("type");
    } else if (m_explicit_type) {
        reactionNode["type"] = type();
    } else if (ba::ends_with(rtype, "Blowers-Masel")) {
        reactionNode["type"] = "Blowers-Masel";
    }

    if (m_third_body) {
        m_third_body->getParameters(reactionNode);
    }
}

void Reaction::setParameters(const AnyMap& node, const Kinetics& kin)
{
    if (node.empty()) {
        throw InputFileError("Reaction::setParameters", input,
            "Cannot set reaction parameters from empty node.");
    }

    input = node;
    input.copyMetadata(node);
    setEquation(node["equation"].asString(), &kin);
    // Non-stoichiometric reaction orders
    if (node.hasKey("orders")) {
        for (const auto& [name, order] : node["orders"].asMap<double>()) {
            orders[name] = order;
            if (kin.kineticsSpeciesIndex(name) == npos) {
                setValid(false);
            }
        }
    }

    // Flags
    id = node.getString("id", "");
    duplicate = node.getBool("duplicate", false);
    allow_negative_orders = node.getBool("negative-orders", false);
    allow_nonreactant_orders = node.getBool("nonreactant-orders", false);

    if (m_third_body) {
        m_third_body->setParameters(node);
        if (m_third_body->name() == "M" && m_third_body->efficiencies.size() == 1) {
            m_third_body->explicit_3rd = true;
        }
    } else if (node.hasKey("default-efficiency") || node.hasKey("efficiencies")) {
        throw InputFileError("Reaction::setParameters", input,
            "Reaction '{}' specifies efficiency parameters\n"
            "but does not involve third body colliders.", equation());
    }
}

void Reaction::setRate(shared_ptr<ReactionRate> rate)
{
    if (!rate) {
        throw InputFileError("Reaction::setRate", input,
            "Reaction rate for reaction '{}' must not be empty.", equation());
    }
    m_rate = rate;
}

string Reaction::reactantString() const
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
    if (m_third_body) {
        result << m_third_body->collider();
    }
    return result.str();
}

string Reaction::productString() const
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
    if (m_third_body) {
        result << m_third_body->collider();
    }
    return result.str();
}

string Reaction::equation() const
{
    if (reversible) {
        return reactantString() + " <=> " + productString();
    } else {
        return reactantString() + " => " + productString();
    }
}

void Reaction::setEquation(const string& equation, const Kinetics* kin)
{
    parseReactionEquation(*this, equation, input, kin);
    string rate_type = (m_rate) ? m_rate->type() : input.getString("type", "");
    if (ba::starts_with(rate_type, "three-body")) {
        // state type when serializing
        m_explicit_type = true;
    } else if (rate_type == "elementary") {
        // user override
        m_explicit_type = true;
        return;
    } else if (kin && kin->thermo(0).nDim() != 3) {
        // interface reactions
        return;
    } else if (rate_type == "electron-collision-plasma") {
        // does not support third body
        return;
    }

    string third_body;
    size_t count = 0;
    size_t countM = 0;
    for (const auto& [name, stoich] : reactants) {
        // detect explicitly specified collision partner
        if (products.count(name)) {
            third_body = name;
            size_t generic = third_body == "M"
                || third_body == "(+M)"  || third_body == "(+ M)";
            count++;
            countM += generic;
            if (stoich > 1 && products[third_body] > 1) {
                count++;
                countM += generic;
            }
        }
    }

    if (count == 0) {
        if (ba::starts_with(rate_type, "three-body")) {
            throw InputFileError("Reaction::setEquation", input,
                "Reactants for reaction '{}'\n"
                "do not contain a valid third body collider", equation);
        }
        return;

    } else if (countM > 1) {
        throw InputFileError("Reaction::setEquation", input,
            "Multiple generic third body colliders 'M' are not supported", equation);

    } else if (count > 1) {
        // equations with more than one explicit third-body collider are handled as a
        // regular elementary reaction unless the equation contains a generic third body
        if (countM) {
            // generic collider 'M' is selected as third body
            third_body = "M";
        } else if (m_third_body) {
            // third body is defined as explicit object
            auto& effs = m_third_body->efficiencies;
            if (effs.size() != 1 || !reactants.count(effs.begin()->first)) {
                throw InputFileError("Reaction::setEquation", input,
                    "Detected ambiguous third body colliders in reaction '{}'\n"
                    "ThirdBody object needs to specify a single species", equation);
            }
            third_body = effs.begin()->first;
            m_third_body->explicit_3rd = true;
        } else if (input.hasKey("efficiencies")) {
            // third body is implicitly defined by efficiency
            auto effs = input["efficiencies"].asMap<double>();
            if (effs.size() != 1 || !reactants.count(effs.begin()->first)) {
                throw InputFileError("Reaction::setEquation", input,
                    "Detected ambiguous third body colliders in reaction '{}'\n"
                    "Collision efficiencies need to specify single species", equation);
            }
            third_body = effs.begin()->first;
            m_third_body = make_shared<ThirdBody>(third_body);
            m_third_body->explicit_3rd = true;
        } else if (input.hasKey("default-efficiency")) {
            // insufficient disambiguation of third bodies
            throw InputFileError("Reaction::setEquation", input,
                "Detected ambiguous third body colliders in reaction '{}'\n"
                "Third-body definition requires specification of efficiencies",
                equation);
        } else if (ba::starts_with(rate_type, "three-body")) {
            // no disambiguation of third bodies
            throw InputFileError("Reaction::setEquation", input,
                "Detected ambiguous third body colliders in reaction '{}'\n"
                "A valid ThirdBody or collision efficiency definition is required",
                equation);
        } else {
            return;
        }

    } else if (third_body != "M" && !ba::starts_with(rate_type, "three-body")
            && !ba::starts_with(third_body, "(+"))
    {
        // check for conditions of three-body reactions:
        // - integer stoichiometric conditions
        // - either reactant or product side involves exactly three species
        size_t nreac = 0;
        size_t nprod = 0;

        // ensure that all reactants have integer stoichiometric coefficients
        for (const auto& [name, stoich] : reactants) {
            if (trunc(stoich) != stoich) {
                return;
            }
            nreac += static_cast<size_t>(stoich);
        }

        // ensure that all products have integer stoichiometric coefficients
        for (const auto& [name, stoich] : products) {
            if (trunc(stoich) != stoich) {
                return;
            }
            nprod += static_cast<size_t>(stoich);
        }

        // either reactant or product side involves exactly three species
        if (nreac != 3 && nprod != 3) {
            return;
        }
    }

    if (m_third_body) {
        string tName = m_third_body->name();
        if (tName != third_body && third_body != "M" && tName != "M") {
            throw InputFileError("Reaction::setEquation", input,
                "Detected incompatible third body colliders in reaction '{}'\n"
                "ThirdBody definition does not match equation", equation);
        }
        m_third_body->setName(third_body);
    } else {
        m_third_body = make_shared<ThirdBody>(third_body);
    }

    // adjust reactant coefficients
    auto reac = reactants.find(third_body);
    if (trunc(reac->second) != 1) {
        reac->second -= 1.;
    } else {
        reactants.erase(reac);
    }

    // adjust product coefficients
    auto prod = products.find(third_body);
    if (trunc(prod->second) != 1) {
        prod->second -= 1.;
    } else {
        products.erase(prod);
    }
}

string Reaction::type() const
{
    if (!m_rate) {
        throw CanteraError("Reaction::type", "Empty Reaction does not have a type");
    }

    string rate_type = m_rate->type();
    string sub_type = m_rate->subType();
    if (sub_type != "") {
        return rate_type + "-" + sub_type;
    }

    if (m_third_body) {
        return "three-body-" + rate_type;
    }

    return rate_type;
}

UnitStack Reaction::calculateRateCoeffUnits(const Kinetics& kin)
{
    if (!valid()) {
        // If a reaction is invalid because of missing species in the Kinetics
        // object, determining the units of the rate coefficient is impossible.
        return UnitStack({});
    }

    // Determine the units of the rate coefficient
    UnitStack rate_units(kin.thermo(0).standardConcentrationUnits());

    // Set output units to standardConcentrationUnits per second
    rate_units.join(1.);
    rate_units.update(Units(1.0, 0, 0, -1), 1.);

    for (const auto& [name, order] : orders) {
        const auto& phase = kin.speciesPhase(name);
        // Account for specified reaction orders
        rate_units.update(phase.standardConcentrationUnits(), -order);
    }
    for (const auto& [name, stoich] : reactants) {
        // Order for each reactant is the reactant stoichiometric coefficient,
        // unless already overridden by user-specified orders
        if (name == "M" || ba::starts_with(name, "(+")) {
            // calculateRateCoeffUnits may be called before these pseudo-species
            // have been stripped from the reactants
            continue;
        } else if (orders.find(name) == orders.end()) {
            const auto& phase = kin.speciesPhase(name);
            // Account for each reactant species
            rate_units.update(phase.standardConcentrationUnits(), -stoich);
        }
    }

    if (m_third_body && m_third_body->mass_action) {
        // Account for third-body collision partner as the last entry
        rate_units.join(-1);
    }

    Reaction::rate_units = rate_units.product();
    return rate_units;
}

void updateUndeclared(vector<string>& undeclared,
                      const Composition& comp, const Kinetics& kin)
{
    for (const auto& [name, stoich]: comp) {
        if (kin.kineticsSpeciesIndex(name) == npos) {
            undeclared.emplace_back(name);
        }
    }
}

void Reaction::checkBalance(const Kinetics& kin) const
{
    Composition balr, balp;

    // iterate over products and reactants
    for (const auto& [name, stoich] : products) {
        const ThermoPhase& ph = kin.speciesPhase(name);
        size_t k = ph.speciesIndex(name);
        for (size_t m = 0; m < ph.nElements(); m++) {
            balr[ph.elementName(m)] = 0.0; // so that balr contains all species
            balp[ph.elementName(m)] += stoich * ph.nAtoms(k, m);
        }
    }
    for (const auto& [name, stoich] : reactants) {
        const ThermoPhase& ph = kin.speciesPhase(name);
        size_t k = ph.speciesIndex(name);
        for (size_t m = 0; m < ph.nElements(); m++) {
            balr[ph.elementName(m)] += stoich * ph.nAtoms(k, m);
        }
    }

    string msg;
    bool ok = true;
    for (const auto& [elem, balance] : balr) {
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

    if (kin.thermo(0).nDim() == 3) {
        return;
    }

    // Check that the number of surface sites is balanced
    double reac_sites = 0.0;
    double prod_sites = 0.0;
    auto& surf = dynamic_cast<const SurfPhase&>(kin.thermo(0));
    for (const auto& [name, stoich] : reactants) {
        size_t k = surf.speciesIndex(name);
        if (k != npos) {
            reac_sites += stoich * surf.size(k);
        }
    }
    for (const auto& [name, stoich] : products) {
        size_t k = surf.speciesIndex(name);
        if (k != npos) {
            prod_sites += stoich * surf.size(k);
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
    vector<string> undeclared;
    updateUndeclared(undeclared, reactants, kin);
    updateUndeclared(undeclared, products, kin);
    if (!undeclared.empty()) {
        if (kin.skipUndeclaredSpecies()) {
            return false;
        } else {
            throw InputFileError("Reaction::checkSpecies", input, "Reaction '{}'\n"
                "contains undeclared species: '{}'",
                equation(), ba::join(undeclared, "', '"));
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
                    equation(), ba::join(undeclared, "', '"));
            }
            // Error for empty input AnyMap (that is, XML)
            throw InputFileError("Reaction::checkSpecies", input, "Reaction '{}'\n"
                "defines reaction orders for undeclared species: '{}'",
                equation(), ba::join(undeclared, "', '"));
        }
    }

    if (m_third_body) {
        return m_third_body->checkSpecies(*this, kin);
    }

    checkBalance(kin);

    return true;
}

bool Reaction::usesElectrochemistry(const Kinetics& kin) const
{
    // Check electrochemistry
    vector<double> e_counter(kin.nPhases(), 0.0);

    // Find the number of electrons in the products for each phase
    for (const auto& [name, stoich] : products) {
        size_t kkin = kin.kineticsSpeciesIndex(name);
        size_t i = kin.speciesPhaseIndex(kkin);
        size_t kphase = kin.thermo(i).speciesIndex(name);
        e_counter[i] += stoich * kin.thermo(i).charge(kphase);
    }

    // Subtract the number of electrons in the reactants for each phase
    for (const auto& [name, stoich] : reactants) {
        size_t kkin = kin.kineticsSpeciesIndex(name);
        size_t i = kin.speciesPhaseIndex(kkin);
        size_t kphase = kin.thermo(i).speciesIndex(name);
        e_counter[i] -= stoich * kin.thermo(i).charge(kphase);
    }

    // If the electrons change phases then the reaction is electrochemical
    for (double delta_e : e_counter) {
        if (std::abs(delta_e) > 1e-4) {
            return true;
        }
    }

    return false;
}


ThirdBody::ThirdBody(const string& third_body)
{
    setName(third_body);
}

void ThirdBody::setName(const string& third_body)
{
    string name = third_body;
    if (ba::starts_with(third_body, "(+ ")) {
        mass_action = false;
        name = third_body.substr(3, third_body.size() - 4);
    } else if (ba::starts_with(third_body, "(+")) {
        mass_action = false;
        name = third_body.substr(2, third_body.size() - 3);
    }

    if (name == m_name) {
        return;
    }
    if (name == "M" && efficiencies.size() == 1) {
        // revert from explicit name to generic collider
        m_name = name;
        return;
    }
    if (efficiencies.size()) {
        throw CanteraError("ThirdBody::setName",
            "Conflicting efficiency definition for explicit third body '{}'", name);
    }
    m_name = name;
    default_efficiency = 0.;
    efficiencies[m_name] = 1.;
}

ThirdBody::ThirdBody(const AnyMap& node)
{
    setParameters(node);
}

void ThirdBody::setParameters(const AnyMap& node)
{
    if (node.hasKey("default-efficiency")) {
        double value = node["default-efficiency"].asDouble();
        if (m_name != "M" && value != 0.) {
            throw InputFileError("ThirdBody::setParameters", node["default-efficiency"],
                "Invalid default efficiency for explicit collider {};\n"
                "value is optional and/or needs to be zero", m_name);
        }
        default_efficiency = value;
    }
    if (node.hasKey("efficiencies")) {
        efficiencies = node["efficiencies"].asMap<double>();
    }
    if (m_name != "M"
        && (efficiencies.size() != 1 || efficiencies.begin()->first != m_name))
    {
        throw InputFileError("ThirdBody::setParameters", node,
            "Detected incompatible third body colliders definitions");
    }
}

void ThirdBody::getParameters(AnyMap& node) const
{
    if (m_name == "M" || explicit_3rd) {
        if (efficiencies.size()) {
            node["efficiencies"] = efficiencies;
            node["efficiencies"].setFlowStyle();
        }
        if (default_efficiency != 1.0 && !explicit_3rd) {
            node["default-efficiency"] = default_efficiency;
        }
    }
}

double ThirdBody::efficiency(const string& k) const
{
    return getValue(efficiencies, k, default_efficiency);
}

string ThirdBody::collider() const
{
    if (mass_action) {
        return " + " + m_name;
    }
    return " (+" + m_name + ")";
}

bool ThirdBody::checkSpecies(const Reaction& rxn, const Kinetics& kin) const
{
    vector<string> undeclared;
    updateUndeclared(undeclared, efficiencies, kin);

    if (!undeclared.empty()) {
        if (!kin.skipUndeclaredThirdBodies()) {
            if (rxn.input.hasKey("efficiencies")) {
                throw InputFileError("ThirdBody::checkSpecies",
                    rxn.input["efficiencies"], "Reaction '{}'\n"
                    "defines third-body efficiencies for undeclared species: '{}'",
                    rxn.equation(), ba::join(undeclared, "', '"));
            }
            // Error for specified ThirdBody or empty input AnyMap
            throw InputFileError("ThirdBody::checkSpecies", rxn.input, "Reaction '{}'\n"
                "is a three-body reaction with undeclared species: '{}'",
                rxn.equation(), ba::join(undeclared, "', '"));
        } else if (kin.skipUndeclaredThirdBodies() && m_name != "M") {
            // Prevent addition of reaction silently as "skip-undeclared-third-bodies"
            // is set to true
            return false;
        }
    }
    return true;
}


unique_ptr<Reaction> newReaction(const string& type)
{
    return make_unique<Reaction>();
}

unique_ptr<Reaction> newReaction(const AnyMap& rxn_node, const Kinetics& kin)
{
    return make_unique<Reaction>(rxn_node, kin);
}

void parseReactionEquation(Reaction& R, const string& equation,
                           const AnyBase& reactionNode, const Kinetics* kin)
{
    // Parse the reaction equation to determine participating species and
    // stoichiometric coefficients
    vector<string> tokens;
    tokenizeString(equation, tokens);
    tokens.push_back("+"); // makes parsing last species not a special case

    size_t last_used = npos; // index of last-used token
    bool reactants = true;
    for (size_t i = 1; i < tokens.size(); i++) {
        if (tokens[i] == "+" || ba::starts_with(tokens[i], "(+") ||
            tokens[i] == "<=>" || tokens[i] == "=" || tokens[i] == "=>") {
            string species = tokens[i-1];

            double stoich = 1.0;
            bool mass_action = true;
            if (last_used != npos && tokens[last_used] == "(+"
                    && ba::ends_with(species, ")")) {
                // Falloff third body with space, such as "(+ M)"
                mass_action = false;
                species = "(+" + species;
            } else if (last_used == i - 1 && ba::starts_with(species, "(+")
                    && ba::ends_with(species, ")")) {
                // Falloff 3rd body written without space, such as "(+M)"
                mass_action = false;
            } else if (last_used == i - 2) {
                // Species with no stoich. coefficient
            } else if (last_used == i - 3) {
                // Stoich. coefficient and species
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
                         && mass_action && species != "M"))
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

vector<shared_ptr<Reaction>> getReactions(const AnyValue& items, Kinetics& kinetics)
{
    vector<shared_ptr<Reaction>> all_reactions;
    for (const auto& node : items.asVector<AnyMap>()) {
        auto R = make_shared<Reaction>(node, kinetics);
        R->validate(kinetics);
        if (R->valid() && R->checkSpecies(kinetics)) {
            all_reactions.emplace_back(R);
        }
    }
    return all_reactions;
}

}
