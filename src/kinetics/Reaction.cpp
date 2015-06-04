/**
 *  @file Reaction.cpp
 */

#include "cantera/kinetics/Reaction.h"
#include "cantera/kinetics/FalloffFactory.h"
#include "cantera/base/ctml.h"
#include "cantera/base/Array.h"
#include <sstream>

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
        for (Composition::iterator iter = orders.begin();
             iter != orders.end();
             ++iter) {
            if (reactants.find(iter->first) == reactants.end()) {
                throw CanteraError("Reaction::validate", "Reaction order "
                    "specified for non-reactant species '" + iter->first + "'");
            }
        }
    }

    if (!allow_negative_orders) {
        for (Composition::iterator iter = orders.begin();
             iter != orders.end();
             ++iter) {
            if (iter->second < 0.0) {
                throw CanteraError("Reaction::validate", "Negative reaction "
                    "order specified for species '" + iter->first + "'");
            }
        }
    }
}

std::string Reaction::reactantString() const
{
    std::ostringstream result;
    for (Composition::const_iterator iter = reactants.begin();
         iter != reactants.end();
         ++iter) {
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
    for (Composition::const_iterator iter = products.begin();
         iter != products.end();
         ++iter) {
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
    if (low_rate.preExponentialFactor() < 0 ||
        high_rate.preExponentialFactor() < 0) {
        throw CanteraError("FalloffReaction::validate", "Negative "
            "pre-exponential factor found for reaction '" + equation() + "'");
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
    : FalloffReaction(reactants_, products_, low_rate, high_rate, tbody)
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
{
    reaction_type = INTERFACE_RXN;
}

InterfaceReaction::InterfaceReaction(const Composition& reactants_,
                                     const Composition& products_,
                                     const Arrhenius& rate_,
                                     bool isStick)
    : ElementaryReaction(reactants_, products_, rate_)
    , is_sticking_coefficient(isStick)
{
    reaction_type = INTERFACE_RXN;
}

ElectrochemicalReaction::ElectrochemicalReaction()
    : film_resistivity(0.0)
    , beta(0.0)
    , exchange_current_density_formulation(false)
{
}

ElectrochemicalReaction::ElectrochemicalReaction(const Composition& reactants_,
                                                 const Composition& products_,
                                                 const Arrhenius& rate_)
    : InterfaceReaction(reactants_, products_, rate_)
    , film_resistivity(0.0)
    , beta(0.0)
    , exchange_current_density_formulation(false)
{
}


Arrhenius readArrhenius(const XML_Node& arrhenius_node)
{
    return Arrhenius(getFloat(arrhenius_node, "A", "toSI"),
                     getFloat(arrhenius_node, "b"),
                     getFloat(arrhenius_node, "E", "actEnergy") / GasConstant);
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

    int falloff_type = 0;
    if (lowercase(falloff["type"]) == "lindemann") {
        falloff_type = SIMPLE_FALLOFF;
        if (np != 0) {
            throw CanteraError("readFalloff", "Lindemann parameterization "
                "takes no parameters, but " + int2str(np) + "were given");
        }
    } else if (lowercase(falloff["type"]) == "troe") {
        falloff_type = TROE_FALLOFF;
        if (np != 3 && np != 4) {
            throw CanteraError("readFalloff", "Troe parameterization takes "
                "3 or 4 parameters, but " + int2str(np) + "were given");
        }
    } else if (lowercase(falloff["type"]) == "sri") {
        falloff_type = SRI_FALLOFF;
        if (np != 3 && np != 5) {
            throw CanteraError("readFalloff", "SRI parameterization takes "
                "3 or 5 parameters, but " + int2str(np) + "were given");
        }
    } else {
        throw CanteraError("readFalloff", "Unrecognized falloff type: '" +
            falloff["type"] + "'");
    }
    R.falloff = newFalloff(falloff_type, falloff_parameters);
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
    setupReaction(R, rxn_node);
}

void setupThreeBodyReaction(ThreeBodyReaction& R, const XML_Node& rxn_node)
{
    readEfficiencies(R.third_body, rxn_node.child("rateCoeff"));
    setupElementaryReaction(R, rxn_node);
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
    readFalloff(R, rc_node);
    readEfficiencies(R.third_body, rc_node);
    setupReaction(R, rxn_node);
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
        rates.insert(std::make_pair(getFloat(node, "P", "toSI"),
                                    readArrhenius(node)));
    }
    R.rate = Plog(rates);
    setupReaction(R, rxn_node);
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

void setupInterfaceReaction(InterfaceReaction& R, const XML_Node& rxn_node)
{
    if (lowercase(rxn_node["type"]) == "global") {
        R.reaction_type = GLOBAL_RXN;
    }
    XML_Node& arr = rxn_node.child("rateCoeff").child("Arrhenius");
    if (lowercase(arr["type"]) == "stick") {
        R.is_sticking_coefficient = true;
        R.sticking_species = arr["species"];
    }
    std::vector<XML_Node*> cov = arr.getChildren("coverage");
    for (std::vector<XML_Node*>::iterator iter = cov.begin();
         iter != cov.end();
         ++iter)
    {
        CoverageDependency& cdep = R.coverage_deps[(*iter)->attrib("species")];
        cdep.a = getFloat(**iter, "a", "toSI");
        cdep.m = getFloat(**iter, "m");
        cdep.E = getFloat(**iter, "e", "actEnergy") / GasConstant;
    }
    setupElementaryReaction(R, rxn_node);
}

void setupElectrochemicalReaction(ElectrochemicalReaction& R,
                                  const XML_Node& rxn_node)
{
    // Fix reaction_type for some specialized reaction types
    std::string type = lowercase(rxn_node["type"]);
    if (type == "butlervolmer") {
        R.reaction_type = BUTLERVOLMER_RXN;
    } else if (type == "butlervolmer_noactivitycoeffs") {
        R.reaction_type = BUTLERVOLMER_NOACTIVITYCOEFFS_RXN;
    } else if (type == "surfaceaffinity") {
        R.reaction_type = SURFACEAFFINITY_RXN;
    } else if (type == "global") {
        R.reaction_type = GLOBAL_RXN;
    }

    XML_Node& rc = rxn_node.child("rateCoeff");
    std::string rc_type = lowercase(rc["type"]);
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
        for (Composition::const_iterator iter = R.reactants.begin();
             iter != R.reactants.end();
             ++iter) {
            R.orders[iter->first] += iter->second * (1.0 - R.beta);
        }
        for (Composition::const_iterator iter = R.products.begin();
             iter != R.products.end();
             ++iter) {
            R.orders[iter->first] += iter->second * R.beta;
        }
    }

    // For affinity reactions, fill in the global reaction formulation terms
    if (rxn_node.hasChild("reactionOrderFormulation")) {
        Composition initial_orders = R.orders;
        R.orders.clear();
        R.allow_nonreactant_orders = true;
        const XML_Node& rof_node = rxn_node.child("reactionOrderFormulation");
        if (lowercase(rof_node["model"]) == "reactantorders") {
            R.orders = initial_orders;
        } else if (lowercase(rof_node["model"]) == "zeroorders") {
            for (Composition::const_iterator iter = R.reactants.begin();
                 iter != R.reactants.end();
                 ++iter) {
                R.orders[iter->first] = 0.0;
            }
        } else if (lowercase(rof_node["model"]) == "butlervolmerorders") {
            // Reaction orders based on provided reaction orders
            for (Composition::const_iterator iter = R.reactants.begin();
                 iter != R.reactants.end();
                 ++iter) {
                double c = getValue(initial_orders, iter->first, iter->second);
                R.orders[iter->first] += c * (1.0 - R.beta);
            }
            for (Composition::const_iterator iter = R.products.begin();
                 iter != R.products.end();
                 ++iter) {
                double c = getValue(initial_orders, iter->first, iter->second);
                R.orders[iter->first] += c * R.beta;
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
        for (Composition::iterator iter = orders.begin();
             iter != orders.end();
             ++iter) {
            R.orders[iter->first] = iter->second;
        }
    }
}

shared_ptr<Reaction> newReaction(const XML_Node& rxn_node)
{
    std::string type = lowercase(rxn_node["type"]);

    // Modify the reaction type for edge reactions which contain electrochemical
    // reaction data
    if (rxn_node.child("rateCoeff").hasChild("electrochem") && type == "edge") {
        type = "electrochemical";
    }

    // Create a new Reaction object of the appropriate type
    if (type == "elementary" || type == "arrhenius" || type == "") {
        shared_ptr<ElementaryReaction> R(new ElementaryReaction());
        setupElementaryReaction(*R, rxn_node);
        return R;

    } else if (type == "threebody" || type == "three_body") {
        shared_ptr<ThreeBodyReaction> R(new ThreeBodyReaction());
        setupThreeBodyReaction(*R, rxn_node);
        return R;

    } else if (type == "falloff") {
        shared_ptr<FalloffReaction> R(new FalloffReaction());
        setupFalloffReaction(*R, rxn_node);
        return R;

    } else if (type == "chemact" || type == "chemically_activated") {
        shared_ptr<ChemicallyActivatedReaction> R(new ChemicallyActivatedReaction());
        setupChemicallyActivatedReaction(*R, rxn_node);
        return R;

    } else if (type == "plog" || type == "pdep_arrhenius") {
        shared_ptr<PlogReaction> R(new PlogReaction());
        setupPlogReaction(*R, rxn_node);
        return R;

    } else if (type == "chebyshev") {
        shared_ptr<ChebyshevReaction> R(new ChebyshevReaction());
        setupChebyshevReaction(*R, rxn_node);
        return R;

    } else if (type == "interface" || type == "surface" || type == "edge" ||
               type == "global") {
        shared_ptr<InterfaceReaction> R(new InterfaceReaction());
        setupInterfaceReaction(*R, rxn_node);
        return R;

    } else if (type == "electrochemical" ||
               type == "butlervolmer_noactivitycoeffs" ||
               type == "butlervolmer" ||
               type == "surfaceaffinity") {
        shared_ptr<ElectrochemicalReaction> R(new ElectrochemicalReaction());
        setupElectrochemicalReaction(*R, rxn_node);
        return R;

    } else {
        throw CanteraError("newReaction",
            "Unknown reaction type '" + rxn_node["type"] + "'");
    }
}

std::vector<shared_ptr<Reaction> > getReactions(const XML_Node& node)
{
    std::vector<shared_ptr<Reaction> > all_reactions;
    std::vector<XML_Node*> reaction_nodes =
        node.child("reactionData").getChildren("reaction");

    for (std::vector<XML_Node*>::iterator iter = reaction_nodes.begin();
         iter != reaction_nodes.end();
         ++iter)
    {
        all_reactions.push_back(newReaction(**iter));
    }
    return all_reactions;
}

}
