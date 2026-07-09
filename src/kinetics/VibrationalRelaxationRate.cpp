//! @file VibrationalRelaxationRate.cpp
//! @since New in %Cantera 4.0

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/VibrationalRelaxationRate.h"
#include "cantera/kinetics/Reaction.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/base/global.h"

#include <cmath>
#include <map>
#include <set>
#include <sstream>
#include <utility>
#include <vector>

namespace Cantera
{

namespace
{

constexpr double Tiny = 1e-12;

const string WhereSetParameters = "VibrationalRelaxationRate::setParameters";
const string WhereGetParameters = "VibrationalRelaxationRate::getParameters";
const string WhereSetContext = "VibrationalRelaxationRate::setContext";

const string ModelConstant = "constant";
const string ModelMultiState = "multi-state-resolved";
const string ModelStarikovskiy = "starikovskiy";
const string ModelCastela = "castela";


// helpers to check the correctness of a chosen vibrational model input data.
void requireNoKey(const AnyMap& node, const string& key,
                  const string& model, const string& where)
{
    if (node.hasKey(key)) {
        throw InputFileError(where, node,
            "Key '{}' is not allowed for vibration-model '{}'.", key, model);
    }
}

void requireKey(const AnyMap& node, const string& key,
                const string& model, const string& where)
{
    if (!node.hasKey(key)) {
        throw InputFileError(where, node,
            "Missing required key '{}' for vibration-model '{}'.", key, model);
    }
}

void requireKeys(const AnyMap& node, const string& model,
                 const string& where, std::initializer_list<string> keys)
{
    for (const auto& key : keys) {
        requireKey(node, key, model, where);
    }
}

void forbidKeys(const AnyMap& node, const string& model,
                const string& where, std::initializer_list<string> keys)
{
    for (const auto& key : keys) {
        requireNoKey(node, key, model, where);
    }
}

// check whether a species is a vibrational excited species 
// or a vibrational reservoir from its name
bool isVibrationalSpecies(const string& name)
{
    // Supported names formats:
    //
    //   N2(v)
    //   O2(v0)
    //   O2(v1)
    //   O2(v12)
    //
    // Style like O2(v=1) should not be used.
    const auto pos = name.find("(v");
    return pos != string::npos && !name.empty() && name.back() == ')';
}

// find the name of the ground state of a vibrationally excited species.
string groundStateName(const string& name)
{
    const auto pos = name.find("(v");
    if (pos == string::npos) {
        return name;
    }
    return name.substr(0, pos);
}

string vibrationalFamilyName(const string& name)
{
    // Collapse O2(v1), O2(v2), O2(v12), and O2(v) into O2(v) as a "family
    // name" to identify which species have the same ground state "parent". 
    return groundStateName(name) + "(v)";
}

double compositionSum(const Composition& comp)
{
    double sum = 0.0;
    for (const auto& item : comp) {
        sum += item.second;
    }
    return sum;
}

Composition replaceVibrationalSpeciesByGroundState(const Composition& comp)
{
    Composition out;
    for (const auto& item : comp) {
        const string& name = item.first;
        const double value = item.second;

        if (isVibrationalSpecies(name)) {
            out[groundStateName(name)] += value;
        } else {
            out[name] += value;
        }
    }
    return out;
}

bool sameComposition(const Composition& a, const Composition& b,
                     double tol = 1e-12)
{
    std::set<string> names;

    for (const auto& item : a) {
        names.insert(item.first);
    }
    for (const auto& item : b) {
        names.insert(item.first);
    }

    for (const auto& name : names) {
        double av = 0.0;
        double bv = 0.0;

        const auto ait = a.find(name);
        if (ait != a.end()) {
            av = ait->second;
        }

        const auto bit = b.find(name);
        if (bit != b.end()) {
            bv = bit->second;
        }

        if (std::abs(av - bv) > tol) {
            return false;
        }
    }

    return true;
}

std::vector<string> vibrationalSpeciesInComposition(const Composition& comp)
{
    std::vector<string> out;

    for (const auto& item : comp) {
        const string& name = item.first;
        const double value = item.second;

        if (value != 0.0 && isVibrationalSpecies(name)) {
            out.push_back(name);
        }
    }

    return out;
}

// Registry used to check that a given vibrational family uses exactly one
// relaxation model inside one Kinetics object.
//
// It is allowed for each "vibrational family" 
// to have its own vibrational relaxation model:
//
//   N2(v)  -> constant
//   O2(v)  -> multi-state-resolved
//   NH3(v) -> starikovskiy
//
// But it is forbidden for each reaction of the same "vibrational family"
// to have different relaxation models:
//
//   N2(v) + O  -> castela
//   N2(v) + N2 -> starikovskiy
//
std::map<const Kinetics*, std::map<string, string>> s_modelByFamily;
std::set<const Kinetics*> s_warnedMixedModels;

void registerVibrationalModelConsistency(const Kinetics& kin,
                                          const string& family,
                                          const string& model,
                                          const AnyMap& input)
{
    auto& modelByFamily = s_modelByFamily[&kin];

    const auto existing = modelByFamily.find(family);
    if (existing != modelByFamily.end() && existing->second != model) {
        throw InputFileError(WhereSetContext, input,
            "Inconsistent vibration-model for vibrational family '{}'. "
            "This family was already registered with model '{}', but the "
            "current reaction uses model '{}'. A given vibrational family "
            "must use exactly one relaxation model.",
            family, existing->second, model);
    }

    modelByFamily[family] = model;

    std::set<string> models;
    for (const auto& item : modelByFamily) {
        models.insert(item.second);
    }

    if (models.size() > 1 && !s_warnedMixedModels.count(&kin)) {
        std::ostringstream msg;
        msg << "Multiple vibrational relaxation models were detected in the "
            << "same kinetics object. This is allowed only if each "
            << "vibrational family is internally consistent:\n";

        for (const auto& item : modelByFamily) {
            msg << "  - " << item.first << ": " << item.second << "\n";
        }

        warn_user(WhereSetContext, msg.str());
        s_warnedMixedModels.insert(&kin);
    }
}

string inferRelaxingFamily(const Reaction& rxn)
{
    const auto vibReactants = vibrationalSpeciesInComposition(rxn.reactants);

    if (vibReactants.empty()) {
        throw InputFileError(WhereSetContext, rxn.input,
            "A vibrational-relaxation reaction must contain at least one "
            "vibrational reactant, for example O2(v1), N2(v), or NH3(v1).");
    }

    return vibrationalFamilyName(vibReactants.front());
}

void validateSimpleRelaxationToGroundState(const Reaction& rxn,
                                           const string& model)
{
    const auto vibReactants = vibrationalSpeciesInComposition(rxn.reactants);
    const auto vibProducts = vibrationalSpeciesInComposition(rxn.products);

    if (vibReactants.size() != 1) {
        throw InputFileError(WhereSetContext, rxn.input,
            "vibration-model '{}' expects exactly one vibrational reactant.",
            model);
    }

    if (!vibProducts.empty()) {
        throw InputFileError(WhereSetContext, rxn.input,
            "vibration-model '{}' describes relaxation to the ground state "
            "and therefore does not allow vibrational products.", model);
    }

    const Composition relaxedReactants =
        replaceVibrationalSpeciesByGroundState(rxn.reactants);

    if (!sameComposition(relaxedReactants, rxn.products)) {
        throw InputFileError(WhereSetContext, rxn.input,
            "Invalid vibrational relaxation stoichiometry for model '{}'. "
            "Expected a reaction equivalent to X(v) + M => X + M, where the "
            "collider M is unchanged.", model);
    }

    if (std::abs(compositionSum(rxn.reactants) - 2.0) > 1e-12
        || std::abs(compositionSum(rxn.products) - 2.0) > 1e-12)
    {
        throw InputFileError(WhereSetContext, rxn.input,
            "vibration-model '{}' expects a bimolecular relaxation reaction "
            "of the form X(v) + M => X + M.", model);
    }
}

void validateCastelaReaction(const Reaction& rxn)
{
    validateSimpleRelaxationToGroundState(rxn, ModelCastela);

    const auto vibReactants = vibrationalSpeciesInComposition(rxn.reactants);
    const string family = vibrationalFamilyName(vibReactants.front());

    if (family != "N2(v)") {
        throw InputFileError(WhereSetContext, rxn.input,
            "The Castela relaxation model is only valid for N2 vibrational "
            "relaxation. Found vibrational family '{}'.", family);
    }

    // Castela parameters used here are only intended for the following
    // colliders: N2, O2, and O.
    const Composition collapsedReactants =
        replaceVibrationalSpeciesByGroundState(rxn.reactants);

    string collider = "";
    for (const auto& item : collapsedReactants) {
        const string& name = item.first;
        const double value = item.second;

        if (name == "N2") {
            if (std::abs(value - 2.0) < 1e-12) {
                collider = "N2";
            }
        } else if (std::abs(value - 1.0) < 1e-12) {
            collider = name;
        }
    }

    if (collider != "N2" && collider != "O2" && collider != "O") {
        throw InputFileError(WhereSetContext, rxn.input,
            "The Castela relaxation model only supports colliders 'N2', "
            "'O2', and 'O'. Found collider '{}'.", collider);
    }
}

void validateDetailedRelaxationReaction(const Reaction& rxn)
{
    const auto vibReactants = vibrationalSpeciesInComposition(rxn.reactants);
    const auto vibProducts = vibrationalSpeciesInComposition(rxn.products);

    if (vibReactants.empty()) {
        throw InputFileError(WhereSetContext, rxn.input,
            "A detailed vibrational relaxation reaction must contain at least "
            "one vibrational reactant.");
    }

    // All vibrational species in a detailed relaxation reaction are required
    // to belong to the same "vibrational family".
    const string family = vibrationalFamilyName(vibReactants.front());

    for (const auto& sp : vibReactants) {
        const string spFamily = vibrationalFamilyName(sp);
        if (spFamily != family) {
            throw InputFileError(WhereSetContext, rxn.input,
                "Invalid detailed vibrational relaxation reaction: all "
                "vibrational reactants must belong to the same vibrational "
                "family. Found '{}' and '{}'.", family, spFamily);
        }
    }

    for (const auto& sp : vibProducts) {
        const string spFamily = vibrationalFamilyName(sp);
        if (spFamily != family) {
            throw InputFileError(WhereSetContext, rxn.input,
                "Invalid detailed vibrational relaxation reaction: all "
                "vibrational products must belong to the same vibrational "
                "family. Found '{}' and '{}'.", family, spFamily);
        }
    }

    // The reaction must conserve the ground-state composition once all
    // vibrational labels are removed.
    const Composition collapsedReactants =
        replaceVibrationalSpeciesByGroundState(rxn.reactants);

    const Composition collapsedProducts =
        replaceVibrationalSpeciesByGroundState(rxn.products);

    if (!sameComposition(collapsedReactants, collapsedProducts)) {
        throw InputFileError(WhereSetContext, rxn.input,
            "Invalid detailed vibrational relaxation reaction: replacing all "
            "vibrational species by their ground-state species does not "
            "conserve stoichiometry.");
    }

    // The current detailed VV/VT formulation is bimolecular.
    if (std::abs(compositionSum(rxn.reactants) - 2.0) > 1e-12
        || std::abs(compositionSum(rxn.products) - 2.0) > 1e-12)
    {
        throw InputFileError(WhereSetContext, rxn.input,
            "Invalid detailed vibrational relaxation reaction: expected a "
            "bimolecular reaction with two reactant molecules and two product "
            "molecules.");
    }
}

} // end of namespace where all the helpers and safety functions are defined.

bool DetailedVibData::update(const ThermoPhase& phase, const Kinetics& kin)
{
    const double T = phase.temperature();

    if (T == temperature) {
        return false;
    }

    ReactionData::update(T);
    return true;
}

// Default constructor
VibrationalRelaxationRate::VibrationalRelaxationRate()
{
}

// Constructor
VibrationalRelaxationRate::VibrationalRelaxationRate(
    double A, double B, double C, double D,
    double b, double scaling, double m, double E, double z)
    : ArrheniusBase(A, b, 0.0)
    , m_B(B)
    , m_C(C)
    , m_D(D)
    , m_scaling(scaling)
    , m_m(m)
    , m_E(E)
    , m_z(z)
{
}

// AnyMap constructor
VibrationalRelaxationRate::VibrationalRelaxationRate(
    const AnyMap& node, const UnitStack& rate_units)
    : VibrationalRelaxationRate()
{
    setParameters(node, rate_units);
}

void VibrationalRelaxationRate::configureBaseFromInternalA(
    const AnyMap& node, const UnitStack& rate_units, double A, double b)
{
    // Store the original input and configure reaction-rate units.
    //
    // We intentionally do not call ArrheniusBase::setParameters here because
    // some vibration models do not expose a standard YAML rate-constant with
    // both A and b. Castela is the main example.
    // This is also done for future class extension compatibility.
    ReactionRate::setParameters(node, rate_units);
    setRateUnits(rate_units);

    m_negativeA_ok = node.getBool("negative-A", false);

    m_A = A;
    m_b = b;

    if (m_A != 0.0) {
        m_logA = std::log(std::abs(m_A));
    } else {
        m_logA = NAN;
    }

    // ArrheniusBase::validate checks this flag.
    m_valid = true;
}

void VibrationalRelaxationRate::configureBaseFromYamlA(
    const AnyMap& node, const UnitStack& rate_units,
    const AnyValue& A, double b)
{
    // Store the original input and configure reaction-rate units first, so
    // conversionUnits() is available for A.
    ReactionRate::setParameters(node, rate_units);
    setRateUnits(rate_units);

    m_negativeA_ok = node.getBool("negative-A", false);

    // Convert the user-facing YAML A value with Cantera's standard
    // rate-coefficient unit conversion.
    m_A = node.units().convertRateCoeff(A, conversionUnits());
    m_b = b;

    if (m_A != 0.0) {
        m_logA = std::log(std::abs(m_A));
    } else {
        m_logA = NAN;
    }

    // ArrheniusBase::validate checks this flag.
    m_valid = true;
}

void VibrationalRelaxationRate::setParameters(const AnyMap& node,
                                              const UnitStack& rate_units)
{
    // The reaction rate type is always:
    //
    //   type: vibrational-relaxation
    //
    // The physical model is selected separately by:
    //
    //   vibration-model: constant
    //   vibration-model: multi-state-resolved
    //   vibration-model: starikovskiy
    //   vibration-model: castela
    //
    // This is done in a spriti of class adaptability:
    // Any user wishing to develop more advanced vibration models
    // is welcome to add them in this class.
    //
    // The default model is multi-state-resolved.

    m_vibration_model = node.getString("vibration-model", ModelMultiState);

    if (!node.hasKey("rate-constant")) {
        throw InputFileError(WhereSetParameters, node,
            "A vibrational-relaxation reaction requires a 'rate-constant' "
            "mapping.");
    }

    const auto& rate = node["rate-constant"];

    if (!rate.is<AnyMap>()) {
        throw InputFileError(WhereSetParameters, node,
            "The 'rate-constant' field must be a mapping.");
    }

    const auto& rateMap = rate.as<AnyMap>();

    if (m_vibration_model == ModelConstant) {
        // Constant model:
        //
        //   k(T) = A
        requireKeys(rateMap, m_vibration_model, WhereSetParameters, {m_A_str});

        forbidKeys(rateMap, m_vibration_model, WhereSetParameters,
            {m_b_str, "n", m_B_str, m_C_str, m_D_str, m_m_str,
            m_E_str, m_z_str, m_scaling_str});

        configureBaseFromYamlA(node, rate_units, rateMap[m_A_str], 0.0);

        m_B = 0.0;
        m_C = 0.0;
        m_D = 0.0;
        m_m = 2.0 / 3.0;
        m_E = 0.0;
        m_z = 1.0;
        m_scaling = 1.0;
    }
    else if (m_vibration_model == ModelMultiState) {
        // Detailed VV/VT model:
        //
        //   k(T) = scaling * A * exp(
        //       b * log(T)
        //       + B
        //       + C * T^(-1/3)
        //       + D * T^(-2/3)
        //   )
        requireKeys(rateMap, m_vibration_model, WhereSetParameters, {m_A_str});
        forbidKeys(rateMap, m_vibration_model, WhereSetParameters, {"n", m_m_str, m_E_str, m_z_str});

        ArrheniusBase::setParameters(node, rate_units);

        m_B = rateMap.getDouble(m_B_str, 0.0);
        m_C = rateMap.getDouble(m_C_str, 0.0);
        m_D = rateMap.getDouble(m_D_str, 0.0);
        m_m = 2.0 / 3.0;
        m_E = 0.0;
        m_z = 1.0;
        m_scaling = rateMap.getDouble(m_scaling_str, 1.0);
    }
    else if (m_vibration_model == ModelStarikovskiy) {
        // User-facing formula:
        //
        //   k(T) = A * T^n * exp(
        //       K
        //       + B * T^(-1/3)
        //       + C * T^(-m)
        //       + D * T^(-z)
        //   )
        //
        // Internal mapping:
        //
        //   m_b = n
        //   m_B = K
        //   m_C = B
        //   m_D = C
        //   m_m = m
        //   m_E = D
        //   m_z = z

        requireKeys(rateMap, m_vibration_model, WhereSetParameters, {m_A_str});
        forbidKeys(rateMap, m_vibration_model, WhereSetParameters, {m_b_str, m_scaling_str});

        const double n = rateMap.hasKey("n") ? rateMap["n"].asDouble() : 0.0;
        configureBaseFromYamlA(node, rate_units, rateMap[m_A_str], n);

        m_B = rateMap.getDouble("K", 0.0);
        m_C = rateMap.getDouble("B", 0.0);
        m_D = rateMap.getDouble("C", 0.0);
        m_m = rateMap.getDouble("m", 1.0);
        m_E = rateMap.getDouble("D", 0.0);
        m_z = rateMap.getDouble("z", 1.0);
        m_scaling = 1.0;

        if (m_m <= 0.0 || m_z <= 0.0) {
            throw InputFileError(WhereSetParameters, node,
                "The Starikovskiy exponents 'm' and 'z' must be positive.");
        }
    }
    else if (m_vibration_model == ModelCastela) {
        // Castela model:
        //
        // Original relaxation time:
        //
        //   tau_k = p0 / p_k
        //           * exp[a_k * (T^(-1/3) - b_k) - 18.42]
        //
        // Equivalent bimolecular rate coefficient:
        //
        //   k_k(T) = R T / p0
        //            * exp[18.42 + a_k b_k - a_k T^(-1/3)]
        //
        // Internal mapping:
        //
        //   A = R / p0
        //   b = 1
        //   B = 18.42 + a_k b_k
        //   C = -a_k
        //   D = 0
        //   E = 0
        //   scaling = 1

        requireKeys(rateMap, m_vibration_model, WhereSetParameters, {"a", "b"});
        forbidKeys(rateMap, m_vibration_model, WhereSetParameters, {m_A_str, "n", "K",
                   m_B_str, m_C_str, m_D_str, m_m_str, m_E_str, m_z_str, m_scaling_str});

        m_castela_a = rateMap["a"].asDouble();
        m_castela_b = rateMap["b"].asDouble();

        if (rateMap.hasKey(m_reference_pressure_str)) {
            m_referencePressure = rateMap.convert(m_reference_pressure_str, "Pa");
        } else {
            m_referencePressure = OneAtm;
        }

        if (m_referencePressure <= 0.0) {
            throw InputFileError(WhereSetParameters, node,
                "Castela reference-pressure must be positive.");
        }

        configureBaseFromInternalA(
            node, rate_units, GasConstant / m_referencePressure, 1.0);

        m_B = 18.42 + m_castela_a * m_castela_b;
        m_C = -m_castela_a;
        m_D = 0.0;
        m_m = 2.0 / 3.0;
        m_E = 0.0;
        m_z = 1.0;
        m_scaling = 1.0;
    }
    else {
        throw InputFileError(WhereSetParameters, node,
            "Unrecognized vibration-model '{}'. Expected 'multi-state-resolved', "
            "'starikovskiy', 'castela', or 'constant'.",
            m_vibration_model);
    }
}

void VibrationalRelaxationRate::getParameters(AnyMap& node) const
{
    if (!valid()) {
        return;
    }

    if (allowNegativePreExponentialFactor()) {
        node["negative-A"] = true;
    }

    node["vibration-model"] = m_vibration_model;

    AnyMap rateNode;

    auto storePreExponentialFactor = [&](AnyMap& target, double A) {
        if (conversionUnits().factor() != 0.0) {
            target[m_A_str].setQuantity(A, conversionUnits());
        } else {
            target[m_A_str] = A;
            target["__unconvertible__"] = true;
        }
    };

    if (m_vibration_model == ModelConstant) {
        const double tol = 1e-12;

        if (std::abs(m_b) > tol
            || std::abs(m_B) > tol
            || std::abs(m_C) > tol
            || std::abs(m_D) > tol
            || std::abs(m_E) > tol)
        {
            throw InputFileError(WhereGetParameters, node,
                "Cannot serialize this rate as 'constant': the internal "
                "parameters contain temperature-dependent terms.");
        }

        storePreExponentialFactor(rateNode, m_scaling * m_A);
    }
    else if (m_vibration_model == ModelMultiState) {
        storePreExponentialFactor(rateNode, m_A);

        rateNode[m_b_str] = m_b;
        rateNode[m_B_str] = m_B;
        rateNode[m_C_str] = m_C;
        rateNode[m_D_str] = m_D;
        rateNode[m_scaling_str] = m_scaling;
    }
    else if (m_vibration_model == ModelStarikovskiy) {
        storePreExponentialFactor(rateNode, m_A);

        rateNode["n"] = m_b;
        rateNode["K"] = m_B;
        rateNode["B"] = m_C;
        rateNode["C"] = m_D;
        rateNode["m"] = m_m;
        rateNode["D"] = m_E;
        rateNode["z"] = m_z;
    }
    else if (m_vibration_model == ModelCastela) {
        const double tol = 1e-12;

        if (std::abs(m_b - 1.0) > tol
            || std::abs(m_D) > tol
            || std::abs(m_E) > tol
            || std::abs(m_scaling - 1.0) > tol)
        {
            throw InputFileError(WhereGetParameters, node,
                "Cannot serialize this rate as 'castela': the internal "
                "parameters are not consistent with the Castela form. "
                "Expected b = 1, D = 0, E = 0, and scaling = 1.");
        }

        if (m_referencePressure <= 0.0) {
            throw InputFileError(WhereGetParameters, node,
                "Cannot serialize this rate as 'castela': "
                "reference-pressure must be positive.");
        }

        rateNode["a"] = m_castela_a;
        rateNode["b"] = m_castela_b;
        rateNode[m_reference_pressure_str].setQuantity(m_referencePressure, "Pa");
    }
    else {
        throw InputFileError(WhereGetParameters, node,
            "Unrecognized vibration-model '{}'. Expected 'multi-state-resolved', "
            "'starikovskiy', 'castela', or 'constant'.",
            m_vibration_model);
    }

    rateNode.setFlowStyle();
    node["rate-constant"] = std::move(rateNode);
}

double VibrationalRelaxationRate::ddTScaledFromStruct(
    const DetailedVibData& shared_data) const
{
    const double invT = shared_data.recipT;
    const double invT13 = std::cbrt(invT);

    return m_b * invT
           - (m_C / 3.0) * invT13 * invT
           - m_m * m_D * std::pow(invT, m_m) * invT
           - m_z * m_E * std::pow(invT, m_z) * invT;
}

void VibrationalRelaxationRate::setContext(const Reaction& rxn, const Kinetics& kin)
{
    if (rxn.reversible) {
        throw InputFileError(WhereSetContext, rxn.input,
            "Vibrational relaxation rates do not support reversible "
            "reactions.");
    }

    const string family = inferRelaxingFamily(rxn);

    if (m_vibration_model == ModelConstant) {
        validateSimpleRelaxationToGroundState(rxn, ModelConstant);
    } else if (m_vibration_model == ModelCastela) {
        validateCastelaReaction(rxn);
    } else if (m_vibration_model == ModelStarikovskiy) {
        validateSimpleRelaxationToGroundState(rxn, ModelStarikovskiy);
    } else if (m_vibration_model == ModelMultiState) {
        validateDetailedRelaxationReaction(rxn);
    } else {
        throw InputFileError(WhereSetContext, rxn.input,
            "Unrecognized vibration-model '{}'.", m_vibration_model);
    }

    registerVibrationalModelConsistency(
        kin, family, m_vibration_model, rxn.input);
}

} // namespace Cantera