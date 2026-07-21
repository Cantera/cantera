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

constexpr double VibTolerance = 1e-12;

const string WhereSetParameters = "VibrationalRelaxationRate::setParameters";
const string WhereGetParameters = "VibrationalRelaxationRate::getParameters";
const string WhereSetContext = "VibrationalRelaxationRate::setContext";

const string ModelConstant = "constant";
const string ModelMultiState = "multi-state-resolved";
const string ModelStarikovskiy = "Starikovskiy";
const string ModelCastela = "Castela";

enum class VibModel {
    Constant,
    MultiStateResolved,
    Starikovskiy,
    Castela
};

VibModel parseVibrationModel(const string& model, const AnyMap& input,
                             const string& where)
{
    if (model == ModelConstant) {
        return VibModel::Constant;
    } else if (model == ModelMultiState) {
        return VibModel::MultiStateResolved;
    } else if (model == ModelStarikovskiy) {
        return VibModel::Starikovskiy;
    } else if (model == ModelCastela) {
        return VibModel::Castela;
    }

    throw InputFileError(where, input,
        "Unrecognized vibration-model '{}'. Expected 'multi-state-resolved', "
        "'Starikovskiy', 'Castela', or 'constant'.", model);
}

const AnyMap& getRateConstantMap(const AnyMap& node)
{
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

    return rate.as<AnyMap>();
}

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

bool sameComposition(Composition diff, const Composition& b,
                     double tol = VibTolerance)
{
    for (const auto& [name, value] : b) {
        diff[name] -= value;
    }

    for (const auto& [name, value] : diff) {
        if (std::abs(value) > tol) {
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
//   NH3(v) -> Starikovskiy
//
// But it is forbidden for each reaction of the same "vibrational family"
// to have different relaxation models:
//
//   N2(v) + O  -> Castela
//   N2(v) + N2 -> Starikovskiy
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

    if (std::abs(compositionSum(rxn.reactants) - 2.0) > VibTolerance
        || std::abs(compositionSum(rxn.products) - 2.0) > VibTolerance)
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
            if (std::abs(value - 2.0) < VibTolerance) {
                collider = "N2";
            }
        } else if (std::abs(value - 1.0) < VibTolerance) {
            collider = name;
        }
    }

    if (collider != "N2" && collider != "O2" && collider != "O") {
        throw InputFileError(WhereSetContext, rxn.input,
            "The Castela relaxation model only supports colliders 'N2', "
            "'O2', and 'O'. Found collider '{}'.", collider);
    }
}

void validateSameVibrationalFamily(const std::vector<string>& species,
                                    const string& family,
                                    const AnyMap& input,
                                    const string& role)
{
    for (const auto& sp : species) {
        const string spFamily = vibrationalFamilyName(sp);
        if (spFamily != family) {
            throw InputFileError(WhereSetContext, input,
                "Invalid detailed vibrational relaxation reaction: all "
                "vibrational {} must belong to the same vibrational "
                "family. Found '{}' and '{}'.", role, family, spFamily);
        }
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

    validateSameVibrationalFamily(vibReactants, family, rxn.input, "reactants");
    validateSameVibrationalFamily(vibProducts, family, rxn.input, "products");

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
    if (std::abs(compositionSum(rxn.reactants) - 2.0) > VibTolerance
        || std::abs(compositionSum(rxn.products) - 2.0) > VibTolerance)
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
    //   vibration-model: Starikovskiy
    //   vibration-model: Castela
    //
    // The default model is multi-state-resolved.

    m_vibration_model = node.getString("vibration-model", ModelMultiState);
    const auto& rateMap = getRateConstantMap(node);

    switch (parseVibrationModel(m_vibration_model, node, WhereSetParameters)) {
    case VibModel::Constant:
        setConstantParameters(node, rateMap, rate_units);
        break;
    case VibModel::MultiStateResolved:
        setMultiStateParameters(node, rateMap, rate_units);
        break;
    case VibModel::Starikovskiy:
        setStarikovskiyParameters(node, rateMap, rate_units);
        break;
    case VibModel::Castela:
        setCastelaParameters(node, rateMap, rate_units);
        break;
    }
}

void VibrationalRelaxationRate::setConstantParameters(
    const AnyMap& node, const AnyMap& rateMap, const UnitStack& rate_units)
{
    // Constant model:
    //
    //   k(T) = A
    requireKeys(rateMap, m_vibration_model, WhereSetParameters, {m_A_str});

    forbidKeys(rateMap, m_vibration_model, WhereSetParameters,
        {m_b_str, "n", m_B_str, m_C_str, m_D_str, m_m_str,
         m_E_str, m_z_str, m_scaling_str});

    configureBaseFromYamlA(node, rate_units, rateMap[m_A_str], 0.0);
    setGenericParameters(0.0, 0.0, 0.0, 2.0 / 3.0, 0.0, 1.0, 1.0);
}

void VibrationalRelaxationRate::setMultiStateParameters(
    const AnyMap& node, const AnyMap& rateMap, const UnitStack& rate_units)
{
    // Detailed VV/VT model:
    //
    //   k(T) = scaling * A * exp(
    //       b * log(T)
    //       + B
    //       + C * T^(-1/3)
    //       + D * T^(-2/3)
    //   )
    requireKeys(rateMap, m_vibration_model, WhereSetParameters, {m_A_str});

    forbidKeys(rateMap, m_vibration_model, WhereSetParameters,
        {"n", m_m_str, m_E_str, m_z_str});

    ArrheniusBase::setParameters(node, rate_units);

    setGenericParameters(
        rateMap.getDouble(m_B_str, 0.0),
        rateMap.getDouble(m_C_str, 0.0),
        rateMap.getDouble(m_D_str, 0.0),
        2.0 / 3.0,
        0.0,
        1.0,
        rateMap.getDouble(m_scaling_str, 1.0));
}

void VibrationalRelaxationRate::setStarikovskiyParameters(
    const AnyMap& node, const AnyMap& rateMap, const UnitStack& rate_units)
{
    // User-facing formula:
    //
    //   k(T) = A * T^n * exp(
    //       K
    //       + B * T^(-1/3)
    //       + C * T^(-m)
    //       + D * T^(-z)
    //   )
    //
    // B, C, and D are signed coefficients read directly from YAML.
    requireKeys(rateMap, m_vibration_model, WhereSetParameters, {m_A_str});

    forbidKeys(rateMap, m_vibration_model, WhereSetParameters,
        {m_b_str, m_scaling_str});

    const double m = rateMap.getDouble("m", 1.0);
    const double z = rateMap.getDouble("z", 1.0);

    if (m <= 0.0 || z <= 0.0) {
        throw InputFileError(WhereSetParameters, node,
            "The Starikovskiy exponents 'm' and 'z' must be positive.");
    }

    configureBaseFromYamlA(
        node, rate_units, rateMap[m_A_str], rateMap.getDouble("n", 0.0));

    setGenericParameters(
        rateMap.getDouble("K", 0.0),
        rateMap.getDouble("B", 0.0),
        rateMap.getDouble("C", 0.0),
        m,
        rateMap.getDouble("D", 0.0),
        z,
        1.0);
}

void VibrationalRelaxationRate::setCastelaParameters(
    const AnyMap& node, const AnyMap& rateMap, const UnitStack& rate_units)
{
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

    forbidKeys(rateMap, m_vibration_model, WhereSetParameters,
        {m_A_str, "n", "K", m_B_str, m_C_str, m_D_str, m_m_str,
         m_E_str, m_z_str, m_scaling_str});

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

    setGenericParameters(
        18.42 + m_castela_a * m_castela_b,
        -m_castela_a,
        0.0,
        2.0 / 3.0,
        0.0,
        1.0,
        1.0);
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

    switch (parseVibrationModel(m_vibration_model, node, WhereGetParameters)) {
    case VibModel::Constant:
        getConstantParameters(node, rateNode);
        break;
    case VibModel::MultiStateResolved:
        getMultiStateParameters(rateNode);
        break;
    case VibModel::Starikovskiy:
        getStarikovskiyParameters(rateNode);
        break;
    case VibModel::Castela:
        getCastelaParameters(node, rateNode);
        break;
    }

    rateNode.setFlowStyle();
    node["rate-constant"] = std::move(rateNode);
}

void VibrationalRelaxationRate::storePreExponentialFactor(
    AnyMap& target, double A) const
{
    if (conversionUnits().factor() != 0.0) {
        target[m_A_str].setQuantity(A, conversionUnits());
    } else {
        target[m_A_str] = A;
        target["__unconvertible__"] = true;
    }
}

void VibrationalRelaxationRate::getConstantParameters(
    AnyMap& node, AnyMap& rateNode) const
{
    if (std::abs(m_b) > VibTolerance
        || std::abs(m_B) > VibTolerance
        || std::abs(m_C) > VibTolerance
        || std::abs(m_D) > VibTolerance
        || std::abs(m_E) > VibTolerance)
    {
        throw InputFileError(WhereGetParameters, node,
            "Cannot serialize this rate as 'constant': the internal "
            "parameters contain temperature-dependent terms.");
    }

    storePreExponentialFactor(rateNode, m_scaling * m_A);
}

void VibrationalRelaxationRate::getMultiStateParameters(
    AnyMap& rateNode) const
{
    storePreExponentialFactor(rateNode, m_A);

    rateNode[m_b_str] = m_b;
    rateNode[m_B_str] = m_B;
    rateNode[m_C_str] = m_C;
    rateNode[m_D_str] = m_D;
    rateNode[m_scaling_str] = m_scaling;
}

void VibrationalRelaxationRate::getStarikovskiyParameters(
    AnyMap& rateNode) const
{
    storePreExponentialFactor(rateNode, m_A);

    rateNode["n"] = m_b;
    rateNode["K"] = m_B;
    rateNode["B"] = m_C;
    rateNode["C"] = m_D;
    rateNode["m"] = m_m;
    rateNode["D"] = m_E;
    rateNode["z"] = m_z;
}

void VibrationalRelaxationRate::getCastelaParameters(
    AnyMap& node, AnyMap& rateNode) const
{
    if (std::abs(m_b - 1.0) > VibTolerance
        || std::abs(m_D) > VibTolerance
        || std::abs(m_E) > VibTolerance
        || std::abs(m_scaling - 1.0) > VibTolerance)
    {
        throw InputFileError(WhereGetParameters, node,
            "Cannot serialize this rate as 'Castela': the internal "
            "parameters are not consistent with the Castela form. "
            "Expected b = 1, D = 0, E = 0, and scaling = 1.");
    }

    if (m_referencePressure <= 0.0) {
        throw InputFileError(WhereGetParameters, node,
            "Cannot serialize this rate as 'Castela': "
            "reference-pressure must be positive.");
    }

    rateNode["a"] = m_castela_a;
    rateNode["b"] = m_castela_b;
    rateNode[m_reference_pressure_str].setQuantity(m_referencePressure, "Pa");
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

    switch (parseVibrationModel(m_vibration_model, rxn.input, WhereSetContext)) {
    case VibModel::Constant:
        validateSimpleRelaxationToGroundState(rxn, ModelConstant);
        break;
    case VibModel::Castela:
        validateCastelaReaction(rxn);
        break;
    case VibModel::Starikovskiy:
        validateSimpleRelaxationToGroundState(rxn, ModelStarikovskiy);
        break;
    case VibModel::MultiStateResolved:
        validateDetailedRelaxationReaction(rxn);
        break;
    }

    registerVibrationalModelConsistency(
        kin, family, m_vibration_model, rxn.input);
}

void VibrationalRelaxationRate::setGenericParameters(
    double B, double C, double D, double m, double E, double z, double scaling)
{
    m_B = B;
    m_C = C;
    m_D = D;
    m_m = m;
    m_E = E;
    m_z = z;
    m_scaling = scaling;
}

} // namespace Cantera