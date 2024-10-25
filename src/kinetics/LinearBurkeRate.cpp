//! @file LinearBurkeRate.cpp
// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/LinearBurkeRate.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/kinetics/Kinetics.h"
#include <boost/variant.hpp>

namespace Cantera
{

LinearBurkeData::LinearBurkeData()
{
    moleFractions.resize(1, NAN);
}

bool LinearBurkeData::update(const ThermoPhase& phase, const Kinetics& kin)
{
    int X = phase.stateMFNumber();
    double T = phase.temperature();
    double P = phase.pressure();
    if (moleFractions.empty()) {
        moleFractions.resize(kin.nTotalSpecies());
    }
    if (P != pressure || T != temperature || X != mf_number) {
        ReactionData::update(T);
        pressure = P;
        logP = std::log(P);
        mf_number = X;
        phase.getMoleFractions(moleFractions.data());
        return true;
    }
    return false;
}

void LinearBurkeData::perturbPressure(double deltaP)
{
    if (m_pressure_buf > 0.) {
        throw CanteraError("LinearBurkeData::perturbPressure",
            "Cannot apply another perturbation as state is already perturbed.");
    }
    m_pressure_buf = pressure;
    update(temperature, pressure * (1. + deltaP));
}

void LinearBurkeData::restore()
{
    ReactionData::restore();
    if (m_pressure_buf < 0.) {
        return;
    }
    update(temperature, m_pressure_buf);
    m_pressure_buf = -1.;
}

LinearBurkeRate::LinearBurkeRate(const AnyMap& node, const UnitStack& rate_units)
{
    setParameters(node, rate_units);
}

void LinearBurkeRate::setParameters(const AnyMap& node, const UnitStack& rate_units)
{
    UnitStack eps_units{{Units(1.0), 1.0}};
    ReactionRate::setParameters(node, rate_units);
    auto eqn = node["equation"].asString();
    if (!node.hasKey("colliders")) {
        throw InputFileError("LinearBurkeRate::setParameters", m_input,
            "'colliders' key missing from reaction '{}'.",eqn);
    }
    auto colliders = node["colliders"].asVector<AnyMap>();
    if (!colliders[0].hasKey("name")) {
        throw InputFileError("LinearBurkeRate::setParameters", m_input,
            "'colliders' key missing from reaction '{}'.",eqn);
    } else if (colliders[0]["name"].as<string>() != "M") {
        throw InputFileError("LinearBurkeRate::setParameters", m_input,
            "The first collider defined in reaction '{}' must be 'M'.",eqn);
    } else if (!colliders[0].hasKey("type")) {
        throw InputFileError("LinearBurkeRate::setParameters", m_input,
            "'type' key missing for 'M' from reaction '{}'. Must be either 'falloff'"
            " (Troe format), 'pressure-dependent-Arrhenius', or 'Chebyshev'.",eqn);
    } else if (colliders[0].hasKey("efficiency")) { //
        if (colliders[0]["efficiency"]["A"] != 1 ||
            colliders[0]["efficiency"]["b"] != 0 ||
            colliders[0]["efficiency"]["Ea"] != 0) {
            throw InputFileError("LinearBurkeRate::setParameters", m_input,
                "It is not necessary to provide an 'efficiency' for 'M' in reaction "
                "'{}', as it is always equal to 1, by definition. However, if it is "
                "entered, then it must be: 'efficiency: {A: 1, b: 0, Ea: 0}'.", eqn);
        }
    }
    if (colliders[0]["type"] == "pressure-dependent-Arrhenius") {
        m_rateObj_M = PlogRate(colliders[0], rate_units);
        m_dataObj_M = PlogData();
    } else if (colliders[0]["type"] == "falloff" && colliders[0].hasKey("Troe")) {
        TroeRate troeRateObj = TroeRate(colliders[0], rate_units);
        troeRateObj.setRateIndex(0);
        m_rateObj_M = troeRateObj;
        m_dataObj_M = FalloffData();
    } else if (colliders[0]["type"] == "Chebyshev") {
        m_rateObj_M = ChebyshevRate(colliders[0], rate_units);
        m_dataObj_M = ChebyshevData();
    } else {
        throw InputFileError("LinearBurkeRate::setParameters", m_input,
            "Collider 'M' for reaction '{}' must be specified in a"
            " pressure-dependent-Arrhenius (PLOG), falloff (Troe form), or Chebyshev"
            " format.", eqn);
    }
    AnyMap params;
    params["A"] = 1.0;
    params["b"] = 0.0;
    params["Ea"] = 0.0;
    m_epsObj_M = ArrheniusRate(AnyValue(params), colliders[0].units(), eps_units);
    string eig_eps_key;
    // If using eig0 for all colliders, then it is mandatory to specify an eig0 value
    // for the reference collider
    if (colliders[0].hasKey("eig0")) {
        eig_eps_key = "eig0";
    } else {
        eig_eps_key = "efficiency";
        colliders[0][eig_eps_key] = params;
    }
    // Start at 1 because index 0 is for "M"
    for (size_t i = 1; i < colliders.size(); i++) {
        if (!colliders[i].hasKey("name")) {
            throw InputFileError("LinearBurkeRate::setParameters", m_input,
                "The collider located at index '{}' in 'colliders' of reaction '{}'"
                " has no 'name' defined.", i, eqn);
        }
        auto nm = colliders[i]["name"].asString();
        if (eig_eps_key == "eig0") { // 'M' has 'eig0'
            if (!colliders[i].hasKey("eig0") && !colliders[i].hasKey("efficiency")) {
                throw InputFileError("LinearBurkeRate::setParameters", m_input,
                    "Collider '{}' in reaction '{}' lacks an 'eig0' key.",
                    nm, eqn);
            } else if (!colliders[i].hasKey("eig0") && colliders[i].hasKey("efficiency")) {
                throw InputFileError("LinearBurkeRate::setParameters", m_input,
                    "Since 'M' has been explicitly given 'eig0', all other "
                    "colliders must receive 'eig0' as well. No mixing and matching "
                    "of 'eig0' and 'efficiency' is allowed.", eqn);
            } else if (colliders[i].hasKey("eig0") && colliders[i].hasKey("efficiency")) {
                throw InputFileError("LinearBurkeRate::setParameters", m_input,
                    "Collider '{}' in reaction '{}' cannot contain both 'efficiency'"
                    " and 'eig0'. All additional colliders must also make the"
                    " same choice as that of 'M'.", nm, eqn);
            }
        } else { // 'M' has 'efficiency'
            if (!colliders[i].hasKey("efficiency") && !colliders[i].hasKey("eig0")) {
                throw InputFileError("LinearBurkeRate::setParameters", m_input,
                    "Collider '{}' in reaction '{}' lacks an 'efficiency' key.",
                    nm, eqn);
            } else if (!colliders[i].hasKey("efficiency") && colliders[i].hasKey("eig0")) {
                throw InputFileError("LinearBurkeRate::setParameters", m_input,
                    "Since 'M' has been (implicitly) given 'efficiency', all other "
                    "colliders must also receive an 'efficiency'. No mixing and matching "
                    "of 'efficiency' and 'eig0' is allowed.\nIf you wish to define all "
                    "colliders according to 'eig0' instead, then an 'eig0' value must "
                    "be explicitly provided for 'M' as well.", eqn);
            } else if (colliders[i].hasKey("efficiency") && colliders[i].hasKey("eig0")) {
                throw InputFileError("LinearBurkeRate::setParameters", m_input,
                    "Collider '{}' in reaction '{}' cannot contain both 'eig0'"
                    " and 'efficiency'. All additional colliders must also make the"
                    " same choice as that of 'M'.", nm, eqn);
            }
        }
        m_colliderNames.push_back(colliders[i]["name"].as<string>());

        ArrheniusRate epsObj_i;
        // eig0 and eps are ONLY interchangeable due to the requirement that eps_M have
        // parameters {A: 1, b: 0, Ea: 0}
        params["A"] = colliders[i][eig_eps_key]["A"].asDouble() /
            colliders[0][eig_eps_key]["A"].asDouble();
        params["b"] = colliders[i][eig_eps_key]["b"].asDouble() -
            colliders[0][eig_eps_key]["b"].asDouble();
        params["Ea"] = colliders[i][eig_eps_key]["Ea"].asDouble() -
            colliders[0][eig_eps_key]["Ea"].asDouble();

        if (params["A"].asDouble() < 0) {
            throw InputFileError("LinearBurkeRate::setParameters", m_input,
                "Invalid 'eig0' or 'efficiency' entry for one or more colliders.");
        }

        epsObj_i = ArrheniusRate(AnyValue(params), colliders[i].units(), eps_units);
        if (colliders[i].hasKey("type")) {
            if (colliders[i]["type"] == "pressure-dependent-Arrhenius") {
                m_rateObjs.push_back(PlogRate(colliders[i], rate_units));
                m_dataObjs.push_back(PlogData());
                m_epsObjs1.push_back(epsObj_i);
                m_epsObjs2.push_back(epsObj_i);
            } else if (colliders[i]["type"] == "falloff"
                       && colliders[i].hasKey("Troe"))
            {
                TroeRate troeRateObj = TroeRate(colliders[i], rate_units);
                troeRateObj.setRateIndex(0);
                m_rateObjs.push_back(troeRateObj);
                m_dataObjs.push_back(FalloffData());
                m_epsObjs1.push_back(epsObj_i);
                m_epsObjs2.push_back(epsObj_i);
            } else if (colliders[i]["type"] == "Chebyshev") {
                m_rateObjs.push_back(ChebyshevRate(colliders[i], rate_units));
                m_dataObjs.push_back(ChebyshevData());
                m_epsObjs1.push_back(epsObj_i);
                m_epsObjs2.push_back(epsObj_i);
            }
            m_hasRateConstant.push_back(true);
        } else {
            // Collider has an 'efficiency' specified, but no other info is provided.
            // Assign it the same rate and data objects as "M"
            m_rateObjs.push_back(m_rateObj_M);
            m_dataObjs.push_back(m_dataObj_M);
            m_epsObjs1.push_back(epsObj_i);
            m_epsObjs2.push_back(m_epsObj_M);
            m_hasRateConstant.push_back(false);
        }
    }
}

void LinearBurkeRate::validate(const string& equation, const Kinetics& kin) {}

void LinearBurkeRate::setContext(const Reaction& rxn, const Kinetics& kin)
{
    for (size_t i = 0; i<m_colliderNames.size(); i++) {
        m_colliderIndices.push_back(kin.kineticsSpeciesIndex(m_colliderNames[i]));
    }
    m_nSpecies = kin.nTotalSpecies();
}

double LinearBurkeRate::evalPlogRate(const LinearBurkeData& shared_data,
    DataTypes& dataObj, RateTypes& rateObj, double logPeff)
{
    PlogData& data = boost::get<PlogData>(dataObj);
    PlogRate& rate = boost::get<PlogRate>(rateObj);
    // Replace logP with log of the effective pressure with respect to eps
    data.logP = logPeff;
    data.logT = shared_data.logT;
    data.recipT = shared_data.recipT;
    rate.updateFromStruct(data);
    return rate.evalFromStruct(data);
}

double LinearBurkeRate::evalTroeRate(const LinearBurkeData& shared_data,
    DataTypes& dataObj, RateTypes& rateObj, double logPeff)
{
    FalloffData& data = boost::get<FalloffData>(dataObj);
    TroeRate& rate = boost::get<TroeRate>(rateObj);
    data.conc_3b = {exp(logPeff) / GasConstant / shared_data.temperature};
    data.logT = shared_data.logT;
    data.recipT = shared_data.recipT;
    data.temperature = shared_data.temperature;
    return rate.evalFromStruct(data);
}

double LinearBurkeRate::evalChebyshevRate(const LinearBurkeData& shared_data,
    DataTypes& dataObj, RateTypes& rateObj, double logPeff)
{
    ChebyshevData& data = boost::get<ChebyshevData>(dataObj);
    ChebyshevRate& rate = boost::get<ChebyshevRate>(rateObj);
    data.log10P = log10(exp(logPeff));
    data.logT = shared_data.logT;
    data.recipT = shared_data.recipT;
    rate.updateFromStruct(data);
    return rate.evalFromStruct(data);
}

double LinearBurkeRate::evalFromStruct(const LinearBurkeData& shared_data)
{
    double sigmaX_M = 0.0;
    // Test each species listed at the top of the YAML file
    for (size_t i = 0; i<m_nSpecies; i++) {
        // Total sum will be essentially 1, but perhaps not exactly due to Cantera's
        // rounding conventions
        sigmaX_M += shared_data.moleFractions[i];
    }
    double eps_mix = 0.0; // mole-fraction-weighted overall eps value of the mixtures
    for (size_t j = 0; j < m_colliderIndices.size(); j++) {
        size_t i = m_colliderIndices[j];
        eps_mix += shared_data.moleFractions[i]
                   * m_epsObjs1[j].evalRate(shared_data.logT, shared_data.recipT);
        sigmaX_M -= shared_data.moleFractions[i];
    }
    // Add all M colliders to eps_mix in a single step
    eps_mix += sigmaX_M; // eps_mix += sigmaX_M * eps_M, but eps_M = 1 always
    if (eps_mix == 0) {
        throw InputFileError("LinearBurkeRate::evalFromStruct", m_input,
            "eps_mix == 0 for some reason");
    }
    double k_LMR_ = 0.0;
    double logPeff;
    for (size_t j = 0; j < m_colliderIndices.size(); j++) {
        size_t i = m_colliderIndices[j];
        double eps1 = m_epsObjs1[j].evalRate(shared_data.logT, shared_data.recipT);
        double eps2 = m_epsObjs2[j].evalRate(shared_data.logT, shared_data.recipT);
        // eps2 equals either eps_M or eps_i, depending on the scenario
        // effective pressure as a function of eps
        logPeff = shared_data.logP + log(eps_mix) - log(eps2);
        if (m_rateObjs[j].which() == 0) { // 0 means PlogRate
            k_LMR_ += evalPlogRate(shared_data, m_dataObjs[j], m_rateObjs[j], logPeff)
                 * eps1 * shared_data.moleFractions[i] / eps_mix;
        } else if (m_rateObjs[j].which() == 1) { // 1 means TroeRate
            k_LMR_ += evalTroeRate(shared_data, m_dataObjs[j], m_rateObjs[j], logPeff)
                 * eps1 * shared_data.moleFractions[i] / eps_mix;
        } else if (m_rateObjs[j].which() == 2) { // 2 means ChebyshevRate
            k_LMR_ += evalChebyshevRate(shared_data, m_dataObjs[j], m_rateObjs[j], logPeff)
                 * eps1 * shared_data.moleFractions[i] / eps_mix;
        } else {
            throw InputFileError("LinearBurkeRate::evalFromStruct", m_input,
                "Something went wrong...");
        }
    }
    // We actually have
    // logPeff = shared_data.logP + log(eps_mix) - log(eps_M)
    // but log(eps_M)=0 always
    logPeff = shared_data.logP + log(eps_mix);
    // We actually have
    // k_LMR_+=evalPlogRate(shared_data,m_dataObj_M,m_rateObj_M)*eps_M*sigmaX_M/eps_mix
    // k_LMR_+=evalTroeRate(shared_data,m_dataObj_M,m_rateObj_M)*eps_M*sigmaX_M/eps_mix
    // etc., but eps_M = 1 always
    if (m_rateObj_M.which() == 0) { // 0 means PlogRate
        k_LMR_ += evalPlogRate(shared_data, m_dataObj_M, m_rateObj_M, logPeff) *
             sigmaX_M / eps_mix;
    } else if (m_rateObj_M.which() == 1) { // 1 means TroeRate
        k_LMR_ += evalTroeRate(shared_data, m_dataObj_M, m_rateObj_M, logPeff) *
             sigmaX_M / eps_mix;
    } else if (m_rateObj_M.which() == 2) { // 2 means ChebyshevRate
        k_LMR_ += evalChebyshevRate(shared_data, m_dataObj_M, m_rateObj_M, logPeff) *
             sigmaX_M / eps_mix;
    }
    return k_LMR_;
}

void LinearBurkeRate::getParameters(AnyMap& rateNode) const
{
    vector<AnyMap> topLevelList;
    AnyMap M_node, M_params;
    M_node["name"] = "M";
    if (m_rateObj_M.which() == 0) {
        auto& rate = boost::get<PlogRate>(m_rateObj_M);
        M_params = rate.parameters();
    } else if (m_rateObj_M.which() == 1) {
        auto& rate = boost::get<TroeRate>(m_rateObj_M);
        M_params = rate.parameters();
    } else if (m_rateObj_M.which() == 2) {
        auto& rate = boost::get<ChebyshevRate>(m_rateObj_M);
        M_params = rate.parameters();
    }
    M_node.update(M_params);
    topLevelList.push_back(std::move(M_node));

    for (size_t i = 0; i < m_colliderNames.size(); i++) {
        AnyMap collider, efficiency, params;
        collider["name"] = m_colliderNames[i];
        m_epsObjs1[i].getRateParameters(efficiency);
        collider["efficiency"] = std::move(efficiency);
        if (m_hasRateConstant[i]) {
            const auto& var_rate = m_rateObjs[i];
            if (var_rate.which() == 0) {
                auto& rate = boost::get<PlogRate>(var_rate);
                params = rate.parameters();
            } else if (var_rate.which() == 1) {
                auto& rate = boost::get<TroeRate>(var_rate);
                params = rate.parameters();
            } else if (var_rate.which() == 2) {
                auto& rate = boost::get<ChebyshevRate>(var_rate);
                params = rate.parameters();
            }
            collider.update(params);
        }
        topLevelList.push_back(std::move(collider));
    }
    rateNode["colliders"] = std::move(topLevelList);
}

}
