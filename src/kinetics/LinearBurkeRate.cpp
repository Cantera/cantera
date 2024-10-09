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
    }
    else if (!colliders[0].hasKey("eps") && !colliders[0].hasKey("eig0")) {
        throw InputFileError("LinearBurkeRate::setParameters", m_input,
                "Collider 'M' in reaction '{}' is missing an 'eps' or 'eig0' key.", eqn);
    }
    else if (colliders[0]["name"].as<string>() != "M") {
        throw InputFileError("LinearBurkeRate::setParameters", m_input,
            "The first collider defined in reaction '{}' must be 'M'.",eqn);
    }
    else if (colliders[0].hasKey("eps")) {
        if (colliders[0]["eps"]["A"] != 1 || colliders[0]["eps"]["b"] != 0 || colliders[0]["eps"]["Ea"] != 0) {
            throw InputFileError("LinearBurkeRate::setParameters", m_input,
                "The third-body efficiency (eps) must be entered for M as 'eps: {A: 1, b: 0, Ea: 0}' in reaction '{}'.", eqn);
        }
    }
    if (colliders[0].hasKey("rate-constants")) {
        rateObj_M = PlogRate(colliders[0], rate_units);
        dataObj_M = PlogData();
    }
    else if (colliders[0].hasKey("Troe")) {
        // Value of "type" is unimportant; just needed to make falloff.cpp run
        colliders[0]["type"];
        rateObj_M = TroeRate(colliders[0], rate_units);
        dataObj_M = FalloffData();
    }
    else if (colliders[0].hasKey("pressure-range")) {
        rateObj_M = ChebyshevRate(colliders[0], rate_units);
        dataObj_M = ChebyshevData();
    }
    else {
        throw InputFileError("LinearBurkeRate::setParameters", m_input,
            "Collider 'M' for reaction '{}' must be specified in a pressure-dependent-Arrhenius (PLOG), falloff (Troe form), or Chebyshev format.", eqn);
    }
    string eig_eps_key;
    if (colliders[0].hasKey("eig0") && !colliders[0].hasKey("eps")) {
        eig_eps_key = "eig0";
    }
    else if (colliders[0].hasKey("eps") && !colliders[0].hasKey("eig0")) {
        eig_eps_key = "eps";
    }
    else{
        throw InputFileError("LinearBurkeRate::setParameters", m_input,
            "Collider 'M' in reaction '{}' cannot contain both eps and eig0. Any additional colliders must also make same choice as that of M.", eqn);
    }
    AnyMap params;
    params["A"] = 1.0;
    params["b"] = 0.0;
    params["Ea"] = 0.0;
    epsObj_M = ArrheniusRate(AnyValue(params), colliders[0].units(), eps_units);
    // Start at 1 because index 0 is for "M"
    for (size_t i = 1; i < colliders.size(); i++){
        if (!colliders[i].hasKey("name")) {
                throw InputFileError("LinearBurkeRate::setParameters", m_input,
                    "The collider located at index '{}' in 'colliders' of reaction '{}' has no 'name' defined.", std::to_string(i), eqn);
        }
        auto nm = colliders[i]["name"].asString();
        if (!colliders[i].hasKey("eps") && !colliders[0].hasKey("eig0")) {
            throw InputFileError("LinearBurkeRate::setParameters", m_input,
                "Collider '{}' in reaction '{}' is missing an 'eps' or 'eig0' key.", nm, eqn);
        }
        if (!colliders[i].hasKey(eig_eps_key)) {
            throw InputFileError("LinearBurkeRate::setParameters", m_input,
                "All collision efficiencies in reaction '{}' must be defined uniformly as either eig0 or eps. No mixing and matching is allowed.", eqn);
        }
        // Save data to colliderInfo, which will make it accessible by getParameters
        colliderInfo.insert({colliders[i]["name"].as<string>(), colliders[i]});
        colliderNames.push_back(colliders[i]["name"].as<string>());
        ArrheniusRate epsObj_i;
        // eig0 and eps are ONLY interchangeable due to the requirement that eps_M have parameters {A: 1, b: 0, Ea: 0}
        params["A"] = colliders[i][eig_eps_key]["A"].asDouble() / colliders[0][eig_eps_key]["A"].asDouble();
        params["b"] = colliders[i][eig_eps_key]["b"].asDouble() - colliders[0][eig_eps_key]["b"].asDouble();
        params["Ea"] = colliders[i][eig_eps_key]["Ea"].asDouble() - colliders[0][eig_eps_key]["Ea"].asDouble();

        if (params["A"].asDouble() < 0) {
            throw InputFileError("LinearBurkeRate::setParameters", m_input,
                "Invalid eig0 or eps entry for one or more colliders.");
        }

        epsObj_i = ArrheniusRate(AnyValue(params), colliders[i].units(), eps_units);
        if (colliders[i].hasKey("rate-constants")) {
            rateObjs.push_back(PlogRate(colliders[i], rate_units));
            dataObjs.push_back(PlogData());
            epsObjs1.push_back(epsObj_i);
            epsObjs2.push_back(epsObj_i);
        }
        else if (colliders[i].hasKey("Troe")) {
            // Value of "type" is unimportant; just needed to make falloff.cpp run
            colliders[i]["type"];
            rateObjs.push_back(TroeRate(colliders[i], rate_units));
            dataObjs.push_back(FalloffData());
            epsObjs1.push_back(epsObj_i);
            epsObjs2.push_back(epsObj_i);
        }
        else if (colliders[i].hasKey("pressure-range")) {
            rateObjs.push_back(ChebyshevRate(colliders[i], rate_units));
            dataObjs.push_back(ChebyshevData());
            epsObjs1.push_back(epsObj_i);
            epsObjs2.push_back(epsObj_i);
        }
        // Collider has an eps specified, but no other info is provided. Assign it the same rate and data objects as "M"
        else {
            rateObjs.push_back(rateObj_M);
            dataObjs.push_back(dataObj_M);
            epsObjs1.push_back(epsObj_i);
            epsObjs2.push_back(epsObj_M);
        }
    }
}

void LinearBurkeRate::validate(const string& equation, const Kinetics& kin){}

void LinearBurkeRate::setContext(const Reaction& rxn, const Kinetics& kin)
{
    for (size_t i = 0; i<colliderNames.size(); i++){
        colliderIndices.push_back(kin.kineticsSpeciesIndex(colliderNames[i]));
    }
    nSpecies = kin.nTotalSpecies();
}

double LinearBurkeRate::evalPlogRate(const LinearBurkeData& shared_data, DataTypes& dataObj, RateTypes& rateObj)
{
    PlogData& data = boost::get<PlogData>(dataObj);
    PlogRate& rate = boost::get<PlogRate>(rateObj);
    // Replace logP with log of the effective pressure with respect to eps
    data.logP = logPeff_;
    data.logT = shared_data.logT;
    data.recipT = shared_data.recipT;
    rate.updateFromStruct(data);
    return rate.evalFromStruct(data);
}

double LinearBurkeRate::evalTroeRate(const LinearBurkeData& shared_data, DataTypes& dataObj, RateTypes& rateObj)
{
    FalloffData& data = boost::get<FalloffData>(dataObj);
    TroeRate& rate = boost::get<TroeRate>(rateObj);
    data.conc_3b = {exp(logPeff_) / GasConstant / shared_data.temperature};
    data.logT = shared_data.logT;
    data.recipT = shared_data.recipT;
    data.temperature = shared_data.temperature;
    rate.setRateIndex(0);
    return rate.evalFromStruct(data);
}

double LinearBurkeRate::evalChebyshevRate(const LinearBurkeData& shared_data, DataTypes& dataObj, RateTypes& rateObj)
{
    ChebyshevData& data = boost::get<ChebyshevData>(dataObj);
    ChebyshevRate& rate = boost::get<ChebyshevRate>(rateObj);
    data.log10P = log10(exp(logPeff_));
    data.logT = shared_data.logT;
    data.recipT = shared_data.recipT;
    rate.updateFromStruct(data);
    return rate.evalFromStruct(data);
}

double LinearBurkeRate::evalFromStruct(const LinearBurkeData& shared_data)
{
    double sigmaX_M = 0.0;
    // Test each species listed at the top of the YAML file
    for (size_t i = 0; i<nSpecies; i++){
        // Total sum will be essentially 1, but perhaps not exactly due to Cantera's rounding conventions
        sigmaX_M += shared_data.moleFractions[i];
    }
    eps_mix = 0.0;
    for (size_t j = 0; j < colliderIndices.size(); j++) {
        size_t i = colliderIndices[j];
        eps_mix += shared_data.moleFractions[i] * epsObjs1[j].evalRate(shared_data.logT, shared_data.recipT);
        sigmaX_M -= shared_data.moleFractions[i];
    }
    // Add all M colliders to eps_mix in a single step
    eps_mix += sigmaX_M; // eps_mix += sigmaX_M * eps_M, but eps_M = 1 always
    if (eps_mix == 0) {
        throw InputFileError("LinearBurkeRate::evalFromStruct", m_input, "eps_mix == 0 for some reason");
    }
    double k_LMR_ = 0.0;
    for (size_t j = 0; j < colliderIndices.size(); j++) {
        size_t i = colliderIndices[j];
        double eps1 = epsObjs1[j].evalRate(shared_data.logT, shared_data.recipT);
        double eps2 = epsObjs2[j].evalRate(shared_data.logT, shared_data.recipT);
        // eps2 equals either eps_M or eps_i, depending on the scenario
        logPeff_ = shared_data.logP + log(eps_mix) - log(eps2);
        // 0 means PlogRate
        if (rateObjs[j].which() == 0) {
            k_LMR_ += evalPlogRate(shared_data, dataObjs[j], rateObjs[j]) * eps1 * shared_data.moleFractions[i] / eps_mix;
        }
        // 1 means TroeRate
        else if (rateObjs[j].which() == 1) {
            k_LMR_ += evalTroeRate(shared_data, dataObjs[j], rateObjs[j]) * eps1 * shared_data.moleFractions[i] / eps_mix;
        }
        // 2 means ChebyshevRate
        else if (rateObjs[j].which() == 2) {
            k_LMR_ += evalChebyshevRate(shared_data, dataObjs[j], rateObjs[j]) * eps1 * shared_data.moleFractions[i] / eps_mix;
        }
        else {
            throw InputFileError("LinearBurkeRate::evalFromStruct", m_input, "Something went wrong...");
        }
    }
    // logPeff_ = shared_data.logP + log(eps_mix) - log(eps_M), but log(eps_M)=0 always
    logPeff_ = shared_data.logP + log(eps_mix);
    if (rateObj_M.which() == 0) { // 0 means PlogRate
        // k_LMR_ += evalPlogRate(shared_data, dataObj_M, rateObj_M) * eps_M * sigmaX_M / eps_mix, but eps_M = 1 always
        k_LMR_ += evalPlogRate(shared_data, dataObj_M, rateObj_M) * sigmaX_M / eps_mix;
    }
    else if (rateObj_M.which() == 1) { // 1 means TroeRate
        k_LMR_ += evalTroeRate(shared_data, dataObj_M, rateObj_M) * sigmaX_M / eps_mix;
    }
    else if (rateObj_M.which() == 2) { // 2 means ChebyshevRate
        k_LMR_ += evalChebyshevRate(shared_data, dataObj_M, rateObj_M) * sigmaX_M / eps_mix;
    }
    return k_LMR_;
}

void LinearBurkeRate::getParameters(AnyMap& rateNode) const
{
    vector<AnyMap> topLevelList;
    for (const auto& entry : colliderInfo) {
        string name = entry.first;
        auto colliders_i = entry.second;
        AnyMap colliderNode;
        if(colliders_i.hasKey("rate-constants")) {
            colliderNode["name"] = name;
            colliderNode["eps"] = colliders_i["eps"];
            colliderNode["rate-constants"] = colliders_i["rate-constants"];
        }
        else if(colliders_i.hasKey("Troe")) {
            colliderNode["name"] = name;
            colliderNode["eps"] = colliders_i["eps"];
            colliderNode["low-P-rate-constant"] = colliders_i["low-P-rate-constant"];
            colliderNode["high-P-rate-constant"] = colliders_i["high-P-rate-constant"];
            colliderNode["Troe"] = colliders_i["Troe"];
        }
        else if(colliders_i.hasKey("data") && colliders_i.hasKey("pressure-range") && colliders_i.hasKey("temperature-range")) {
            colliderNode["name"] = name;
            colliderNode["eps"] = colliders_i["eps"];
            colliderNode["temperature-range"] = colliders_i["temperature-range"];
            colliderNode["pressure-range"] = colliders_i["pressure-range"];
            colliderNode["data"] = colliders_i["data"];
        }
        else {
            colliderNode["name"] = name;
            colliderNode["eps"] = colliders_i["eps"];
        }
        topLevelList.push_back(std::move(colliderNode));
    }
    rateNode["colliders"] = std::move(topLevelList);
}

}
