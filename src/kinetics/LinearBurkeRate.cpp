//! @file LinearBurkeRate.cpp
// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/LinearBurkeRate.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/kinetics/Kinetics.h"
#include <boost/variant.hpp>

namespace Cantera
{

LmrData::LmrData()
{
    moleFractions.resize(1, NAN);
}

bool LmrData::update(const ThermoPhase& phase, const Kinetics& kin)
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

void LmrData::perturbPressure(double deltaP)
{
    if (m_pressure_buf > 0.) {
        throw CanteraError("LmrData::perturbPressure",
            "Cannot apply another perturbation as state is already perturbed.");
    }
    m_pressure_buf = pressure;
    update(temperature, pressure * (1. + deltaP));
}

void LmrData::restore()
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
    if (!node.hasKey("collider-list")) {
        throw InputFileError("LinearBurkeRate::setParameters", m_input,
            "Incorrect YAML input for LMR-R reaction. Please review implementation guide.");
    }
    auto colliders = node["collider-list"].asVector<AnyMap>();
    if (!colliders[0].hasKey("collider")) {
        throw InputFileError("LinearBurkeRate::setParameters", m_input,
            "Incorrect YAML input for LMR-R reaction. Please review implementation guide.");
    }
    else if (colliders[0]["collider"].as<string>() != "M") {
        throw InputFileError("LinearBurkeRate::setParameters", m_input,
            "The first species defined in LMR-R YAML input must be 'M'. Please review implementation guide.");
    }
    else if (!colliders[0].hasKey("eig0") && !colliders[0].hasKey("eps")) {
        throw InputFileError("LinearBurkeRate::setParameters", m_input,
            "A third-body efficiency (eps) or ME eigenvalue (eig0) must be provided for collider M. Please review implementation guide.");
    }
    else if (colliders[0].hasKey("eps")) {
        if (colliders[0]["eps"]["A"] != 1 || colliders[0]["eps"]["b"] != 0 || colliders[0]["eps"]["Ea"] != 0) {
            throw InputFileError("LinearBurkeRate::setParameters", m_input,
                "The third-body efficiency (eps) must be entered for M as 'eps: {A: 1, b: 0, Ea: 0}'. Please review implementation guide.");
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
            "The M-collider must be specified in a pressure-dependent-Arrhenius (PLOG), falloff (Troe form), or Chebyshev format.");
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
            "Cannot have both eig0 and eps chosen for M. Any additional colliders must make same choice as that of M. Please review implementation guide.");
    }
    AnyMap params;
    params["A"] = 1.0;
    params["b"] = 0.0;
    params["Ea"] = 0.0;
    epsObj_M = ArrheniusRate(AnyValue(params), colliders[0].units(), eps_units);
    // Start at 1 because index 0 is for "M"
    for (size_t i = 1; i < colliders.size(); i++){
        if (!colliders[i].hasKey("collider")) {
                throw InputFileError("LinearBurkeRate::setParameters", m_input,
                    "Incorrect YAML input for LMR-R reaction. Please review implementation guide.");
        }
        if (!colliders[i].hasKey(eig_eps_key)) {
            throw InputFileError("LinearBurkeRate::setParameters", m_input,
                "All collider strengths must be defined uniformly as either eig0 or eps. No mixing and matching is allowed.");
        }
        // Save data to colliderInfo, which will make it accessible by getParameters
        colliderInfo.insert({colliders[i]["collider"].as<string>(), colliders[i]});
        colliderNames.push_back(colliders[i]["collider"].as<string>());
        ArrheniusRate epsObj_i;
        // eig0 and eps are ONLY interchangeable due to the requirement that eps_M have parameters {A: 1, b: 0, Ea: 0}
        params["A"] = colliders[i][eig_eps_key]["A"].asDouble() / colliders[0][eig_eps_key]["A"].asDouble();
        params["b"] = colliders[i][eig_eps_key]["b"].asDouble() - colliders[0][eig_eps_key]["b"].asDouble();
        params["Ea"] = colliders[i][eig_eps_key]["Ea"].asDouble() - colliders[0][eig_eps_key]["Ea"].asDouble();
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

void LinearBurkeRate::validate(const string& equation, const Kinetics& kin)
{
    vector<double> T = {500, 1000, 1500, 2000};
    for (size_t i = 0; i < T.size(); i++){
        if (epsObj_M.evalRate(log(T[i]), 1 / T[i]) < 0) {
            throw InputFileError("LinearBurkeRate::validate", m_input,
                "Invalid eig0 or eps entry for M.");
        }
        for (size_t j = 0; j<colliderIndices.size(); j++) {
            if (epsObjs1[j].evalRate(log(T[i]), 1 / T[i]) < 0){
                throw InputFileError("LinearBurkeRate::validate", m_input,
                    "Invalid eig0 or eps entry for one of the specified colliders.");
            }
            if (rateObjs[j].which() != 0 && rateObjs[j].which() != 1 && rateObjs[j].which() != 2) {
                throw InputFileError("LinearBurkeRate::validate", m_input,
                    "Something went wrong... Please review implementation guide and check your k(T,P) definitions.");
            }
        }
    }
}

void LinearBurkeRate::setContext(const Reaction& rxn, const Kinetics& kin)
{
    for (size_t i = 0; i<colliderNames.size(); i++){
        colliderIndices.push_back(kin.kineticsSpeciesIndex(colliderNames[i]));
    }
    nSpecies = kin.nTotalSpecies();
}

double LinearBurkeRate::evalPlogRate(const LmrData& shared_data, DataTypes& dataObj, RateTypes& rateObj)
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

double LinearBurkeRate::evalTroeRate(const LmrData& shared_data, DataTypes& dataObj, RateTypes& rateObj)
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

double LinearBurkeRate::evalChebyshevRate(const LmrData& shared_data, DataTypes& dataObj, RateTypes& rateObj)
{
    ChebyshevData& data = boost::get<ChebyshevData>(dataObj);
    ChebyshevRate& rate = boost::get<ChebyshevRate>(rateObj);
    data.log10P = log10(exp(logPeff_));
    data.logT = shared_data.logT;
    data.recipT = shared_data.recipT;
    rate.updateFromStruct(data);
    return rate.evalFromStruct(data);
}

double LinearBurkeRate::evalFromStruct(const LmrData& shared_data)
{
    double sigmaX_M = 0.0;
    // Test each species listed at the top of the YAML file
    for (size_t i = 0; i<nSpecies; i++){
        // Total sum will be essentially 1, but perhaps not exactly due to Cantera's rounding conventions
        sigmaX_M += shared_data.moleFractions[i];
    }
    eps_mix = 0.0;
    size_t counter = 0;
    // Break loop when all colliders have been located to avoid unnecessary iterations
    while(counter<colliderIndices.size()){
        // Test each species listed at the top of the YAML file
        for (size_t i = 0; i<nSpecies; i++){
            for (size_t j = 0; j<colliderIndices.size(); j++){
                // Check if the species is specified as an LMR-R collider
                if (i == colliderIndices[j]) {
                    eps_mix += shared_data.moleFractions[i] * epsObjs1[j].evalRate(shared_data.logT, shared_data.recipT);
                    sigmaX_M -= shared_data.moleFractions[i];
                    counter += 1;
                     // Break after collider has been located to avoid unnecessary iterations
                    break;
                }
            }
        }
    }
    double eps_M = epsObj_M.evalRate(shared_data.logT, shared_data.recipT);
    // Add all M colliders to eps_mix in a single step
    eps_mix += sigmaX_M * eps_M;
    if (eps_mix == 0) {
        throw InputFileError("LinearBurkeRate::evalFromStruct", m_input, "eps_mix == 0 for some reason");
    }
    double k_LMR_ = 0.0;
    counter = 0;
    while(counter<colliderIndices.size()){
        for (size_t i = 0; i < nSpecies; i++){
            for (size_t j = 0; j < colliderIndices.size(); j++){
                if (i == colliderIndices[j]) {
                    double eps1 = epsObjs1[j].evalRate(shared_data.logT, shared_data.recipT);
                    double eps2 = epsObjs2[j].evalRate(shared_data.logT, shared_data.recipT);
                    // eps2 equals either eps_M or eps_i, depending on the scenario
                    logPeff_ = shared_data.logP + log(eps_mix) - log(eps2);
                    // 0 means PlogRate
                    if (rateObjs[j].which() == 0) {
                        k_LMR_ += evalPlogRate(shared_data, dataObjs[j], rateObjs[j]) * eps1 * shared_data.moleFractions[i] / eps_mix;
                        counter += 1;
                        break;
                    }
                    // 1 means TroeRate
                    else if (rateObjs[j].which() == 1) {
                        k_LMR_ += evalTroeRate(shared_data, dataObjs[j], rateObjs[j]) * eps1 * shared_data.moleFractions[i] / eps_mix;
                        counter += 1;
                        break;
                    }
                    // 2 means ChebyshevRate
                    else if (rateObjs[j].which() == 2) {
                        k_LMR_ += evalChebyshevRate(shared_data, dataObjs[j], rateObjs[j]) * eps1 * shared_data.moleFractions[i] / eps_mix;
                        counter += 1;
                        break;
                    }
                    else {
                        throw InputFileError("LinearBurkeRate::evalFromStruct", m_input, "Something went wrong...");
                    }
                }
            }
        }
    }
    logPeff_ = shared_data.logP + log(eps_mix) - log(eps_M);
    // 0 means PlogRate
    if (rateObj_M.which() == 0) {
        k_LMR_ += evalPlogRate(shared_data, dataObj_M, rateObj_M) * eps_M * sigmaX_M / eps_mix;
    }
    // 1 means TroeRate
    else if (rateObj_M.which() == 1) {
        k_LMR_ += evalTroeRate(shared_data, dataObj_M, rateObj_M) * eps_M * sigmaX_M / eps_mix;
    }
    // 2 means ChebyshevRate
    else if (rateObj_M.which() == 2) {
        k_LMR_ += evalChebyshevRate(shared_data, dataObj_M, rateObj_M) * eps_M * sigmaX_M / eps_mix;
    }
    return k_LMR_;
}

void LinearBurkeRate::getParameters(AnyMap& rateNode, const Units& rate_units) const
{
    vector<AnyMap> topLevelList;
    for (const auto& entry : colliderInfo) {
        string name = entry.first;
        auto colliders_i = entry.second;
        AnyMap colliderNode;
        if(colliders_i.hasKey("rate-constants")) {
            colliderNode["collider"] = name;
            colliderNode["eps"] = colliders_i["eps"];
            colliderNode["rate-constants"] = colliders_i["rate-constants"];
        }
        else if(colliders_i.hasKey("Troe")) {
            colliderNode["collider"] = name;
            colliderNode["eps"] = colliders_i["eps"];
            colliderNode["low-P-rate-constant"] = colliders_i["low-P-rate-constant"];
            colliderNode["high-P-rate-constant"] = colliders_i["high-P-rate-constant"];
            colliderNode["Troe"] = colliders_i["Troe"];
        }
        else if(colliders_i.hasKey("data") && colliders_i.hasKey("pressure-range") && colliders_i.hasKey("temperature-range")) {
            colliderNode["collider"] = name;
            colliderNode["eps"] = colliders_i["eps"];
            colliderNode["temperature-range"] = colliders_i["temperature-range"];
            colliderNode["pressure-range"] = colliders_i["pressure-range"];
            colliderNode["data"] = colliders_i["data"];
        }
        else {
            colliderNode["collider"] = name;
            colliderNode["eps"] = colliders_i["eps"];
        }
        topLevelList.push_back(std::move(colliderNode));
    }
    rateNode["collider-list"] = std::move(topLevelList);
}

}