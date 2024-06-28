//! @file LmrRate.cpp
// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/LmrRate.h"
#include "cantera/thermo/ThermoPhase.h"
// #include "cantera/kinetics/Falloff.h"
#include "cantera/kinetics/Kinetics.h"
// #include "cantera/kinetics/Reaction.h"
// #include "cantera/kinetics/ChebyshevRate.h"
// #include "cantera/kinetics/PlogRate.h"
#include <boost/variant.hpp>
// #include "cantera/kinetics/ReactionRateFactory.h"
// #include "cantera/base/FactoryBase.h"
// #include "cantera/kinetics/Reaction.h"
// #include <sstream>

namespace Cantera{

LmrData::LmrData(){ //THIS METHOD WAS ADAPTED SOMEWHAT BLINDLY FROM FALLOFF.CPP, PLEASE VERIFY IF CORRECT
    moleFractions.resize(1, NAN);
}

bool LmrData::update(const ThermoPhase& phase, const Kinetics& kin){
    // writelog("update::1"); writelog("\n");
    double T = phase.temperature();
    double P = phase.pressure();
    int X = phase.stateMFNumber();
    if (allSpecies.empty()){
        allSpecies = phase.speciesNames(); //Get the list of all species in yaml (not just the ones for which LMRR data exists)
    }
    if (moleFractions.empty()){
        moleFractions.resize(allSpecies.size());
    }
    if (P != pressure || T != temperature || X != mfNumber) {
        ReactionData::update(T);
        pressure = P;
        logP = std::log(P);
        mfNumber=X;
        phase.getMoleFractions(moleFractions.data());
        return true;
    }
    return false;
}


void LmrData::perturbPressure(double deltaP){
    if (m_pressure_buf > 0.) {
        throw CanteraError("LmrData::perturbPressure",
            "Cannot apply another perturbation as state is already perturbed.");
    }
    m_pressure_buf = pressure;
    update(temperature,pressure*(1. + deltaP));
}

void LmrData::restore(){
    ReactionData::restore();
    if (m_pressure_buf < 0.) {
        return;
    }
    update(temperature,m_pressure_buf);
    m_pressure_buf = -1.;
}

LmrRate::LmrRate(const AnyMap& node, const UnitStack& rate_units){
    setParameters(node, rate_units);
}

void LmrRate::setParameters(const AnyMap& node, const UnitStack& rate_units){
    UnitStack eps_units{{Units(1.0), 1.0}};
    ReactionRate::setParameters(node, rate_units);
    if(!node.hasKey("collider-list")){
        throw InputFileError("LmrRate::setParameters", m_input,"Incorrect YAML input for LMR-R reaction. Please review implementation guide.");
    }
    const auto& colliders = node["collider-list"].asVector<AnyMap>();
    if(!colliders[0].hasKey("collider")){
        throw InputFileError("LmrRate::setParameters", m_input,"Incorrect YAML input for LMR-R reaction. Please review implementation guide.");
    }
    if (colliders[0]["collider"].as<std::string>() != "M"){
        throw InputFileError("LmrRate::setParameters", m_input,"The first species defined in LMR-R YAML input must be 'M'. Please review implementation guide.");
    }
    if (colliders[0].hasKey("rate-constants")){
        rateObj_M = PlogRate(colliders[0], rate_units);
        dataObj_M = PlogData();
    } else if (colliders[0].hasKey("Troe")){
        // writelog("setParameters::3"); writelog("\n");
        rateObj_M = TroeRate(colliders[0], rate_units);
        dataObj_M = FalloffData(); 
    } else if (colliders[0].hasKey("pressure-range")){
        // writelog("setParameters::4"); writelog("\n");
        rateObj_M = ChebyshevRate(colliders[0], rate_units);
        dataObj_M = ChebyshevData();
    } else {
        throw InputFileError("LmrRate::setParameters", m_input,"The M-collider must be specified in a PLOG, Troe, or Chebyshev format.");    
    }
    if (!colliders[0].hasKey("eig0") && !colliders[0].hasKey("eps")){
        throw InputFileError("LmrRate::setParameters", m_input,"A third-body efficiency (eps) or ME eigenvalue (eig0) must be provided for collider M. Please review implementation guide.");
    }
    if (colliders[0].hasKey("eps")){
        if (colliders[0]["eps"]["A"]!=1 || colliders[0]["eps"]["b"]!=0 || colliders[0]["eps"]["Ea"]!=0){ 
            throw InputFileError("LmrRate::setParameters", m_input,"The third-body efficiency (eps) must be entered for M as 'eps: {A: 1, b: 0, Ea: 0}'. Please review implementation guide.");
        }
    }
    string eig_eps_key;
    if (colliders[0].hasKey("eig0") && !colliders[0].hasKey("eps")){
        eig_eps_key="eig0";
    }
    else if (colliders[0].hasKey("eps") && !colliders[0].hasKey("eig0")){
        eig_eps_key="eps";
    }
    else{
        throw InputFileError("LmrRate::setParameters", m_input,"Cannot have both eig0 and eps provided for M. Only one is allowed, and the same choice must be maintained across any additional colliders. Please review implementation guide.");
    }
    AnyMap params;
    params["A"]=1.0;
    params["b"]=0.0;
    params["Ea"]=0.0;
    epsObj_M=ArrheniusRate(AnyValue(params),colliders[0].units(),eps_units);
    for (size_t i = 1; i < colliders.size(); i++){ //Starts at 1 because idx 0 is for "M"
        if(!colliders[i].hasKey("collider")){
                throw InputFileError("LmrRate::setParameters", m_input,"Incorrect YAML input for LMR-R reaction. Please review implementation guide.");
        }
        if (!colliders[i].hasKey(eig_eps_key)){
            throw InputFileError("LmrRate::setParameters",m_input,"All collider strengths must be defined uniformly as either eig0 or eps. No mixing and matching is allowed.");
        }
        colliderInfo.insert({colliders[i]["collider"].as<std::string>(), colliders[i]}); //Legacy parameter, used b.c. getParameters would have to be rewritten otherwise
        colliderNames.push_back(colliders[i]["collider"].as<std::string>());
        ArrheniusRate epsObj_i;
        // note: eig0 and eps are ONLY interchangeable here due to the requirement that eps_M has parameters {A: 1, b: 0, Ea: 0}
        params["A"]=colliders[i][eig_eps_key]["A"].asDouble() / colliders[0][eig_eps_key]["A"].asDouble();
        params["b"]=colliders[i][eig_eps_key]["b"].asDouble() - colliders[0][eig_eps_key]["b"].asDouble();
        params["Ea"]=colliders[i][eig_eps_key]["Ea"].asDouble() - colliders[0][eig_eps_key]["Ea"].asDouble();
        epsObj_i = ArrheniusRate(AnyValue(params),colliders[i].units(),eps_units);
        if (colliders[i].hasKey("rate-constants")){
            rateObjs.push_back(PlogRate(colliders[i], rate_units));
            dataObjs.push_back(PlogData());
            epsObjs1.push_back(epsObj_i);
            epsObjs2.push_back(epsObj_i);
        } 
        else if (colliders[i].hasKey("Troe")){
            rateObjs.push_back(TroeRate(colliders[i], rate_units));
            dataObjs.push_back(FalloffData());
            epsObjs1.push_back(epsObj_i);
            epsObjs2.push_back(epsObj_i);
        } 
        else if (colliders[i].hasKey("pressure-range")){
            rateObjs.push_back(ChebyshevRate(colliders[i], rate_units));
            dataObjs.push_back(ChebyshevData());
            epsObjs1.push_back(epsObj_i);
            epsObjs2.push_back(epsObj_i);
        } 
        else { //Collider has an eps specified, but no other info is provided. Assign it the same rate and data objects as "M"
            rateObjs.push_back(rateObj_M);
            dataObjs.push_back(dataObj_M);
            epsObjs1.push_back(epsObj_i);
            epsObjs2.push_back(epsObj_M); //NOTE: epsObj_M because logPeff is taken w.r.t. M only in this scenario
        }
    }
}

void LmrRate::validate(const string& equation, const Kinetics& kin){
    vector<double> T = {500,1000,1500,2000};
    for (size_t i=0; i<T.size(); i++){
        if (epsObj_M.evalRate(log(T[i]), 1/T[i])<0){
            throw InputFileError("LmrRate::validate", m_input,"Invalid eig0 or eps entry for M.");
        }
        for (size_t j=0; j<colliderIndices.size(); j++){
            if (epsObjs1[j].evalRate(log(T[i]), 1/T[i])<0){
                throw InputFileError("LmrRate::validate", m_input,"Invalid eig0 or eps entry for one of the specified colliders.");
            }
            if (rateObjs[j].which()!=0 && rateObjs[j].which()!=1 && rateObjs[j].which()!=2){
                throw InputFileError("LmrRate::validate", m_input,"Something went wrong... Please review implementation guide and check your k(T,P) definitions.");
            }
        }
    }
}

void LmrRate::setContext(const Reaction& rxn, const Kinetics& kin){
    for (size_t i=0; i<colliderNames.size();i++){
        colliderIndices.push_back(kin.kineticsSpeciesIndex(colliderNames[i]));
    }
    nSpecies = kin.nTotalSpecies();
}

double LmrRate::evalPlogRate(const LmrData& shared_data, DataTypes& dataObj, RateTypes& rateObj){
    PlogData& data = boost::get<PlogData>(dataObj);
    PlogRate& rate = boost::get<PlogRate>(rateObj);
    data.logP = logPeff_; //replaces logP with log of the effective pressure w.r.t. eps
    data.logT = shared_data.logT;
    // dataObj.pressure = shared_data.pressure; // THIS IS NOT CORRECT
    data.recipT = shared_data.recipT;
    // dataObj.temperature = shared_data.temperature;
    // rate.convert("P", "Pa");
    rate.updateFromStruct(data);
    return rate.evalFromStruct(data);
}

double LmrRate::evalTroeRate(const LmrData& shared_data, DataTypes& dataObj, RateTypes& rateObj){
    FalloffData& data = boost::get<FalloffData>(dataObj);
    TroeRate& rate = boost::get<TroeRate>(rateObj);
    data.conc_3b = shared_data.moleFractions;
    data.logT = shared_data.logT;
    // dataObj.molar_density = shared_data.pressure; //
    data.ready = shared_data.ready;
    data.recipT = shared_data.recipT;
    data.temperature = shared_data.temperature;
    return rate.evalFromStruct(data);
}

double LmrRate::evalChebyshevRate(const LmrData& shared_data, DataTypes& dataObj, RateTypes& rateObj){
    ChebyshevData& data = boost::get<ChebyshevData>(dataObj);
    ChebyshevRate& rate = boost::get<ChebyshevRate>(rateObj);
    data.log10P=log10(exp(logPeff_));
    // data.pressure=shared_data.pressure;
    data.recipT=shared_data.recipT;
    // data.temperature=shared_data.temperature;
    rate.updateFromStruct(data);
    return rate.evalFromStruct(data);  
}

double LmrRate::evalFromStruct(const LmrData& shared_data){
    double sigmaX_M=0.0;
    for (size_t i=0; i<nSpecies; i++){ //testing each species listed at the top of yaml file
        sigmaX_M += shared_data.moleFractions[i]; //total sum will be essentially 1, but perhaps not exactly due to Cantera's rounding conventions
    }
    eps_mix=0.0;
    size_t counter=0;
    while(counter<colliderIndices.size()){//breaks after all colliders have been located to prevent unnecessary iterations
        for (size_t i=0; i<nSpecies; i++){ //testing each species listed at the top of yaml file
            for (size_t j=0; j<colliderIndices.size(); j++){
                if (i==colliderIndices[j]){ // Species is in collider list
                    eps_mix += shared_data.moleFractions[i]*epsObjs1[j].evalRate(shared_data.logT, shared_data.recipT);
                    sigmaX_M -= shared_data.moleFractions[i];
                    counter+=1;
                    break; //breaks after collider has been located to prevent unnecessary iterations
                }
            }
        }
    }
    double eps_M = epsObj_M.evalRate(shared_data.logT, shared_data.recipT);
    eps_mix += sigmaX_M*eps_M; // add all M colliders to eps_mix in a single step
    if (eps_mix==0){
        throw InputFileError("LmrRate::evalFromStruct", m_input,"eps_mix==0 for some reason");
    }
    double k_LMR_=0.0;
    counter=0;
    while(counter<colliderIndices.size()){//breaks after all colliders have been located to prevent unnecessary iterations
        for (size_t i=0; i<nSpecies; i++){ //testing each species listed at the top of yaml file
            for (size_t j=0; j<colliderIndices.size(); j++){
                if (i==colliderIndices[j]){ // Species is in collider list
                    double eps1 = epsObjs1[j].evalRate(shared_data.logT, shared_data.recipT); 
                    double eps2 = epsObjs2[j].evalRate(shared_data.logT, shared_data.recipT); 
                    logPeff_ = shared_data.logP+log(eps_mix)-log(eps2); //NOTE: eps2 equals either eps_M or eps_i, depending on the scenario
                    if (rateObjs[j].which()==0){ // 0 means PlogRate    
                        k_LMR_ += evalPlogRate(shared_data,dataObjs[j],rateObjs[j])*eps1*shared_data.moleFractions[i]/eps_mix;
                        counter+=1;
                        break; //breaks after collider has been located to prevent unnecessary iterations
                    }
                    else if (rateObjs[j].which()==1){ // 1 means TroeRate  
                        k_LMR_ += evalTroeRate(shared_data,dataObjs[j],rateObjs[j])*eps1*shared_data.moleFractions[i]/eps_mix;
                        counter+=1;
                        break;
                    }
                    else if (rateObjs[j].which()==2){ // 2 means ChebyshevRate  
                        k_LMR_ += evalChebyshevRate(shared_data,dataObjs[j],rateObjs[j])*eps1*shared_data.moleFractions[i]/eps_mix;
                        counter+=1;
                        break;
                    }
                    else {
                        throw InputFileError("LmrRate::evalFromStruct", m_input,"Something went wrong...");
                    }
                }
            }
        }
    }
    logPeff_ = shared_data.logP+log(eps_mix)-log(eps_M);
    if (rateObj_M.which()==0){ // 0 means PlogRate
        k_LMR_ += evalPlogRate(shared_data,dataObj_M,rateObj_M)*eps_M*sigmaX_M/eps_mix;
    }
    else if (rateObj_M.which()==1){ // 1 means TroeRate 
        k_LMR_ += evalTroeRate(shared_data,dataObj_M,rateObj_M)*eps_M*sigmaX_M/eps_mix;
    }
    else if (rateObj_M.which()==2){ // 2 means ChebyshevRate
        k_LMR_ += evalChebyshevRate(shared_data,dataObj_M,rateObj_M)*eps_M*sigmaX_M/eps_mix;
    }
    return k_LMR_;
}

void LmrRate::getParameters(AnyMap& rateNode, const Units& rate_units) const{ //STILL NEED TO ACCOUNT FOR TROE AND PLOG DATA ENTRY
    vector<AnyMap> topLevelList;
    for (const auto& entry : colliderInfo) {
        string name = entry.first;
        auto colliders_i = entry.second;
        AnyMap colliderNode;
        if(colliders_i.hasKey("rate-constants")){
            colliderNode["collider"]=name;
            colliderNode["eps"]=colliders_i["eps"];
            colliderNode["rate-constants"]=colliders_i["rate-constants"];
        } 
        else if(colliders_i.hasKey("Troe")){
            colliderNode["collider"]=name;
            colliderNode["eps"]=colliders_i["eps"];
            colliderNode["low-P-rate-constant"]=colliders_i["low-P-rate-constant"];
            colliderNode["high-P-rate-constant"]=colliders_i["high-P-rate-constant"];
            colliderNode["Troe"]=colliders_i["Troe"];
        } 
        else if(colliders_i.hasKey("data")&&colliders_i.hasKey("pressure-range")&&colliders_i.hasKey("temperature-range")){ //"M" is of Chebyshev type
            colliderNode["collider"]=name;
            colliderNode["eps"]=colliders_i["eps"];
            colliderNode["temperature-range"]=colliders_i["temperature-range"];
            colliderNode["pressure-range"]=colliders_i["pressure-range"];
            colliderNode["data"]=colliders_i["data"];
        } 
        else {
            colliderNode["collider"]=name;
            colliderNode["eps"]=colliders_i["eps"];
        }
        topLevelList.push_back(std::move(colliderNode));
    }
    rateNode["collider-list"]=std::move(topLevelList);
}
}

