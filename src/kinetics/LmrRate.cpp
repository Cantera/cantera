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
    
    //Determine whether collider strength is represented by ME eigenvalues (eig0), or third-body efficiency relative to M (eps) 
    if (colliders[0].hasKey("eig0")){
        AnyMap params;
        params["A"]=1.0;
        params["b"]=0.0;
        params["Ea"]=0.0;
        epsObj_M=ArrheniusRate(AnyValue(params),colliders[0].units(),eps_units);
        for (size_t i = 1; i < colliders.size(); i++){
            if(!colliders[i].hasKey("collider")){
                throw InputFileError("LmrRate::setParameters", m_input,"Incorrect YAML input for LMR-R reaction. Please review implementation guide.");
            }
            if (!colliders[i].hasKey("eig0")){
                throw InputFileError("LmrRate::setParameters",m_input,"All collider strengths must be defined uniformly as either eig0 or eps. No mixing and matching is allowed.");
            } 
            AnyMap params;
            params["A"]=colliders[i]["eig0"]["A"].asDouble() / colliders[0]["eig0"]["A"].asDouble();
            params["b"]=colliders[i]["eig0"]["b"].asDouble() - colliders[0]["eig0"]["b"].asDouble();
            params["Ea"]=colliders[i]["eig0"]["Ea"].asDouble() - colliders[0]["eig0"]["Ea"].asDouble();
            // AnyMap params = {{"A", colliders[i]["eig0"]["A"].asDouble() / colliders[0]["eig0"]["A"].asDouble()},
            //                  {"b", colliders[i]["eig0"]["b"].asDouble() - colliders[0]["eig0"]["b"].asDouble()},
            //                  {"Ea", colliders[i]["eig0"]["Ea"].asDouble() - colliders[0]["eig0"]["Ea"].asDouble()}};
            epsObjs.push_back(ArrheniusRate(AnyValue(params),colliders[i].units(),eps_units));
        }
    }
    else if (colliders[0].hasKey("eps") && colliders[0].hasKey("collider")){ //Already in relative terms, no need to convert
        epsObj_M=ArrheniusRate(AnyValue(colliders[0]["eps"]),colliders[0].units(),eps_units);
        for (size_t i = 1; i < colliders.size(); i++){
            if(!colliders[i].hasKey("collider")){
                throw InputFileError("LmrRate::setParameters", m_input,"Incorrect YAML input for LMR-R reaction. Please review implementation guide.");
            }
            if (!colliders[i].hasKey("eps")){
                throw InputFileError("LmrRate::setParameters",m_input,"All collider strengths must be defined uniformly as either eig0 or eps. No mixing and matching is allowed.");
            }
            epsObjs.push_back(ArrheniusRate(AnyValue(colliders[i]["eps"]),colliders[i].units(),eps_units));
        }
    }
    else {
        throw InputFileError("LmrRate::setParameters", m_input,"A third-body efficiency (eps) or ME eigenvalue (eig0) must be provided for all explicitly declared colliders in LMRR yaml entry. Please review implementation guide.");
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
    for (size_t i = 1; i < colliders.size(); i++){ //Starts at 1 because idx 0 is for "M"
        colliderInfo.insert({colliders[i]["collider"].as<std::string>(), colliders[i]}); //Legacy parameter, used b.c. getParameters would have to be rewritten otherwise
        colliderNames.push_back(colliders[i]["collider"].as<std::string>());
        if (colliders[i].hasKey("rate-constants")){
            rateObjs.push_back(PlogRate(colliders[i], rate_units));
            dataObjs.push_back(PlogData());
        } else if (colliders[i].hasKey("Troe")){
            rateObjs.push_back(TroeRate(colliders[i], rate_units));
            dataObjs.push_back(FalloffData());
        } else if (colliders[i].hasKey("pressure-range")){
            rateObjs.push_back(ChebyshevRate(colliders[i], rate_units));
            dataObjs.push_back(ChebyshevData());
        } else { //Collider has an eps specified, but no other info is provided. Assign it the same rate and data objects as "M"
            rateObjs.push_back(rateObj_M);
            dataObjs.push_back(dataObj_M);
        }
    }
}

void LmrRate::validate(const string& equation, const Kinetics& kin){
    // writelog("validate::1"); writelog("\n");
}

void LmrRate::setContext(const Reaction& rxn, const Kinetics& kin){
    // colliderIndices.reserve(colliderNames.size());
    for (size_t i=0; i<colliderNames.size();i++){
        // writelog("setContext::2"); writelog("\n");
        colliderIndices.push_back(kin.kineticsSpeciesIndex(colliderNames[i]));
        // // writelog("coll idx = ",val); writelog("\n");
        // writelog("coll idx = " + std::to_string(kin.kineticsSpeciesIndex(colliderNames[i])) + "\n");
    }
    nSpecies = kin.nTotalSpecies();
    for (size_t i=0; i<colliderIndices.size();i++){
        // writelog("colliderIndices["+std::to_string(i)+"] = " + std::to_string(colliderIndices[i]) + "\n");
    }
    // writelog("nSpecies = " + std::to_string(nSpecies) + "\n");
    // if (colliderIndices.empty()){
    //     throw InputFileError("LmrRate::setContext", m_input,"collidersIndices is empty for some reason");
    // }
}

double LmrRate::evalPlogRate(const LmrData& shared_data, DataTypes& dataObj, RateTypes& rateObj, double& eps){
    PlogData& data = boost::get<PlogData>(dataObj);
    PlogRate& rate = boost::get<PlogRate>(rateObj);
    //writelog("eps (3) = "+std::to_string(eps)+"\n");
    data.logP = shared_data.logP+log(eps_mix)-log(eps); //replaces logP with log of the effective pressure w.r.t. eps
    data.logT = shared_data.logT;
    // dataObj.pressure = shared_data.pressure; // THIS IS NOT CORRECT
    data.recipT = shared_data.recipT;
    // dataObj.temperature = shared_data.temperature;
    // rate.convert("P", "Pa");
    rate.updateFromStruct(data);
    //writelog("eps (4) = "+std::to_string(eps)+"\n");
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

double LmrRate::evalChebyshevRate(const LmrData& shared_data, DataTypes& dataObj, RateTypes& rateObj, double& eps){
    ChebyshevData& data = boost::get<ChebyshevData>(dataObj);
    ChebyshevRate& rate = boost::get<ChebyshevRate>(rateObj);
    data.log10P=log10(exp(shared_data.logP+log(eps_mix)-log(eps)));
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
    // writelog("sigmaX_M (0) = " + std::to_string(sigmaX_M) + "\n"); 
    
    eps_mix=0.0;
    size_t counter=0;
    for (size_t i=0; i<nSpecies; i++){ //testing each species listed at the top of yaml file
        for (size_t j=0; j<colliderIndices.size(); j++){
            if (i==colliderIndices[j]){ // Species is in collider list
                eps_mix += shared_data.moleFractions[i]*epsObjs[j].evalRate(shared_data.logT, shared_data.recipT);
                // writelog("eps[j] = " + std::to_string(epsObjs[j].evalRate(logT_, recipT_)) + " logT = " + std::to_string(logT_) + " recipT = " + std::to_string(recipT_) + "\n"); 
                sigmaX_M -= shared_data.moleFractions[i];
                // writelog("sigmaX_M (1) = " + std::to_string(sigmaX_M) + "\n"); 
                counter+=1;
                break; //breaks after collider has been located to prevent unnecessary iterations
            }
        }
        if (counter==colliderIndices.size()){
            break; //breaks after all colliders have been located to prevent unnecessary iterations
        }
    }

    // writelog("eps_mix (0) = " + std::to_string(eps_mix) + "\n"); 
    // writelog("sigmaX_M (2) = " + std::to_string(sigmaX_M) + "\n"); 
    double eps_M = epsObj_M.evalRate(shared_data.logT, shared_data.recipT);
    //writelog("eps_M (1) = "+std::to_string(eps_M)+"\n");
    // writelog("eps_M = " + std::to_string(eps_M) + "\n"); 
    eps_mix += sigmaX_M*eps_M; // add all M colliders to eps_mix in a single step
    // writelog("eps_mix (1) = " + std::to_string(eps_mix) + "\n"); 
    if (eps_mix==0){
        throw InputFileError("LmrRate::evalFromStruct", m_input,"eps_mix==0 for some reason");
    }

    double k_LMR_=0.0;
    // writelog("k_LMR (1) = " + std::to_string(k_LMR_) + "\n"); 
    counter=0;
    for (size_t i=0; i<nSpecies; i++){ //testing each species listed at the top of yaml file
        // writelog("i = " + std::to_string(i) + ", X_i = " + std::to_string(shared_data.moleFractions[i]) + "\n"); 
        for (size_t j=0; j<colliderIndices.size(); j++){
            //writelog("colliderIndices[" + std::to_string(j) + "] = " + std::to_string(colliderIndices[j]) + "\n"); 
            // writelog("j = " + std::to_string(j) + "\n"); 
            if (i==colliderIndices[j]){ // Species is in collider list
                // writelog("i = colliderIndices[j] = " + std::to_string(colliderIndices[j]) + ", X_i = " + std::to_string(shared_data.moleFractions[i]) + "\n");
                
                // writelog("eps_i = " + std::to_string(epsObjs[j].evalRate(logT_, recipT_))+"\n");
                // if (eps==0){
                //     throw InputFileError("LmrRate::evalFromStruct", m_input,"eps is 0 for some reason");
                // }
                double eps = epsObjs[j].evalRate(shared_data.logT, shared_data.recipT); 
                //writelog("eps (1) = "+std::to_string(eps)+"\n");
                if (rateObjs[j].which()==0){ // 0 means PlogRate    
                    //writelog("eps (2) = "+std::to_string(eps)+"\n");
                    k_LMR_ += evalPlogRate(shared_data,dataObjs[j],rateObjs[j],eps_M)*eps*shared_data.moleFractions[i]/eps_mix;
                    //writelog("eps (5) = "+std::to_string(eps)+"\n");
                    counter+=1;
                    break; //breaks after collider has been located to prevent unnecessary iterations
                }
                else if (rateObjs[j].which()==1){ // 1 means TroeRate  
                    double eps = epsObjs[j].evalRate(shared_data.logT, shared_data.recipT);
                    k_LMR_ += evalTroeRate(shared_data,dataObjs[j],rateObjs[j])*eps*shared_data.moleFractions[i]/eps_mix;
                }
                else if (rateObjs[j].which()==2){ // 2 means ChebyshevRate  
                    double eps = epsObjs[j].evalRate(shared_data.logT, shared_data.recipT);
                    k_LMR_ += evalChebyshevRate(shared_data,dataObjs[j],rateObjs[j],eps)*eps*shared_data.moleFractions[i]/eps_mix;
                }
                else {
                    throw InputFileError("LmrRate::evalFromStruct", m_input,"Something went wrong...");
                }
            }
        }
        if (counter==colliderIndices.size()){
            break; //breaks after all colliders have been located to prevent unnecessary iterations
        }
    }
    // writelog("eps (6) = "+std::to_string(eps)+"\n");
    // writelog("k_LMR (2) = " + std::to_string(k_LMR_) + "\n"); 

    if (rateObj_M.which()==0){ // 0 means PlogRate
        //writelog("eps_M (2) = "+std::to_string(eps_M)+"\n");
        k_LMR_ += evalPlogRate(shared_data,dataObj_M,rateObj_M,eps_M)*eps_M*sigmaX_M/eps_mix;
        //writelog("eps_M (3) = "+std::to_string(eps_M)+"\n");
        // writelog("k_M = " + std::to_string(rate.evalFromStruct(data)) + "\n"); 
    }
    else if (rateObj_M.which()==1){ // 1 means TroeRate 
        k_LMR_ += evalTroeRate(shared_data,dataObj_M,rateObj_M)*eps_M*sigmaX_M/eps_mix;
    }
    else if (rateObj_M.which()==2){ // 2 means ChebyshevRate
        k_LMR_ += evalChebyshevRate(shared_data,dataObj_M,rateObj_M,eps_M)*eps_M*sigmaX_M/eps_mix;
    }

    // writelog("k_LMR (3) = " + std::to_string(k_LMR_) + "\n\n\n"); 
    // writelog("k_LMR = "+std::to_string(k_LMR_)+"\n");
    // writelog("T = "+std::to_string(shared_data.temperature)+"\n");
    // writelog("\n\n\nNEW REACTION\n");
    //writelog("\n\n\n");
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

