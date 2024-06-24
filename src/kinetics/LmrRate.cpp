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
    ReactionRate::setParameters(node, rate_units);
    rate_units_ = rate_units;
    UnitStack eig0_units{{Units(1.0), 1.0}};

    // writelog("setParameters::1"); writelog("\n");
    
    if(!node.hasKey("collider-list")){
        throw InputFileError("LmrRate::setParameters", m_input,"Yaml input for LMR-R does not follow the necessary structure.");
    }
    const auto& colliders = node["collider-list"].asVector<AnyMap>();
    if (colliders[0]["collider"].as<std::string>() != "M"){
        throw InputFileError("LmrRate::setParameters", m_input,"The first species defined in yaml input must be 'M'.");
    }
    // rate_units_.join(1);
    // setRateUnits(rate_units_);
    eigObj_M=ArrheniusRate(AnyValue(colliders[0]["eig0"]),colliders[0].units(),eig0_units);
    // eigObj_M=ArrheniusRate(AnyValue(colliders[0]["eig0"]),Units(1.0));
    if (colliders[0].hasKey("rate-constants")){
        // writelog("setParameters::2"); writelog("\n");
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
        if (!colliders[i].hasKey("collider")) {
            throw InputFileError("LmrRate::setParameters", m_input,"Yaml input for LMR-R does not follow the necessary structure.");
        } else if (!colliders[i].hasKey("eig0")) {
            throw InputFileError("LmrRate::setParameters", m_input,"An eig0 value must be provided for all explicitly declared colliders in LMRR yaml entry.");
        }
        colliderInfo.insert({colliders[i]["collider"].as<std::string>(), colliders[i]}); //Legacy parameter, used b.c. getParameters would have to be rewritten otherwise
        colliderNames.push_back(colliders[i]["collider"].as<std::string>());
        eigObjs.push_back(ArrheniusRate(AnyValue(colliders[i]["eig0"]),colliders[i].units(),eig0_units));
        if (colliders[i].hasKey("rate-constants")){
            // writelog("setParameters::5"); writelog("\n");
            rateObjs.push_back(PlogRate(colliders[i], rate_units));
            dataObjs.push_back(PlogData());
        } else if (colliders[i].hasKey("Troe")){
            // writelog("setParameters::6"); writelog("\n");
            rateObjs.push_back(TroeRate(colliders[i], rate_units));
            dataObjs.push_back(FalloffData());
        } else if (colliders[i].hasKey("pressure-range")){
            // writelog("setParameters::7"); writelog("\n");
            rateObjs.push_back(ChebyshevRate(colliders[i], rate_units));
            dataObjs.push_back(ChebyshevData());
        } else { //Collider has an eig0 specified, but no other info is provided. Assign it the same rate and data objects as "M"
            // writelog("setParameters::8"); writelog("\n");
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

double LmrRate::evalPlogRate(const LmrData& shared_data, DataTypes& dataObj, RateTypes& rateObj, double& eig0){
    PlogData& data = boost::get<PlogData>(dataObj);
    PlogRate& rate = boost::get<PlogRate>(rateObj);
    data.logP = shared_data.logP+log(eig0_mix)-log(eig0); //replaces logP with log of the effective pressure w.r.t. eig0
    data.logT = shared_data.logT;
    // dataObj.pressure = shared_data.pressure; // THIS IS NOT CORRECT
    data.recipT = shared_data.recipT;
    // dataObj.temperature = shared_data.temperature;
    rate.updateFromStruct(data);
    return rate.evalFromStruct(data)*1;  
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

double LmrRate::evalChebyshevRate(const LmrData& shared_data, DataTypes& dataObj, RateTypes& rateObj, double& eig0){
    ChebyshevData& data = boost::get<ChebyshevData>(dataObj);
    ChebyshevRate& rate = boost::get<ChebyshevRate>(rateObj);
    data.log10P=log10(exp(shared_data.logP+log(eig0_mix)-log(eig0)));
    // data.pressure=shared_data.pressure;
    data.recipT=shared_data.recipT;
    // data.temperature=shared_data.temperature;
    rate.updateFromStruct(data);
    return rate.evalFromStruct(data);  
}

double LmrRate::evalFromStruct(const LmrData& shared_data){
    // writelog("logT (0) = " + std::to_string(shared_data.logT)+"\n");
    // logP_=shared_data.logP;
    // logT_=shared_data.logT;
    // pressure_=shared_data.pressure;
    // recipT_=shared_data.recipT;
    // temperature_=shared_data.temperature;
    // ready_=shared_data.ready;
    // vector<double> moleFractions_=shared_data.moleFractions;
    
    double sigmaX_M=0.0;
    for (size_t i=0; i<nSpecies; i++){ //testing each species listed at the top of yaml file
        sigmaX_M += shared_data.moleFractions[i]; //total sum will be essentially 1, but perhaps not exactly due to Cantera's rounding conventions
    }
    // writelog("sigmaX_M (0) = " + std::to_string(sigmaX_M) + "\n"); 
    
    eig0_mix=0.0;
    size_t counter=0;
    for (size_t i=0; i<nSpecies; i++){ //testing each species listed at the top of yaml file
        for (size_t j=0; j<colliderIndices.size(); j++){
            if (i==colliderIndices[j]){ // Species is in collider list
                eig0_mix += shared_data.moleFractions[i]*eigObjs[j].evalRate(shared_data.logT, shared_data.recipT);
                // writelog("eig[j] = " + std::to_string(eigObjs[j].evalRate(logT_, recipT_)) + " logT = " + std::to_string(logT_) + " recipT = " + std::to_string(recipT_) + "\n"); 
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

    // writelog("eig0_mix (0) = " + std::to_string(eig0_mix) + "\n"); 
    // writelog("sigmaX_M (2) = " + std::to_string(sigmaX_M) + "\n"); 
    double eig0_M = eigObj_M.evalRate(shared_data.logT, shared_data.recipT);
    // writelog("eig0_M = " + std::to_string(eig0_M) + "\n"); 
    eig0_mix += sigmaX_M*eig0_M; // add all M colliders to eig0_mix in a single step
    // writelog("eig0_mix (1) = " + std::to_string(eig0_mix) + "\n"); 
    if (eig0_mix==0){
        throw InputFileError("LmrRate::evalFromStruct", m_input,"eig0_mix==0 for some reason");
    }

    double k_LMR_=0.0;
    writelog("k_LMR (1) = " + std::to_string(k_LMR_) + "\n"); 
    counter=0;
    for (size_t i=0; i<nSpecies; i++){ //testing each species listed at the top of yaml file
        // writelog("i = " + std::to_string(i) + ", X_i = " + std::to_string(shared_data.moleFractions[i]) + "\n"); 
        for (size_t j=0; j<colliderIndices.size(); j++){
            // writelog("j = " + std::to_string(j) + "\n"); 
            if (i==colliderIndices[j]){ // Species is in collider list
                // writelog("i = colliderIndices[j] = " + std::to_string(colliderIndices[j]) + ", X_i = " + std::to_string(shared_data.moleFractions[i]) + "\n");
                
                // writelog("eig0_i = " + std::to_string(eigObjs[j].evalRate(logT_, recipT_))+"\n");
                // if (eig0==0){
                //     throw InputFileError("LmrRate::evalFromStruct", m_input,"eig0 is 0 for some reason");
                // }
                double eig0 = eigObjs[j].evalRate(shared_data.logT, shared_data.recipT); 
                if (rateObjs[j].which()==0){ // 0 means PlogRate    
                    k_LMR_ += evalPlogRate(shared_data,dataObjs[j],rateObjs[j],eig0)*eig0*shared_data.moleFractions[i]/eig0_mix;
                    counter+=1;
                    break; //breaks after collider has been located to prevent unnecessary iterations
                }
                else if (rateObjs[j].which()==1){ // 1 means TroeRate  
                    double eig0 = eigObjs[j].evalRate(shared_data.logT, shared_data.recipT);
                    k_LMR_ += evalTroeRate(shared_data,dataObjs[j],rateObjs[j])*eig0*shared_data.moleFractions[i]/eig0_mix;
                }
                else if (rateObjs[j].which()==2){ // 2 means ChebyshevRate  
                    double eig0 = eigObjs[j].evalRate(shared_data.logT, shared_data.recipT);
                    k_LMR_ += evalChebyshevRate(shared_data,dataObjs[j],rateObjs[j],eig0)*eig0*shared_data.moleFractions[i]/eig0_mix;
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
    writelog("k_LMR (2) = " + std::to_string(k_LMR_) + "\n"); 

    if (rateObj_M.which()==0){ // 0 means PlogRate
        k_LMR_ += evalPlogRate(shared_data,dataObj_M,rateObj_M,eig0_M)*eig0_M*sigmaX_M/eig0_mix;
        // writelog("k_M = " + std::to_string(rate.evalFromStruct(data)) + "\n"); 
    }
    else if (rateObj_M.which()==1){ // 1 means TroeRate 
        k_LMR_ += evalTroeRate(shared_data,dataObj_M,rateObj_M)*eig0_M*sigmaX_M/eig0_mix;
    }
    else if (rateObj_M.which()==2){ // 2 means ChebyshevRate
        k_LMR_ += evalChebyshevRate(shared_data,dataObj_M,rateObj_M,eig0_M)*eig0_M*sigmaX_M/eig0_mix;
    }

    writelog("k_LMR (3) = " + std::to_string(k_LMR_) + "\n\n\n"); 
    // writelog("k_LMR = "+std::to_string(k_LMR_)+"\n");
    // writelog("T = "+std::to_string(shared_data.temperature)+"\n");
    // writelog("\n\n\nNEW REACTION\n");
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
            colliderNode["eig0"]=colliders_i["eig0"];
            colliderNode["rate-constants"]=colliders_i["rate-constants"];
        } 
        else if(colliders_i.hasKey("Troe")){
            colliderNode["collider"]=name;
            colliderNode["eig0"]=colliders_i["eig0"];
            colliderNode["low-P-rate-constant"]=colliders_i["low-P-rate-constant"];
            colliderNode["high-P-rate-constant"]=colliders_i["high-P-rate-constant"];
            colliderNode["Troe"]=colliders_i["Troe"];
        } 
        else if(colliders_i.hasKey("data")&&colliders_i.hasKey("pressure-range")&&colliders_i.hasKey("temperature-range")){ //"M" is of Chebyshev type
            colliderNode["collider"]=name;
            colliderNode["eig0"]=colliders_i["eig0"];
            colliderNode["temperature-range"]=colliders_i["temperature-range"];
            colliderNode["pressure-range"]=colliders_i["pressure-range"];
            colliderNode["data"]=colliders_i["data"];
        } 
        else {
            colliderNode["collider"]=name;
            colliderNode["eig0"]=colliders_i["eig0"];
        }
        topLevelList.push_back(std::move(colliderNode));
    }
    rateNode["collider-list"]=std::move(topLevelList);
}
}

