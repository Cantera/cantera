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
    colliderIndices.reserve(colliderNames.size());
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

void LmrRate::updatePlogData(const LmrData& shared_data, PlogData& dataObj, double eig0){
    dataObj.logP = shared_data.logP+log(eig0_mix)-log(eig0); //replaces logP with log of the effective pressure w.r.t. eig0
    dataObj.logT = shared_data.logT;
    dataObj.pressure = shared_data.pressure;
    dataObj.recipT = shared_data.recipT;
    dataObj.temperature = shared_data.temperature;
}

void LmrRate::updateTroeData(const LmrData& shared_data, FalloffData& dataObj){
    dataObj.conc_3b = shared_data.moleFractions;
    dataObj.logT = shared_data.logT;
    // dataObj.molar_density = shared_data.pressure; //
    dataObj.ready = shared_data.ready;
    dataObj.recipT = shared_data.recipT;
    dataObj.temperature = shared_data.temperature;
}

void LmrRate::updateChebyshevData(const LmrData& shared_data, ChebyshevData& dataObj){
    dataObj.log10P=log10(exp(shared_data.logP));
    dataObj.pressure=shared_data.pressure;
    dataObj.recipT=shared_data.recipT;
    dataObj.temperature=shared_data.temperature;
}

double LmrRate::evalFromStruct(const LmrData& shared_data){
    // writelog("evalFromStruct::1"); writelog("\n");
    // logP_=shared_data.logP;
    logT_=shared_data.logT;
    // pressure_=shared_data.pressure;
    recipT_=shared_data.recipT;
    // temperature_=shared_data.temperature;
    // ready_=shared_data.ready;
    moleFractions_=shared_data.moleFractions;
    double sigmaX_M=0.0;
    for (size_t i=0; i<nSpecies; i++){ //testing each species listed at the top of yaml file
        sigmaX_M += moleFractions_[i]; //total sum will be essentially 1, but perhaps not exactly due to Cantera's rounding conventions
    }
    writelog("sigmaX_M (1) = " + std::to_string(sigmaX_M) + "\n"); 
    eig0_mix=0.0;
    for (size_t i=0; i<nSpecies; i++){ //testing each species listed at the top of yaml file
        // writelog("evalFromStruct::2"); writelog("\n");
        for (size_t j=0; j<colliderIndices.size(); j++){
            // writelog("evalFromStruct::3"); writelog("\n");
            if (i==colliderIndices[j]){ // Species is in collider list
                // writelog("evalFromStruct::4"); writelog("\n");
                eig0_mix += moleFractions_[i]*eigObjs[j].evalRate(logT_, recipT_);
                sigmaX_M -= moleFractions_[i];
            }
        }
    }
    writelog("sigmaX_M (2) = " + std::to_string(sigmaX_M) + "\n"); 
    double eig0_M = eigObj_M.evalRate(logT_, recipT_);
    // writelog("eig0_M = " + std::to_string(eig0_M) + "\n"); 
    // double eig0_M = 10;
    eig0_mix += sigmaX_M*eig0_M; // add all M colliders to eig0_mix in a single step
    
    // writelog("eig0_mix = " + std::to_string(eig0_mix) + "\n"); 
    if (eig0_mix==0){
        throw InputFileError("LmrRate::evalFromStruct", m_input,"eig0_mix==0 for some reason");
    }
    k_LMR_=0;
    for (size_t i=0; i<nSpecies; i++){ //testing each species listed at the top of yaml file
        // writelog("evalFromStruct::7"); writelog("\n");
        for (size_t j=0; j<colliderIndices.size(); j++){
            // writelog("evalFromStruct::8"); writelog("\n");
            if (i==colliderIndices[j]){ // Species is in collider list
                double eig0 = eigObjs[j].evalRate(logT_, recipT_);
                // double eig0 = 1;
                if (eig0==0){
                    throw InputFileError("LmrRate::evalFromStruct", m_input,"eig0 is 0 for some reason");
                }
                // // writelog("evalFromStruct::9"); writeMsg(" eig0 = ",eig0*1e33);
                // writelog("evalFromStruct::9"); writeMsg(" eig0 = ",eig0);
                if (rateObjs[j].which()==0){ // 0 means PlogRate     
                    PlogData& data = boost::get<PlogData>(dataObjs[j]);
                    writelog("data.logP (1) = " + std::to_string(data.logP) + "\n"); 
                    writelog("data.logT (1) = " + std::to_string(data.logT) + "\n"); 
                    writelog("data.pressure (1) = " + std::to_string(data.pressure) + "\n"); 
                    writelog("data.temperature (1) = " + std::to_string(data.temperature) + "\n"); 
                    PlogRate& rate = boost::get<PlogRate>(rateObjs[j]);
                    updatePlogData(shared_data,data,eig0);
                    writelog("data.logP (2) = " + std::to_string(data.logP) + "\n"); 
                    writelog("data.logT (2) = " + std::to_string(data.logT) + "\n"); 
                    writelog("data.pressure (2) = " + std::to_string(data.pressure) + "\n"); 
                    writelog("data.temperature (2) = " + std::to_string(data.temperature) + "\n"); 
                    writelog("logP (2) = " + std::to_string(shared_data.logP) + "\n"); 
                    writelog("logT (2) = " + std::to_string(shared_data.logT) + "\n"); 
                    writelog("pressure (2) = " + std::to_string(shared_data.pressure) + "\n"); 
                    writelog("temperature (2) = " + std::to_string(shared_data.temperature) + "\n"); 
                    rate.updateFromStruct(data);
                    k_LMR_ += rate.evalFromStruct(data)*eig0*moleFractions_[i]/eig0_mix;
                    writelog(" eig0 = " + std::to_string(eig0) + "\n"); 
                    writelog(" eig0_mix = " + std::to_string(eig0_mix) + "\n"); 
                    writelog(" molefractions_[i] = " + std::to_string(moleFractions_[i]) + "\n"); 
                    // writelog("evalFromStruct::10"); writelog(" k_i_plog = " + std::to_string(evalPlogRate(rate,data,colliderNodes[j])) + "\n"); 
                    // writelog("evalFromStruct::10"); writelog(" k_i_plog = " + std::to_string(rate.evalFromStruct(data)) + "\n"); 
                }
                else if (rateObjs[j].which()==1){ // 1 means TroeRate  
                    FalloffData& data = boost::get<FalloffData>(dataObjs[j]);
                    TroeRate& rate = boost::get<TroeRate>(rateObjs[j]);
                    updateTroeData(shared_data,data);
                    k_LMR_ += rate.evalFromStruct(data)*eig0*moleFractions_[i]/eig0_mix;
                    // k_LMR_ += evalTroeRate(rate,data)*eig0*moleFractions_[i]/eig0_mix;
                    // writelog("evalFromStruct::11"); writelog(" k_i_troe = " + std::to_string(evalTroeRate(rate,data,colliderNodes[j])) + "\n");
                }
                else if (rateObjs[j].which()==2){ // 2 means ChebyshevRate  
                    ChebyshevData& data = boost::get<ChebyshevData>(dataObjs[j]);
                    ChebyshevRate& rate = boost::get<ChebyshevRate>(rateObjs[j]);
                    updateChebyshevData(shared_data,data);
                    rate.updateFromStruct(data);
                    k_LMR_ += rate.evalFromStruct(data)*eig0*moleFractions_[i]/eig0_mix;
                    // k_LMR_ += evalChebyshevRate(rate,data)*eig0*moleFractions_[i]/eig0_mix;
                    // writelog("evalFromStruct::12"); writelog(" k_i_cheb = " + std::to_string(evalChebyshevRate(rate,data,colliderNodes[j])) + "\n");
                }
                else {
                    throw InputFileError("LmrRate::evalFromStruct", m_input,"Something went wrong...");
                }
            }
        }
    }
    // writelog("evalFromStruct::14"); writelog("\n");
    if (rateObj_M.which()==0){ // 0 means PlogRate
        PlogData& data = boost::get<PlogData>(dataObj_M);
        PlogRate& rate = boost::get<PlogRate>(rateObj_M);
        updatePlogData(shared_data,data,eig0_M);
        rate.updateFromStruct(data);
        k_LMR_ += rate.evalFromStruct(data)*eig0_M*sigmaX_M/eig0_mix;
        // writelog("evalFromStruct::15"); writelog(" k_M_plog = " + std::to_string(evalPlogRate(boost::get<PlogRate>(rateObj_M),boost::get<PlogData>(dataObj_M),node_M)*1e13) + "\n"); 
        // writelog("evalFromStruct::15"); writelog(" k_M_plog = " + std::to_string(rate.evalFromStruct(data)) + "\n");
    }
    else if (rateObj_M.which()==1){ // 1 means TroeRate  
        // writelog("evalFromStruct::16"); writelog(" k_M_troe = " + std::to_string(evalTroeRate(boost::get<TroeRate>(rateObj_M),boost::get<FalloffData>(dataObj_M),node_M)*1e13) + "\n"); 
        // writelog("evalFromStruct::16"); writelog(" k_M_troe = " + std::to_string(evalTroeRate(boost::get<TroeRate>(rateObj_M),boost::get<FalloffData>(dataObj_M),node_M)) + "\n"); 
        FalloffData& data = boost::get<FalloffData>(dataObj_M);
        TroeRate& rate = boost::get<TroeRate>(rateObj_M);
        updateTroeData(shared_data,data);
        k_LMR_ += rate.evalFromStruct(data)*eig0_M*sigmaX_M/eig0_mix;
    }
    else if (rateObj_M.which()==2){ // 2 means ChebyshevRate
        ChebyshevData& data = boost::get<ChebyshevData>(dataObj_M);
        ChebyshevRate& rate = boost::get<ChebyshevRate>(rateObj_M);
        // writelog("evalFromStruct::17"); writelog(" k_M_cheb = " + std::to_string(evalChebyshevRate(boost::get<ChebyshevRate>(rateObj_M),boost::get<ChebyshevData>(dataObj_M),node_M)*1e13) + "\n"); 
        // writelog("evalFromStruct::17"); writelog(" k_M_cheb = " + std::to_string(evalChebyshevRate(boost::get<ChebyshevRate>(rateObj_M),boost::get<ChebyshevData>(dataObj_M),node_M)) + "\n"); 
        updateChebyshevData(shared_data,data);
        rate.updateFromStruct(data);
        k_LMR_ += rate.evalFromStruct(data)*eig0_M*sigmaX_M/eig0_mix;
        // k_LMR_ += evalChebyshevRate(boost::get<ChebyshevRate>(rateObj_M),boost::get<ChebyshevData>(dataObj_M))*eig0_M*sigmaX_M/eig0_mix;
    }
    // writelog("evalFromStruct::18");writelog("\n T = "+std::to_string(shared_data.temperature)+"\n");writelog("\n P = "+std::to_string(shared_data.pressure)+"\n"); writelog("\n\n\nNEW REACTION\n");
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

