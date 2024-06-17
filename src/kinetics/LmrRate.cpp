//! @file LmrRate.cpp
// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/LmrRate.h"
#include "cantera/thermo/ThermoPhase.h"
// #include "cantera/kinetics/Falloff.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/kinetics/Reaction.h"
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
    std::string rxn = node["equation"].as<std::string>();
    // writelog("setParameters::1"); writelog("\n");
    ReactionRate::setParameters(node, rate_units);
    rate_units_=rate_units;
    // for (int i = 0; i < rxn.length(); i++) {
    //     if (i==0 && rxn[i]=='2' && rxn[i+1]==' ') {
    //         rate_units_.join(1);
    //         break;
    //     } else if (rxn[i]=='<'){
    //         break;
    //     }
    // } 

    if(!node.hasKey("collider-list")){
        throw InputFileError("LmrRate::setParameters", m_input,"Yaml input for LMR-R does not follow the necessary structure.");
    }
    auto& colliders = node["collider-list"].asVector<AnyMap>();
    if (colliders[0]["collider"].as<std::string>() != "M"){
        throw InputFileError("LmrRate::setParameters", m_input,"The first species defined in yaml input must be 'M'.");
    }
    // node_M=colliders[0];
    eigObj_M=ArrheniusRate(AnyValue(colliders[0]["eig0"]),colliders[0].units(), rate_units);
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
        // eigObjs.push_back(ArrheniusRate(colliders[i]["eig0"], colliders[i].units(), rate_units));
        eigObjs.push_back(ArrheniusRate(AnyValue(colliders[i]["eig0"]),colliders[i].units(), rate_units));
        if (colliders[i].hasKey("rate-constants")){
            // writelog("setParameters::5"); writelog("\n");
            rateObjs.push_back(PlogRate(colliders[i], rate_units));
            dataObjs.push_back(PlogData());
            // colliderNodes.push_back(colliders[i]);
        } else if (colliders[i].hasKey("Troe")){
            // writelog("setParameters::6"); writelog("\n");
            rateObjs.push_back(TroeRate(colliders[i], rate_units));
            dataObjs.push_back(FalloffData());   
            // colliderNodes.push_back(colliders[i]);
        } else if (colliders[i].hasKey("pressure-range")){
            // writelog("setParameters::7"); writelog("\n");
            rateObjs.push_back(ChebyshevRate(colliders[i], rate_units));
            dataObjs.push_back(ChebyshevData());
        } else { //Collider has an eig0 specified, but no other info is provided. Assign it the same rate and data objects as "M"
            // writelog("setParameters::8"); writelog("\n");
            rateObjs.push_back(rateObj_M);
            dataObjs.push_back(dataObj_M);
            // colliderNodes.push_back(node_M);
        }
    }
}

void LmrRate::validate(const string& equation, const Kinetics& kin){
    // writelog("validate::1"); writelog("\n");
}

void LmrRate::setContext(const Reaction& rxn, const Kinetics& kin){   
    // m_stoichCoeffs.clear();
    // for (const auto& [name, stoich] : rxn.reactants) {
    //     m_stoichCoeffs.emplace_back(kin.kineticsSpeciesIndex(name), -stoich);
    // }
    // for (const auto& [name, stoich] : rxn.products) {
    //     m_stoichCoeffs.emplace_back(kin.kineticsSpeciesIndex(name), stoich);
    // }


    // writelog("setContext::1"); writelog("\n");
    // for (size_t i=0; i<colliderNames.size();i++){ //Starts at 1, because colliderNames[0] == "M" //THIS LOOP IS NOT WORKING
    //     writelog("setContext::2"); writelog("\n");
    //     for (const auto& [name, stoich] : rxn.reactants){
    //         writelog("setContext::3"); writelog("\n");
    //         if (name == colliderNames[i]){
    //             writelog("setContext::4"); writelog("\n");
    //             colliderIndices.push_back(kin.kineticsSpeciesIndex(name));
    //         }
    //     }
    // }

    for (size_t i=0; i<colliderNames.size();i++){ //Starts at 1, because names[0] == "M" //THIS LOOP IS NOT WORKING
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
    if (colliderIndices.empty()){
        throw InputFileError("LmrRate::setContext", m_input,"collidersIndices is empty for some reason");
    }
}

double LmrRate::evalPlogRate(PlogRate& rate, PlogData& data){
    // writelog("evalPlogRate::1"); writelog("\n");
    data.logP = logP_; //replaces logP with log of the effective pressure w.r.t. eig0_M
    data.logT = logT_;
    data.pressure = pressure_;
    data.recipT = recipT_;
    data.temperature = temperature_;
    // rate.setParameters(node,rate_units_); //it->second refers to the yaml data for the ith collider
    rate.updateFromStruct(data);
    return rate.evalFromStruct(data);
}
double LmrRate::evalTroeRate(TroeRate& rate, FalloffData& data){
    // writelog("evalTroeRate::1"); writelog("\n");
    data.conc_3b = moleFractions_;
    data.logT = logT_;
    // data.molar_density = pressure_; //
    data.ready = ready_;
    data.recipT = recipT_;
    data.temperature = temperature_;
    // rate.setParameters(node,rate_units_); //it->second refers to the yaml data for the ith collider
    return rate.evalFromStruct(data);
}

double LmrRate::evalChebyshevRate(ChebyshevRate& rate, ChebyshevData& data){
    // writelog("evalChebyshevRate::1"); writelog("\n");
    data.log10P=logP_; //THIS IS INCORRECT. logP_ is natural log, not base 10!!!
    data.logT=logT_;
    data.pressure=pressure_;
    data.recipT=recipT_;
    data.temperature=temperature_;
    // rate.setParameters(node,rate_units_); //it->second refers to the yaml data for the ith collider
    rate.updateFromStruct(data);
    return rate.evalFromStruct(data);
}

double LmrRate::evalFromStruct(const LmrData& shared_data){
    // writelog("evalFromStruct::1"); writelog("\n");
    logP_=shared_data.logP;
    logT_=shared_data.logT;
    pressure_=shared_data.pressure;
    recipT_=shared_data.recipT;
    temperature_=shared_data.temperature;
    ready_=shared_data.ready;
    moleFractions_=shared_data.moleFractions;

    double eig0_mix=0;
    double sigmaX_M=0.0;
    for (size_t i=0; i<nSpecies; i++){ //testing each species listed at the top of yaml file
        // writelog("evalFromStruct::2"); writelog("\n");
        for (size_t j=0; j<colliderIndices.size(); j++){
            // writelog("evalFromStruct::3"); writelog("\n");
            if (i==colliderIndices[j]){ // Species is in collider list
                // writelog("evalFromStruct::4"); writelog("\n");
                eig0_mix += moleFractions_[i]*eigObjs[j].evalRate(logT_, recipT_);
            } else { // Species not in collider list so treat as "M"
                // writelog("evalFromStruct::5"); writelog("\n");
                sigmaX_M += moleFractions_[i];
            }
        }
    }
    double eig0_M = eigObj_M.evalRate(logT_, recipT_);
    // writelog(" eig0_M = " + std::to_string(eig0_M*1e33) + "\n"); 
    writelog(" eig0_M = " + std::to_string(eig0_M) + "\n"); 
    // double eig0_M = 10;
    eig0_mix += sigmaX_M*eig0_M; // add all M colliders to eig0_mix in a single step
    // writelog("evalFromStruct::6"); writelog(" eig0_mix = " + std::to_string(eig0_mix) + "\n"); 
    if (eig0_mix==0){
        throw InputFileError("LmrRate::evalFromStruct", m_input,"eig0_mix==0 for some reason");
    }
    k_LMR_=0;
    sigmaX_M=0.0;
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
                    logP_= shared_data.logP+log(eig0_mix)-log(eig0); //replaces logP with log of the effective pressure w.r.t. eig0
                    // PlogData& data = boost::get<PlogData>(dataObjs.at(j));
                    // PlogRate& rate = boost::get<PlogRate>(rateObjs.at(j));
                    PlogData& data = boost::get<PlogData>(dataObjs[j]);
                    PlogRate& rate = boost::get<PlogRate>(rateObjs[j]);
                    // k_LMR_ += evalPlogRate(rate,data,colliderNodes[j])*eig0*moleFractions_[i]/eig0_mix;
                    
                    k_LMR_ += evalPlogRate(rate,data)*eig0*moleFractions_[i]/eig0_mix;
                    // evalPlogRate(boost::get<PlogRate>(rateObj_M),boost::get<PlogData>(dataObj_M),node_M)*eig0_M*sigmaX_M/eig0_mix;
                    logP_ = shared_data.logP; //return to the "normal" logP value to avoid messing up other calcs
                    // writelog("evalFromStruct::10"); writelog(" k_i_plog = " + std::to_string(evalPlogRate(rate,data,colliderNodes[j])*1e13) + "\n"); 
                    writelog("evalFromStruct::10"); writelog(" k_i_plog = " + std::to_string(evalPlogRate(rate,data)) + "\n"); 
                }
                else if (rateObjs[j].which()==1){ // 1 means TroeRate  
                    FalloffData& data = boost::get<FalloffData>(dataObjs[j]);
                    TroeRate& rate = boost::get<TroeRate>(rateObjs[j]);
                    k_LMR_ += evalTroeRate(rate,data)*eig0*moleFractions_[i]/eig0_mix;
                    // writelog("evalFromStruct::11"); writelog(" k_i_troe = " + std::to_string(evalTroeRate(rate,data,colliderNodes[j])*1e13) + "\n"); 
                    // writelog("evalFromStruct::11"); writelog(" k_i_troe = " + std::to_string(evalTroeRate(rate,data,colliderNodes[j])) + "\n"); 
                }
                else if (rateObjs[j].which()==2){ // 2 means ChebyshevRate  
                    ChebyshevData& data = boost::get<ChebyshevData>(dataObjs[j]);
                    ChebyshevRate& rate = boost::get<ChebyshevRate>(rateObjs[j]);
                    k_LMR_ += evalChebyshevRate(rate,data)*eig0*moleFractions_[i]/eig0_mix;
                    // writelog("evalFromStruct::12"); writelog(" k_i_cheb = " + std::to_string(evalChebyshevRate(rate,data,colliderNodes[j])*1e13) + "\n"); 
                    // writelog("evalFromStruct::12"); writelog(" k_i_cheb = " + std::to_string(evalChebyshevRate(rate,data,colliderNodes[j])) + "\n"); 
                }
                else {
                    throw InputFileError("LmrRate::evalFromStruct", m_input,"Something went wrong...");
                }
            }
            else {
                // writelog("evalFromStruct::13"); writelog("\n");
                sigmaX_M += moleFractions_[i];
            }
        }
    }
    // writelog("evalFromStruct::14"); writelog("\n");
    if (rateObj_M.which()==0){ // 0 means PlogRate
        logP_= shared_data.logP+log(eig0_mix)-log(eig0_M); //replaces logP with log of the effective pressure w.r.t. eig0_M
        k_LMR_ += evalPlogRate(boost::get<PlogRate>(rateObj_M),boost::get<PlogData>(dataObj_M))*eig0_M*sigmaX_M/eig0_mix;
        logP_ = shared_data.logP; //return to the "normal" logP value to avoid messing up other calcs
        // writelog("evalFromStruct::15"); writelog(" k_M_plog = " + std::to_string(evalPlogRate(boost::get<PlogRate>(rateObj_M),boost::get<PlogData>(dataObj_M),node_M)*1e13) + "\n"); 
        writelog("evalFromStruct::15"); writelog(" k_M_plog = " + std::to_string(evalPlogRate(boost::get<PlogRate>(rateObj_M),boost::get<PlogData>(dataObj_M))) + "\n"); 
    }
    else if (rateObj_M.which()==1){ // 1 means TroeRate  
        // writelog("evalFromStruct::16"); writelog(" k_M_troe = " + std::to_string(evalTroeRate(boost::get<TroeRate>(rateObj_M),boost::get<FalloffData>(dataObj_M),node_M)*1e13) + "\n"); 
        // writelog("evalFromStruct::16"); writelog(" k_M_troe = " + std::to_string(evalTroeRate(boost::get<TroeRate>(rateObj_M),boost::get<FalloffData>(dataObj_M),node_M)) + "\n"); 
        k_LMR_ += evalTroeRate(boost::get<TroeRate>(rateObj_M),boost::get<FalloffData>(dataObj_M))*eig0_M*sigmaX_M/eig0_mix;
    }
    else if (rateObj_M.which()==2){ // 2 means ChebyshevRate
        // writelog("evalFromStruct::17"); writelog(" k_M_cheb = " + std::to_string(evalChebyshevRate(boost::get<ChebyshevRate>(rateObj_M),boost::get<ChebyshevData>(dataObj_M),node_M)*1e13) + "\n"); 
        // writelog("evalFromStruct::17"); writelog(" k_M_cheb = " + std::to_string(evalChebyshevRate(boost::get<ChebyshevRate>(rateObj_M),boost::get<ChebyshevData>(dataObj_M),node_M)) + "\n"); 
        k_LMR_ += evalChebyshevRate(boost::get<ChebyshevRate>(rateObj_M),boost::get<ChebyshevData>(dataObj_M))*eig0_M*sigmaX_M/eig0_mix;
    }
    // writelog("evalFromStruct::18");writelog("\n T = "+std::to_string(shared_data.temperature)+"\n");writelog("\n P = "+std::to_string(shared_data.pressure)+"\n"); writelog("\n\n\nNEW REACTION\n");
    // writemsg("k_LMR = ",k_LMR_);
    // writelog("NEW REACTION\n");
    // writelog("T = "+std::to_string(shared_data.temperature)+"\n");
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

