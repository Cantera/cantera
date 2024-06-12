//! @file LmrRate.cpp
// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/LmrRate.h"
#include "cantera/thermo/ThermoPhase.h"
// #include "cantera/kinetics/Falloff.h"
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
    ReactionRate::setParameters(node, rate_units);
    if(!node.hasKey("collider-list")){
        throw InputFileError("LmrRate::setParameters", m_input,"Yaml input for LMR-R does not follow the necessary structure.");
    }
    
    auto& colliders = node["collider-list"].asVector<AnyMap>();        
    for (size_t i = 0; i < colliders.size(); i++){
        if (!colliders[i].hasKey("collider")) {
            throw InputFileError("LmrRate::setParameters", m_input,"Yaml input for LMR-R does not follow the necessary structure.");
        } else if (!colliders[i].hasKey("eig0")) {
            throw InputFileError("LmrRate::setParameters", m_input,"An eig0 value must be provided for all explicitly declared colliders in LMRR yaml entry.");
        }
        colliderInfo.insert({colliders[i]["collider"].as<std::string>(), colliders[i]});
    }

    ArrheniusRate eig0_M;
    RateTypes rateObj_M;
    DataTypes dataObj_M;
    AnyMap node_M;

    auto it1 = colliderInfo.find("M");
    if (it1 != colliderInfo.end() && it1->second.hasKey("rate-constants")){ 
        node_M=it1->second;
        eig0_M=ArrheniusRate(AnyValue(it1->second["eig0"]), it1->second.units(), rate_units);
        rateObj_M = PlogRate();
        dataObj_M = PlogData();
    } 
    else if(it1 != colliderInfo.end() && it1->second.hasKey("Troe")){ 
        node_M=it1->second;
        eig0_M=ArrheniusRate(AnyValue(it1->second["eig0"]), it1->second.units(), rate_units);
        rateObj_M = TroeRate();
        dataObj_M = FalloffData();
    }
    else if(it1 != colliderInfo.end() && it1->second.hasKey("pressure-range")){
        node_M=it1->second;
        eig0_M=ArrheniusRate(AnyValue(it1->second["eig0"]), it1->second.units(), rate_units);
        rateObj_M = ChebyshevRate();
        dataObj_M = ChebyshevData();
    }
    else {
        throw InputFileError("LmrRate::setParameters", m_input,"Not enough data provided for species 'M'.");
    }
    
    vector<string> allSpecies=LmrData().allSpecies;

    for (size_t i=0; i<allSpecies.size(); i++){ //testing each species listed at the top of yaml file
        auto it2 = colliderInfo.find(allSpecies[i]);
        if (it2 != colliderInfo.end() && it2->second.hasKey("rate-constants")){ 
            node_M=it1->second;
            eigObjs.push_back(ArrheniusRate(AnyValue(it2->second["eig0"]), it2->second.units(), rate_units));
            rateObjs.push_back(PlogRate());
            dataObjs.push_back(PlogData());
        } 
        else if(it2 != colliderInfo.end() && it2->second.hasKey("Troe")){ 
            colliderNodes.push_back(it2->second);
            eigObjs.push_back(ArrheniusRate(AnyValue(it2->second["eig0"]), it2->second.units(), rate_units));
            rateObjs.push_back(TroeRate());
            dataObjs.push_back(FalloffData());
        }
        else if(it2 != colliderInfo.end() && it2->second.hasKey("pressure-range")){ 
            colliderNodes.push_back(it2->second);
            eigObjs.push_back(ArrheniusRate(AnyValue(it2->second["eig0"]), it2->second.units(), rate_units));
            rateObjs.push_back(ChebyshevRate());
            dataObjs.push_back(ChebyshevData());   
        }
        else if(it2 != colliderInfo.end() && !(it2->second.hasKey("pressure-range")) && !(it2->second.hasKey("Troe")) && !(it2->second.hasKey("rate-constants"))){ //yaml species has an eig0 but no additional LMRR data, so treat its rate as same as "M"
            colliderNodes.push_back(node_M);
            eigObjs.push_back(ArrheniusRate(AnyValue(it2->second["eig0"]), it2->second.units(), rate_units));
            rateObjs.push_back(rateObj_M);
            dataObjs.push_back(dataObj_M);
        }
        else if(it2 == colliderInfo.end()){ //yaml species has no LMRR data, so treat its rate and eig0 as same as "M"
            colliderNodes.push_back(node_M);
            eigObjs.push_back(eig0_M);
            rateObjs.push_back(rateObj_M);
            dataObjs.push_back(dataObj_M);
        }
        else{
            throw InputFileError("LmrRate::setParameters", m_input,"LMRR reaction has invalid yaml input.");
        }
    }
}

void LmrRate::validate(const string& equation, const Kinetics& kin){}

double LmrRate::evalFromStruct(const LmrData& shared_data){
    logP_=shared_data.logP;
    logT_=shared_data.logT;
    pressure_=shared_data.pressure;
    recipT_=shared_data.recipT;
    temperature_=shared_data.temperature;
    ready_=shared_data.ready;
    moleFractions_=shared_data.moleFractions;

    double eig0_mix=0;
    for (size_t i=0; i<shared_data.allSpecies.size(); i++){ //testing each species listed at the top of yaml file
        eig0_mix += moleFractions_[i]*eigObjs[i].evalRate(logT_, recipT_);
    }

    k_LMR_=0;
    for (size_t i=0; i<shared_data.allSpecies.size(); i++){ //testing each species listed at the top of yaml file
        double Xi=moleFractions_[i];
        double eig0 = eigObjs[i].evalRate(logT_, recipT_);
        if (rateObjs.at(i).which()==0){ // 0 means PlogRate   
            logPeff_= logP_+log(eig0_mix)-log(eig0); //need to update this every time before evalPlogRate
            PlogData& data = boost::get<PlogData>(dataObjs.at(i));
            PlogRate& rate = boost::get<PlogRate>(rateObjs.at(i));
            data.logP = logP_; //replaces logP with log of the effective pressure w.r.t. eig0_M
            data.logT = logT_;
            data.pressure = pressure_;
            data.recipT = recipT_;
            data.temperature = temperature_;
            rate.setParameters(colliderNodes[i],rate_units_); //it->second refers to the yaml data for the ith collider
            rate.updateFromStruct(data);
            double k_i = rate.evalFromStruct(data);
            double X_i = eig0*Xi/eig0_mix;
            k_LMR_ += k_i*X_i;
            writelog("species = ",shared_data.allSpecies[i]); 
            writeMsg("k_i_plog = ",k_i);
        }
        else if (rateObjs.at(i).which()==1){ // 1 means TroeRate   
            FalloffData& data = boost::get<FalloffData>(dataObjs.at(i));
            TroeRate& rate = boost::get<TroeRate>(rateObjs.at(i));
            data.conc_3b = moleFractions_;
            data.logT = logT_;
            // data.molar_density = pressure_; //
            data.ready = ready_;
            data.recipT = recipT_;
            data.temperature = temperature_;
            rate.setParameters(colliderNodes[i],rate_units_); //it->second refers to the yaml data for the ith collider
            double k_i = rate.evalFromStruct(data);
            double X_i = eig0*Xi/eig0_mix;
            k_LMR_ += k_i*X_i;
            writelog("species = ",shared_data.allSpecies[i]); 
            writeMsg("k_i_troe = ",k_i);
        }
        else if (rateObjs.at(i).which()==2){ // 2 means ChebyshevRate  
            ChebyshevData& data = boost::get<ChebyshevData>(dataObjs.at(i));
            ChebyshevRate& rate = boost::get<ChebyshevRate>(rateObjs.at(i));
            data.log10P=logP_; //THIS IS INCORRECT. logP_ is natural log, not base 10!!!
            data.logT=logT_;
            data.pressure=pressure_;
            data.recipT=recipT_;
            data.temperature=temperature_;
            rate.setParameters(colliderNodes[i],rate_units_); //it->second refers to the yaml data for the ith collider
            rate.updateFromStruct(data);
            double k_i = rate.evalFromStruct(data);
            double X_i = eig0*Xi/eig0_mix;
            k_LMR_ += k_i*X_i;
            writelog("species = ",shared_data.allSpecies[i]); 
            writeMsg("k_i_cheb = ",k_i);
        }
        else {
            throw InputFileError("LmrRate::evalFromStruct", m_input,"Something went wrong...");
        }
    }
    writeMsg("k_LMR = ",k_LMR_);
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

