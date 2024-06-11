//! @file LmrRate.cpp
// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/LmrRate.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/kinetics/Falloff.h"
#include "cantera/kinetics/ChebyshevRate.h"
#include "cantera/kinetics/PlogRate.h"
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
    rate_units_=rate_units;
    if(node.hasKey("collider-list")){
        auto& colliders = node["collider-list"].asVector<AnyMap>();
        for (int i = 0; i < colliders.size(); i++){
            if (colliders[i].hasKey("collider") && colliders[i].hasKey("eig0")) {
                colliderInfo.insert({colliders[i]["collider"].as<std::string>(), colliders[i]});
            } else {
                throw InputFileError("LmrRate::setParameters", m_input,"An eig0 value must be provided for all explicitly declared colliders in LMRR yaml entry.");
            }
        }
    } else {
        throw InputFileError("LmrRate::setParameters", m_input,"Yaml input for LMR-R does not follow the necessary structure.");
    }
}

void LmrRate::validate(const string& equation, const Kinetics& kin){}

double LmrRate::evalPlogRate(map<string,AnyMap>::iterator it){
    PlogData plog_data;
    PlogRate plog_rate;
    // plog_data.logP = logPeff_; //replaces logP with log of the effective pressure w.r.t. eig0_M
    plog_data.logP = logP_; //replaces logP with log of the effective pressure w.r.t. eig0_M
    plog_data.logT = logT_;
    plog_data.pressure = pressure_;
    plog_data.recipT = recipT_;
    plog_data.temperature = temperature_;
    plog_rate.setParameters(it->second,rate_units_); //it->second refers to the yaml data for the ith collider
    plog_rate.updateFromStruct(plog_data);
    return plog_rate.evalFromStruct(plog_data);
}

double LmrRate::evalTroeRate(map<string,AnyMap>::iterator it){
    FalloffData troe_data;
    TroeRate troe_rate;
    troe_data.conc_3b = moleFractions_;
    troe_data.logT = logT_;
    // troe_data.molar_density = pressure_; //
    troe_data.ready = ready_;
    troe_data.recipT = recipT_;
    troe_data.temperature = temperature_;
    // colliders_i["type"]; 
    troe_rate.setParameters(it->second,rate_units_); //it->second refers to the yaml data for the ith collider
    return troe_rate.evalFromStruct(troe_data);
}

double LmrRate::evalChebyshevRate(map<string,AnyMap>::iterator it){
    ChebyshevData cheb_data;
    ChebyshevRate cheb_rate;
    cheb_data.log10P=logP_; //THIS IS INCORRECT. logP_ is natural log, not base 10!!!
    cheb_data.logT=logT_;
    cheb_data.pressure=pressure_;
    cheb_data.recipT=recipT_;
    cheb_data.temperature=temperature_;
    cheb_rate.setParameters(it->second,rate_units_); //it->second refers to the yaml data for the ith collider
    cheb_rate.updateFromStruct(cheb_data);
    return cheb_rate.evalFromStruct(cheb_data);
}

// double LmrRate::geteig0(map<string,AnyMap>::iterator it){
//     ArrheniusRate eig = ArrheniusRate(AnyValue(it->second["eig0"]), it->second.units(), rate_units_);
//     return eig.evalRate(logT_, recipT_);
// }

// vector<double> LmrRate::get_eig0M_kM(){
//     double kM;
//     double eig0M;
//     auto it = colliderInfo.find("M");
//     if (it != colliderInfo.end() && it->second.hasKey("rate-constants")){ 
//         eig0M=ArrheniusRate(AnyValue(it->second["eig0"]), it->second.units(), rate_units_).evalRate(logT_, recipT_);
//         kM = evalPlogRate(eig0M,it);
//         writeMsg("eig0_M = ",eig0M);
//         writeMsg("kM_plog = ",kM);
//     } 
//     else if(it != colliderInfo.end() && it->second.hasKey("Troe")){ 
//         eig0M=ArrheniusRate(AnyValue(it->second["eig0"]), it->second.units(), rate_units_).evalRate(logT_, recipT_);
//         kM = evalTroeRate(it);
//         writeMsg("eig0_M = ",eig0M);
//         writeMsg("kM_troe = ",kM);
//         // colliders_i["type"]="LMR_R"; //this dummy key needs to be defined because falloff.cpp requires it
//     }
//     else if(it != colliderInfo.end() && it->second.hasKey("pressure-range")){
//         eig0M=ArrheniusRate(AnyValue(it->second["eig0"]), it->second.units(), rate_units_).evalRate(logT_, recipT_);
//         kM = evalChebyshevRate(it);
//         writeMsg("eig0_M = ",eig0M);
//         writeMsg("kM_cheb = ",kM);
//     }
//     else {
//         throw InputFileError("LmrRate::getkM", m_input,"Not enough data provided for species 'M'.");
//     }
//     return {eig0M,kM};
// }

// double LmrRate::geteig0mix(const LmrData& shared_data){
//     double eig0mix=0.0;
//     for (size_t i=0; i<speciesList_.size(); i++){ //testing each species listed at the top of yaml file
//         double Xi=moleFractions_[i];
//         auto it = colliderInfo.find(speciesList_[i]);
//         if (it != colliderInfo.end()) { //yaml species has at least an eig0 value provided for LMRR 
//             eig0mix += Xi*geteig0(it);
//         }
//         else if (it == colliderInfo.end()) { //yaml species has no data provided for LMRR (treat as "M")
//             eig0mix += Xi*eig0_M;
//         }
//         else{
//             throw InputFileError("LmrRate::geteig0mix", m_input,"Cannot compute eig0mix due to invalid LMRR yaml input.");
//         }
//     }
//     writeMsg("eig0mix = ",eig0mix);
//     return eig0mix;
// }

double LmrRate::evalFromStruct(const LmrData& shared_data){
    k_LMR_=0;
    double eig0_mix=0;
    moleFractions_=shared_data.moleFractions;

    // GET EIG0_M
    double eig0_M;
    auto it1 = colliderInfo.find("M");
    if (it1 != colliderInfo.end() && it1->second.hasKey("rate-constants")){ 
        eig0_M=ArrheniusRate(AnyValue(it1->second["eig0"]), it1->second.units(), rate_units_).evalRate(logT_, recipT_);
        writeMsg("eig0_M = ",eig0_M);
    } 
    else if(it1 != colliderInfo.end() && it1->second.hasKey("Troe")){ 
        eig0_M=ArrheniusRate(AnyValue(it1->second["eig0"]), it1->second.units(), rate_units_).evalRate(logT_, recipT_);
        writeMsg("eig0_M = ",eig0_M);
        // colliders_i["type"]="LMR_R"; //this dummy key needs to be defined because falloff.cpp requires it
    }
    else if(it1 != colliderInfo.end() && it1->second.hasKey("pressure-range")){
        eig0_M=ArrheniusRate(AnyValue(it1->second["eig0"]), it1->second.units(), rate_units_).evalRate(logT_, recipT_);
        writeMsg("eig0_M = ",eig0_M);
    }
    else {
        throw InputFileError("LmrRate::getkM", m_input,"Not enough data provided for species 'M'.");
    }

    // GET EIG0_mix

    for (size_t i=0; i<shared_data.allSpecies.size(); i++){ //testing each species listed at the top of yaml file
        double Xi=moleFractions_[i];
        auto it2 = colliderInfo.find(shared_data.allSpecies[i]);
        if (it2 != colliderInfo.end()) { //yaml species has at least an eig0 value provided for LMRR 
            eig0_mix += Xi*ArrheniusRate(AnyValue(it2->second["eig0"]), it2->second.units(), rate_units_).evalRate(logT_, recipT_);
        }
        else if (it2 == colliderInfo.end()) { //yaml species has no data provided for LMRR (treat as "M")
            eig0_mix += Xi*eig0_M;
        }
        else{
            throw InputFileError("LmrRate::geteig0mix", m_input,"Cannot compute eig0mix due to invalid LMRR yaml input.");
        }
    }

    double eig0; //eig0 val of a single species

    logP_=shared_data.logP;
    logT_=shared_data.logT;
    pressure_=shared_data.pressure;
    recipT_=shared_data.recipT;
    temperature_=shared_data.temperature;
    ready_=shared_data.ready;
    

    // GET K_M
    double k_M;
    auto it3 = colliderInfo.find("M");
    if (it3 != colliderInfo.end() && it3->second.hasKey("rate-constants")){ 
        logPeff_= logP_+log(eig0_mix)-log(eig0_M); //need to update this every time before evalPlogRate
        k_M = evalPlogRate(it3);
        writeMsg("kM_plog = ",k_M);
    } 
    else if(it3 != colliderInfo.end() && it3->second.hasKey("Troe")){ 
        k_M = evalTroeRate(it3);
        writeMsg("kM_troe = ",k_M);
        // colliders_i["type"]="LMR_R"; //this dummy key needs to be defined because falloff.cpp requires it
    }
    else if(it3 != colliderInfo.end() && it3->second.hasKey("pressure-range")){
        k_M = evalChebyshevRate(it3);
        writeMsg("kM_cheb = ",k_M);
    }
    else {
        throw InputFileError("LmrRate::getkM", m_input,"Not enough data provided for species 'M'.");
    }

    //GET K_LMRR
    for (size_t i=0; i<shared_data.allSpecies.size(); i++){ //testing each species listed at the top of yaml file
        auto it4 = colliderInfo.find(shared_data.allSpecies[i]);
        double Xi=moleFractions_[i];
        // double k_i;
        // double eig0_i;
        if (it4 != colliderInfo.end() && it4->second.hasKey("rate-constants")){ 
            writelog(it4->first);writelog("\n"); //speciesID
            double eig0 = ArrheniusRate(AnyValue(it4->second["eig0"]), it4->second.units(), rate_units_).evalRate(logT_, recipT_);
            logPeff_= logP_+log(eig0_mix)-log(eig0); //need to update this every time before evalPlogRate
            double k_i = evalPlogRate(it4);
            double Xtilde = eig0*Xi/eig0_mix;
            k_LMR_ += k_i*Xtilde;
            // k_LMR_ += evalPlogRate(it4)*eig0*Xi/eig0_mix;
            writeMsg("eig0_i_plog = ",eig0); writeMsg("k_i_plog = ",k_i);
        } 
        else if(it4 != colliderInfo.end() && it4->second.hasKey("Troe")){ 
            writelog(it4->first);writelog("\n"); //speciesID
            // colliders_i["type"]="LMR_R"; //this dummy key needs to be defined because falloff.cpp requires it
            double eig0 = ArrheniusRate(AnyValue(it4->second["eig0"]), it4->second.units(), rate_units_).evalRate(logT_, recipT_);
            double k_i = evalTroeRate(it4);



            double Xtilde = eig0*Xi/eig0_mix;
            k_LMR_ += k_i*Xtilde;
            // k_LMR_ += evalTroeRate(it4)*eig0*Xi/eig0_mix;
            writeMsg("eig0_i_troe = ",eig0); writeMsg("k_i_troe = ",k_i);
        }
        else if(it4 != colliderInfo.end() && it4->second.hasKey("pressure-range")){ 
            writelog(it4->first);writelog("\n"); //speciesID
            double eig0 = ArrheniusRate(AnyValue(it4->second["eig0"]), it4->second.units(), rate_units_).evalRate(logT_, recipT_);
            double k_i = evalChebyshevRate(it4);
            double Xtilde = eig0*Xi/eig0_mix;
            k_LMR_ += k_i*Xtilde;    
            // k_LMR_ += evalChebyshevRate(it4)*eig0*Xi/eig0_mix;        
            writeMsg("eig0_i_cheb = ",eig0); writeMsg("k_i_cheb = ",k_i);
        }
        else if(it4 != colliderInfo.end() && !(it4->second.hasKey("pressure-range")) && !(it4->second.hasKey("Troe")) && !(it4->second.hasKey("rate-constants"))){ //yaml species has an eig0 but no additional LMRR data, so treat its rate as same as "M"
            writelog(it4->first);writelog("\n"); //speciesID
            double eig0 = ArrheniusRate(AnyValue(it4->second["eig0"]), it4->second.units(), rate_units_).evalRate(logT_, recipT_);
            double Xtilde = eig0*Xi/eig0_mix;
            k_LMR_ += k_M*Xtilde;  
            // k_LMR_ += k_M*eig0*Xi/eig0_mix;
            writeMsg("eig0_i_eigOnly = ",eig0); writeMsg("k_i_eigOnly = ",k_M);
        }
        else if(it4 == colliderInfo.end()){ //yaml species has no LMRR data, so treat its rate and eig0 as same as "M"
            double Xtilde = eig0_M*Xi/eig0_mix;
            k_LMR_ += k_M*Xtilde; 
            // k_LMR_ += k_M*eig0_M*Xi/eig0_mix;
            writeMsg("eig0_i_noLMR = ",eig0_M); writeMsg("k_i_noLMR = ",k_M);
        }
        else{
            throw InputFileError("LmrRate::evalFromStruct", m_input,"LMRR reaction has invalid yaml input.");
        }
        // writeMsg("k_i_final = ",k_i);
        // k_LMR += k_i*eig0_i*moleFractions_[i]/eig0_mix; //Note: Xtilde = eig0_i*moleFractions_[i]/eig0_mix;
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