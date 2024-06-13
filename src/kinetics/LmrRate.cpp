//! @file LmrRate.cpp
// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/LmrRate.h"
#include "cantera/thermo/ThermoPhase.h"
// #include "cantera/kinetics/Falloff.h"
#include "cantera/kinetics/Kinetics.h"
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
    writelog("1"); writelog("\n");
    ReactionRate::setParameters(node, rate_units);
    rate_units_=rate_units;
    if(!node.hasKey("collider-list")){
        throw InputFileError("LmrRate::setParameters", m_input,"Yaml input for LMR-R does not follow the necessary structure.");
    }
    auto& colliders = node["collider-list"].asVector<AnyMap>();
    if (colliders[0]["collider"].as<std::string>() != "M"){
        throw InputFileError("LmrRate::setParameters", m_input,"The first species defined in yaml input must be 'M'.");
    }
    node_M=colliders[0];
    eig0_M=ArrheniusRate(AnyValue(colliders[0]["eig0"]), colliders[0].units(), rate_units_);
    if (colliders[0].hasKey("rate-constants")){
        rateObj_M = PlogRate();
        dataObj_M = PlogData();
    } else if (colliders[0].hasKey("Troe")){
        rateObj_M = TroeRate();
        dataObj_M = FalloffData(); 
    } else if (colliders[0].hasKey("pressure-range")){
        rateObj_M = ChebyshevRate();
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
        names.push_back(colliders[i]["collider"].as<std::string>());
        colliderNodes.push_back(colliders[i]);
        eigObjs.push_back(ArrheniusRate(AnyValue(colliders[i]["eig0"]), colliders[i].units(), rate_units_));
        if (colliders[i].hasKey("rate-constants")){
            rateObjs.push_back(PlogRate());
            dataObjs.push_back(PlogData());
        } else if (colliders[i].hasKey("Troe")){
            rateObjs.push_back(TroeRate());
            dataObjs.push_back(FalloffData());   
        } else if (colliders[i].hasKey("pressure-range")){
            rateObjs.push_back(ChebyshevRate());
            dataObjs.push_back(ChebyshevData());
        } else { //Collider has an eig0 specified, but no other info is provided. Assign it the same rate and data objects as "M"
            rateObjs.push_back(rateObj_M);
            dataObjs.push_back(dataObj_M);   
        }
    }
}

void LmrRate::setContext(const Reaction& rxn, const Kinetics& kin){   
    for (int i=1; i<names.size();i++){ //Starts at 1, because names[0] == "M"
        // colliderIncides.push_back(kin.kineticsSpeciesIndex(names[i]));
        colliderIncides.push_back(1);
    }
    nSpecies = kin.nTotalSpecies();
}

void LmrRate::validate(const string& equation, const Kinetics& kin){}

double LmrRate::evalPlogRate(PlogRate& rate, PlogData& data, AnyMap node){
    data.logP = logP_; //replaces logP with log of the effective pressure w.r.t. eig0_M
    data.logT = logT_;
    data.pressure = pressure_;
    data.recipT = recipT_;
    data.temperature = temperature_;
    rate.setParameters(node,rate_units_); //it->second refers to the yaml data for the ith collider
    rate.updateFromStruct(data);
    return rate.evalFromStruct(data);
}
double LmrRate::evalTroeRate(TroeRate& rate, FalloffData& data, AnyMap node){
    data.conc_3b = moleFractions_;
    data.logT = logT_;
    // data.molar_density = pressure_; //
    data.ready = ready_;
    data.recipT = recipT_;
    data.temperature = temperature_;
    rate.setParameters(node,rate_units_); //it->second refers to the yaml data for the ith collider
    return rate.evalFromStruct(data);
}

double LmrRate::evalChebyshevRate(ChebyshevRate& rate, ChebyshevData& data, AnyMap node){
    data.log10P=logP_; //THIS IS INCORRECT. logP_ is natural log, not base 10!!!
    data.logT=logT_;
    data.pressure=pressure_;
    data.recipT=recipT_;
    data.temperature=temperature_;
    rate.setParameters(node,rate_units_); //it->second refers to the yaml data for the ith collider
    rate.updateFromStruct(data);
    return rate.evalFromStruct(data);
}

double LmrRate::evalFromStruct(const LmrData& shared_data){
    logP_=shared_data.logP;
    logT_=shared_data.logT;
    pressure_=shared_data.pressure;
    recipT_=shared_data.recipT;
    temperature_=shared_data.temperature;
    ready_=shared_data.ready;
    moleFractions_=shared_data.moleFractions;

    double eig0_mix=0;
    double sigmaX_M=0.0;
    for (int i=0; i<nSpecies; i++){ //testing each species listed at the top of yaml file
        for (int j=0; j<colliderIndices.size(); j++){
            if (i==colliderIncides[j]){ // Species is in collider list
                eig0_mix += moleFractions_[i]*eigObjs[j].evalRate(logT_, recipT_);
            } else { // Species not in collider list so treat as "M"
                sigmaX_M += moleFractions_[i];
            }
        }
    }
    eig0_mix += sigmaX_M*eig0_M.evalRate(logT_, recipT_); // add all M colliders to eig0_mix in a single step

    k_LMR_=0;
    sigmaX_M=0.0;
    for (int i=0; i<nSpecies; i++){ //testing each species listed at the top of yaml file
        // double Xi=moleFractions_[i];
        for (int j=0; j<colliderIndices.size(); j++){
            if (i==colliderIncides[j]){ // Species is in collider list
                double eig0 = eigObjs[j].evalRate(logT_, recipT_);
                if (rateObjs.at(j).which()==0){ // 0 means PlogRate 
                    logP_= logP_+log(eig0_mix)-log(eig0); //replaces logP with log of the effective pressure w.r.t. eig0_M
                    PlogData& data = boost::get<PlogData>(dataObjs.at(j));
                    PlogRate& rate = boost::get<PlogRate>(rateObjs.at(j));
                    k_LMR_ += evalPlogRate(rate,data,colliderNodes[j])*eig0*moleFractions_[i]/eig0_mix;
                    logP_ = shared_data.logP; //return to the "normal" logP value to avoid messing up other calcs
                    writeMsg("species: ",i); writeMsg("k_i_plog = ",evalPlogRate(rate,data,colliderNodes[j]));
                }
                else if (rateObjs.at(i).which()==1){ // 1 means TroeRate   
                    FalloffData& data = boost::get<FalloffData>(dataObjs.at(i));
                    TroeRate& rate = boost::get<TroeRate>(rateObjs.at(i));
                    k_LMR_ += evalTroeRate(rate,data,colliderNodes[j])*eig0*moleFractions_[i]/eig0_mix;
                    writeMsg("species: ",i); writeMsg("k_i_troe = ",evalTroeRate(rate,data,colliderNodes[j]));
                }
                else if (rateObjs.at(i).which()==2){ // 2 means ChebyshevRate  
                    ChebyshevData& data = boost::get<ChebyshevData>(dataObjs.at(i));
                    ChebyshevRate& rate = boost::get<ChebyshevRate>(rateObjs.at(i));
                    k_LMR_ += evalChebyshevRate(rate,data,colliderNodes[j])*eig0*moleFractions_[i]/eig0_mix;
                    writeMsg("species: ",i); writeMsg("k_i_cheb = ",evalChebyshevRate(rate,data,colliderNodes[j]));
                }
                else {
                    throw InputFileError("LmrRate::evalFromStruct", m_input,"Something went wrong...");
                }
            }
            else {
                sigmaX_M += moleFractions_[i];
            }
        }
    }

    if (rateObj_M.which()==0){ // 0 means PlogRate 
        k_LMR_ += evalPlogRate(boost::get<PlogRate>(rateObj_M),boost::get<PlogData>(dataObj_M),node_M)*eig0_M.evalRate(logT_, recipT_)*sigmaX_M/eig0_mix;
    }
    else if (rateObj_M.which()==1){ // 1 means TroeRate  
        k_LMR_ += evalTroeRate(boost::get<TroeRate>(rateObj_M),boost::get<FalloffData>(dataObj_M),node_M)*eig0_M.evalRate(logT_, recipT_)*sigmaX_M/eig0_mix;
    }
    else if (rateObj_M.which()==2){ // 2 means ChebyshevRate 
        k_LMR_ += evalChebyshevRate(boost::get<ChebyshevRate>(rateObj_M),boost::get<ChebyshevData>(dataObj_M),node_M)*eig0_M.evalRate(logT_, recipT_)*sigmaX_M/eig0_mix;
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

