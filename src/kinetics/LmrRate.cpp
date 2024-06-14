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
    writelog("setParameters::1"); writelog("\n");
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
    eigObj_M=ArrheniusRate(AnyValue(colliders[0]["eig0"]), colliders[0].units(), rate_units_);
    if (colliders[0].hasKey("rate-constants")){
        writelog("setParameters::2"); writelog("\n");
        rateObj_M = PlogRate();
        dataObj_M = PlogData();
    } else if (colliders[0].hasKey("Troe")){
        writelog("setParameters::3"); writelog("\n");
        rateObj_M = TroeRate();
        dataObj_M = FalloffData(); 
    } else if (colliders[0].hasKey("pressure-range")){
        writelog("setParameters::4"); writelog("\n");
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
            writelog("setParameters::5"); writelog("\n");
            rateObjs.push_back(PlogRate());
            dataObjs.push_back(PlogData());
        } else if (colliders[i].hasKey("Troe")){
            writelog("setParameters::6"); writelog("\n");
            rateObjs.push_back(TroeRate());
            dataObjs.push_back(FalloffData());   
        } else if (colliders[i].hasKey("pressure-range")){
            writelog("setParameters::7"); writelog("\n");
            rateObjs.push_back(ChebyshevRate());
            dataObjs.push_back(ChebyshevData());
        } else { //Collider has an eig0 specified, but no other info is provided. Assign it the same rate and data objects as "M"
            writelog("setParameters::8"); writelog("\n");
            rateObjs.push_back(rateObj_M);
            dataObjs.push_back(dataObj_M);   
        }
    }
}

void LmrRate::validate(const string& equation, const Kinetics& kin){
    writelog("validate::1"); writelog("\n");
}

void LmrRate::setContext(const Reaction& rxn, const Kinetics& kin){   
    writelog("setContext::1"); writelog("\n");
    for (size_t i=0; i<names.size();i++){ //Starts at 1, because names[0] == "M" //THIS LOOP IS NOT WORKING
        writelog("setContext::2"); writelog("\n");
        colliderIndices.push_back(kin.kineticsSpeciesIndex(names[i]));
        // double val = kin.kineticsSpeciesIndex(names[i]);
        // writelog("coll idx = ", val); writelog("\n");
    }
    nSpecies = kin.nTotalSpecies();
    // writelog("collisionIndices[0] = ", colliderIndices[0]);

    if (colliderIndices.empty()){
        throw InputFileError("LmrRate::setContext", m_input,"collidersIndices is empty for some reason");
    }
    // if (nSpecies!=9){
    //     throw InputFileError("LmrRate::setContext", m_input,"nSpecies is incorrect.");
    // }
}



double LmrRate::evalPlogRate(PlogRate& rate, PlogData& data, AnyMap node){
    writelog("evalPlogRate::1"); writelog("\n");
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
    writelog("evalTroeRate::1"); writelog("\n");
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
    writelog("evalChebyshevRate::1"); writelog("\n");
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
    writelog("evalFromStruct::1"); writelog("\n");
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
        writelog("evalFromStruct::2"); writelog("\n");
        for (size_t j=0; j<colliderIndices.size(); j++){
            writelog("evalFromStruct::3"); writelog("\n");
            if (i==colliderIndices[j]){ // Species is in collider list
                writelog("evalFromStruct::4"); writelog("\n");
                eig0_mix += moleFractions_[i]*eigObjs[j].evalRate(logT_, recipT_);
            } else { // Species not in collider list so treat as "M"
                writelog("evalFromStruct::5"); writelog("\n");
                sigmaX_M += moleFractions_[i];
            }
        }
    }
    if (eig0_mix==0){
        throw InputFileError("LmrRate::evalFromStruct", m_input,"eig0_mix==0 for some reason");
    }
    double eig0_M = eigObj_M.evalRate(logT_, recipT_);
    eig0_mix += sigmaX_M*eig0_M; // add all M colliders to eig0_mix in a single step
    writelog("evalFromStruct::6"); writeMsg(" eig0_mix = ",eig0_mix);
    k_LMR_=0;
    sigmaX_M=0.0;
    for (size_t i=0; i<nSpecies; i++){ //testing each species listed at the top of yaml file
        writelog("evalFromStruct::7"); writelog("\n");
        for (size_t j=0; j<colliderIndices.size(); j++){
            writelog("evalFromStruct::8"); writelog("\n");
            if (i==colliderIndices[j]){ // Species is in collider list
                double eig0 = eigObjs[j].evalRate(logT_, recipT_);
                if (eig0==0){
                    throw InputFileError("LmrRate::evalFromStruct", m_input,"eig0 is 0 for some reason");
                }
                writelog("evalFromStruct::9"); writeMsg(" eig0 = ",eig0);
                if (rateObjs[j].which()==0){ // 0 means PlogRate     
                    // logP_= shared_data.logP+log(eig0_mix)-log(eig0); //replaces logP with log of the effective pressure w.r.t. eig0
                    // PlogData& data = boost::get<PlogData>(dataObjs.at(j));
                    // PlogRate& rate = boost::get<PlogRate>(rateObjs.at(j));
                    PlogData& data = boost::get<PlogData>(dataObjs[j]);
                    PlogRate& rate = boost::get<PlogRate>(rateObjs[j]);
                    k_LMR_ += evalPlogRate(rate,data,colliderNodes[j])*eig0*moleFractions_[i]/eig0_mix;
                    evalPlogRate(boost::get<PlogRate>(rateObj_M),boost::get<PlogData>(dataObj_M),node_M)*eig0_M*sigmaX_M/eig0_mix;
                    // logP_ = shared_data.logP; //return to the "normal" logP value to avoid messing up other calcs
                    writelog("evalFromStruct::10"); writeMsg(" k_i_plog = ",evalPlogRate(rate,data,colliderNodes[j]));
                }
                else if (rateObjs[j].which()==1){ // 1 means TroeRate  
                    FalloffData& data = boost::get<FalloffData>(dataObjs[j]);
                    TroeRate& rate = boost::get<TroeRate>(rateObjs[j]);
                    k_LMR_ += evalTroeRate(rate,data,colliderNodes[j])*eig0*moleFractions_[i]/eig0_mix;
                    writelog("evalFromStruct::11"); writeMsg(" k_i_troe = ",evalTroeRate(rate,data,colliderNodes[j]));
                }
                else if (rateObjs[j].which()==2){ // 2 means ChebyshevRate  
                    ChebyshevData& data = boost::get<ChebyshevData>(dataObjs[j]);
                    ChebyshevRate& rate = boost::get<ChebyshevRate>(rateObjs[j]);
                    k_LMR_ += evalChebyshevRate(rate,data,colliderNodes[j])*eig0*moleFractions_[i]/eig0_mix;
                    writelog("evalFromStruct::12"); writeMsg(" k_i_cheb = ",evalChebyshevRate(rate,data,colliderNodes[j]));
                }
                else {
                    throw InputFileError("LmrRate::evalFromStruct", m_input,"Something went wrong...");
                }
            }
            else {
                writelog("evalFromStruct::13"); writelog("\n");
                sigmaX_M += moleFractions_[i];
            }
        }
    }
    writelog("evalFromStruct::14"); writelog("\n");
    if (rateObj_M.which()==0){ // 0 means PlogRate
        // logP_= shared_data.logP+log(eig0_mix)-log(eig0_M); //replaces logP with log of the effective pressure w.r.t. eig0_M
        k_LMR_ += evalPlogRate(boost::get<PlogRate>(rateObj_M),boost::get<PlogData>(dataObj_M),node_M)*eig0_M*sigmaX_M/eig0_mix;
        // logP_ = shared_data.logP; //return to the "normal" logP value to avoid messing up other calcs
        writelog("evalFromStruct::15"); writeMsg(" k_M_plog = ",evalPlogRate(boost::get<PlogRate>(rateObj_M),boost::get<PlogData>(dataObj_M),node_M));
    }
    else if (rateObj_M.which()==1){ // 1 means TroeRate  
        writelog("evalFromStruct::16"); writeMsg(" k_M_troe = ",evalTroeRate(boost::get<TroeRate>(rateObj_M),boost::get<FalloffData>(dataObj_M),node_M));
        k_LMR_ += evalTroeRate(boost::get<TroeRate>(rateObj_M),boost::get<FalloffData>(dataObj_M),node_M)*eig0_M*sigmaX_M/eig0_mix;
    }
    else if (rateObj_M.which()==2){ // 2 means ChebyshevRate 
        writelog("evalFromStruct::17"); writeMsg(" k_M_cheb = ",evalChebyshevRate(boost::get<ChebyshevRate>(rateObj_M),boost::get<ChebyshevData>(dataObj_M),node_M));
        k_LMR_ += evalChebyshevRate(boost::get<ChebyshevRate>(rateObj_M),boost::get<ChebyshevData>(dataObj_M),node_M)*eig0_M*sigmaX_M/eig0_mix;
    }
    writelog("evalFromStruct::18"); writelog("\n");
    // writemsg("k_LMR = ",k_LMR_);
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

