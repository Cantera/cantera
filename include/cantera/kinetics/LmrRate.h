//! @file LmrRate.h
// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_LMRRATE_H
#define CT_LMRRATE_H
#include "cantera/kinetics/Arrhenius.h"
#include <boost/variant.hpp>
#include "cantera/kinetics/Falloff.h"
#include "cantera/kinetics/ChebyshevRate.h"
#include "cantera/kinetics/PlogRate.h"
// #include "cantera/base/Array.h"
// #include "cantera/base/FactoryBase.h"
// #include "cantera/kinetics/Reaction.h"

namespace Cantera{

struct LmrData : public ReactionData{
    // LmrData() = default;
    LmrData();
    void update(double T, double P) override {
        ReactionData::update(T);
        pressure = P;
        logP = std::log(P);
    }
    bool update(const ThermoPhase& phase, const Kinetics& kin) override;
    using ReactionData::update;
    void perturbPressure(double deltaP);
    void restore() override;
    virtual void resize(size_t nSpecies, size_t nReactions, size_t nPhases) override {
        moleFractions.resize(nSpecies, NAN);
        ready = true;
    }
    void invalidateCache() override {
        ReactionData::invalidateCache();
        pressure = NAN;
    }
    double pressure = NAN; //!< pressure
    double logP = 0.0; //!< logarithm of pressure
    bool ready = false; //!< boolean indicating whether vectors are accessible
    vector<double> moleFractions;
    int mfNumber; 
    vector<string> allSpecies; //list of all yaml species (not just those for which LMRR data exists)  
// protected:
    double m_pressure_buf = -1.0; //!< buffered pressure
};

class LmrRate final : public ReactionRate
{
public:
    LmrRate() = default;//! Default constructor.
    vector<size_t> colliderIndices;
    map<string, AnyMap> colliderInfo;
    using RateTypes = boost::variant<PlogRate, TroeRate, ChebyshevRate>;
    using DataTypes = boost::variant<PlogData, FalloffData, ChebyshevData>;
    vector<RateTypes> rateObjs;
    vector<DataTypes> dataObjs;
    // vector<AnyMap> colliderNodes;
    
    vector<ArrheniusRate> eigObjs;
    vector<string> colliderNames;
    RateTypes rateObj_M;
    DataTypes dataObj_M;
    ArrheniusRate eigObj_M;
    // AnyMap node_M;
    size_t nSpecies;
    UnitStack rate_units_;
    // AnyMap node_;
    // double logP_;
    // double logT_;
    // double pressure_;
    // double recipT_;
    // double temperature_;
    // bool ready_;
    // vector<double> moleFractions_;
    // double logPeff_;
    // double k_LMR_;
    double eig0_mix;

    // PlogData plog_data;
    // FalloffData troe_data;
    // ChebyshevData cheb_data;

    explicit LmrRate(const std::multimap<double, ArrheniusRate>& rates);
    LmrRate(const AnyMap& node, const UnitStack& rate_units={});
    unique_ptr<MultiRateBase> newMultiRate() const override {
        return make_unique<MultiRate<LmrRate, LmrData>>();
    }
    const string type() const override { return "LMR_R"; } //! Identifier of reaction rate type
    void setParameters(const AnyMap& node, const UnitStack& rate_units) override;
    void getParameters(AnyMap& rateNode, const Units& rate_units) const;
    void getParameters(AnyMap& rateNode) const override {
        return getParameters(rateNode, Units(0));
    }
    double evalPlogRate(const LmrData& shared_data, DataTypes& dataObj, RateTypes& rateObj, double& eig0);
    double evalTroeRate(const LmrData& shared_data, DataTypes& dataObj, RateTypes& rateObj);
    double evalChebyshevRate(const LmrData& shared_data, DataTypes& dataObj, RateTypes& rateObj, double& eig0);

    void setContext(const Reaction& rxn, const Kinetics& kin) override;

    double evalFromStruct(const LmrData& shared_data);
    void validate(const string& equation, const Kinetics& kin) override; //removed from cpp, but re-insert later
    void writeMsg(string str, double dbl){writelog(str); writelog(std::to_string(dbl)); writelog("\n");}

protected:
    // double logP_ = -1000;


};


}
#endif