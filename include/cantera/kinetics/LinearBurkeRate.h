//! @file LINEARBURKERATE.h
// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_LINEARBURKERATE_H
#define CT_LINEARBURKERATE_H
#include "cantera/kinetics/Arrhenius.h"
#include <boost/variant.hpp>
#include "cantera/kinetics/Falloff.h"
#include "cantera/kinetics/ChebyshevRate.h"
#include "cantera/kinetics/PlogRate.h"

namespace Cantera
{

//! Data container holding shared data specific to LinearBurkeRate
/**
 * The data container `LmrData` holds precalculated data common to
 * all `LinearBurkeRate` objects.
 */
struct LmrData : public ReactionData
{
    LmrData();

    void update(double T, double P) override
    {
        ReactionData::update(T);
        pressure = P;
        logP = std::log(P);
    }

    bool update(const ThermoPhase& phase, const Kinetics& kin) override;

    using ReactionData::update;

    //! Perturb pressure of data container
    /**
     * The method is used for the evaluation of numerical derivatives.
     * @param  deltaP  relative pressure perturbation
     */
    void perturbPressure(double deltaP);

    void restore() override;

    virtual void resize(size_t nSpecies, size_t nReactions, size_t nPhases) override
    {
        moleFractions.resize(nSpecies, NAN);
        ready = true;
    }

    void invalidateCache() override
    {
        ReactionData::invalidateCache();
        pressure = NAN;
    }

    double pressure = NAN; //!< Pressure
    double logP = 0.0; //!< Logarithm of pressure
    bool ready = false; //!< Boolean indicating whether vectors are accessible
    vector<double> moleFractions;
    int mf_number;

protected:
    double m_pressure_buf = -1.0;
};


//! Pressure-dependent and composition-dependent reaction rate calculated
//! according to the reduced-pressure linear mixture rule (LMR-R) developed
//! at Columbia University. @cite singal2025 [CITATION NOT YET ACTIVE]
/*!
 * [ADD IN THE MATHEMATICAL FORMULA]
 */
class LinearBurkeRate final : public ReactionRate
{
public:
    //! Default constructor.
    LinearBurkeRate() = default;

    //! Constructor from Arrhenius rate expressions at a set of pressures
    explicit LinearBurkeRate(const std::multimap<double, ArrheniusRate>& rates);

    LinearBurkeRate(const AnyMap& node, const UnitStack& rate_units={});

    unique_ptr<MultiRateBase> newMultiRate() const override {
        return make_unique<MultiRate<LinearBurkeRate, LmrData>>();
    }

    //! Identifier of reaction rate type
    const string type() const override { return "linear-burke"; }

    //! Perform object setup based on AnyMap node information
    /*!
     *  @param node  AnyMap containing rate information
     *  @param rate_units  Unit definitions specific to rate information
     */
    void setParameters(const AnyMap& node, const UnitStack& rate_units) override;

    void getParameters(AnyMap& rateNode, const Units& rate_units) const;
    void getParameters(AnyMap& rateNode) const override {
        return getParameters(rateNode, Units(0));
    }

    //! Evaluate overall reaction rate, using Troe/PLOG/Chebyshev to evaluate
    //! pressure-dependent aspect of the reaction
    /*!
     *  @param shared_data  data shared by all reactions of a given type
     */
    using RateTypes = boost::variant<PlogRate, TroeRate, ChebyshevRate>;
    using DataTypes = boost::variant<PlogData, FalloffData, ChebyshevData>;
    double evalPlogRate(const LmrData& shared_data, DataTypes& dataObj, RateTypes& rateObj);
    double evalTroeRate(const LmrData& shared_data, DataTypes& dataObj, RateTypes& rateObj);
    double evalChebyshevRate(const LmrData& shared_data, DataTypes& dataObj, RateTypes& rateObj);
    double evalFromStruct(const LmrData& shared_data);

    void setContext(const Reaction& rxn, const Kinetics& kin) override;

    void validate(const string& equation, const Kinetics& kin) override;

    vector<size_t> colliderIndices;
    map<string, AnyMap> colliderInfo;
    // Third-body collision efficiency objects (eps = eig0_i/eig0_M)
    // epsObjs1 used for k(T,P,X) and eig0_mix calculation
    // epsObjs2 used for logPeff calculation
    vector<ArrheniusRate> epsObjs1;
    vector<ArrheniusRate> epsObjs2;
    vector<string> colliderNames;
    vector<RateTypes> rateObjs;
    vector<DataTypes> dataObjs;
    RateTypes rateObj_M;
    DataTypes dataObj_M;
    // Third-body collision efficiency object for M (eig0_M/eig0_M = 1)
    ArrheniusRate epsObj_M;
    size_t nSpecies;
    double logPeff_;
    double eps_mix;
};

}
#endif