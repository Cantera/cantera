//! @file LinearBurkeRate.h
// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_LINEARBURKERATE_H
#define CT_LINEARBURKERATE_H

#include "cantera/kinetics/Arrhenius.h"
#include "cantera/kinetics/Falloff.h"
#include "cantera/kinetics/ChebyshevRate.h"
#include "cantera/kinetics/PlogRate.h"

#include <variant>

namespace Cantera
{

//! Data container holding shared data specific to LinearBurkeRate
/**
 * The data container `LinearBurkeData` holds precalculated data common to
 * all `LinearBurkeRate` objects.
 */
struct LinearBurkeData : public ReactionData
{
    LinearBurkeData();

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

//! Pressure-dependent and composition-dependent reaction rate calculated according to
//! the reduced-pressure linear mixture rule (LMR-R).
//!
//! This parameterization is described by Singal et al. @cite singal2024 and in the
//! [science reference](../reference/kinetics/rate-constants.html#linear-burke-rate-expressions)
//! documentation.
//! @ingroup otherRateGroup
class LinearBurkeRate final : public ReactionRate
{
public:
    //! Default constructor.
    LinearBurkeRate() = default;

    LinearBurkeRate(const AnyMap& node, const UnitStack& rate_units={});

    unique_ptr<MultiRateBase> newMultiRate() const override {
        return make_unique<MultiRate<LinearBurkeRate, LinearBurkeData>>();
    }

    //! Identifier of reaction rate type
    const string type() const override { return "linear-Burke"; }

    //! Perform object setup based on AnyMap node information
    /*!
     *  @param node  AnyMap containing rate information
     *  @param rate_units  Unit definitions specific to rate information
     */
    void setParameters(const AnyMap& node, const UnitStack& rate_units) override;

    void getParameters(AnyMap& rateNode) const override;

    //! Type alias that refers to PlogRate, TroeRate, and ChebyshevRate
    using RateTypes = std::variant<PlogRate, TroeRate, ChebyshevRate>;
    //! Type alias that refers to PlogData, FalloffData, and ChebyshevData
    using DataTypes = std::variant<PlogData, FalloffData, ChebyshevData>;

    double evalFromStruct(const LinearBurkeData& shared_data);

    void setContext(const Reaction& rxn, const Kinetics& kin) override;

    void validate(const string& equation, const Kinetics& kin) override;

protected:
    //! Evaluate overall reaction rate using PLOG to evaluate pressure-dependent aspect
    //! of the reaction
    double evalPlogRate(const LinearBurkeData& shared_data, DataTypes& dataObj,
                        RateTypes& rateObj, double logPeff);

    //! Evaluate overall reaction rate using Troe to evaluate pressure-dependent aspect
    //! of the reaction
    double evalTroeRate(const LinearBurkeData& shared_data, DataTypes& dataObj,
                        RateTypes& rateObj, double logPeff);

    //! Evaluate overall reaction rate using Chebyshev to evaluate pressure-dependent
    //! aspect of the reaction
    double evalChebyshevRate(const LinearBurkeData& shared_data, DataTypes& dataObj,
                             RateTypes& rateObj, double logPeff);

    //! String name of each collider, appearing in the same order as that of the
    //! original reaction input.
    vector<string> m_colliderNames;

    //! Index of each collider in the kinetics object species list where the vector
    //! elements appear in the same order as that of the original reaction input.
    vector<size_t> m_colliderIndices;
    //! Indicates which colliders have a distinct k(T,P) versus only an efficiency
    vector<bool> m_hasRateConstant;

    //! Third-body collision efficiency object for k(T,P,X) and eig0_mix calculation
    vector<ArrheniusRate> m_epsObjs1;
    //! Third-body collision efficiency object for logPeff calculation
    vector<ArrheniusRate> m_epsObjs2;
    //! Third-body collision efficiency object for the reference collider M
    //! (eig0_M/eig0_M = 1 always)
    ArrheniusRate m_epsObj_M;

    //! Stores rate objects corresponding the reference collider M, which can be
    //! either PlogRate, TroeRate, or ChebyshevRate
    RateTypes m_rateObj_M;
    //! Stores rate objects corresponding to each non-M collider, which can be either
    //! PlogRate, TroeRate, or ChebyshevRate
    vector<RateTypes> m_rateObjs;

    //! Stores data objects corresponding to the reference collider M, which can be
    //! either PlogData, TroeData, or ChebyshevData
    DataTypes m_dataObj_M; //!< collider M
    //! Stores data objects corresponding to each non-M collider, which can be either
    //! PlogData, TroeData, or ChebyshevData
    vector<DataTypes> m_dataObjs; //!< list for non-M colliders
};

}
#endif
