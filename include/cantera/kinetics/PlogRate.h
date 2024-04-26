//! @file PlogRate.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_PLOGRATE_H
#define CT_PLOGRATE_H

#include "cantera/kinetics/Arrhenius.h"

namespace Cantera
{

//! Data container holding shared data specific to PlogRate
/**
 * The data container `PlogData` holds precalculated data common to
 * all `PlogRate` objects.
 */
struct PlogData : public ReactionData
{
    PlogData() = default;

    void update(double T) override;

    void update(double T, double P) override {
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

    void invalidateCache() override {
        ReactionData::invalidateCache();
        pressure = NAN;
    }

    double pressure = NAN; //!< pressure
    double logP = 0.0; //!< logarithm of pressure

protected:
    double m_pressure_buf = -1.0; //!< buffered pressure
};


//! Pressure-dependent reaction rate expressed by logarithmically interpolating
//! between Arrhenius rate expressions at various pressures.
/*!
 * Given two rate expressions at two specific pressures:
 *
 *   * @f$ P_1: k_1(T) = A_1 T^{b_1} e^{-E_1 / RT} @f$
 *   * @f$ P_2: k_2(T) = A_2 T^{b_2} e^{-E_2 / RT} @f$
 *
 * The rate at an intermediate pressure @f$ P_1 < P < P_2 @f$ is computed as
 * @f[
 *  \ln k(T,P) = \ln k_1(T) + \bigl(\ln k_2(T) - \ln k_1(T)\bigr)
 *      \frac{\ln P - \ln P_1}{\ln P_2 - \ln P_1}
 * @f]
 * Multiple rate expressions may be given at the same pressure, in which case
 * the rate used in the interpolation formula is the sum of all the rates given
 * at that pressure @cite gou2011. For pressures outside the given range, the
 * rate expression at the nearest pressure is used.
 * @ingroup otherRateGroup
 */
class PlogRate final : public ReactionRate
{
public:
    //! Default constructor.
    PlogRate() = default;

    //! Constructor from Arrhenius rate expressions at a set of pressures
    explicit PlogRate(const std::multimap<double, ArrheniusRate>& rates);

    PlogRate(const AnyMap& node, const UnitStack& rate_units={});

    unique_ptr<MultiRateBase> newMultiRate() const override {
        return make_unique<MultiRate<PlogRate, PlogData>>();
    }

    //! Identifier of reaction rate type
    const string type() const override { return "pressure-dependent-Arrhenius"; }

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

    //! Update information specific to reaction
    /*!
     *  @param shared_data  data shared by all reactions of a given type
     */
    void updateFromStruct(const PlogData& shared_data) {
        if (shared_data.logP != logP_) {
            logP_ = shared_data.logP;
            if (logP_ > logP1_ && logP_ < logP2_) {
                return;
            }

            auto iter = pressures_.upper_bound(logP_);
            AssertThrowMsg(iter != pressures_.end(), "PlogRate::updateFromStruct",
                           "Pressure out of range: {}", logP_);
            AssertThrowMsg(iter != pressures_.begin(), "PlogRate::updateFromStruct",
                           "Pressure out of range: {}", logP_);

            // upper interpolation pressure
            logP2_ = iter->first;
            ihigh1_ = iter->second.first;
            ihigh2_ = iter->second.second;

            // lower interpolation pressure
            logP1_ = (--iter)->first;
            ilow1_ = iter->second.first;
            ilow2_ = iter->second.second;

            rDeltaP_ = 1.0 / (logP2_ - logP1_);
        }
    }

    //! Evaluate reaction rate
    /*!
     *  @param shared_data  data shared by all reactions of a given type
     */
    double evalFromStruct(const PlogData& shared_data) {
        double log_k1, log_k2;
        if (ilow1_ == ilow2_) {
            log_k1 = rates_[ilow1_].evalLog(shared_data.logT, shared_data.recipT);
        } else {
            double k = 1e-300; // non-zero to make log(k) finite
            for (size_t i = ilow1_; i < ilow2_; i++) {
                k += rates_[i].evalRate(shared_data.logT, shared_data.recipT);
            }
            log_k1 = std::log(k);
        }

        if (ihigh1_ == ihigh2_) {
            log_k2 = rates_[ihigh1_].evalLog(shared_data.logT, shared_data.recipT);
        } else {
            double k = 1e-300; // non-zero to make log(k) finite
            for (size_t i = ihigh1_; i < ihigh2_; i++) {
                k += rates_[i].evalRate(shared_data.logT, shared_data.recipT);
            }
            log_k2 = std::log(k);
        }

        return std::exp(log_k1 + (log_k2 - log_k1) * (logP_ - logP1_) * rDeltaP_);
    }

    //! Set up Plog object
    void setRates(const std::multimap<double, ArrheniusRate>& rates);

    //! Check to make sure that the rate expression is finite over a range of
    //! temperatures at each interpolation pressure. This is potentially an
    //! issue when one of the Arrhenius expressions at a particular pressure
    //! has a negative pre-exponential factor.
    void validate(const string& equation, const Kinetics& kin) override;

    //! Return the pressures and Arrhenius expressions which comprise this
    //! reaction.
    std::multimap<double, ArrheniusRate> getRates() const;

protected:
    //! log(p) to (index range) in the rates_ vector
    map<double, pair<size_t, size_t>> pressures_;

    // Rate expressions which are referenced by the indices stored in pressures_
    vector<ArrheniusRate> rates_;

    double logP_ = -1000; //!< log(p) at the current state
    double logP1_ = 1000; //!< log(p) at the lower pressure reference
    double logP2_ = -1000; //!< log(p) at the upper pressure reference

    //! Indices to the ranges within rates_ for the lower / upper pressure, such
    //! that rates_[ilow1_] through rates_[ilow2_] (inclusive) are the rates
    //! expressions which are combined to form the rate at the lower reference
    //! pressure.
    size_t ilow1_, ilow2_, ihigh1_, ihigh2_;

    double rDeltaP_ = -1.0; //!< reciprocal of (logP2 - logP1)
};

}
#endif
