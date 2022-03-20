//! @file PlogRate.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_PLOGRATE_H
#define CT_PLOGRATE_H

#include "cantera/kinetics/Arrhenius.h"

namespace Cantera
{

class Arrhenius2;

//! Data container holding shared data specific to PlogRate
/**
 * The data container `PlogData` holds precalculated data common to
 * all `PlogRate` objects.
 */
struct PlogData : public ReactionData
{
    PlogData() : pressure(NAN), logP(0.), m_pressure_buf(-1.) {}

    virtual void update(double T) override;

    virtual void update(double T, double P) override {
        ReactionData::update(T);
        pressure = P;
        logP = std::log(P);
    }

    virtual bool update(const ThermoPhase& phase, const Kinetics& kin) override;

    using ReactionData::update;

    //! Perturb pressure of data container
    /**
     * The method is used for the evaluation of numerical derivatives.
     * @param  deltaP  relative pressure perturbation
     */
    void perturbPressure(double deltaP);

    virtual void restore() override;

    virtual void invalidateCache() override {
        ReactionData::invalidateCache();
        pressure = NAN;
    }

    double pressure; //!< pressure
    double logP; //!< logarithm of pressure

protected:
    double m_pressure_buf; //!< buffered pressure
};


//! Pressure-dependent reaction rate expressed by logarithmically interpolating
//! between Arrhenius rate expressions at various pressures.
/*!
 * Given two rate expressions at two specific pressures:
 *
 *   * \f$ P_1: k_1(T) = A_1 T^{b_1} e^{-E_1 / RT} \f$
 *   * \f$ P_2: k_2(T) = A_2 T^{b_2} e^{-E_2 / RT} \f$
 *
 * The rate at an intermediate pressure \f$ P_1 < P < P_2 \f$ is computed as
 * \f[
 *  \log k(T,P) = \log k_1(T) + \bigl(\log k_2(T) - \log k_1(T)\bigr)
 *      \frac{\log P - \log P_1}{\log P_2 - \log P_1}
 * \f]
 * Multiple rate expressions may be given at the same pressure, in which case
 * the rate used in the interpolation formula is the sum of all the rates given
 * at that pressure. For pressures outside the given range, the rate expression
 * at the nearest pressure is used.
 */
class PlogRate final : public ReactionRate
{
public:
    //! Default constructor.
    PlogRate();

    //! Constructor from Arrhenius rate expressions at a set of pressures
    explicit PlogRate(const std::multimap<double, ArrheniusRate>& rates);

    //! Constructor using legacy Arrhenius2 framework
    explicit PlogRate(const std::multimap<double, Arrhenius2>& rates);

    PlogRate(const AnyMap& node, const UnitStack& rate_units={}) : PlogRate() {
        setParameters(node, rate_units);
    }

    unique_ptr<MultiRateBase> newMultiRate() const {
        return unique_ptr<MultiRateBase>(new MultiRate<PlogRate, PlogData>);
    }

    //! Identifier of reaction rate type
    const std::string type() const { return "pressure-dependent-Arrhenius"; }

    //! Perform object setup based on AnyMap node information
    /*!
     *  @param node  AnyMap containing rate information
     *  @param units  Unit definitions specific to rate information
     */
    void setParameters(const AnyMap& node, const UnitStack& units);

    void getParameters(AnyMap& rateNode, const Units& rate_units) const;
    void getParameters(AnyMap& rateNode) const {
        return getParameters(rateNode, Units(0));
    }

    //! Update information specific to reaction
    /*!
     *  @param shared_data  data shared by all reactions of a given type
     */
    void updateFromStruct(const PlogData& shared_data) {
        if (shared_data.logP != logP_) {
            update_C(&shared_data.logP);
        }
    }

    //! Evaluate reaction rate
    /*!
     *  @param shared_data  data shared by all reactions of a given type
     */
    double evalFromStruct(const PlogData& shared_data) {
        return updateRC(shared_data.logT, shared_data.recipT);
    }

    //! Set up Plog object
    /*!
     * @deprecated   Deprecated in Cantera 2.6. Replaced by setRates.
     */
    void setup(const std::multimap<double, Arrhenius2>& rates);

    //! Set up Plog object
    void setRates(const std::multimap<double, ArrheniusRate>& rates);

    //! Update concentration-dependent parts of the rate coefficient.
    //! @param c natural log of the pressure in Pa
    //! @deprecated To be removed after Cantera 2.6. Implementation will be moved to
    //! the updateFromStruct() method.
    void update_C(const double* c) {
        logP_ = c[0];
        if (logP_ > logP1_ && logP_ < logP2_) {
            return;
        }

        auto iter = pressures_.upper_bound(c[0]);
        AssertThrowMsg(iter != pressures_.end(), "PlogRate::update_C",
                       "Pressure out of range: {}", logP_);
        AssertThrowMsg(iter != pressures_.begin(), "PlogRate::update_C",
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

    /**
     * Update the value the rate constant.
     *
     * This function returns the actual value of the rate constant.
     * @deprecated To be removed after Cantera 2.6. Implementation will be moved to
     * the evalFromStruct() method.
     */
    double updateRC(double logT, double recipT) const {
        double log_k1, log_k2;
        if (ilow1_ == ilow2_) {
            log_k1 = rates_[ilow1_].evalLog(logT, recipT);
        } else {
            double k = 1e-300; // non-zero to make log(k) finite
            for (size_t i = ilow1_; i < ilow2_; i++) {
                k += rates_[i].evalRate(logT, recipT);
            }
            log_k1 = std::log(k);
        }

        if (ihigh1_ == ihigh2_) {
            log_k2 = rates_[ihigh1_].evalLog(logT, recipT);
        } else {
            double k = 1e-300; // non-zero to make log(k) finite
            for (size_t i = ihigh1_; i < ihigh2_; i++) {
                k += rates_[i].evalRate(logT, recipT);
            }
            log_k2 = std::log(k);
        }

        return std::exp(log_k1 + (log_k2-log_k1) * (logP_-logP1_) * rDeltaP_);
    }

    //! Check to make sure that the rate expression is finite over a range of
    //! temperatures at each interpolation pressure. This is potentially an
    //! issue when one of the Arrhenius expressions at a particular pressure
    //! has a negative pre-exponential factor.
    void validate(const std::string& equation, const Kinetics& kin) {
        validate(equation);
    }

    void validate(const std::string& equation);

    //! Return the pressures and Arrhenius expressions which comprise this
    //! reaction.
    /*!
     * @deprecated  Behavior to change after Cantera 2.6.
     *              @see getRates for new behavior.
     */
    std::vector<std::pair<double, Arrhenius2>> rates() const;

    //! Return the pressures and Arrhenius expressions which comprise this
    //! reaction.
    std::multimap<double, ArrheniusRate> getRates() const;

protected:
    //! log(p) to (index range) in the rates_ vector
    std::map<double, std::pair<size_t, size_t>> pressures_;

    // Rate expressions which are referenced by the indices stored in pressures_
    std::vector<ArrheniusRate> rates_;

    double logP_; //!< log(p) at the current state
    double logP1_, logP2_; //!< log(p) at the lower / upper pressure reference

    //! Indices to the ranges within rates_ for the lower / upper pressure, such
    //! that rates_[ilow1_] through rates_[ilow2_] (inclusive) are the rates
    //! expressions which are combined to form the rate at the lower reference
    //! pressure.
    size_t ilow1_, ilow2_, ihigh1_, ihigh2_;

    double rDeltaP_; //!< reciprocal of (logP2 - logP1)
};

}
#endif