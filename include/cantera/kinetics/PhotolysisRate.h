/**
 * @file PhotolysisRate.h
 * This header file includes definitions to implement Photolysis reactions.
 * The photolysis rate, J, can be parameterized in two ways. A constant J value or
 * with the equation: \f$ J(x) = l(\cos x)^m e^{-n \sec x } \f$ described by
 * Saunders et al., 2003. This parameterization is chosen as the
 * Master Chemical Mechanism (MCM) is one of the most extensive atmospheric models and
 * uses this parameterization. So it makes model translation into Cantera format
 * simpler.
 *
 *  Saunders, S. M., Jenkin, M. E., Derwent, R. G., & Pilling, M. J. (2003). Protocol
 * for the development of the Master Chemical Mechanism, MCM v3 (Part A): tropospheric
 * degradation of non-aromatic volatile organic compounds. Atmospheric Chemistry and
 * Physics, 3(1), 161-180.
 *
 * @since  New in Cantera 3.0.
 **/

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_PHOTOLYSISRATE_H
#define CT_PHOTOLYSISRATE_H

#include "cantera/base/ct_defs.h"
#include "cantera/base/Units.h"
#include "cantera/kinetics/ReactionData.h"
#include "ReactionRate.h"
#include "MultiRate.h"

namespace Cantera
{

class AnyValue;
class AnyMap;

//! Data container holding shared data specific to PhotolysisRate
/**
 * The data container `PhotolysisData` holds precalculated data common to
 * all `PhotolysisRate` objects.
 */
struct PhotolysisData : public ReactionData
{
    double X = NAN; //!< Zenith angle initially set to NAN
    virtual bool update(const ThermoPhase& phase, const Kinetics& kin);
    using ReactionData::update;
};

//! Base class for Photolysis parameterizations
/*!
 * This base class provides a minimally functional interface that allows for parameter
 * access from derived classes as well as classes that use Photolysis expressions.
 *
 */
class PhotolysisRate : public ReactionRate
{
public:
    //! Default constructor.
    PhotolysisRate() {}

    //! Constructor with constant parameterization.
    /*!
     *  @param J constant photolysis rate to be used
     */
    PhotolysisRate(double J);

    //! Constructor with MCM parameterization.
    /*!
     *  @param l  First parameter in photolysis parameterization.
     *  @param m  Second parameter in photolysis parameterization
     *  @param n  Third parameter in photolysis parameterization
     */
    PhotolysisRate(double l, double m, double n);

    //! Constructor based on AnyValue content
    PhotolysisRate(const AnyValue& rate, const UnitSystem& units,
                   const UnitStack& rate_units);

    explicit PhotolysisRate(const AnyMap& node, const UnitStack& rate_units={});

    //! Perform object setup based on AnyValue node information
    /*!
     *  Used to set parameters from a child of the reaction node, which may have
     *  different names for different rate parameterizations, such as falloff rates.
     *
     *  @param rate  Child of the reaction node containing rate parameters.
     *      For example, the `rate-coefficient` node for a standard Photolysis reaction.
     *  @param units  Unit system
     *  @param rate_units  Unit definitions specific to rate information
     */
    void setRateParameters(const AnyValue& rate,
                           const UnitSystem& units,
                           const UnitStack& rate_units);

    //! Get parameters used to populate the `rate-coefficient` or
    //! equivalent field
    void getRateParameters(AnyMap& node) const;

    virtual void setParameters(
        const AnyMap& node, const UnitStack& rate_units) override;

    virtual void getParameters(AnyMap& node) const override;

    //! Check rate expression
    virtual void check(const std::string& equation) override;

    virtual void validate(const std::string& equation, const Kinetics& kin) override;

    //! Return the parameter l in the photolysis parameterization
    virtual double l_param() const {
        return m_l;
    }

    //! Return the parameter m in the photolysis parameterization
    virtual double m_param() const {
        return m_m;
    }

    //! Return the parameter n in the photolysis parameterization
    virtual double n_param() const {
        return m_n;
    }

    unique_ptr<MultiRateBase> newMultiRate() const override {
        return make_unique<MultiRate<PhotolysisRate, PhotolysisData>>();
    }

    virtual const std::string type() const override {
        return "photolysis";
    }

    //! Evaluate reaction rate
    double evalRate(double X) const {
        if (!m_rate_evaluated && !std::isnan(X)) {
            double J = m_l * std::pow(std::cos(X), m_m);
            J *= std::exp(-m_n * 1 / std::cos(X));
            return J;
        } else {
            return m_J;
        }
    }

    //! Evaluate reaction rate
    /*!
     *  @param shared_data  data shared by all reactions of a given type
     */
    double evalFromStruct(const PhotolysisData& shared_data) const {
        double X = shared_data.X;
        if (!m_rate_evaluated && !std::isnan(X)) {
            double J = m_l * std::pow(std::cos(X), m_m);
            J *= std::exp(-m_n * 1 / std::cos(X));
            return J;
        } else {
            return m_J;
        }
    }

    //! Evaluate derivative of reaction rate with respect to temperature
    //! divided by reaction rate
    /*!
     *  @param shared_data  data shared by all reactions of a given type
     */
    double ddTScaledFromStruct(const PhotolysisData& shared_data) const {
        return 0;
    }

protected:
    double m_rate_evaluated = false; //!< a flag indicating if the rate needs evaluated
    double m_J = 0; //!< Photolysis rate
    double m_l = NAN; //!< First parameter in photolysis parameterization
    double m_m = NAN; //!< Second parameter in photolysis parameterization
    double m_n = NAN; //!< Third parameter in photolysis parameterization
    std::string m_J_str = "J"; //!< Parameter string for J
    std::string m_l_str = "l"; //!< Parameter string for l
    std::string m_m_str = "m"; //!< Parameter string for m
    std::string m_n_str = "n"; //!< Parameter string for n
};

}

#endif
