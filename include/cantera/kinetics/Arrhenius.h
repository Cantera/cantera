/**
 * @file Arrhenius.h
 * Header for reaction rates that involve Arrhenius-type kinetics.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_ARRHENIUS_H
#define CT_ARRHENIUS_H

#include "cantera/base/ct_defs.h"
#include "cantera/base/Units.h"
#include "cantera/kinetics/ReactionData.h"
#include "ReactionRate.h"
#include "MultiRate.h"

namespace Cantera
{

class AnyValue;
class AnyMap;

//! Data container holding shared data specific to ArrheniusRate
/**
 * The data container `ArrheniusData` holds precalculated data common to
 * all `ArrheniusRate` objects.
 */
struct ArrheniusData : public ReactionData
{
    bool update(const ThermoPhase& phase, const Kinetics& kin) override;
    using ReactionData::update;
};


//! Base class for Arrhenius-type Parameterizations
/*!
 * This base class provides a minimally functional interface that allows for parameter
 * access from derived classes as well as classes that use Arrhenius-type expressions
 * internally, for example FalloffRate and PlogRate.
 * @ingroup arrheniusGroup
 */
class ArrheniusBase : public ReactionRate
{
public:
    //! Default constructor.
    ArrheniusBase() {}

    //! Constructor.
    /*!
     *  @param A  Pre-exponential factor. The unit system is (kmol, m, s); actual units
     *      depend on the reaction order and the dimensionality (surface or bulk).
     *  @param b  Temperature exponent (non-dimensional)
     *  @param Ea  Activation energy in energy units [J/kmol]
     */
    ArrheniusBase(double A, double b, double Ea);

    //! Constructor based on AnyValue content
    ArrheniusBase(const AnyValue& rate, const UnitSystem& units,
                  const UnitStack& rate_units);

    explicit ArrheniusBase(const AnyMap& node, const UnitStack& rate_units={});

    //! Perform object setup based on AnyValue node information
    /*!
     *  Used to set parameters from a child of the reaction node, which may have
     *  different names for different rate parameterizations, such as falloff rates.
     *
     *  @param rate  Child of the reaction node containing Arrhenius rate parameters.
     *      For example, the `rate-coefficient` node for a standard Arrhenius reaction.
     *  @param units  Unit system
     *  @param rate_units  Unit definitions specific to rate information
     */
    void setRateParameters(const AnyValue& rate,
                           const UnitSystem& units,
                           const UnitStack& rate_units);

    //! Get Arrhenius parameters used to populate the `rate-coefficient` or
    //! equivalent field
    void getRateParameters(AnyMap& node) const;

    void setParameters(const AnyMap& node, const UnitStack& rate_units) override;

    void getParameters(AnyMap& node) const override;

    //! Check rate expression
    void check(const string& equation) override;

    void validate(const string& equation, const Kinetics& kin) override;

    //! Return the pre-exponential factor *A* (in m, kmol, s to powers depending
    //! on the reaction order)
    /*!
     * Class specializations may provide alternate definitions that describe
     * an effective pre-exponential factor that depends on the thermodynamic state.
     */
    virtual double preExponentialFactor() const {
        return m_A;
    }

    //! Return the temperature exponent *b*
    /*!
     * Class specializations may provide alternate definitions that describe
     * an effective temperature exponent that depends on the thermodynamic state.
     */
    virtual double temperatureExponent() const {
        return m_b;
    }

    //! Return the activation energy *Ea* [J/kmol]
    //! The value corresponds to the constant specified by input parameters;
    /*!
     * Class specializations may provide alternate definitions that describe
     * an effective activation energy that depends on the thermodynamic state.
     */
    virtual double activationEnergy() const {
        return m_Ea_R * GasConstant;
    }

    //! Return reaction order associated with the reaction rate
    double order() const {
        return m_order;
    }

    //! Set units of the reaction rate expression
    void setRateUnits(const UnitStack& rate_units) override {
        ReactionRate::setRateUnits(rate_units);
        if (rate_units.size() > 1) {
            m_order = 1 - rate_units.product().dimension("quantity");
        } else {
            m_order = NAN;
        }
    }

    //! Get flag indicating whether negative A values are permitted
    bool allowNegativePreExponentialFactor() const {
        return m_negativeA_ok;
    }

    //! Set flag indicating whether negative A values are permitted
    void setAllowNegativePreExponentialFactor(bool value) {
        m_negativeA_ok = value;
    }

protected:
    bool m_negativeA_ok = false; //!< Permissible negative A values
    double m_A = NAN; //!< Pre-exponential factor
    double m_b = NAN; //!< Temperature exponent
    double m_Ea_R = 0.; //!< Activation energy (in temperature units)
    double m_E4_R = 0.; //!< Optional 4th energy parameter (in temperature units)
    double m_logA = NAN; //!< Logarithm of pre-exponential factor
    double m_order = NAN; //!< Reaction order
    string m_A_str = "A"; //!< The string for the pre-exponential factor
    string m_b_str = "b"; //!< The string for temperature exponent
    string m_Ea_str = "Ea"; //!< The string for activation energy
    string m_E4_str = ""; //!< The string for an optional 4th parameter
};

//! Arrhenius reaction rate type depends only on temperature
/*!
 * A reaction rate coefficient of the following form.
 *
 *   @f[
 *        k_f =  A T^b \exp (-Ea/RT)
 *   @f]
 *
 * @ingroup arrheniusGroup
 */
class ArrheniusRate : public ArrheniusBase
{
public:
    using ArrheniusBase::ArrheniusBase; // inherit constructors

    unique_ptr<MultiRateBase> newMultiRate() const override {
        return make_unique<MultiRate<ArrheniusRate, ArrheniusData>>();
    }

    const string type() const override {
        return "Arrhenius";
    }

    //! Evaluate reaction rate
    double evalRate(double logT, double recipT) const {
        return m_A * std::exp(m_b * logT - m_Ea_R * recipT);
    }

    //! Evaluate natural logarithm of the rate constant.
    double evalLog(double logT, double recipT) const {
        return m_logA + m_b * logT - m_Ea_R * recipT;
    }

    //! Evaluate reaction rate
    /*!
     *  @param shared_data  data shared by all reactions of a given type
     */
    double evalFromStruct(const ArrheniusData& shared_data) const {
        return m_A * std::exp(m_b * shared_data.logT - m_Ea_R * shared_data.recipT);
    }

    //! Evaluate derivative of reaction rate with respect to temperature
    //! divided by reaction rate
    /*!
     *  @param shared_data  data shared by all reactions of a given type
     */
    double ddTScaledFromStruct(const ArrheniusData& shared_data) const {
        return (m_Ea_R * shared_data.recipT + m_b) * shared_data.recipT;
    }
};

}

#endif
