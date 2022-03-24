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
    virtual bool update(const ThermoPhase& phase, const Kinetics& kin);
    using ReactionData::update;
};


/**
 *  @defgroup arrheniusGroup  Arrhenius-type Parameterizations
 *
 *  This section describes the parameterizations used to describe the standard
 *  Arrhenius rate parameterization and derived models.
 *
 *  @ingroup chemkinetics
 */

//! Base class for Arrhenius-type Parameterizations
/*!
 * This base class provides a minimally functional interface that allows for parameter
 * access from derived classes as well as classes that use Arrhenius-type expressions
 * internally, for example FalloffRate and PlogRate.
 */
class ArrheniusBase : public ReactionRate
{
public:
    //! Default constructor.
    ArrheniusBase();

    //! Constructor.
    /*!
     *  @param A  Pre-exponential factor. The unit system is (kmol, m, s); actual units
     *      depend on the reaction order and the dimensionality (surface or bulk).
     *  @param b  Temperature exponent (non-dimensional)
     *  @param Ea  Activation energy in energy units [J/kmol]
     */
    ArrheniusBase(double A, double b, double Ea);

    //! Constructor based on AnyValue content
    ArrheniusBase(const AnyValue& rate,
                  const UnitSystem& units, const UnitStack& rate_units)
    {
        setRateParameters(rate, units, rate_units);
    }

    explicit ArrheniusBase(const AnyMap& node, const UnitStack& rate_units={})
        : ArrheniusBase()
    {
        setParameters(node, rate_units);
    }

    //! Perform object setup based on AnyValue node information
    /*!
     *  @param rate  AnyValue containing rate information
     *  @param units  Unit system
     *  @param rate_units  Unit definitions specific to rate information
     */
    void setRateParameters(const AnyValue& rate,
                           const UnitSystem& units,
                           const UnitStack& rate_units);

    //! Return parameters
    void getRateParameters(AnyMap& node) const;

    virtual void setParameters(
        const AnyMap& node, const UnitStack& rate_units) override;

    virtual void getParameters(AnyMap& node) const override;

    //! Check rate expression
    virtual void check(const std::string& equation, const AnyMap& node) override;

    virtual void validate(const std::string& equation, const Kinetics& kin) override;

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

    // Return units of the reaction rate expression
    const Units& rateUnits() const {
        return m_rate_units;
    }

    //! Return reaction order associated with the reaction rate
    double order() const {
        return m_order;
    }

    //! Set units of the reaction rate expression
    void setRateUnits(const UnitStack& rate_units) {
        if (rate_units.size() > 1) {
            m_rate_units = rate_units.product();
            m_order = 1 - m_rate_units.dimension("quantity");
        } else {
            m_order = NAN;
            m_rate_units = rate_units.standardUnits();
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
    bool m_negativeA_ok; //!< Flag indicating whether negative A values are permitted
    double m_A; //!< Pre-exponential factor
    double m_b; //!< Temperature exponent
    double m_Ea_R; //!< Activation energy (in temperature units)
    double m_E4_R; //!< Optional 4th energy parameter (in temperature units)
    double m_logA; //!< Logarithm of pre-exponential factor
    double m_order; //!< Reaction order
    std::string m_A_str = "A"; //!< The string for temperature exponent
    std::string m_b_str = "b"; //!< The string for temperature exponent
    std::string m_Ea_str = "Ea"; //!< The string for activation energy
    std::string m_E4_str = ""; //!< The string for an optional 4th parameter
    Units m_rate_units; //!< Reaction rate units
};

//! Arrhenius reaction rate type depends only on temperature
/*!
 * A reaction rate coefficient of the following form.
 *
 *   \f[
 *        k_f =  A T^b \exp (-Ea/RT)
 *   \f]
 *
 * @ingroup arrheniusGroup
 *
 * @todo supersedes Arrhenius2 and will replace Arrhenius after Cantera 2.6. The new
 *      behavior can be forced in self-compiled Cantera installations by defining
 *      CT_NO_LEGACY_REACTIONS_26 via the 'no_legacy_reactions' option in SCons.
 */
class ArrheniusRate : public ArrheniusBase
{
public:
    using ArrheniusBase::ArrheniusBase; // inherit constructors

    unique_ptr<MultiRateBase> newMultiRate() const override {
        return unique_ptr<MultiRateBase>(new MultiRate<ArrheniusRate, ArrheniusData>);
    }

    virtual const std::string type() const override {
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
