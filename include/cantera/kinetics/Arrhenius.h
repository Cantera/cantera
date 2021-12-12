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

/**
 *  @defgroup arrheniusGroup  Arrhenius-type Parameterizations
 *
 *  This section describes the parameterizations used to describe the standard
 *  Arrhenius rate parameterization and derived models.
 *
 *  @ingroup chemkinetics
 */

//! Base class for Arrhenius-type Parameterizations
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

    //! Perform object setup based on AnyValue node information
    /*!
     *  @param rate  AnyValue containing rate information
     *  @param units  Unit system
     *  @param rate_units  Unit definitions specific to rate information
     */
    virtual void setRateParameters(const AnyValue& rate,
                                   const UnitSystem& units,
                                   const UnitStack& rate_units);

    virtual void getParameters(AnyMap& node) const;

    //! Return parameters - required for legacy framework
    //! @todo: merge with single-parameter version after removal of old framework
    virtual void getParameters(AnyMap& node, const Units& rate_units) const;

    void check(const std::string& equation, const AnyMap& node);

    //! Return the pre-exponential factor *A* (in m, kmol, s to powers depending
    //! on the reaction order)
    double preExponentialFactor() const {
        return m_A;
    }

    //! Return the temperature exponent *b*
    double temperatureExponent() const {
        return m_b;
    }

    //! Return the intrinsic activation energy *Ea* [J/kmol]
    //! The value corresponds to the constant specified by input parameters;
    //! in the simplest case, this is equal to the activation energy itself
    double intrinsicActivationEnergy() const {
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
    double m_order; //!< Reaction order
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
 * @todo supersedes Arrhenius: replace existing instances with this class, rename,
 *      and deprecate Arrhenius3.
 *
 * @ingroup arrheniusGroup
 */
class ArrheniusRate final : public ArrheniusBase
{
public:
    ArrheniusRate() = default;
    using ArrheniusBase::ArrheniusBase; // inherit constructors

    //! Constructor based on AnyMap content
    ArrheniusRate(const AnyMap& node, const UnitStack& rate_units={}) {
        setParameters(node, rate_units);
    }

    unique_ptr<MultiRateBase> newMultiRate() const override {
        return unique_ptr<MultiRateBase>(new MultiBulkRate<ArrheniusRate, ArrheniusData>);
    }

    //! Identifier of reaction rate type
    virtual const std::string type() const override {
        return "Arrhenius";
    }

    //! Perform object setup based on AnyMap node information
    /*!
     *  @param node  AnyMap containing rate information
     *  @param rate_units  Unit definitions specific to rate information
     */
    virtual void setParameters(const AnyMap& node, const UnitStack& rate_units) override;

    virtual void getParameters(AnyMap& node) const override;

    //! Evaluate reaction rate
    //! @param shared_data  data shared by all reactions of a given type
    double evalFromStruct(const ReactionData& shared_data) const {
        return m_A * std::exp(m_b * shared_data.logT - m_Ea_R * shared_data.recipT);
    }

    //! Evaluate derivative of reaction rate with respect to temperature
    //! divided by reaction rate
    /*!
     *  @param shared_data  data shared by all reactions of a given type
     */
    virtual double ddTScaledFromStruct(const ReactionData& shared_data) const {
        return (m_Ea_R * shared_data.recipT + m_b) * shared_data.recipT;
    }

    //! Return the activation energy *Ea* [J/kmol]
    double activationEnergy() const {
        return m_Ea_R * GasConstant;
    }

    //! Return the activation energy divided by the gas constant (i.e. the
    //! activation temperature) [K]
    double activationEnergy_R() const {
        return m_Ea_R;
    }
};


//! Blowers Masel reaction rate type depends on the enthalpy of reaction
/**
 * The Blowers Masel approximation is written by Paul Blowers,
 * Rich Masel (DOI: https://doi.org/10.1002/aic.690461015) to
 * adjust the activation energy based on enthalpy change of a reaction:
 *
 *   \f{eqnarray*}{
 *        E_a &=& 0\; \text{if }\Delta H < -4E_0 \\
 *        E_a &=& \Delta H\; \text{if }\Delta H > 4E_0 \\
 *        E_a &=& \frac{(w + \Delta H / 2)(V_P - 2w +
 *               \Delta H)^2}{(V_P^2 - 4w^2 + (\Delta H)^2)}\; \text{Otherwise}
 *   \f}
 * where
 *   \f[
 *        V_P = \frac{2w (w + E_0)}{w - E_0},
 *   \f]
 * \f$ w \f$ is the average bond dissociation energy of the bond breaking
 * and that being formed in the reaction. Since the expression is
 * very insensitive to \f$ w \f$ for \f$ w >= 2 E_0 \f$, \f$ w \f$
 * can be approximated to an arbitrary high value like 1000 kJ/mol.
 *
 * After the activation energy is determined by Blowers-Masel approximation,
 * it can be plugged into Arrhenius function to calculate the rate constant.
 *   \f[
 *        k_f =  A T^b \exp (-E_a/RT)
 *   \f]
 *
 * @ingroup arrheniusGroup
 */
class BlowersMaselRate : public ArrheniusBase
{
public:
    //! Default constructor.
    BlowersMaselRate();

    //! Constructor.
    /*!
     *  @param A  Pre-exponential factor. The unit system is (kmol, m, s); actual units
     *      depend on the reaction order and the dimensionality (surface or bulk).
     *  @param b  Temperature exponent (non-dimensional)
     *  @param Ea0  Intrinsic activation energy in energy units [J/kmol]
     *  @param w  Average bond dissociation energy of the bond being formed and
     *      broken in the reaction, in energy units [J/kmol]
     */
    BlowersMaselRate(double A, double b, double Ea0, double w);

    unique_ptr<MultiRateBase> newMultiRate() const {
        return unique_ptr<MultiRateBase>(
            new MultiBulkRate<BlowersMaselRate, BlowersMaselData>);
    }

    //! Constructor based on AnyMap content
    BlowersMaselRate(const AnyMap& node, const UnitStack& rate_units={})
        : BlowersMaselRate()
    {
        setParameters(node, rate_units);
    }

    //! Identifier of reaction rate type
    virtual const std::string type() const {
        return "Blowers-Masel";
    }

    virtual void setRateParameters(const AnyValue& rate,
                                   const UnitSystem& units,
                                   const UnitStack& rate_units);

    //! Perform object setup based on AnyMap node information
    /*!
     *  @param node  AnyMap containing rate information
     *  @param rate_units  Unit definitions specific to rate information
     */
    virtual void setParameters(const AnyMap& node, const UnitStack& rate_units);

    virtual void getParameters(AnyMap& node) const;

    double evalFromStruct(const BlowersMaselData& shared_data) {
        double deltaH_R;
        if (shared_data.ready) {
            deltaH_R = shared_data.dH[m_rate_index] / GasConstant;
        } else {
            deltaH_R = shared_data.dH[0] / GasConstant;
        }
        double Ea_R = activationEnergy_R(deltaH_R);
        return m_A * std::exp(m_b * shared_data.logT - Ea_R * shared_data.recipT);
    }

    //! Return the actual activation energy (a function of the delta H of reaction)
    //! divided by the gas constant (i.e. the activation temperature) [K]
    double activationEnergy_R(double deltaH_R) const {
        if (deltaH_R < -4 * m_Ea_R) {
            return 0.;
        }
        if (deltaH_R > 4 * m_Ea_R) {
            return deltaH_R;
        }
        double vp = 2 * m_w_R * ((m_w_R + m_Ea_R) / (m_w_R - m_Ea_R)); // in Kelvin
        double vp_2w_dH = (vp - 2 * m_w_R + deltaH_R); // (Vp - 2 w + dH)
        return (m_w_R + deltaH_R / 2) * (vp_2w_dH * vp_2w_dH) /
            (vp * vp - 4 * m_w_R * m_w_R + deltaH_R * deltaH_R); // in Kelvin
    }

    //! Return the actual activation energy [J/kmol]
    /*!
     *  @param deltaH  Enthalpy change of reaction [J/kmol]
     */
    double activationEnergy(double deltaH) const {
        return activationEnergy_R(deltaH / GasConstant) * GasConstant;
    }

    //! Return the intrinsic activation energy [J/kmol]
    double activationEnergy0() const {
        return m_Ea_R * GasConstant;
    }

    //! Return the bond dissociation energy *w* [J/kmol]
    double bondEnergy() const {
        return m_w_R * GasConstant;
    }

protected:
    double m_w_R; //!< Bond dissociation energy (in temperature units)
};

}


#endif
