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
/*!
 * This base class provides a minimally functional interface that allows for parameter
 * access from derived classes as well as classes that use Arrhenius-type expressions
 * internally, for example FalloffRate and PlogRate.
 *
 * @todo supersedes Arrhenius2 and will replace Arrhenius(2) after Cantera 2.6,
 *      The class should be renamed to Arrhenius after removal of Arrhenius2. The new
 *      behavior can be forced in self-compiled Cantera installations by defining
 *      CT_NO_LEGACY_REACTIONS_26 via the 'no_legacy_reactions' option in SCons.
 */
class ArrheniusBase
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
    void setRateParameters(const AnyValue& rate,
                           const UnitSystem& units,
                           const UnitStack& rate_units);

    //! Return parameters
    void getRateParameters(AnyMap& node) const;

    //! Check rate expression
    void checkRate(const std::string& equation, const AnyMap& node);

    //! Evaluate reaction rate
    /*!
     *  @internal  Non-virtual method that should not be overloaded
     */
    double evalRate(double logT, double recipT) const {
        return m_A * std::exp(m_b * logT - m_Ea_R * recipT);
    }

    //! Evaluate natural logarithm of the rate constant.
    /*!
     *  @internal  Non-virtual method that should not be overloaded
     */
    double evalLog(double logT, double recipT) const {
        return m_logA + m_b * logT - m_Ea_R * recipT;
    }

    //! Return the pre-exponential factor *A* (in m, kmol, s to powers depending
    //! on the reaction order)
    double preExponentialFactor() const {
        return m_A;
    }

    //! Return the temperature exponent *b*
    double temperatureExponent() const {
        return m_b;
    }

    //! Return the activation energy *Ea* [J/kmol]
    //! The value corresponds to the constant specified by input parameters;
    //! class specializations may provide alternate definitions that describe
    //! an effective activation energy that depends on the thermodynamic state.
    double activationEnergy() const {
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
 */
class ArrheniusRate final : public ArrheniusBase, public ReactionRate
{
public:
    ArrheniusRate() = default;
    using ArrheniusBase::ArrheniusBase; // inherit constructors

    //! Constructor based on AnyMap content
    ArrheniusRate(const AnyMap& node, const UnitStack& rate_units={}) {
        setParameters(node, rate_units);
    }

    unique_ptr<MultiRateBase> newMultiRate() const override {
        return unique_ptr<MultiRateBase>(new MultiRate<ArrheniusRate, ArrheniusData>);
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

    void check(const std::string& equation, const AnyMap& node) override {
        checkRate(equation, node);
    }

    //! Evaluate reaction rate
    /*!
     *  @param shared_data  data shared by all reactions of a given type
     */
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

    //! Return the activation energy divided by the gas constant (i.e. the
    //! activation temperature) [K]
    double activationEnergy_R() const {
        return m_Ea_R;
    }
};


//! Two temperature plasma reaction rate type depends on both
//! gas temperature and electron temperature.
/*!
 * The form of the two temperature plasma reaction rate coefficient is similar to an
 * Arrhenius reaction rate coefficient. The temperature exponent (b) is applied to
 * the electron temperature instead. In addition, the exponential term with
 * activation energy for electron is included.
 *
 *   \f[
 *        k_f =  A T_e^b \exp (-E_a/RT) \exp (-E_{a,e}/RT_e)
 *   \f]
 *
 * Ref.: Kossyi, I. A., Kostinsky, A. Y., Matveyev, A. A., & Silakov, V. P. (1992).
 * Kinetic scheme of the non-equilibrium discharge in nitrogen-oxygen mixtures.
 * Plasma Sources Science and Technology, 1(3), 207.
 * doi: 10.1088/0963-0252/1/3/011
 * @ingroup arrheniusGroup
 */
class TwoTempPlasmaRate final : public ArrheniusBase, public ReactionRate
{
public:
    TwoTempPlasmaRate();

    //! Constructor.
    /*!
     *  @param A  Pre-exponential factor. The unit system is (kmol, m, s); actual units
     *      depend on the reaction order and the dimensionality (surface or bulk).
     *  @param b  Temperature exponent (non-dimensional)
     *  @param Ea  Activation energy in energy units [J/kmol]
     *  @param EE  Activation electron energy in energy units [J/kmol]
     */
    TwoTempPlasmaRate(double A, double b, double Ea=0.0, double EE=0.0);

    unique_ptr<MultiRateBase> newMultiRate() const override {
        return unique_ptr<MultiRateBase>(
            new MultiRate<TwoTempPlasmaRate, TwoTempPlasmaData>);
    }

    //! Constructor based on AnyMap content
    TwoTempPlasmaRate(const AnyMap& node, const UnitStack& rate_units={})
        : TwoTempPlasmaRate()
    {
        setParameters(node, rate_units);
    }

    //! Identifier of reaction rate type
    virtual const std::string type() const override {
        return "two-temperature-plasma";
    }

    //! Perform object setup based on AnyMap node information
    /*!
     *  @param node  AnyMap containing rate information
     *  @param rate_units  Unit definitions specific to rate information
     */
    virtual void setParameters(const AnyMap& node, const UnitStack& rate_units) override;

    virtual void getParameters(AnyMap& node) const override;

    void check(const std::string& equation, const AnyMap& node) override {
        checkRate(equation, node);
    }

    //! Evaluate reaction rate
    /*!
     *  @param shared_data  data shared by all reactions of a given type
     */
    double evalFromStruct(const TwoTempPlasmaData& shared_data) {
        // m_E4_R is the electron activation (in temperature units)
        return m_A * std::exp(m_b * shared_data.logTe -
                              m_Ea_R * shared_data.recipT +
                              m_E4_R * (shared_data.electronTemp - shared_data.temperature)
                              * shared_data.recipTe * shared_data.recipT);
    }

    //! Return the activation energy divided by the gas constant (i.e. the
    //! activation temperature) [K]
    double activationEnergy_R() const {
        return m_Ea_R;
    }

    //! Return the electron activation energy *Ea* [J/kmol]
    double activationElectronEnergy() const {
        return m_E4_R * GasConstant;
    }

    //! Return the electron activation energy divided by the gas constant (i.e. the
    //! activation temperature) [K]
    double activationElectronEnergy_R() const {
        return m_E4_R;
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
class BlowersMaselRate final : public ArrheniusBase, public ReactionRate
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

    unique_ptr<MultiRateBase> newMultiRate() const override {
        return unique_ptr<MultiRateBase>(
            new MultiRate<BlowersMaselRate, BlowersMaselData>);
    }

    //! Constructor based on AnyMap content
    BlowersMaselRate(const AnyMap& node, const UnitStack& rate_units={})
        : BlowersMaselRate()
    {
        setParameters(node, rate_units);
    }

    //! Identifier of reaction rate type
    virtual const std::string type() const override {
        return "Blowers-Masel";
    }

    //! Perform object setup based on AnyMap node information
    /*!
     *  @param node  AnyMap containing rate information
     *  @param rate_units  Unit definitions specific to rate information
     */
    virtual void setParameters(
        const AnyMap& node, const UnitStack& rate_units) override;

    virtual void getParameters(AnyMap& node) const override;

    void check(const std::string& equation, const AnyMap& node) override {
        checkRate(equation, node);
    }

    virtual void setContext(const Reaction& rxn, const Kinetics& kin) override;

    //! Update information specific to reaction
    /*!
     *  @param shared_data  data shared by all reactions of a given type
     */
    void updateFromStruct(const BlowersMaselData& shared_data) {
        if (shared_data.ready) {
            m_deltaH_R = 0.;
            for (const auto& item : m_stoich_coeffs) {
                m_deltaH_R += shared_data.grt[item.first] * item.second;
            }
            m_deltaH_R /= GasConstant;
        } else {
            m_deltaH_R = shared_data.dH_direct / GasConstant;
        }
    }

    //! Evaluate reaction rate
    /*!
     *  @param shared_data  data shared by all reactions of a given type
     */
    double evalFromStruct(const BlowersMaselData& shared_data) {
        double Ea_R = activationEnergy_R(m_deltaH_R);
        return m_A * std::exp(m_b * shared_data.logT - Ea_R * shared_data.recipT);
    }

    //! Evaluate derivative of reaction rate with respect to temperature
    //! divided by reaction rate
    /*!
     *  This method is used to override the numerical derivative, which does not
     *  consider potential changes due to a changed reaction enthalpy. A corresponding
     *  warning is raised.
     *  @param shared_data  data shared by all reactions of a given type
     */
    virtual double ddTScaledFromStruct(const BlowersMaselData& shared_data) const;

    //! Return the effective activation energy (a function of the delta H of reaction)
    //! divided by the gas constant (i.e. the activation temperature) [K]
    double activationEnergy_R(double deltaH_R) const {
        if (deltaH_R < -4 * m_Ea_R) {
            return 0.;
        }
        if (deltaH_R > 4 * m_Ea_R) {
            return deltaH_R;
        }
        // m_E4_R is the bond dissociation energy "w" (in temperature units)
        double vp = 2 * m_E4_R * ((m_E4_R + m_Ea_R) / (m_E4_R - m_Ea_R)); // in Kelvin
        double vp_2w_dH = (vp - 2 * m_E4_R + deltaH_R); // (Vp - 2 w + dH)
        return (m_E4_R + deltaH_R / 2) * (vp_2w_dH * vp_2w_dH) /
            (vp * vp - 4 * m_E4_R * m_E4_R + deltaH_R * deltaH_R); // in Kelvin
    }

    //! Return the effective activation energy [J/kmol]
    /*!
     *  @param deltaH  Enthalpy change of reaction [J/kmol]
     */
    double effectiveActivationEnergy(double deltaH) const {
        return activationEnergy_R(deltaH / GasConstant) * GasConstant;
    }

    //! Return the bond dissociation energy *w* [J/kmol]
    double bondEnergy() const {
        return m_E4_R * GasConstant;
    }

protected:
    //! Pairs of species indices and multiplers to calculate enthalpy change
    std::vector<std::pair<size_t, double>> m_stoich_coeffs;

    double m_deltaH_R; //!< enthalpy change of reaction (in temperature units)
};

}


#endif
