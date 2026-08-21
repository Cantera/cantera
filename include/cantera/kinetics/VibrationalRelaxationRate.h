//! @file VibrationalRelaxationRate.h
//! Header for vibrational relaxation reaction rates in plasma kinetics.
//! @since New in %Cantera 4.0
//!
//! This file is part of Cantera. See License.txt in the top-level directory or
//! at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_VIBRATIONALRELAXATIONRATE_H
#define CT_VIBRATIONALRELAXATIONRATE_H

#include "cantera/kinetics/Reaction.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/base/global.h"
#include "cantera/kinetics/ReactionRate.h"
#include "cantera/kinetics/ReactionData.h"
#include "cantera/kinetics/MultiRate.h"
#include "cantera/base/Units.h"

namespace Cantera
{
//! Shared temperature data for vibrational relaxation rates.
struct VibrationalRelaxationData : public ReactionData
{
    //! Update cached temperature-dependent data directly from temperature.
    void update(double T) override;

    //! Update cached temperature-dependent data.
    /**
     * @param phase Thermodynamic phase used to retrieve the gas temperature.
     * @param kin   Kinetics object. Not used here, but required by the
     *              ReactionData interface.
     * @return `true` if the temperature has changed and rates need to be
     *         recomputed; `false` otherwise.
     */
    bool update(const ThermoPhase& phase, const Kinetics& kin) override;

    using ReactionData::update;

    double recipT13 = 1.0; //!< @f$T^(-1/3)@f$
};

//! Vibrational relaxation reaction rate class definition.
/**
 * This class provides a common implementation for several vibrational
 * relaxation models:
 * - `constant`
 * - `multi-state-resolved`
 * - `Starikovskiy`
 * - `Castela`
 * 
 * Internally, all models are mapped to the following generic expression:
 * @f[
 * k_f =
 * A \,
 * \exp \left(
 *     b \ln T
 *     + B
 *     + C T^{-1/3}
 *     + D T^{-m}
 *     + E T^{-z}
 * \right)
 * @f]
 * where `T` is the gas temperature in Kelvin (K).
 *
 * The `constant` model relaxes the vibrational species with a constant rate 
 * coefficient. It could just as well be an Arrhenius rate, but the constant 
 * model is provided for convenience and to avoid confusion with conventional 
 * Arrhenius rates for the YAML file user.
 *
 * The `multi-state-resolved` model fully resolves vibrational relaxation by
 * taking into account all vibrational species in the phase, for example
 * `N₂(v=1-8)`, and solves for their vibrational-translational (V-T) and 
 * vibrational-vibrational (V-V) relaxation. The scaling of the
 * rates is based on the Schwartz–Slawsky–Herzfeld (SSH) theory detailed 
 * in Chapter 7 of Capitelli et al. @cite capitelli2000. The simplified SSH 
 * theory implemented here is based on the harmonic oscillator approximation 
 * and can be found in equations 18 and 19 of Guerra et al. @cite guerra2019. 
 * The @f$ k_{10} @f$ rates are taken from Zhong et al. @cite zhong2023, 
 * Capitelli et al. @cite capitelli2000, and Starikovskiy and 
 * Aleksandrov @cite starikovskiy2013.
 *
 * The `Castela` model is meant to be used only for N₂ vibrational relaxation,
 * by collisions with N₂, O₂, and O exclusively. It implements the mean
 * vibrational energy relaxation model using a fictitious Cantera species and
 * is based on Castela et al. @cite castela2016.
 *
 * The `Starikovskiy` model is an extension of the Castela model to several
 * vibrational species and additional colliders. Many vibrational relaxation
 * rates can be found in Table 1 of Starikovskiy and Aleksandrov @cite starikovskiy2013, 
 * hence the model's name. The rates for the vibrational relaxation of NH₃ can be 
 * found in the reaction mechanism provided in the supplementary material of
 * Zhong et al. @cite zhong2023. More rates for the vibrational relaxation of
 * CH₄ can be found in Popov @cite popov2016.
 *
 * Unit conventions:
 *
 * - `A` uses standard Cantera rate coefficient units. Its units depend on the
 *   reaction order and are converted using the rate-coefficient units
 *   configured by `ReactionRate`.
 * - `b`, `B`, `m`, `z` are dimensionless.
 * - `C` is interpreted as @f$ K^{1/3} @f$, assuming `T` is in K.
 * - `D` is interpreted as @f$ K^{m} @f$, assuming `T` is in K.
 * - `E` is interpreted as @f$ K^{z} @f$, assuming `T` is in K.
 *
 * The coefficients `B`, `C`, `D`, `E`, `m`, `z` are read as
 * raw floating-point values. They are not converted by Cantera's unit system.
 * 
 * IMPORTANT: No internal checks are conducted to verify the physical compatibility of
 * the modelling of vibrational relaxation rates chosen by the user.
 * The `Castela` model should be used alone, exception being made for the `constant`
 * model for other vibrational species than N₂ provided that they are relaxed fast
 * enough to be considered as fast gas heating.
 * Generally, it should be avoided to mix different models together for a same ground 
 * species, except in the case of conventional reactions producing vibration, for example:
 * 
 * ```yaml
 * - equation: O(1D) + O2 => O + 0.0098321587 O2(v) + 0.9901678413 O2(a1)
 *   rate-constant: {A: 6.022e11, b: 0.0, Ea: 0.0}
 * ```
 * in which case it is impossible to know which vibrational state is generated, thus
 * requiring a description with either one of `constant`, `Castela` or `Starikovskiy`
 * for O₂(v). It is in this case possible to combine for example this mean vibrational
 * energy equation treatment for O₂ with a detailed vibration treament via `multi-state-resolved`
 * since the electron impact reactions themselves are able to discriminate between the
 * different vibrational excitations of O₂.
 * To avoid any mixup, the user is invite to thoroughly read the present documentation
 * and to consult the YAML examples. 
 * 
 * For further information on the YAML implementation of this class, please 
 * refer to [the corresponding YAML documentation section]
 * (../YAML/reactions.html#sec-yaml-vibrational-relaxation).
 *
 * @ingroup otherRateGroup
 */
class VibrationalRelaxationRate : public ReactionRate
{
public:
    
    void setParameters(const AnyMap& node, const UnitStack& rate_units) override;

    void getParameters(AnyMap& node) const override;

    //! Set context of reaction rate evaluation.
    /**
     * Vibrational relaxation rates are intended for irreversible
     * non-equilibrium plasma reactions. Reversible reactions are rejected
     * because the reverse rate cannot be obtained from conventional
     * thermochemistry for these models.
     */
    void setContext(const Reaction& rxn, const Kinetics& kin) override;

    double evalFromStruct(const VibrationalRelaxationData& shared_data) const {

        return m_A * std::exp(
            m_b * shared_data.logT
            + m_B
            + m_C * shared_data.recipT13
            + m_D * std::pow(shared_data.recipT, m_m)
            + m_E * std::pow(shared_data.recipT, m_z)
        );
    }

    // Evaluate the scaled temperature derivative.
    //
    // This returns:
    //
    // @f[
    // \frac{1}{k_f} \frac{d k_f}{dT}
    // =
    // \frac{d \ln k_f}{dT}
    // @f]
    // 
    // For the internal generic expression, this is:
    //
    // @f[
    // \frac{b}{T}
    // - \frac{C}{3} T^{-4/3}
    // - m D T^{-m-1}
    // - z E T^{-z-1}
    // @f]
    double ddTScaledFromStruct(const VibrationalRelaxationData& shared_data) const;

    void check(const string& equation) override;

    void validate(
        const string& equation, const Kinetics& kin) override;

protected:
    //! Default constructor.
    VibrationalRelaxationRate();

    //! Constructor using the internal representation for the default 
    //! model multi-state-resolved.
    /**
     * @param A       Pre-exponential factor.
     * @param B       Dimensionless constant in the exponential.
     * @param C       Coefficient multiplying @f$ T^{-1/3} @f$.
     * @param D       Coefficient multiplying @f$ T^{-m} @f$.
     * @param b       Dimensionless temperature exponent.
     */
    VibrationalRelaxationRate(double A, double B, double C, double D, double b);

    static void requireKeys(
        const AnyMap& node,
        const string& rateType,
        const string& where,
        std::initializer_list<string> keys);

    static void forbidKeys(
        const AnyMap& node,
        const string& rateType,
        const string& where,
        std::initializer_list<string> keys);

    void getPreExponentialFactor(AnyMap& rateNode) const;

    //! Pre-exponential constant for the reaction rate
    double m_A = NAN;

    //! Temperature exponent for the reaction rate.
    double m_b = NAN;

    //! Dimensionless constant in the exponential.
    double m_B = 0.0;

    //! Coefficient multiplying @f$ T^{-1/3} @f$.
    double m_C = 0.0;

    //! Coefficient multiplying @f$ T^{-m} @f$.
    double m_D = 0.0;

    //! Temperature exponent used by the D term.
    double m_m = 2.0 / 3.0;

    //! Coefficient multiplying @f$ T^{-z} @f$.
    double m_E = 0.0;

    //! Temperature exponent used by the E term.
    double m_z = 1.0;

private:

    //! Whether a negative leading coefficient is explicitly allowed.
    bool m_negativeA_ok = false;

};

}
#endif