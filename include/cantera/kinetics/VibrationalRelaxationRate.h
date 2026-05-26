//! @file VibrationalRelaxationRate.h
//! Header for vibrational relaxation reaction rates in plasma kinetics.
//!
//! This file is part of Cantera. See License.txt in the top-level directory or
//! at https://cantera.org/license.txt for license and copyright information.


//! It implements the relaxation rate of vibrationally excited species in plasma kinetics.
//! Four models are currently supported by this class: 'constant', 'multi-state-resolved', 'starikovskiy', and 'castela'.

//! The 'constant' model relaxes the vibrational species with a constant rate coefficient. It could just as well be an
//! Arrhenius rate, but the constant model is provided for convenience and to avoid confusion with conventional Arrhenius rates.

//! The 'multi-state-resolved' model is meant to fully resolve vibrational relaxation by taking into account all vibrational species
//! in the phase (Ex: N2(v=1-8)) and solves for their V-T and V-V relaxation. The scaling of the rates are based on the SSH theory
//! detailed in Chapter 7 of the following book:
//! Capitelli, M., Ferreira, C. M., Gordiets, B. F., & Osipov, A. I. (2013). 
//! Plasma kinetics in atmospheric gases (Vol. 31). Springer Science & Business Media.
//! The simplified version of the SSH theory implemented in this model is based on the harmonic oscillator approximation and can be found
//! in the following paper (equations 18 and 19):
//! Vasco Guerra et al 2019 Plasma Sources Sci. Technol. 28 073001
//! The k10 rates are taken from //! Zhong, H. et al. (2023) “Understanding non-equilibrium N2O/NOx chemistry in plasma-assisted 
//! low-temperature NH3 oxidation,” Combustion and Flame, 256, which itself took those rates from both the Capitelli book and 
//! the article from Strikovskiy and Aleksandrov (2013), see more below.

//! The 'castela' model is meant to be used only for N2 vibrational relaxation, by collisions with N2, O2 and O exclusively.
//! It implements the mean vibrational energy relaxation model using a fictive cantera species.
//! The Castela model is based on the following paper:
//! Castela, M., Fiorina, B., Coussement, A., Gicquel, O., Darabiha, N., & Laux, C. O. (2016). 
//! Modelling the impact of non-equilibrium discharges on reactive mixtures for simulations of 
//! plasma-assisted ignition in turbulent flows. Combustion and flame, 166, 133-147.

//! The so-called 'starikovskiy' model is an extension of the Castela model to several vibrational species and additional colliders.
//! A lot of vibrational relaxation rates can be found in Table 1 of the following paper:
//! Starikovskiy, A., & Aleksandrov, N. (2013). Plasma-assisted ignition and combustion. 
//! Progress in Energy and Combustion Science, 39(1), 61-110.
//! The rates for the vibrational relaxation of NH3 can be found in the reaction mechanism provided in supplementary material
//! of the following paper:
//! Zhong, H. et al. (2023) “Understanding non-equilibrium N2O/NOx chemistry in plasma-assisted low-temperature NH3 oxidation,” 
//! Combustion and Flame, 256. 
//! More rates for the vibrational relaxation of CH$ can be found in the following paper:
//! Popov, N.A. (2016) “Kinetics of plasma-assisted combustion: Effect of non-equilibrium excitation on the ignition 
//! and oxidation of combustible mixtures,” Plasma Sources Science and Technology. Institute of Physics Publishing. 



#ifndef CT_VIBRATIONALRELAXATIONRATE_H
#define CT_VIBRATIONALRELAXATIONRATE_H

#include "Arrhenius.h"

#include <cmath>

namespace Cantera
{

//! Shared temperature data for vibrational relaxation rates.
struct DetailedVibData : public ReactionData
{
    //! Update cached temperature-dependent data.
    /*!
     * @param phase Thermodynamic phase used to retrieve the gas temperature.
     * @param kin   Kinetics object. Not used here, but required by the
     *              ReactionData interface.
     *
     * @return `true` if the temperature has changed and rates need to be
     *         recomputed; `false` otherwise.
     */
    bool update(const ThermoPhase& phase, const Kinetics& kin) override;

    using ReactionData::update;
};


//! Vibrational relaxation reaction rate.
/*!
 * This class provides a common implementation for several vibrational
 * relaxation models:
 *
 * - `constant`
 * - `multi-state-resolved`
 * - `starikovskiy`
 * - `castela`
 *
 * Internally, all models are mapped to the following generic expression:
 *
 * @f[
 * k_f =
 * scaling \, A \,
 * \exp \left(
 *     b \ln T
 *     + B
 *     + C T^{-1/3}
 *     + D T^{-m}
 *     + E T^{-z}
 * \right)
 * @f]
 *
 * where `T` is the gas temperature in K.
 *
 * The YAML reaction type is:
 *
 * @code{.yaml}
 * type: vibrational-relaxation
 * @endcode
 *
 * The selected physical model is specified separately using:
 *
 * @code{.yaml}
 * vibration_model: constant
 * vibration_model: multi-state-resolved
 * vibration_model: starikovskiy
 * vibration_model: castela
 * @endcode
 *
 * Unit conventions:
 *
 * - `A` uses standard Cantera rate coefficient units. Its units depend on the
 *   reaction order and are converted by `ArrheniusBase`.
 * - `b`, `B`, `m`, `z`, and `scaling` are dimensionless.
 * - `C` is interpreted as K^(1/3), assuming `T` is in K.
 * - `D` is interpreted as K^m, assuming `T` is in K.
 * - `E` is interpreted as K^z, assuming `T` is in K.
 *
 * The coefficients `B`, `C`, `D`, `E`, `m`, `z`, and `scaling` are read as
 * raw floating-point values. They are not converted by Cantera's unit system.
 *
 * @ingroup arrheniusGroup
 */
class VibrationalRelaxationRate : public ArrheniusBase
{
public:
    //! Default constructor.
    VibrationalRelaxationRate();

    //! Constructor using the internal generic representation.
    /*!
     * @param A       Pre-exponential factor.
     * @param B       Dimensionless constant in the exponential.
     * @param C       Coefficient multiplying T^(-1/3).
     * @param D       Coefficient multiplying T^(-m).
     * @param b       Dimensionless temperature exponent.
     * @param scaling Dimensionless scaling factor.
     * @param m       Temperature exponent used by the D term.
     * @param E       Coefficient multiplying T^(-z).
     * @param z       Temperature exponent used by the E term.
     */
    VibrationalRelaxationRate(double A, double B, double C, double D,
                     double b, double scaling = 1.0,
                     double m = 2.0 / 3.0,
                     double E = 0.0, double z = 1.0);

    //! Constructor based on AnyMap content.
    explicit VibrationalRelaxationRate(const AnyMap& node,
                              const UnitStack& rate_units = {});

    //! Set rate parameters from an AnyMap.
    void setParameters(const AnyMap& node,
                       const UnitStack& rate_units) override;

    //! Get rate parameters for YAML serialization.
    void getParameters(AnyMap& node) const override;

    //! Create a rate evaluator for this reaction rate type.
    unique_ptr<MultiRateBase> newMultiRate() const override {
        return make_unique<MultiRate<VibrationalRelaxationRate, DetailedVibData>>();
    }

    //! String identifying this reaction rate type.
    const string type() const override {
        return "vibrational-relaxation";
    }

    //! Set context of reaction rate evaluation.
    /*!
     * Vibrational relaxation rates are intended for irreversible
     * non-equilibrium plasma reactions. Reversible reactions are rejected
     * because the reverse rate cannot be obtained from conventional
     * thermochemistry for these models.
     */
    void setContext(const Reaction& rxn, const Kinetics& kin) override;

    //! Evaluate the forward rate coefficient.
    double evalFromStruct(const DetailedVibData& shared_data) const {
        const double invT = shared_data.recipT;
        const double invT13 = std::cbrt(invT);

        return m_scaling * m_A * std::exp(
            m_b * shared_data.logT
            + m_B
            + m_C * invT13
            + m_D * std::pow(invT, m_m)
            + m_E * std::pow(invT, m_z)
        );
    }

    //! Evaluate the scaled temperature derivative.
    /*!
     * This returns:
     *
     * @f[
     * \frac{1}{k_f} \frac{d k_f}{dT}
     * =
     * \frac{d \ln k_f}{dT}
     * @f]
     *
     * For the internal generic expression, this is:
     *
     * @f[
     * \frac{b}{T}
     * - \frac{C}{3} T^{-4/3}
     * - m D T^{-m-1}
     * - z E T^{-z-1}
     * @f]
     */
    double ddTScaledFromStruct(const DetailedVibData& shared_data) const;

private:
    //! Dimensionless constant in the exponential.
    double m_B = 0.0;

    //! Coefficient multiplying T^(-1/3).
    double m_C = 0.0;

    //! Coefficient multiplying T^(-m).
    double m_D = 0.0;

    //! Dimensionless scaling factor.
    double m_scaling = 1.0;

    //! Temperature exponent used by the D term.
    double m_m = 2.0 / 3.0;

    //! Coefficient multiplying T^(-z).
    double m_E = 0.0;

    //! Temperature exponent used by the E term.
    double m_z = 1.0;

    //! Castela coefficient a.
    double m_castela_a = 0.0;

    //! Castela coefficient b.
    double m_castela_b = 0.0;

    //! Castela reference pressure.
    double m_referencePressure = OneAtm;

    //! Selected vibrational relaxation model.
    /*!
     * Accepted values:
     *
     * - `constant`
     * - `multi-state-resolved`
     * - `starikovskiy`
     * - `castela`
     */
    string m_vibration_model = "multi-state-resolved";

    //! YAML key names.
    string m_B_str = "B";
    string m_C_str = "C";
    string m_D_str = "D";
    string m_scaling_str = "scaling";
    string m_m_str = "m";
    string m_E_str = "E";
    string m_z_str = "z";
    string m_reference_pressure_str = "reference-pressure";

    //! Configure the ArrheniusBase part from an already-converted internal A value.
    /*!
    * This is needed for models such as Castela, where the user-facing YAML does
    * not contain a standard Arrhenius A coefficient.
    */
    void configureBaseFromInternalA(const AnyMap& node,
                                    const UnitStack& rate_units,
                                    double A, double b);

    //! Configure the ArrheniusBase part from a YAML A value and an explicit b.
    /*!
    * This is needed for models such as constant and Starikovskiy, where the YAML
    * does not contain the standard Arrhenius pair A / b.
    */
    void configureBaseFromYamlA(const AnyMap& node,
                                const UnitStack& rate_units,
                                const AnyValue& A,
                                double b);

    };

}

#endif