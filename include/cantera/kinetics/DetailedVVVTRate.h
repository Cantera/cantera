//! @file DetailedVVVTRate.h
//! Header for plasma reaction rates corresponding to detailed vibration handling.
//!
//! This file is part of Cantera. See License.txt in the top-level directory or
//! at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_DETAILEDVVVTRATE_H
#define CT_DETAILEDVVVTRATE_H

#include "Arrhenius.h"

#include <cmath>

namespace Cantera
{

//! Data container holding shared data specific to DetailedVVVTRate.
struct DetailedVibData : public ReactionData
{
    //! Update cached temperature-dependent data.
    /*!
     * @param phase Thermodynamic phase used to retrieve the gas temperature.
     * @param kin   Kinetics object. Not used here, but required by the
     *              `ReactionData` interface.
     *
     * @return `true` if the temperature has changed and rates need to be
     *         recomputed; `false` otherwise.
     */
    bool update(const ThermoPhase& phase, const Kinetics& kin) override;

    using ReactionData::update;
};


//! Detailed VV-VT reaction rate.
/*!
 * This reaction rate implements the following temperature-dependent
 * expression:
 *
 * @f[
 * k_f =
 * scaling \, A \,
 * \exp \left(
 *     b \ln T
 *     + B
 *     + C T^{-1/3}
 *     + D T^{-2/3}
 * \right)
 * @f]
 *
 * where `T` is the gas temperature in K.
 *
 * Unit conventions:
 *
 * - `A` uses the standard Cantera rate coefficient units. Its units depend
 *   on the reaction order and dimensionality, and are converted by
 *   `ArrheniusBase`.
 * - `b` is dimensionless.
 * - `B` is dimensionless.
 * - `C` is interpreted as being expressed in K^(1/3), assuming `T` is in K.
 * - `D` is interpreted as being expressed in K^(2/3), assuming `T` is in K.
 * - `scaling` is dimensionless.
 *
 * Important: `B`, `C`, `D`, and `scaling` are read as raw floating-point
 * numbers from YAML. They are not converted by Cantera's unit system.
 *
 * @ingroup arrheniusGroup
 */
class DetailedVVVTRate : public ArrheniusBase
{
public:
    //! Default constructor.
    DetailedVVVTRate();

    //! Constructor.
    /*!
     * @param A       Pre-exponential factor. The unit system is (kmol, m, s);
     *                actual units depend on the reaction order and the
     *                dimensionality, surface or bulk.
     * @param B       Dimensionless constant in the exponential.
     * @param C       Constant multiplying T^(-1/3), interpreted in K^(1/3).
     * @param D       Constant multiplying T^(-2/3), interpreted in K^(2/3).
     * @param b       Dimensionless temperature exponent.
     * @param scaling Dimensionless scaling factor from the harmonic oscillator
     *                theory.
     */
    DetailedVVVTRate(double A, double B, double C, double D,
                     double b, double scaling = 1.0);

    //! Constructor based on AnyMap content.
    /*!
     * Expected YAML form:
     *
     * @code{.yaml}
     * type: detailed-vv-vt
     * rate-constant:
     *   A: ...
     *   b: ...
     *   B: ...
     *   C: ...
     *   D: ...
     *   scaling: ...
     * @endcode
     *
     * The `scaling` entry is optional and defaults to 1.0.
     */
    explicit DetailedVVVTRate(const AnyMap& node,
                              const UnitStack& rate_units = {});

    //! Set rate parameters from an AnyMap.
    /*!
     * `A` and `b` are handled by `ArrheniusBase`.
     * `B`, `C`, `D`, and `scaling` are handled by this class.
     */
    void setParameters(const AnyMap& node,
                       const UnitStack& rate_units) override;

    //! Get rate parameters for YAML serialization.
    /*!
     * This method is needed so that Cantera can reconstruct or export the
     * full custom rate, including `B`, `C`, `D`, and `scaling`.
     */
    void getParameters(AnyMap& node) const override;

    //! Create a rate evaluator for this reaction rate type.
    unique_ptr<MultiRateBase> newMultiRate() const override {
        return make_unique<MultiRate<DetailedVVVTRate, DetailedVibData>>();
    }

    //! String identifying this reaction rate specialization.
    const string type() const override {
        return "detailed-vv-vt";
    }

    //! Set context of reaction rate evaluation.
    /*!
     * This rate is intended for irreversible non-equilibrium plasma reactions.
     * Reversible reactions are rejected because the reverse rate cannot be
     * obtained from conventional thermochemistry for this model.
     */
    void setContext(const Reaction& rxn, const Kinetics& kin) override;

    //! Evaluate reaction rate.
    /*!
     * @param shared_data Data shared by all reactions of this type.
     *
     * @return Forward rate coefficient.
     */
    double evalFromStruct(const DetailedVibData& shared_data) const {
        const double recipT13 = std::cbrt(shared_data.recipT);

        return m_scaling * m_A * std::exp(
            m_b * shared_data.logT
            + m_B
            + m_C * recipT13
            + m_D * recipT13 * recipT13
        );
    }

    //! Evaluate derivative of reaction rate with respect to temperature,
    //! divided by the reaction rate.
    /*!
     * This returns:
     *
     * @f[
     * \frac{1}{k_f} \frac{d k_f}{dT}
     * =
     * \frac{b}{T}
     * - \frac{C}{3} T^{-4/3}
     * - \frac{2D}{3} T^{-5/3}
     * @f]
     *
     * This derivative is consistent with `evalFromStruct` and does not modify
     * the reaction rate expression itself.
     *
     * @param shared_data Data shared by all reactions of this type.
     */
    double ddTScaledFromStruct(const DetailedVibData& shared_data) const;

private:
    //! Dimensionless constant in the exponential.
    double m_B = 0.0;

    //! Constant multiplying T^(-1/3), interpreted in K^(1/3).
    double m_C = 0.0;

    //! Constant multiplying T^(-2/3), interpreted in K^(2/3).
    double m_D = 0.0;

    //! Dimensionless scaling factor.
    double m_scaling = 1.0;

    //! YAML key for `B`.
    string m_B_str = "B";

    //! YAML key for `C`.
    string m_C_str = "C";

    //! YAML key for `D`.
    string m_D_str = "D";

    //! YAML key for `scaling`.
    string m_scaling_str = "scaling";
};

}

#endif