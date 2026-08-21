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

//! Common implementation for vibrational relaxation reaction rates.
/**
 * This intermediate class provides the common data storage and evaluation
 * used by the concrete vibrational relaxation rate implementations.
 *
 * Internally, the supported parameterizations are mapped to
 *
 * @f[
 * k_f =
 * A \exp\left(
 *     b \ln T
 *     + B
 *     + C T^{-1/3}
 *     + D T^{-m}
 *     + E T^{-z}
 * \right).
 * @f]
 *
 * Concrete derived classes are responsible for mapping their user-facing
 * parameters onto this internal representation.
 *
 * Vibrational relaxation rates are restricted to irreversible reactions.
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