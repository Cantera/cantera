//! @file BlowersMaselRate.h   Header for Blowers-Masel reaction rates

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_BLOWERSMASELRATE_H
#define CT_BLOWERSMASELRATE_H

#include "Arrhenius.h"

namespace Cantera
{

//! Data container holding shared data specific to BlowersMaselRate
/**
 * The data container `BlowersMaselData` holds precalculated data common to
 * all `BlowersMaselRate` objects.
 */
struct BlowersMaselData : public ReactionData
{
    BlowersMaselData();

    virtual void update(double T) override;
    virtual bool update(const ThermoPhase& phase, const Kinetics& kin) override;
    using ReactionData::update;

    virtual void resize(size_t nSpecies, size_t nReactions, size_t nPhases) override {
        partialMolarEnthalpies.resize(nSpecies, 0.);
        ready = true;
    }

    bool ready; //!< boolean indicating whether vectors are accessible
    double density; //!< used to determine if updates are needed
    vector_fp partialMolarEnthalpies; //!< partial molar enthalpies

protected:
    int m_state_mf_number; //!< integer that is incremented when composition changes
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

    explicit BlowersMaselRate(const AnyMap& node, const UnitStack& rate_units={})
        : BlowersMaselRate()
    {
        setParameters(node, rate_units);
    }

    unique_ptr<MultiRateBase> newMultiRate() const override {
        return unique_ptr<MultiRateBase>(new MultiRate<BlowersMaselRate, BlowersMaselData>);
    }

    virtual const std::string type() const override {
        return "Blowers-Masel";
    }

    virtual void setContext(const Reaction& rxn, const Kinetics& kin) override;

    //! Evaluate reaction rate
    double evalRate(double logT, double recipT) const {
        double Ea_R = effectiveActivationEnergy_R(m_deltaH_R);
        return m_A * std::exp(m_b * logT - Ea_R * recipT);
    }

    //! Update information specific to reaction
    void updateFromStruct(const BlowersMaselData& shared_data) {
        if (shared_data.ready) {
            m_deltaH_R = 0.;
            for (const auto& item : m_stoich_coeffs) {
                m_deltaH_R += shared_data.partialMolarEnthalpies[item.first] * item.second;
            }
            m_deltaH_R /= GasConstant;
        }
    }

    //! Evaluate reaction rate
    /*!
     *  @param shared_data  data shared by all reactions of a given type
     */
    double evalFromStruct(const BlowersMaselData& shared_data) const {
        double Ea_R = effectiveActivationEnergy_R(m_deltaH_R);
        return m_A * std::exp(m_b * shared_data.logT - Ea_R * shared_data.recipT);
    }

    //! Evaluate derivative of reaction rate with respect to temperature
    //! divided by reaction rate
    /*!
     *  This method does not consider potential changes due to a changed reaction
     *  enthalpy. A corresponding warning is raised.
     *  @param shared_data  data shared by all reactions of a given type
     */
    double ddTScaledFromStruct(const BlowersMaselData& shared_data) const;

protected:
    //! Return the effective activation energy (a function of the delta H of reaction)
    //! divided by the gas constant (that is, the activation temperature) [K]
    //! @internal  The enthalpy change of reaction is not an independent parameter
    double effectiveActivationEnergy_R(double deltaH_R) const {
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

public:
    virtual double activationEnergy() const override {
        return effectiveActivationEnergy_R(m_deltaH_R) * GasConstant;
    }

    //! Return the bond dissociation energy *w* [J/kmol]
    double bondEnergy() const {
        return m_E4_R * GasConstant;
    }

    //! Return current enthalpy change of reaction [J/kmol]
    double deltaH() const {
        return m_deltaH_R * GasConstant;
    }

    //! Set current enthalpy change of reaction [J/kmol]
    /*!
     *  @internal  used for testing purposes only; note that this quantity is not an
     *      independent variable and will be overwritten during an update of the state.
     *
     *  @warning  This method is an experimental part of the %Cantera API and
     *      may be changed or removed without notice.
     */
    void setDeltaH(double deltaH) {
        m_deltaH_R = deltaH / GasConstant;
    }

protected:
    //! Pairs of species indices and multiplers to calculate enthalpy change
    std::vector<std::pair<size_t, double>> m_stoich_coeffs;

    double m_deltaH_R; //!< enthalpy change of reaction (in temperature units)
};

}

#endif
