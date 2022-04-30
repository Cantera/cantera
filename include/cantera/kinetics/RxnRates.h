/**
 *  @file RxnRates.h
 *  @deprecated To be removed after Cantera 2.6. See class Cantera::ReactionRate and
 *      derived classes for new reaction rate handlers.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_RXNRATES_H
#define CT_RXNRATES_H

#include "cantera/kinetics/reaction_defs.h"
#include "cantera/kinetics/Arrhenius.h"
#include "cantera/kinetics/PlogRate.h"
#include "MultiRate.h"
#include "cantera/base/Array.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/global.h"

namespace Cantera
{

class Array2D;
class AnyValue;
class AnyMap;
class UnitSystem;
class Units;
class Func1;

//! Arrhenius reaction rate type depends only on temperature
/**
 * A reaction rate coefficient of the following form.
 *
 *   \f[
 *        k_f =  A T^b \exp (-E/RT)
 *   \f]
 *
 * @deprecated To be removed after Cantera 2.6. See ArrheniusRate.
 */
class Arrhenius2 final : public ArrheniusRate
{
public:
    //! Default constructor.
    Arrhenius2();

    /// Constructor.
    /// @param A pre-exponential. The unit system is
    ///     (kmol, m, s). The actual units depend on the reaction
    ///     order and the dimensionality (surface or bulk).
    /// @param b Temperature exponent. Non-dimensional.
    /// @param E Activation energy in temperature units. Kelvin.
    ///
    /// @todo  Add deprecation warning pointing out change of activation
    ///     energy units.
    Arrhenius2(doublereal A, doublereal b, doublereal E);

    //! Constructor based on AnyMap content
    Arrhenius2(const AnyValue& rate,
               const UnitSystem& units, const Units& rate_units);

    //! Converting constructor (to facilitate back-ward compatibility)
    Arrhenius2(const ArrheniusRate& other);

    void setRateParameters(const AnyValue& rate,
                           const UnitSystem& units, const Units& rate_units);
    using ArrheniusRate::setRateParameters;

    //! Return parameters - two-parameter version
    void getParameters(AnyMap& node, const Units& rate_units) const;
    using ArrheniusRate::getParameters;

    //! Update concentration-dependent parts of the rate coefficient.
    /*!
     *   For this class, there are no concentration-dependent parts, so this
     *   method does nothing.
     */
    void update_C(const doublereal* c) {
    }

    /**
     * Update the value of the natural logarithm of the rate constant.
     */
    doublereal updateLog(doublereal logT, doublereal recipT) const {
        return m_logA + m_b*logT - m_Ea_R*recipT;
    }

    /**
     * Update the value the rate constant.
     *
     * This function returns the actual value of the rate constant. It can be
     * safely called for negative values of the pre-exponential factor.
     */
    doublereal updateRC(doublereal logT, doublereal recipT) const {
        return m_A * std::exp(m_b*logT - m_Ea_R*recipT);
    }

    virtual const std::string type() const override {
        return "Arrhenius2";
    }

    //! Return the activation energy divided by the gas constant (that is, the
    //! activation temperature) [K]
    double activationEnergy_R() const {
        return m_Ea_R;
    }
};


/**
 * An Arrhenius rate with coverage-dependent terms.
 *
 * The rate expression is given by [Kee, R. J., Coltrin, M. E., & Glarborg, P.
 * (2005). Chemically reacting flow: theory and practice. John Wiley & Sons.
 * Eq 11.113]:
 * \f[
 *     k_f = A T^b \exp \left(
 *             \ln 10 \sum a_k \theta_k
 *             - \frac{1}{RT} \left( E_a + \sum E_k\theta_k \right)
 *             + \sum m_k \ln \theta_k
 *             \right)
 *   \f]
 * or, equivalently, and as implemented in Cantera,
 * \f[
 *     k_f = A T^b \exp \left( - \frac{E_a}{RT} \right)
 *             \prod_k 10^{a_k \theta_k} \theta_k^{m_k}
 *             \exp \left( \frac{- E_k \theta_k}{RT} \right)
 *   \f]
 * where the parameters \f$ (a_k, E_k, m_k) \f$ describe the dependency on the
 * surface coverage of species \f$k, \theta_k \f$.
 *
 * @deprecated To be removed after Cantera 2.6. See InterfaceRate and StickingRate.
 */
class SurfaceArrhenius
{

public:
    SurfaceArrhenius();
    explicit SurfaceArrhenius(double A, double b, double Ta);

    //! Add a coverage dependency for species *k*, with exponential dependence
    //! *a*, power-law exponent *m*, and activation energy dependence *e*,
    //! where *e* is in Kelvin, that is, energy divided by the molar gas constant.
    void addCoverageDependence(size_t k, doublereal a,
                               doublereal m, doublereal e);

    void update_C(const doublereal* theta) {
        m_acov = 0.0;
        m_ecov = 0.0;
        m_mcov = 0.0;
        size_t k;
        doublereal th;
        for (size_t n = 0; n < m_ac.size(); n++) {
            k = m_sp[n];
            m_acov += m_ac[n] * theta[k];
            m_ecov += m_ec[n] * theta[k];
        }
        for (size_t n = 0; n < m_mc.size(); n++) {
            k = m_msp[n];
            th = std::max(theta[k], Tiny);
            m_mcov += m_mc[n]*std::log(th);
        }
    }

    /**
     * Update the value of the rate constant.
     *
     * This function returns the actual value of the rate constant. It can be
     * safely called for negative values of the pre-exponential factor.
     */
    doublereal updateRC(doublereal logT, doublereal recipT) const {
        return m_A * std::exp(std::log(10.0)*m_acov + m_b*logT -
                              (m_E + m_ecov)*recipT + m_mcov);
    }

    //! Return the pre-exponential factor *A* (in m, kmol, s to powers depending
    //! on the reaction order) accounting coverage dependence.
    /*!
     *  Returns reaction pre-exponent accounting for both *a* and *m*.
     */
    doublereal preExponentialFactor() const {
        return m_A * std::exp(std::log(10.0)*m_acov + m_mcov);
    }

    //! Return effective temperature exponent
    doublereal temperatureExponent() const {
        return m_b;
    }

    //! Return the activation energy divided by the gas constant (that is, the
    //! activation temperature) [K], accounting coverage dependence.
    doublereal activationEnergy_R() const {
        return m_E + m_ecov;
    }

protected:
    doublereal m_b, m_E, m_A;
    doublereal m_acov, m_ecov, m_mcov;
    std::vector<size_t> m_sp, m_msp;
    vector_fp m_ac, m_ec, m_mc;
};


#ifdef CT_NO_LEGACY_REACTIONS_26
typedef ArrheniusRate Arrhenius;
#else
typedef Arrhenius2 Arrhenius;
#endif

typedef PlogRate Plog;

}

#endif
