/**
 *  @file RxnRates.h
 *
 * @todo  at least part of the content of this file should be transferred
 *     to ReactionRate.h once the old XML interface is removed after Cantera 2.6
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_RXNRATES_H
#define CT_RXNRATES_H

#include "cantera/kinetics/reaction_defs.h"
#include "cantera/kinetics/Arrhenius.h"
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
 */
class Arrhenius2 : public ArrheniusBase
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
    Arrhenius2(doublereal A, doublereal b, doublereal E);

    //! Constructor based on AnyMap content
    Arrhenius2(const AnyValue& rate,
               const UnitSystem& units, const Units& rate_units);

    void setRateParameters(const AnyValue& rate,
                           const UnitSystem& units, const Units& rate_units);
    using ArrheniusBase::setRateParameters;

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

    //! Return the activation energy divided by the gas constant (i.e. the
    //! activation temperature) [K]
    doublereal activationEnergy_R() const {
        return m_Ea_R;
    }

    const std::string type() const {
        return "Arrhenius";
    }

    // The following methods are unused but required by the base class
    unique_ptr<MultiRateBase> newMultiRate() const {
        throw NotImplementedError("Arrhenius2::newMultiRate");
    }

protected:
    double m_logA;
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
 */
class BlowersMasel2
{
public:
    //! Default constructor.
    BlowersMasel2();

    /// Constructor.
    /// @param A pre-exponential. The unit system is
    ///     (kmol, m, s). The actual units depend on the reaction
    ///     order and the dimensionality (surface or bulk).
    /// @param b Temperature exponent. Non-dimensional.
    /// @param E0 Intrinsic activation energy in temperature units. Kelvin.
    /// @param w average bond dissociation energy of the bond being formed and
    ///     broken in the reaction, in temperature units. Kelvin.

    BlowersMasel2(double A, double b, double E0, double w);

    unique_ptr<MultiRateBase> newMultiRate() const {
        throw NotImplementedError("BlowersMasel2::newMultiRate");
    }

    void getParameters(AnyMap& rateNode, const Units& rate_units) const;

    //! Update concentration-dependent parts of the rate coefficient.
    /*!
     *   For this class, there are no concentration-dependent parts, so this
     *   method does nothing.
     */
    void update_C(const doublereal* c) {
    }

    /**
     * Update the value the rate constant.
     *
     * This function returns the actual value of the rate constant. It can be
     * safely called for negative values of the pre-exponential factor.
     */
    doublereal updateRC(doublereal logT, doublereal recipT, doublereal deltaH) {
        double m_E = activationEnergy_R(deltaH);
        return m_A * std::exp(m_b * logT - m_E * recipT);
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

    //! Return the actual activation energy (a function of the delta H of reaction)
    //! divided by the gas constant (i.e. the activation temperature) [K]
    doublereal activationEnergy_R(doublereal deltaH) {
        double Ea; // will be in temperature units (Kelvin)
        double deltaH_R = deltaH / GasConstant; // deltaH in temperature units (Kelvin)
        if (deltaH_R < -4 * m_E0) {
            Ea = 0;
        } else if (deltaH_R > 4 * m_E0) {
            Ea = deltaH_R;
        } else {
            // m_w is in Kelvin
            // vp is in Kelvin
            double vp = 2 * m_w * ((m_w + m_E0) / (m_w - m_E0));
            double vp_2w_dH = (vp - 2 * m_w + deltaH_R); // (Vp - 2 w + dH)
            Ea = (m_w + deltaH_R / 2) * (vp_2w_dH * vp_2w_dH) /
                 (vp * vp - 4 * m_w * m_w + deltaH_R * deltaH_R); // in Kelvin
        }
        return Ea;
    }

    //! Return the intrinsic activation energy divided by the gas constant (i.e. the
    //! intrinsic activation temperature) [K]
    doublereal activationEnergy_R0() const {
        return m_E0;
    }
    //! Return the bond dissociation energy *w* divided by the gas constant[K]
    doublereal bondEnergy() const {
        return m_w;
    }

protected:
    doublereal m_logA, m_b, m_A, m_w, m_E0;
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
 */
class SurfaceArrhenius
{

public:
    SurfaceArrhenius();
    explicit SurfaceArrhenius(double A, double b, double Ta);

    //! Add a coverage dependency for species *k*, with exponential dependence
    //! *a*, power-law exponent *m*, and activation energy dependence *e*,
    //! where *e* is in Kelvin, i.e. energy divided by the molar gas constant.
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

    //! Return the activation energy divided by the gas constant (i.e. the
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
typedef BlowersMaselRate BlowersMasel;
#else
typedef Arrhenius2 Arrhenius;
typedef BlowersMasel2 BlowersMasel;
#endif


//! Pressure-dependent reaction rate expressed by logarithmically interpolating
//! between Arrhenius rate expressions at various pressures.
/*!
 * Given two rate expressions at two specific pressures:
 *
 *   * \f$ P_1: k_1(T) = A_1 T^{b_1} e^{-E_1 / RT} \f$
 *   * \f$ P_2: k_2(T) = A_2 T^{b_2} e^{-E_2 / RT} \f$
 *
 * The rate at an intermediate pressure \f$ P_1 < P < P_2 \f$ is computed as
 * \f[
 *  \log k(T,P) = \log k_1(T) + \bigl(\log k_2(T) - \log k_1(T)\bigr)
 *      \frac{\log P - \log P_1}{\log P_2 - \log P_1}
 * \f]
 * Multiple rate expressions may be given at the same pressure, in which case
 * the rate used in the interpolation formula is the sum of all the rates given
 * at that pressure. For pressures outside the given range, the rate expression
 * at the nearest pressure is used.
 */
class PlogRate final : public ReactionRate
{
public:
    //! Default constructor.
    PlogRate();

    //! Constructor from Arrhenius rate expressions at a set of pressures
    explicit PlogRate(const std::multimap<double, Arrhenius>& rates);

    PlogRate(const AnyMap& node, const UnitStack& rate_units={}) : PlogRate() {
        setParameters(node, rate_units);
    }

    unique_ptr<MultiRateBase> newMultiRate() const {
        return unique_ptr<MultiRateBase>(new MultiBulkRate<PlogRate, PlogData>);
    }

    //! Identifier of reaction rate type
    const std::string type() const { return "pressure-dependent-Arrhenius"; }

    //! Perform object setup based on AnyMap node information
    /*!
     *  @param node  AnyMap containing rate information
     *  @param units  Unit definitions specific to rate information
     */
    void setParameters(const AnyMap& node, const UnitStack& units);

    //! Perform object setup based on reaction rate information
    /*!
     *  @param rates  vector of AnyMap containing rate information
     *  @param units  unit system
     *  @param rate_units  unit definitions specific to rate information
     */
    void setRateParameters(const std::vector<AnyMap>& rates,
                           const UnitSystem& units, const Units& rate_units);

    void getParameters(AnyMap& rateNode, const Units& rate_units) const;
    void getParameters(AnyMap& rateNode) const {
        return getParameters(rateNode, Units(0));
    }

    //! Evaluate reaction rate
    /*!
     *  @param shared_data  data shared by all reactions of a given type
     */
    double evalFromStruct(const PlogData& shared_data) {
        if (shared_data.logP != logP_) {
            update_C(&shared_data.logP);
        }
        return updateRC(shared_data.logT, shared_data.recipT);
    }

    //! Set up Plog object
    /*!
     * @deprecated   Deprecated in Cantera 2.6. Renamed to setRates.
     */
    void setup(const std::multimap<double, Arrhenius>& rates);

    //! Set up Plog object
    void setRates(const std::multimap<double, Arrhenius>& rates);

    //! Update concentration-dependent parts of the rate coefficient.
    //! @param c natural log of the pressure in Pa
    void update_C(const doublereal* c) {
        logP_ = c[0];
        if (logP_ > logP1_ && logP_ < logP2_) {
            return;
        }

        auto iter = pressures_.upper_bound(c[0]);
        AssertThrowMsg(iter != pressures_.end(), "PlogRate::update_C",
                       "Pressure out of range: {}", logP_);
        AssertThrowMsg(iter != pressures_.begin(), "PlogRate::update_C",
                       "Pressure out of range: {}", logP_);

        // upper interpolation pressure
        logP2_ = iter->first;
        ihigh1_ = iter->second.first;
        ihigh2_ = iter->second.second;

        // lower interpolation pressure
        logP1_ = (--iter)->first;
        ilow1_ = iter->second.first;
        ilow2_ = iter->second.second;

        rDeltaP_ = 1.0 / (logP2_ - logP1_);
    }

    /**
     * Update the value the rate constant.
     *
     * This function returns the actual value of the rate constant.
     */
    doublereal updateRC(doublereal logT, doublereal recipT) const {
        double log_k1, log_k2;
        if (ilow1_ == ilow2_) {
            log_k1 = rates_[ilow1_].updateLog(logT, recipT);
        } else {
            double k = 1e-300; // non-zero to make log(k) finite
            for (size_t i = ilow1_; i < ilow2_; i++) {
                k += rates_[i].updateRC(logT, recipT);
            }
            log_k1 = std::log(k);
        }

        if (ihigh1_ == ihigh2_) {
            log_k2 = rates_[ihigh1_].updateLog(logT, recipT);
        } else {
            double k = 1e-300; // non-zero to make log(k) finite
            for (size_t i = ihigh1_; i < ihigh2_; i++) {
                k += rates_[i].updateRC(logT, recipT);
            }
            log_k2 = std::log(k);
        }

        return std::exp(log_k1 + (log_k2-log_k1) * (logP_-logP1_) * rDeltaP_);
    }

    //! Check to make sure that the rate expression is finite over a range of
    //! temperatures at each interpolation pressure. This is potentially an
    //! issue when one of the Arrhenius expressions at a particular pressure
    //! has a negative pre-exponential factor.
    void validate(const std::string& equation, const Kinetics& kin) {
        validate(equation);
    }

    void validate(const std::string& equation);

    //! Return the pressures and Arrhenius expressions which comprise this
    //! reaction.
    /*!
     * @deprecated  Behavior to change after Cantera 2.6.
     *              @see getRates for new behavior.
     */
    std::vector<std::pair<double, Arrhenius> > rates() const;

    //! Return the pressures and Arrhenius expressions which comprise this
    //! reaction.
    std::multimap<double, Arrhenius> getRates() const;

protected:
    //! log(p) to (index range) in the rates_ vector
    std::map<double, std::pair<size_t, size_t> > pressures_;

    // Rate expressions which are referenced by the indices stored in pressures_
    std::vector<Arrhenius> rates_;

    double logP_; //!< log(p) at the current state
    double logP1_, logP2_; //!< log(p) at the lower / upper pressure reference

    //! Indices to the ranges within rates_ for the lower / upper pressure, such
    //! that rates_[ilow1_] through rates_[ilow2_] (inclusive) are the rates
    //! expressions which are combined to form the rate at the lower reference
    //! pressure.
    size_t ilow1_, ilow2_, ihigh1_, ihigh2_;

    double rDeltaP_; //!< reciprocal of (logP2 - logP1)
};

typedef PlogRate Plog;

//! Pressure-dependent rate expression where the rate coefficient is expressed
//! as a bivariate Chebyshev polynomial in temperature and pressure.
/*!
 * The rate constant can be written as:
 * \f[
 *     \log k(T,P) = \sum_{t=1}^{N_T} \sum_{p=1}^{N_P} \alpha_{tp}
 *                       \phi_t(\tilde{T}) \phi_p(\tilde{P})
 * \f]
 * where \f$\alpha_{tp}\f$ are the constants defining the rate, \f$\phi_n(x)\f$
 * is the Chebyshev polynomial of the first kind of degree *n* evaluated at
 * *x*, and
 * \f[
 *  \tilde{T} \equiv \frac{2T^{-1} - T_\mathrm{min}^{-1} - T_\mathrm{max}^{-1}}
 *                        {T_\mathrm{max}^{-1} - T_\mathrm{min}^{-1}}
 * \f]
 * \f[
 *  \tilde{P} \equiv \frac{2 \log P - \log P_\mathrm{min} - \log P_\mathrm{max}}
 *                        {\log P_\mathrm{max} - \log P_\mathrm{min}}
 * \f]
 * are reduced temperature and reduced pressures which map the ranges
 * \f$ (T_\mathrm{min}, T_\mathrm{max}) \f$ and
 * \f$ (P_\mathrm{min}, P_\mathrm{max}) \f$ to (-1, 1).
 *
 * A ChebyshevRate3 rate expression is specified in terms of the coefficient matrix
 * \f$ \alpha \f$ and the temperature and pressure ranges. Note that the
 * Chebyshev polynomials are not defined outside the interval (-1,1), and
 * therefore extrapolation of rates outside the range of temperatures and
 * pressures for which they are defined is strongly discouraged.
 */
class ChebyshevRate3 final : public ReactionRate
{
public:
    //! Default constructor.
    ChebyshevRate3() : m_log10P(NAN), m_rate_units(Units(0.)) {}

    //! Constructor directly from coefficient array
    /*!
     *  @param Tmin    Minimum temperature [K]
     *  @param Tmax    Maximum temperature [K]
     *  @param Pmin    Minimum pressure [Pa]
     *  @param Pmax    Maximum pressure [Pa]
     *  @param coeffs  Coefficient array dimensioned `nT` by `nP` where `nT` and
     *      `nP` are the number of temperatures and pressures used in the fit,
     *      respectively.
     */
    ChebyshevRate3(double Tmin, double Tmax, double Pmin, double Pmax,
                   const Array2D& coeffs);

    ChebyshevRate3(const AnyMap& node, const UnitStack& rate_units={})
        : ChebyshevRate3()
    {
        setParameters(node, rate_units);
    }

    unique_ptr<MultiRateBase> newMultiRate() const {
        return unique_ptr<MultiRateBase>(
            new MultiBulkRate<ChebyshevRate3, ChebyshevData>);
    }

    const std::string type() const { return "Chebyshev"; }

    //! Perform object setup based on AnyMap node information
    /*!
     *  @param node  AnyMap containing rate information
     *  @param units  Unit definitions specific to rate information
     */
    void setParameters(const AnyMap& node, const UnitStack& units);
    void getParameters(AnyMap& rateNode, const Units& rate_units) const {
        // @todo: deprecate, as second argument is no longer needed
        return getParameters(rateNode);
    }
    void getParameters(AnyMap& rateNode) const;

    //! Update information specific to reaction
    /*!
     *  @param shared_data  data shared by all reactions of a given type
     */
    double evalFromStruct(const ChebyshevData& shared_data) {
        if (shared_data.log10P != m_log10P) {
            update_C(&shared_data.log10P);
        }
        return updateRC(0., shared_data.recipT);
    }

    //! Set up ChebyshevRate3 object
    /*!
     * @deprecated   Deprecated in Cantera 2.6. Replaceable with
     *               @see setLimits() and @see setCoeffs().
     */
    void setup(double Tmin, double Tmax, double Pmin, double Pmax,
                  const Array2D& coeffs);

    //! Set limits for ChebyshevRate3 object
    /*!
     *  @param Tmin    Minimum temperature [K]
     *  @param Tmax    Maximum temperature [K]
     *  @param Pmin    Minimum pressure [Pa]
     *  @param Pmax    Maximum pressure [Pa]
     */
    void setLimits(double Tmin, double Tmax, double Pmin, double Pmax);

    //! Update concentration-dependent parts of the rate coefficient.
    //! @param c base-10 logarithm of the pressure in Pa
    void update_C(const double* c) {
        m_log10P = c[0];
        double Pr = (2 * c[0] + PrNum_) * PrDen_;
        double Cnm1 = Pr;
        double Cn = 1;
        double Cnp1;
        for (size_t i = 0; i < m_coeffs.nRows(); i++) {
            dotProd_[i] = m_coeffs(i, 0);
        }
        for (size_t j = 1; j < m_coeffs.nColumns(); j++) {
            Cnp1 = 2 * Pr * Cn - Cnm1;
            for (size_t i = 0; i < m_coeffs.nRows(); i++) {
                dotProd_[i] += Cnp1 * m_coeffs(i, j);
            }
            Cnm1 = Cn;
            Cn = Cnp1;
        }
    }

    /**
     * Update the value the rate constant.
     *
     * This function returns the actual value of the rate constant.
     */
    double updateRC(double logT, double recipT) const {
        double Tr = (2 * recipT + TrNum_) * TrDen_;
        double Cnm1 = Tr;
        double Cn = 1;
        double Cnp1;
        double logk = dotProd_[0];
        for (size_t i = 1; i < m_coeffs.nRows(); i++) {
            Cnp1 = 2 * Tr * Cn - Cnm1;
            logk += Cnp1 * dotProd_[i];
            Cnm1 = Cn;
            Cn = Cnp1;
        }
        return std::pow(10, logk);
    }

    //! Minimum valid temperature [K]
    double Tmin() const {
        return Tmin_;
    }

    //! Maximum valid temperature [K]
    double Tmax() const {
        return Tmax_;
    }

    //! Minimum valid pressure [Pa]
    double Pmin() const {
        return Pmin_;
    }

    //! Maximum valid pressure [Pa]
    double Pmax() const {
        return Pmax_;
    }

    //! Number of points in the pressure direction
    size_t nPressure() const {
        return m_coeffs.nColumns();
    }

    //! Number of points in the temperature direction
    size_t nTemperature() const {
        return m_coeffs.nRows();
    }

    //! Access the ChebyshevRate3 coefficients.
    /*!
     *  \f$ \alpha_{t,p} = \mathrm{coeffs}[N_P*t + p] \f$ where
     *  \f$ 0 <= t < N_T \f$ and \f$ 0 <= p < N_P \f$.
     *
     * @deprecated   To be removed after Cantera 2.6. Replaceable by @see data().
     */
    const vector_fp& coeffs() const {
        warn_deprecated("ChebyshevRate3::coeffs", "Deprecated in Cantera 2.6 "
            "and to be removed thereafter; replaceable by data().");
        return chebCoeffs_;
    }

    //! Access Chebyshev coefficients as 2-dimensional array with temperature and
    //! pressure dimensions corresponding to rows and columns, respectively.
    const Array2D& data() const {
        return m_coeffs;
    }

    //! Set the Chebyshev coefficients as 2-dimensional array.
    void setData(const Array2D& coeffs);

protected:
    double m_log10P; //!< value detecting updates
    double Tmin_, Tmax_; //!< valid temperature range
    double Pmin_, Pmax_; //!< valid pressure range
    double TrNum_, TrDen_; //!< terms appearing in the reduced temperature
    double PrNum_, PrDen_; //!< terms appearing in the reduced pressure

    Array2D m_coeffs; //!<< coefficient array
    vector_fp chebCoeffs_; //!< Chebyshev coefficients, length nP * nT
    vector_fp dotProd_; //!< dot product of chebCoeffs with the reduced pressure polynomial

    Units m_rate_units; //!< Reaction rate units
};

typedef ChebyshevRate3 Chebyshev;

/**
 * A Blowers Masel rate with coverage-dependent terms.
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
 * surface coverage of species \f$ k, \theta_k \f$.
 * They are coupled with Blowers Masel approximation written by [Paul Blowers,
 * Rich Masel(2004). Engineering Approximations For Activation Energies in Hydrogen
 * Transfer Reactions. Eq (10)-(11)], which can be used to change the activation
 * based on enthalpy of the reaction:
 *
 *   \f{eqnarray*}{
 *        E_a &=& 0\; \text{if }\Delta H < -4E_0 \\
 *        E_a &=& \Delta H\; \text{if }\Delta H > 4E_0 \\
 *        E_a &=& \frac{(w + \Delta H / 2)(V_P - 2w +
 *               \Delta H)^2}{(V_P^2 - 4w^2 + (\Delta H)^2)}\; \text{Otherwise}
 *   \f}
 * where
 *   \f[
 *        V_P = 2w (w + E_0) / (w - E_0),
 *   \f]
 * \f$ w \f$ is the average bond dissociation energy of the bond breaking
 * and that being formed in the reaction. Since the expression is
 * very insensitive to \f$ w \f$ for \f$ w >= 2 E_0 \f$, \f$ w \f$
 * can be approximated to an arbitrary high value like 1000 kJ/mol.
 */
class BMSurfaceArrhenius
{

public:
    BMSurfaceArrhenius();
    explicit BMSurfaceArrhenius(double A, double b, double Ta, double w);

    //! Add a coverage dependency for species *k*, with exponential dependence
    //! *a*, power-law exponent *m*, and activation energy dependence *e*,
    //! where *e* is in Kelvin, i.e. energy divided by the molar gas constant.
    void addCoverageDependence(size_t k, double a,
                               double m, double e);

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
    doublereal updateRC(doublereal logT, doublereal recipT, doublereal deltaH) {
        double Ea = activationEnergy_R(deltaH);
        return m_A * std::exp(std::log(10.0)*m_acov + m_b*logT -
                              (Ea)*recipT + m_mcov);
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

    //! Return the activation energy divided by the gas constant (i.e. the
    //! activation temperature) [K] based on the input enthalpy
    //! accounting coverage dependence.
    doublereal activationEnergy_R(doublereal deltaH) {
        double deltaH_R = deltaH / GasConstant; // deltaH in temperature units (Kelvin)
        double E;
        if (deltaH_R < -4 * m_E0) {
            E = 0;
        } else if (deltaH_R > 4 * m_E0) {
            E = deltaH_R;
        } else {
            // m_w is in Kelvin
            // vp is in Kelvin
            double vp = 2 * m_w * ((m_w + m_E0) / (m_w - m_E0));
            double vp_2w_dH = (vp - 2 * m_w + deltaH_R); // (Vp - 2 w + dH)
            E = (m_w + deltaH_R / 2) * (vp_2w_dH * vp_2w_dH) /
                (vp * vp - 4 * m_w * m_w + deltaH_R * deltaH_R); // in Kelvin
        }
        return E + m_ecov;
    }

    //! Return the intrinsic activation energy divided by the gas constant (i.e. the
    //! activation temperature) [K], accounting coverage dependence.
    doublereal activationEnergy_R0() const {
        return m_E0 + m_ecov;
    }

    //! Return the bond energy *w* divided by the gas constant[K]
    doublereal bondEnergy() const {
        return m_w;
    }

protected:
    doublereal m_b, m_A, m_E0, m_w;
    doublereal m_acov, m_ecov, m_mcov;
    std::vector<size_t> m_sp, m_msp;
    vector_fp m_ac, m_ec, m_mc;
};

}

#endif
