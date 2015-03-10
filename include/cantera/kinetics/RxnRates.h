/**
 *  @file RxnRates.h
 *
 */
// Copyright 2001  California Institute of Technology


#ifndef CT_RXNRATES_H
#define CT_RXNRATES_H

#include "ReactionData.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/stringUtils.h"

#include <iostream>

namespace Cantera
{

class Array2D;

//! Arrhenius reaction rate type depends only on temperature
/**
 * A reaction rate coefficient of the following form.
 *
 *   \f[
 *        k_f =  A T^b \exp (-E/RT)
 *   \f]
 *
 */
class Arrhenius
{
public:
    //! return the rate coefficient type.
    static int type() {
        return ARRHENIUS_REACTION_RATECOEFF_TYPE;
    }

    //! Default constructor.
    Arrhenius();

    //! Constructor from ReactionData.
    explicit Arrhenius(const ReactionData& rdata);

    /// Constructor.
    /// @param A pre-exponential. The unit system is
    /// (kmol, m, s). The actual units depend on the reaction
    /// order and the dimensionality (surface or bulk).
    /// @param b Temperature exponent. Non-dimensional.
    /// @param E Activation energy in temperature units. Kelvin.
    Arrhenius(doublereal A, doublereal b, doublereal E);

    //! Update concentration-dependent parts of the rate coefficient.
    /*!
     *   For this class, there are no
     *   concentration-dependent parts, so this method does  nothing.
     */
    void update_C(const doublereal* c) {
    }

    /**
     * Update the value of the logarithm of the rate constant.
     *
     * Note, this function should never be called for negative A values.
     * If it does then it will produce a negative overflow result, and
     * a zero net forwards reaction rate, instead of a negative reaction
     * rate constant that is the expected result.
     * @deprecated. To be removed after Cantera 2.2
     */
    doublereal update(doublereal logT, doublereal recipT) const {
        return m_logA + m_b*logT - m_E*recipT;
    }

    /**
     * Update the value the rate constant.
     *
     * This function returns the actual value of the rate constant.
     * It can be safely called for negative values of the pre-exponential
     * factor.
     */
    doublereal updateRC(doublereal logT, doublereal recipT) const {
        return m_A * std::exp(m_b*logT - m_E*recipT);
    }

    //! @deprecated. To be removed after Cantera 2.2
    void writeUpdateRHS(std::ostream& s) const {
        s << " exp(" << m_logA;
        if (m_b != 0.0) {
            s << " + " << m_b << " * tlog";
        }
        if (m_E != 0.0) {
            s << " - " << m_E << " * rt";
        }
        s << ");" << std::endl;
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

    //! Return the activation energy divided by the gas constant (i.e. the
    //! activation temperature) [K]
    doublereal activationEnergy_R() const {
        return m_E;
    }

    //! @deprecated. To be removed after Cantera 2.2
    static bool alwaysComputeRate() {
        return false;
    }

protected:
    doublereal m_logA, m_b, m_E, m_A;
};


/**
 * An Arrhenius rate with coverage-dependent terms.
 *
 * The rate expression is given by:
 * \f[
 *     k_f = A T^b \exp \left(
 *             \ln 10 \sum a_k \theta_k
 *             - \frac{1}{RT} \left( E_a + \sum E_k\theta_k \right)
 *             + \sum m_k \ln \theta_k
 *             \right)
 *   \f]
 * where the parameters \f$ (a_k, E_k, m_k) \f$ describe the dependency on the
 * surface coverage of species \f$k, \theta_k \f$.
 */
class SurfaceArrhenius
{

public:
    static int type() {
        return SURF_ARRHENIUS_REACTION_RATECOEFF_TYPE;
    }

    SurfaceArrhenius();
    explicit SurfaceArrhenius(double A, double b, double Ta);
    explicit SurfaceArrhenius(const ReactionData& rdata);

    void addCoverageDependence(size_t k, doublereal a,
                               doublereal m, doublereal e);

    void update_C(const doublereal* theta) {
        m_acov = 0.0;
        m_ecov = 0.0;
        m_mcov = 0.0;
        size_t k;
        doublereal th;
        for (size_t n = 0; n < m_ncov; n++) {
            k = m_sp[n];
            m_acov += m_ac[n] * theta[k];
            m_ecov += m_ec[n] * theta[k];
        }
        for (size_t n = 0; n < m_nmcov; n++) {
            k = m_msp[n];
            // changed n to k, dgg 1/22/04
            th = std::max(theta[k], Tiny);
            //                th = fmaxx(theta[n], Tiny);
            m_mcov += m_mc[n]*std::log(th);
        }
    }

    /**
     * Update the value of the logarithm of the rate constant.
     *
     * This calculation is not safe for negative values of
     * the preexponential.
     * @deprecated. To be removed after Cantera 2.2
     */
    doublereal update(doublereal logT, doublereal recipT) const {
        return m_logA + m_acov + m_b*logT
               - (m_E + m_ecov)*recipT + m_mcov;
    }

    /**
     * Update the value the rate constant.
     *
     * This function returns the actual value of the rate constant.
     * It can be safely called for negative values of the pre-exponential
     * factor.
     */
    doublereal updateRC(doublereal logT, doublereal recipT) const {
        return m_A * std::exp(std::log(10.0)*m_acov + m_b*logT -
                              (m_E + m_ecov)*recipT + m_mcov);
    }

    doublereal activationEnergy_R() const {
        return m_E + m_ecov;
    }

    //! @deprecated. To be removed after Cantera 2.2
    static bool alwaysComputeRate() {
        return true;
    }

protected:
    doublereal m_logA, m_b, m_E, m_A;
    doublereal m_acov, m_ecov, m_mcov;
    std::vector<size_t> m_sp, m_msp;
    vector_fp m_ac, m_ec, m_mc;
    size_t m_ncov, m_nmcov;
};


//! Arrhenius reaction rate type depends only on temperature
/**
 * A reaction rate coefficient of the following form.
 *
 *   \f[
 *        k_f =  A T^b \exp (-E/RT)
 *   \f]
 *
 * @deprecated Duplicate of class Arrhenius. To be removed after Cantera 2.2.
 */
class ExchangeCurrent
{
public:

    //! return the rate coefficient type.
    static int type() {
        return EXCHANGE_CURRENT_REACTION_RATECOEFF_TYPE;
    }

    //! Default constructor.
    ExchangeCurrent();

    //! Constructor with Arrhenius parameters from a ReactionData struct.
    explicit ExchangeCurrent(const ReactionData& rdata);

    /// Constructor.
    /// @param A pre-exponential. The unit system is
    /// (kmol, m, s). The actual units depend on the reaction
    /// order and the dimensionality (surface or bulk).
    /// @param b Temperature exponent. Non-dimensional.
    /// @param E Activation energy in temperature units. Kelvin.
    ExchangeCurrent(doublereal A, doublereal b, doublereal E);

    //! Update concentration-dependent parts of the rate coefficient.
    /*!
     *   For this class, there are no
     *   concentration-dependent parts, so this method does  nothing.
     */
    void update_C(const doublereal* c) {
    }

    /**
     * Update the value of the logarithm of the rate constant.
     *
     * Note, this function should never be called for negative A values.
     * If it does then it will produce a negative overflow result, and
     * a zero net forwards reaction rate, instead of a negative reaction
     * rate constant that is the expected result.
     * @deprecated. To be removed after Cantera 2.2
     */
    doublereal update(doublereal logT, doublereal recipT) const {
        return m_logA + m_b*logT - m_E*recipT;
    }

    /**
     * Update the value the rate constant.
     *
     * This function returns the actual value of the rate constant.
     * It can be safely called for negative values of the pre-exponential
     * factor.
     */
    doublereal updateRC(doublereal logT, doublereal recipT) const {
        return m_A * std::exp(m_b*logT - m_E*recipT);
    }

    //! @deprecated. To be removed after Cantera 2.2
    void writeUpdateRHS(std::ostream& s) const {
        s << " exp(" << m_logA;
        if (m_b != 0.0) {
            s << " + " << m_b << " * tlog";
        }
        if (m_E != 0.0) {
            s << " - " << m_E << " * rt";
        }
        s << ");" << std::endl;
    }

    //! @deprecated. To be removed after Cantera 2.2
    doublereal activationEnergy_R() const {
        return m_E;
    }

    //! @deprecated. To be removed after Cantera 2.2
    static bool alwaysComputeRate() {
        return false;
    }

protected:
    doublereal m_logA, m_b, m_E, m_A;
};


class Plog
{
public:
    //! return the rate coefficient type.
    static int type() {
        return PLOG_REACTION_RATECOEFF_TYPE;
    }

    //! Default constructor.
    Plog() {}

    //! Constructor from ReactionData.
    explicit Plog(const ReactionData& rdata);

    //! Constructor from Arrhenius rate expressions at a set of pressures
    explicit Plog(const std::multimap<double, Arrhenius>& rates);

    //! Update concentration-dependent parts of the rate coefficient.
    //! @param c natural log of the pressure in Pa
    void update_C(const doublereal* c) {
        logP_ = c[0];
        if (logP_ > logP1_ && logP_ < logP2_) {
            return;
        }

        pressureIter iter = pressures_.upper_bound(c[0]);
        AssertThrowMsg(iter != pressures_.end(), "Plog::update_C",
                       "Pressure out of range: " + fp2str(logP_));
        AssertThrowMsg(iter != pressures_.begin(), "Plog::update_C",
                       "Pressure out of range: " + fp2str(logP_));

        // upper interpolation pressure
        logP2_ = iter->first;
        size_t start = iter->second.first;
        m2_ = iter->second.second - start;
        for (size_t m = 0; m < m2_; m++) {
            A2_[m] = A_[start+m];
            n2_[m] = n_[start+m];
            Ea2_[m] = Ea_[start+m];
        }

        // lower interpolation pressure
        logP1_ = (--iter)->first;
        start = iter->second.first;
        m1_ = iter->second.second - start;
        for (size_t m = 0; m < m1_; m++) {
            A1_[m] = A_[start+m];
            n1_[m] = n_[start+m];
            Ea1_[m] = Ea_[start+m];
        }

        rDeltaP_ = 1.0 / (logP2_ - logP1_);
    }

    /**
     * Update the value of the logarithm of the rate constant.
     * @deprecated. To be removed after Cantera 2.2
     */
    doublereal update(doublereal logT, doublereal recipT) const {
        return std::log(updateRC(logT, recipT));
    }

    /**
     * Update the value the rate constant.
     *
     * This function returns the actual value of the rate constant.
     */
    doublereal updateRC(doublereal logT, doublereal recipT) const {
        double log_k1, log_k2;
        if (m1_ == 1) {
            log_k1 = A1_[0] + n1_[0] * logT - Ea1_[0] * recipT;
        } else {
            double k = 1e-300; // non-zero to make log(k) finite
            for (size_t m = 0; m < m1_; m++) {
                k += A1_[m] * std::exp(n1_[m] * logT - Ea1_[m] * recipT);
            }
            log_k1 = std::log(k);
        }

        if (m2_ == 1) {
            log_k2 = A2_[0] + n2_[0] * logT - Ea2_[0] * recipT;
        } else {
            double k = 1e-300; // non-zero to make log(k) finite
            for (size_t m = 0; m < m2_; m++) {
                k += A2_[m] * std::exp(n2_[m] * logT - Ea2_[m] * recipT);
            }
            log_k2 = std::log(k);
        }

        return std::exp(log_k1 + (log_k2-log_k1) * (logP_-logP1_) * rDeltaP_);
    }

    //! @deprecated. To be removed after Cantera 2.2
    doublereal activationEnergy_R() const {
        throw CanteraError("Plog::activationEnergy_R", "Not implemented");
    }

    //! @deprecated. To be removed after Cantera 2.2
    static bool alwaysComputeRate() {
        return false;
    }

    //! Check to make sure that the rate expression is finite over a range of
    //! temperatures at each interpolation pressure. This is potentially an
    //! issue when one of the Arrhenius expressions at a particular pressure
    //! has a negative pre-exponential factor.
    void validate(const std::string& equation);

protected:
    //! log(p) to (index range) in A_, n, Ea vectors
    std::map<double, std::pair<size_t, size_t> > pressures_;
    typedef std::map<double, std::pair<size_t, size_t> >::iterator pressureIter;

    vector_fp A_; //!< Pre-exponential factor at each pressure (or log(A))
    vector_fp n_; //!< Temperature exponent at each pressure [dimensionless]
    vector_fp Ea_; //!< Activation energy at each pressure [K]

    double logP_; //!< log(p) at the current state
    double logP1_, logP2_; //!< log(p) at the lower / upper pressure reference

    //! Pre-exponential factors at lower / upper pressure reference.
    //! Stored as log(A) when there is only one at the corresponding pressure.
    vector_fp A1_, A2_;
    vector_fp n1_, n2_; //!< n at lower / upper pressure reference
    vector_fp Ea1_, Ea2_; //!< Activation energy at lower / upper pressure reference

    //! Number of Arrhenius expressions at lower / upper pressure references
    size_t m1_, m2_;
    double rDeltaP_; //!< reciprocal of (logP2 - logP1)

    size_t maxRates_; //!< The maximum number of rates at any given pressure
};


class ChebyshevRate
{
public:
    //! return the rate coefficient type.
    static int type() {
        return CHEBYSHEV_REACTION_RATECOEFF_TYPE;
    }

    //! Default constructor.
    ChebyshevRate() {}

    //! Constructor from ReactionData.
    explicit ChebyshevRate(const ReactionData& rdata);

    //! Constructor directly from coefficient array
    /*
     *  @param Pmin    Minimum pressure [Pa]
     *  @param Pmax    Maximum pressure [Pa]
     *  @param Tmin    Minimum temperature [K]
     *  @param Tmax    Maximum temperature [K]
     *  @param coeffs  Coefficient array dimensioned `nT` by `nP` where `nT` and
     *      `nP` are the number of temperatures and pressures used in the fit,
     *      respectively.
     */
    ChebyshevRate(double Pmin, double Pmax, double Tmin, double Tmax,
                  const Array2D& coeffs);

    //! Update concentration-dependent parts of the rate coefficient.
    //! @param c base-10 logarithm of the pressure in Pa
    void update_C(const doublereal* c) {
        double Pr = (2 * c[0] + PrNum_) * PrDen_;
        double Cnm1 = 1;
        double Cn = Pr;
        double Cnp1;
        for (size_t j = 0; j < nT_; j++) {
            dotProd_[j] = chebCoeffs_[nP_*j] + Pr * chebCoeffs_[nP_*j+1];
        }
        for (size_t i = 2; i < nP_; i++) {
            Cnp1 = 2 * Pr * Cn - Cnm1;
            for (size_t j = 0; j < nT_; j++) {
                dotProd_[j] += Cnp1 * chebCoeffs_[nP_*j + i];
            }
            Cnm1 = Cn;
            Cn = Cnp1;
        }
    }

    /**
     * Update the value of the base-10 logarithm of the rate constant.
     * @deprecated. To be removed after Cantera 2.2
     */
    doublereal update(doublereal logT, doublereal recipT) const {
        return std::log10(updateRC(logT, recipT));
    }

    /**
     * Update the value the rate constant.
     *
     * This function returns the actual value of the rate constant.
     */
    doublereal updateRC(doublereal logT, doublereal recipT) const {
        double Tr = (2 * recipT + TrNum_) * TrDen_;
        double Cnm1 = 1;
        double Cn = Tr;
        double Cnp1;
        double logk = dotProd_[0] + Tr * dotProd_[1];
        for (size_t i = 2; i < nT_; i++) {
            Cnp1 = 2 * Tr * Cn - Cnm1;
            logk += Cnp1 * dotProd_[i];
            Cnm1 = Cn;
            Cn = Cnp1;
        }
        return std::pow(10, logk);
    }

    //! @deprecated. To be removed after Cantera 2.2
    doublereal activationEnergy_R() const {
        return 0.0;
    }

    //! @deprecated. To be removed after Cantera 2.2
    static bool alwaysComputeRate() {
        return false;
    }

protected:
    double TrNum_, TrDen_; //!< terms appearing in the reduced temperature
    double PrNum_, PrDen_; //!< terms appearing in the reduced pressure

    size_t nP_; //!< number of points in the pressure direction
    size_t nT_; //!< number of points in the temperature direction
    vector_fp chebCoeffs_; //!< Chebyshev coefficients, length nP * nT
    vector_fp dotProd_; //!< dot product of chebCoeffs with the reduced pressure polynomial
};

}

#endif
