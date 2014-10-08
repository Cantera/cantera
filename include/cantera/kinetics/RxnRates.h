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
    Arrhenius() :
        m_logA(-1.0E300),
        m_b(0.0),
        m_E(0.0),
        m_A(0.0) {}

    //! Constructor from ReactionData.
    explicit Arrhenius(const ReactionData& rdata) :
        m_b(rdata.rateCoeffParameters[1]),
        m_E(rdata.rateCoeffParameters[2]),
        m_A(rdata.rateCoeffParameters[0]) {
        if (m_A  <= 0.0) {
            m_logA = -1.0E300;
        } else {
            m_logA = std::log(m_A);
        }
    }

    /// Constructor.
    /// @param A pre-exponential. The unit system is
    /// (kmol, m, s). The actual units depend on the reaction
    /// order and the dimensionality (surface or bulk).
    /// @param b Temperature exponent. Non-dimensional.
    /// @param E Activation energy in temperature units. Kelvin.
    Arrhenius(doublereal A, doublereal b, doublereal E) :
        m_b(b),
        m_E(E),
        m_A(A) {
        if (m_A  <= 0.0) {
            m_logA = -1.0E300;
        } else {
            m_logA = log(m_A);
        }
    }

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

    SurfaceArrhenius() :
        m_logA(-1.0E300),
        m_b(0.0),
        m_E(0.0),
        m_A(0.0),
        m_acov(0.0),
        m_ecov(0.0),
        m_mcov(0.0),
        m_ncov(0),
        m_nmcov(0) {
    }

    explicit SurfaceArrhenius(const ReactionData& rdata) :
        m_b(rdata.rateCoeffParameters[1]),
        m_E(rdata.rateCoeffParameters[2]),
        m_A(rdata.rateCoeffParameters[0]),
        m_acov(0.0),
        m_ecov(0.0),
        m_mcov(0.0),
        m_ncov(0),
        m_nmcov(0) {
        if (m_A <= 0.0) {
            m_logA = -1.0E300;
        } else {
            m_logA = std::log(m_A);
        }

        const vector_fp& data = rdata.rateCoeffParameters;
        if (data.size() >= 7) {
            for (size_t n = 3; n < data.size()-3; n += 4) {
                addCoverageDependence(size_t(data[n]), data[n+1],
                                      data[n+2], data[n+3]);
            }
        }
    }

    void addCoverageDependence(size_t k, doublereal a,
                               doublereal m, doublereal e) {
        m_ncov++;
        m_sp.push_back(k);
        m_ac.push_back(a);
        m_ec.push_back(e);
        if (m != 0.0) {
            m_msp.push_back(k);
            m_mc.push_back(m);
            m_nmcov++;
        }
    }

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

    //! @deprecated. To be removed after Cantera 2.2
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
 */
class ExchangeCurrent
{
public:

    //! return the rate coefficient type.
    static int type() {
        return EXCHANGE_CURRENT_REACTION_RATECOEFF_TYPE;
    }

    //! Default constructor.
    ExchangeCurrent() :
        m_logA(-1.0E300),
        m_b(0.0),
        m_E(0.0),
        m_A(0.0) {}

    //! Constructor with Arrhenius parameters from a ReactionData struct.
    explicit ExchangeCurrent(const ReactionData& rdata) :
        m_b(rdata.rateCoeffParameters[1]),
        m_E(rdata.rateCoeffParameters[2]),
        m_A(rdata.rateCoeffParameters[0]) {
        if (m_A  <= 0.0) {
            m_logA = -1.0E300;
        } else {
            m_logA = std::log(m_A);
        }
    }

    /// Constructor.
    /// @param A pre-exponential. The unit system is
    /// (kmol, m, s). The actual units depend on the reaction
    /// order and the dimensionality (surface or bulk).
    /// @param b Temperature exponent. Non-dimensional.
    /// @param E Activation energy in temperature units. Kelvin.
    ExchangeCurrent(doublereal A, doublereal b, doublereal E) :
        m_b(b),
        m_E(E),
        m_A(A) {
        if (m_A  <= 0.0) {
            m_logA = -1.0E300;
        } else {
            m_logA = std::log(m_A);
        }
    }

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
    explicit Plog(const ReactionData& rdata) :
        logP_(-1000),
        logP1_(1000),
        logP2_(-1000),
        m1_(npos),
        m2_(npos),
        rDeltaP_(-1.0),
        maxRates_(1) {
        typedef std::multimap<double, vector_fp>::const_iterator iter_t;

        size_t j = 0;
        size_t rateCount = 0;
        // Insert intermediate pressures
        for (iter_t iter = rdata.plogParameters.begin();
                iter != rdata.plogParameters.end();
                iter++) {
            double logp = std::log(iter->first);
            if (pressures_.empty() || pressures_.rbegin()->first != logp) {
                // starting a new group
                pressures_[logp] = std::make_pair(j, j+1);
                rateCount = 1;
            } else {
                // another rate expression at the same pressure
                pressures_[logp].second = j+1;
                rateCount++;
            }
            maxRates_ = std::max(rateCount, maxRates_);

            j++;
            A_.push_back(iter->second[0]);
            n_.push_back(iter->second[1]);
            Ea_.push_back(iter->second[2]);
        }

        // For pressures with only one Arrhenius expression, it is more
        // efficient to work with log(A)
        for (pressureIter iter = pressures_.begin();
                iter != pressures_.end();
                iter++) {
            if (iter->second.first == iter->second.second - 1) {
                A_[iter->second.first] = std::log(A_[iter->second.first]);
            }
        }

        // Duplicate the first and last groups to handle P < P_0 and P > P_N
        pressures_.insert(std::make_pair(-1000.0, pressures_.begin()->second));
        pressures_.insert(std::make_pair(1000.0, pressures_.rbegin()->second));

        // Resize work arrays
        A1_.resize(maxRates_);
        A2_.resize(maxRates_);
        n1_.resize(maxRates_);
        n2_.resize(maxRates_);
        Ea1_.resize(maxRates_);
        Ea2_.resize(maxRates_);

        if (rdata.validate) {
            validate(rdata);
        }
    }

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
    void validate(const ReactionData& rdata) {
        double T[] = {200.0, 500.0, 1000.0, 2000.0, 5000.0, 10000.0};
        for (pressureIter iter = pressures_.begin();
                iter->first < 1000;
                iter++) {
            update_C(&iter->first);
            for (size_t i=0; i < 6; i++) {
                double k = updateRC(log(T[i]), 1.0/T[i]);
                if (!(k >= 0)) {
                    // k is NaN. Increment the iterator so that the error
                    // message will correctly indicate that the problematic rate
                    // expression is at the higher of the adjacent pressures.
                    throw CanteraError("Plog::validate",
                                       "Invalid rate coefficient for reaction #" +
                                       int2str(rdata.number) + ":\n" + rdata.equation + "\n" +
                                       "at P = " + fp2str(std::exp((++iter)->first)) +
                                       ", T = " + fp2str(T[i]));
                }
            }
        }
    }

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
    explicit ChebyshevRate(const ReactionData& rdata) :
        nP_(rdata.chebDegreeP),
        nT_(rdata.chebDegreeT),
        chebCoeffs_(rdata.chebCoeffs),
        dotProd_(rdata.chebDegreeT) {
        double logPmin = std::log10(rdata.chebPmin);
        double logPmax = std::log10(rdata.chebPmax);
        double TminInv = 1.0 / rdata.chebTmin;
        double TmaxInv = 1.0 / rdata.chebTmax;

        TrNum_ = - TminInv - TmaxInv;
        TrDen_ = 1.0 / (TmaxInv - TminInv);
        PrNum_ = - logPmin - logPmax;
        PrDen_ = 1.0 / (logPmax - logPmin);
    }

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
