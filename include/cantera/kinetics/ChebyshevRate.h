//! @file ChebyshevRate.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_CHEBYSHEV_H
#define CT_CHEBYSHEV_H

#include "cantera/kinetics/ReactionRate.h"
#include "cantera/kinetics/ReactionData.h"
#include "cantera/kinetics/MultiRate.h"
#include "cantera/base/Array.h"
#include "cantera/base/global.h"

namespace Cantera
{

//! Data container holding shared data specific to ChebyshevRate
/**
 * The data container `ChebyshevData` holds precalculated data common to
 * all `ChebyshevRate` objects.
 */
struct ChebyshevData : public ReactionData
{
    ChebyshevData() : pressure(NAN), log10P(0.), m_pressure_buf(-1.) {}

    virtual void update(double T) override;

    virtual void update(double T, double P) override {
        ReactionData::update(T);
        pressure = P;
        log10P = std::log10(P);
    }

    virtual bool update(const ThermoPhase& phase, const Kinetics& kin) override;

    using ReactionData::update;

    //! Perturb pressure of data container
    /**
     * The method is used for the evaluation of numerical derivatives.
     * @param  deltaP  relative pressure perturbation
     */
    void perturbPressure(double deltaP);

    virtual void restore() override;

    virtual void invalidateCache() override {
        ReactionData::invalidateCache();
        pressure = NAN;
    }

    double pressure; //!< pressure
    double log10P; //!< base 10 logarithm of pressure

protected:
    double m_pressure_buf; //!< buffered pressure
};

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
 * A ChebyshevRate rate expression is specified in terms of the coefficient matrix
 * \f$ \alpha \f$ and the temperature and pressure ranges. Note that the
 * Chebyshev polynomials are not defined outside the interval (-1,1), and
 * therefore extrapolation of rates outside the range of temperatures and
 * pressures for which they are defined is strongly discouraged.
 */
class ChebyshevRate final : public ReactionRate
{
public:
    //! Default constructor.
    ChebyshevRate() : m_log10P(NAN), m_rate_units(Units(0.)) {}

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
    ChebyshevRate(double Tmin, double Tmax, double Pmin, double Pmax,
                  const Array2D& coeffs);

    ChebyshevRate(const AnyMap& node, const UnitStack& rate_units={})
        : ChebyshevRate()
    {
        setParameters(node, rate_units);
    }

    unique_ptr<MultiRateBase> newMultiRate() const {
        return unique_ptr<MultiRateBase>(
            new MultiRate<ChebyshevRate, ChebyshevData>);
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

    virtual void validate(const std::string& equation, const Kinetics& kin);

    //! Update information specific to reaction
    /*!
     *  @param shared_data  data shared by all reactions of a given type
     */
    void updateFromStruct(const ChebyshevData& shared_data) {
        if (shared_data.log10P != m_log10P) {
            update_C(&shared_data.log10P);
        }
    }

    //! Evaluate reaction rate
    /*!
     *  @param shared_data  data shared by all reactions of a given type
     */
    double evalFromStruct(const ChebyshevData& shared_data) {
        return updateRC(0., shared_data.recipT);
    }

    //! Set up ChebyshevRate object
    /*!
     * @deprecated   Deprecated in Cantera 2.6. Replaceable with
     *               @see setLimits() and @see setCoeffs().
     */
    void setup(double Tmin, double Tmax, double Pmin, double Pmax,
                  const Array2D& coeffs);

    //! Set limits for ChebyshevRate object
    /*!
     *  @param Tmin    Minimum temperature [K]
     *  @param Tmax    Maximum temperature [K]
     *  @param Pmin    Minimum pressure [Pa]
     *  @param Pmax    Maximum pressure [Pa]
     */
    void setLimits(double Tmin, double Tmax, double Pmin, double Pmax);

    //! Update concentration-dependent parts of the rate coefficient.
    //! @param c base-10 logarithm of the pressure in Pa
    //! @deprecated To be removed after Cantera 2.6. Implementation will be moved to
    //! the updateFromStruct() method.
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
     * @deprecated To be removed after Cantera 2.6. Implementation will be moved to
     * the evalFromStruct() method.
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

    //! Access the ChebyshevRate coefficients.
    /*!
     *  \f$ \alpha_{t,p} = \mathrm{coeffs}[N_P*t + p] \f$ where
     *  \f$ 0 <= t < N_T \f$ and \f$ 0 <= p < N_P \f$.
     *
     * @deprecated   To be removed after Cantera 2.6. Replaceable by @see data().
     */
    const vector_fp& coeffs() const {
        warn_deprecated("ChebyshevRate::coeffs", "Deprecated in Cantera 2.6 "
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

}

#endif
