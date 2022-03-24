// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_FALLOFF_H
#define CT_FALLOFF_H

#include "cantera/kinetics/Arrhenius.h"
#include "MultiRate.h"

namespace Cantera
{

class AnyMap;

/**
 *  @defgroup falloffGroup  Falloff Parameterizations
 *
 *  This section describes the parameterizations used to describe
 *  the fall-off in reaction rate constants due to intermolecular
 *  energy transfer.
 *  @ingroup chemkinetics
 */


//! Data container holding shared data specific to Falloff rates
/**
 * The data container `FalloffData` holds precalculated data common to
 * all Falloff related reaction rate classes.
 */
struct FalloffData : public ReactionData
{
    FalloffData();

    virtual bool update(const ThermoPhase& phase, const Kinetics& kin) override;

    virtual void update(double T) override;

    virtual void update(double T, double M) override;

    using ReactionData::update;

    //! Perturb third-body concentration vector of data container
    /**
     * The method is used for the evaluation of numerical derivatives.
     * @param deltaM  relative third-body perturbation
     */
    void perturbThirdBodies(double deltaM);

    virtual void restore() override;

    virtual void resize(size_t nSpecies, size_t nReactions, size_t nPhases) override {
        conc_3b.resize(nReactions, NAN);
        m_conc_3b_buf.resize(nReactions, NAN);
        ready = true;
    }

    virtual void invalidateCache() override {
        ReactionData::invalidateCache();
        molar_density = NAN;
    }

    bool ready; //!< boolean indicating whether vectors are accessible
    double molar_density; //!< used to determine if updates are needed
    vector_fp conc_3b; //!< vector of effective third-body concentrations

protected:
    int m_state_mf_number; //!< integer that is incremented when composition changes
    bool m_perturbed; //!< boolean indicating whether 3-rd body values are perturbed
    vector_fp m_conc_3b_buf; //!< buffered third-body concentrations
};


/**
 * Base class for falloff rate calculators. Each instance of a subclass of FalloffRate
 * calculates the falloff reaction rate based on specific implementations of the
 * falloff function.
 */
class FalloffRate : public ReactionRate
{
public:
    FalloffRate()
        : m_chemicallyActivated(false)
        , m_negativeA_ok(false)
        , m_rc_low(NAN)
        , m_rc_high(NAN)
    {
    }

    FalloffRate(const AnyMap& node, const UnitStack& rate_units={})
        : FalloffRate()
    {
        setParameters(node, rate_units);
    }

    /**
     * Initialize. Must be called before any other method is invoked.
     *
     * @param c Vector of coefficients of the parameterization. The number and
     *     meaning of these coefficients is subclass-dependent.
     *
     * @todo  deprecate; superseded by setFalloffCoeffs
     */
    void init(const vector_fp& c);

    /**
     * Set coefficients of the falloff parameterization.
     *
     * @param c Vector of coefficients of the parameterization. The number and
     *     meaning of these coefficients is subclass-dependent.
     */
    virtual void setFalloffCoeffs(const vector_fp& c);

    /**
     * Retrieve coefficients of the falloff parameterization.
     *
     * @param c Vector of coefficients of the parameterization. The number and
     *     meaning of these coefficients is subclass-dependent.
     */
    virtual void getFalloffCoeffs(vector_fp& c) const;

    /**
     * Update the temperature-dependent portions of the falloff function, if
     * any, and store them in the 'work' array. If not overloaded, the default
     * behavior is to do nothing.
     * @param T Temperature [K].
     * @param work storage space for intermediate results.
     */
    virtual void updateTemp(double T, double* work) const {}

    /**
     * The falloff function. This is defined so that the rate coefficient is
     *
     * \f[  k = F(Pr)\frac{Pr}{1 + Pr}. \f]
     *
     * Here \f$ Pr \f$ is the reduced pressure, defined by
     *
     * \f[
     * Pr = \frac{k_0 [M]}{k_\infty}.
     * \f]
     *
     * @param pr reduced pressure (dimensionless).
     * @param work array of size workSize() containing cached
     *             temperature-dependent intermediate results from a prior call
     *             to updateTemp.
     * @returns the value of the falloff function \f$ F \f$ defined above
     */
    virtual double F(double pr, const double* work) const {
        return 1.0;
    }

    //! Evaluate falloff function at current conditions
    double evalF(double T, double conc3b) {
        updateTemp(T, m_work.data());
        double logT = std::log(T);
        double recipT = 1. / T;
        m_rc_low = m_lowRate.evalRate(logT, recipT);
        m_rc_high = m_highRate.evalRate(logT, recipT);
        double pr = conc3b * m_rc_low / (m_rc_high + SmallNumber);
        return F(pr, m_work.data());
    }

    //! The size of the work array required.
    /**
     * @todo  deprecate; only used by legacy framework
     */
    virtual size_t workSize() const {
        return 0;
    }

    virtual const std::string type() const {
        return "Falloff";
    }

    //! Returns the number of parameters used by this parameterization. The
    //! values of these parameters can be obtained from getParameters().
    virtual size_t nParameters() const {
        return 0;
    }

    virtual void setParameters(const AnyMap& node, const UnitStack& rate_units);

    //! Get the values of the parameters for this object. *params* must be an
    //! array of at least nParameters() elements.
    /**
     * @todo  deprecate; superseded by getFalloffCoeffs
     */
    virtual void getParameters(double* params) const {}

    virtual void getParameters(AnyMap& node) const;

    //! Evaluate reaction rate
    //! @param shared_data  data shared by all reactions of a given type
    virtual double evalFromStruct(const FalloffData& shared_data) {
        updateTemp(shared_data.temperature, m_work.data());
        m_rc_low = m_lowRate.evalRate(shared_data.logT, shared_data.recipT);
        m_rc_high = m_highRate.evalRate(shared_data.logT, shared_data.recipT);
        double thirdBodyConcentration;
        if (shared_data.ready) {
            thirdBodyConcentration = shared_data.conc_3b[m_rate_index];
        } else {
            thirdBodyConcentration = shared_data.conc_3b[0];
        }
        double pr = thirdBodyConcentration * m_rc_low / (m_rc_high + SmallNumber);

        // Apply falloff function
        if (m_chemicallyActivated) {
            // 1 / (1 + Pr) * F
            pr = F(pr, m_work.data()) / (1.0 + pr);
            return pr * m_rc_low;
        }

        // Pr / (1 + Pr) * F
        pr *= F(pr, m_work.data()) / (1.0 + pr);
        return pr * m_rc_high;
    }

    void check(const std::string& equation, const AnyMap& node);
    virtual void validate(const std::string& equation, const Kinetics& kin);

    //! Get flag indicating whether negative A values are permitted
    bool allowNegativePreExponentialFactor() const {
        return m_negativeA_ok;
    }

    //! Set flag indicating whether negative A values are permitted
    void setAllowNegativePreExponentialFactor(bool value) {
        m_negativeA_ok = value;
    }

    //! Get flag indicating whether reaction is chemically activated
    bool chemicallyActivated() const {
        return m_chemicallyActivated;
    }

    //! Set flag indicating whether reaction is chemically activated
    void setChemicallyActivated(bool activated) {
        m_chemicallyActivated = activated;
    }

    //! Get reaction rate in the low-pressure limit
    ArrheniusRate& lowRate() {
        return m_lowRate;
    }

    //! Set reaction rate in the low-pressure limit
    void setLowRate(const ArrheniusRate& low);

    //! Get reaction rate in the high-pressure limit
    ArrheniusRate& highRate() {
        return m_highRate;
    }

    //! Set reaction rate in the high-pressure limit
    void setHighRate(const ArrheniusRate& high);

protected:
    ArrheniusRate m_lowRate; //!< The reaction rate in the low-pressure limit
    ArrheniusRate m_highRate; //!< The reaction rate in the high-pressure limit

    bool m_chemicallyActivated; //!< Flag labeling reaction as chemically activated
    bool m_negativeA_ok; //!< Flag indicating whether negative A values are permitted

    double m_rc_low; //!< Evaluated reaction rate in the low-pressure limit
    double m_rc_high; //!< Evaluated reaction rate in the high-pressure limit
    vector_fp m_work; //!< Work vector
};


//! The Lindemann falloff parameterization.
/**
 * This class implements the trivial falloff function F = 1.0.
 *
 * @ingroup falloffGroup
 */
class LindemannRate final : public FalloffRate
{
public:
    LindemannRate() = default;

    LindemannRate(const AnyMap& node, const UnitStack& rate_units={})
        : LindemannRate()
    {
        setParameters(node, rate_units);
    }

    LindemannRate(
        const ArrheniusRate& low, const ArrheniusRate& high, const vector_fp& c)
        : LindemannRate()
    {
        m_lowRate = low;
        m_highRate = high;
        setFalloffCoeffs(c);
    }

    unique_ptr<MultiRateBase> newMultiRate() const {
        return unique_ptr<MultiRateBase>(
            new MultiRate<LindemannRate, FalloffData>);
    }

    virtual const std::string type() const {
        return "Lindemann";
    }
};


//! The 3- or 4-parameter Troe falloff parameterization.
/*!
 * The falloff function defines the value of \f$ F \f$ in the following
 * rate expression
 *
 *  \f[ k = k_{\infty} \left( \frac{P_r}{1 + P_r} \right) F \f]
 *  where
 *  \f[ P_r = \frac{k_0 [M]}{k_{\infty}} \f]
 *
 * This parameterization is defined by
 *
 * \f[ F = F_{cent}^{1/(1 + f_1^2)} \f]
 *    where
 * \f[ F_{cent} = (1 - A)\exp(-T/T_3) + A \exp(-T/T_1) + \exp(-T_2/T) \f]
 *
 * \f[ f_1 = (\log_{10} P_r + C) /
 *              \left(N - 0.14 (\log_{10} P_r + C)\right) \f]
 *
 * \f[ C = -0.4 - 0.67 \log_{10} F_{cent} \f]
 *
 * \f[ N = 0.75 - 1.27 \log_{10} F_{cent} \f]
 *
 *  - If \f$ T_3 \f$ is zero, then the corresponding term is set to zero.
 *  - If \f$ T_1 \f$ is zero, then the corresponding term is set to zero.
 *  - If \f$ T_2 \f$ is zero, then the corresponding term is set to zero.
 *
 * @ingroup falloffGroup
 */
class TroeRate final : public FalloffRate
{
public:
    //! Constructor
    TroeRate() : m_a(NAN), m_rt3(0.0), m_rt1(0.0), m_t2(0.0) {
        m_work.resize(1);
    }

    TroeRate(const AnyMap& node, const UnitStack& rate_units={})
        : TroeRate()
    {
        setParameters(node, rate_units);
    }

    TroeRate(const ArrheniusRate& low, const ArrheniusRate& high, const vector_fp& c)
        : TroeRate()
    {
        m_lowRate = low;
        m_highRate = high;
        setFalloffCoeffs(c);
    }

    unique_ptr<MultiRateBase> newMultiRate() const {
        return unique_ptr<MultiRateBase>(new MultiRate<TroeRate, FalloffData>);
    }

    //! Set coefficients used by parameterization
    /*!
     * @param c Vector of three or four doubles: The doubles are the parameters,
     *          a, T_3, T_1, and (optionally) T_2 of the Troe parameterization
     */
    virtual void setFalloffCoeffs(const vector_fp& c);

    virtual void getFalloffCoeffs(vector_fp& c) const;

    //! Update the temperature parameters in the representation
    /*!
     *   @param T         Temperature (Kelvin)
     *   @param work      Vector of working space, length 1, representing the
     *                    temperature-dependent part of the parameterization.
     */
    virtual void updateTemp(double T, double* work) const;

    virtual double F(double pr, const double* work) const;

    virtual size_t workSize() const {
        return 1;
    }

    virtual const std::string type() const {
        return "Troe";
    }

    virtual size_t nParameters() const {
        return 4;
    }

    virtual void setParameters(const AnyMap& node, const UnitStack& rate_units);

    //! Sets params to contain, in order, \f[ (A, T_3, T_1, T_2) \f]
    /**
     * @todo  deprecate; superseded by getFalloffCoeffs
     */
    virtual void getParameters(double* params) const;

    virtual void getParameters(AnyMap& node) const;

protected:
    //! parameter a in the 4-parameter Troe falloff function. Dimensionless
    double m_a;

    //! parameter 1/T_3 in the 4-parameter Troe falloff function. [K^-1]
    double m_rt3;

    //! parameter 1/T_1 in the 4-parameter Troe falloff function. [K^-1]
    double m_rt1;

    //! parameter T_2 in the 4-parameter Troe falloff function. [K]
    double m_t2;
};

//! The SRI falloff function
/*!
 * The falloff function defines the value of \f$ F \f$ in the following
 * rate expression
 *
 *  \f[ k = k_{\infty} \left( \frac{P_r}{1 + P_r} \right) F \f]
 *  where
 *  \f[ P_r = \frac{k_0 [M]}{k_{\infty}} \f]
 *
 *  \f[ F = {\left( a \; exp(\frac{-b}{T}) + exp(\frac{-T}{c})\right)}^n
 *              \;  d \; T^e \f]
 *      where
 *  \f[ n = \frac{1.0}{1.0 + (\log_{10} P_r)^2} \f]
 *
 *  \f$ c \f$ s required to greater than or equal to zero. If it is zero, then
 *  the corresponding term is set to zero.
 *
 *  \f$ d \f$ is required to be greater than zero.
 *
 * @ingroup falloffGroup
 */
class SriRate final : public FalloffRate
{
public:
    //! Constructor
    SriRate() : m_a(NAN), m_b(-1.0), m_c(-1.0), m_d(-1.0), m_e(-1.0) {
        m_work.resize(2);
    }

    SriRate(const AnyMap& node, const UnitStack& rate_units={})
        : SriRate()
    {
        setParameters(node, rate_units);
    }

    SriRate(const ArrheniusRate& low, const ArrheniusRate& high, const vector_fp& c)
        : SriRate()
    {
        m_lowRate = low;
        m_highRate = high;
        setFalloffCoeffs(c);
    }

    unique_ptr<MultiRateBase> newMultiRate() const {
        return unique_ptr<MultiRateBase>(new MultiRate<SriRate, FalloffData>);
    }

    //! Set coefficients used by parameterization
    /*!
     * @param c Vector of three or five doubles: The doubles are the parameters,
     *          a, b, c, d (optional; default 1.0), and e (optional; default
     *          0.0) of the SRI parameterization
     */
    virtual void setFalloffCoeffs(const vector_fp& c);

    virtual void getFalloffCoeffs(vector_fp& c) const;

    //! Update the temperature parameters in the representation
    /*!
     *   @param T         Temperature (Kelvin)
     *   @param work      Vector of working space, length 2, representing the
     *                    temperature-dependent part of the parameterization.
     */
    virtual void updateTemp(double T, double* work) const;

    virtual double F(double pr, const double* work) const;

    virtual size_t workSize() const {
        return 2;
    }

    virtual const std::string type() const {
        return "SRI";
    }

    virtual size_t nParameters() const {
        return 5;
    }

    virtual void setParameters(const AnyMap& node, const UnitStack& rate_units);

    //! Sets params to contain, in order, \f[ (a, b, c, d, e) \f]
    /**
     * @todo  deprecate; superseded by getFalloffCoeffs
     */
    virtual void getParameters(double* params) const;

    virtual void getParameters(AnyMap& node) const;

protected:
    //! parameter a in the 5-parameter SRI falloff function. Dimensionless.
    double m_a;

    //! parameter b in the 5-parameter SRI falloff function. [K]
    double m_b;

    //! parameter c in the 5-parameter SRI falloff function. [K]
    double m_c;

    //! parameter d in the 5-parameter SRI falloff function. Dimensionless.
    double m_d;

    //! parameter d in the 5-parameter SRI falloff function. Dimensionless.
    double m_e;
};

//! The 1- or 2-parameter Tsang falloff parameterization.
/*!
 *  The Tsang falloff model is adapted from that of Troe.
 *  It provides a constant or linear in temperature value for \f$ F_{cent} \f$:
 *  \f[ F_{cent} = A + B*T \f]
 *
 *  The value of \f$ F_{cent} \f$ is then applied to Troe's model for the
 *  determination of the value of \f$ F \f$:
 * \f[ F = F_{cent}^{1/(1 + f_1^2)} \f]
 *    where
 * \f[ f_1 = (\log_{10} P_r + C) /
 *              \left(N - 0.14 (\log_{10} P_r + C)\right) \f]
 *
 * \f[ C = -0.4 - 0.67 \log_{10} F_{cent} \f]
 *
 * \f[ N = 0.75 - 1.27 \log_{10} F_{cent} \f]
 *
 *  References:
 *  Example of reaction database developed by Tsang utilizing this format:
 *      https://doi.org/10.1063/1.555890
 *  Example of Chemkin implementation of Tsang format (supplemental materials):
 *      https://doi.org/10.1016/j.combustflame.2011.02.010
 *
 * @ingroup falloffGroup
 */
class TsangRate final : public FalloffRate
{
public:
    //! Constructor
    TsangRate() : m_a(NAN), m_b(0.0) {
        m_work.resize(1);
    }

    TsangRate(const AnyMap& node, const UnitStack& rate_units={})
        : TsangRate()
    {
        setParameters(node, rate_units);
    }

    TsangRate(const ArrheniusRate& low, const ArrheniusRate& high, const vector_fp& c)
        : TsangRate()
    {
        m_lowRate = low;
        m_highRate = high;
        setFalloffCoeffs(c);
    }

    unique_ptr<MultiRateBase> newMultiRate() const {
        return unique_ptr<MultiRateBase>(new MultiRate<TsangRate, FalloffData>);
    }

    //! Set coefficients used by parameterization
    /*!
     * @param c Vector of one or two doubles: The doubles are the parameters,
     *          a and (optionally) b of the Tsang F_cent parameterization
     */
    virtual void setFalloffCoeffs(const vector_fp& c);

    virtual void getFalloffCoeffs(vector_fp& c) const;

    //! Update the temperature parameters in the representation
    /*!
     *   @param T         Temperature (Kelvin)
     *   @param work      Vector of working space, length 1, representing the
     *                    temperature-dependent part of the parameterization.
     */
    virtual void updateTemp(double T, double* work) const;

    virtual double F(double pr, const double* work) const;

    virtual size_t workSize() const {
        return 1;
    }

    virtual const std::string type() const {
        return "Tsang";
    }

    virtual size_t nParameters() const {
        return 2;
    }

    virtual void setParameters(const AnyMap& node, const UnitStack& rate_units);

    //! Sets params to contain, in order, \f[ (A, B) \f]
    /**
     * @todo  deprecate; superseded by getFalloffCoeffs
     */
    virtual void getParameters(double* params) const;

    virtual void getParameters(AnyMap& node) const;

protected:
    //! parameter a in the Tsang F_cent formulation. Dimensionless
    double m_a;

    //! parameter b in the Tsang F_cent formulation. [K^-1]
    double m_b;
};

typedef FalloffRate Falloff;
typedef LindemannRate Lindemann;
typedef TroeRate Troe;
typedef SriRate SRI;
typedef TsangRate Tsang;

}

#endif
