// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_FALLOFF_H
#define CT_FALLOFF_H

#include "cantera/kinetics/Arrhenius.h"
#include "MultiRate.h"

namespace Cantera
{

class AnyMap;

//! Data container holding shared data specific to Falloff rates
/**
 * The data container `FalloffData` holds precalculated data common to
 * all Falloff related reaction rate classes.
 */
struct FalloffData : public ReactionData
{
    FalloffData();

    bool update(const ThermoPhase& phase, const Kinetics& kin) override;

    void update(double T) override;

    void update(double T, double M) override;

    using ReactionData::update;

    //! Perturb third-body concentration vector of data container
    /**
     * The method is used for the evaluation of numerical derivatives.
     * @param deltaM  relative third-body perturbation
     */
    void perturbThirdBodies(double deltaM);

    void restore() override;

    void resize(size_t nSpecies, size_t nReactions, size_t nPhases) override {
        conc_3b.resize(nReactions, NAN);
        m_conc_3b_buf.resize(nReactions, NAN);
        ready = true;
    }

    void invalidateCache() override {
        ReactionData::invalidateCache();
        molar_density = NAN;
    }

    bool ready = false; //!< boolean indicating whether vectors are accessible
    double molar_density = NAN; //!< used to determine if updates are needed
    vector<double> conc_3b; //!< vector of effective third-body concentrations

protected:
    //! integer that is incremented when composition changes
    int m_state_mf_number = -1;
    //! boolean indicating whether 3-rd body values are perturbed
    bool m_perturbed = false;
    vector<double> m_conc_3b_buf; //!< buffered third-body concentrations
};


/**
 * Base class for falloff rate calculators.
 * Each instance of a subclass of FalloffRate calculates the falloff reaction rate
 * based on specific implementations of the falloff function.
 *
 * The falloff function @f$ F(P_r, T) @f$ is implemented by FalloffRate specializations,
 * and is defined so that the rate coefficient is
 * @f[
 *  k = k_\infty \frac{P_r}{1 + P_r} F(P_r,T)
 * @f]
 *
 * Here @f$ P_r @f$ is the reduced pressure, defined by
 * @f[
 *  P_r = \frac{k_0 [M]}{k_\infty}.
 * @f]
 * @ingroup falloffGroup
 */
class FalloffRate : public ReactionRate
{
public:
    FalloffRate() = default;

    FalloffRate(const AnyMap& node, const UnitStack& rate_units={});

    /**
     * Set coefficients of the falloff parameterization.
     *
     * @param c Vector of coefficients of the parameterization. The number and
     *     meaning of these coefficients is subclass-dependent.
     */
    virtual void setFalloffCoeffs(const vector<double>& c);

    /**
     * Retrieve coefficients of the falloff parameterization.
     *
     * @param c Vector of coefficients of the parameterization. The number and
     *     meaning of these coefficients is subclass-dependent.
     */
    virtual void getFalloffCoeffs(vector<double>& c) const;

    /**
     * Update the temperature-dependent portions of the falloff function, if
     * any, and store them in the 'work' array. If not overloaded, the default
     * behavior is to do nothing.
     * @param T Temperature [K].
     * @param work storage space for intermediate results.
     */
    virtual void updateTemp(double T, double* work) const {}

    /**
     * The falloff function.
     *
     * @param pr reduced pressure (dimensionless).
     * @param work array of size workSize() containing cached
     *             temperature-dependent intermediate results from a prior call
     *             to updateTemp.
     * @returns the value of the falloff function @f$ F @f$ defined above
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

    const string type() const override {
        if (m_chemicallyActivated) {
            return "chemically-activated";
        }
        return "falloff";
    }

    //! Returns the number of parameters used by this parameterization. The
    //! values of these parameters can be obtained from getParameters().
    virtual size_t nParameters() const {
        return 0;
    }

    void setParameters(const AnyMap& node, const UnitStack& rate_units) override;

    void getParameters(AnyMap& node) const override;

    //! Evaluate reaction rate
    //! @param shared_data  data shared by all reactions of a given type
    double evalFromStruct(const FalloffData& shared_data) {
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

    void check(const string& equation) override;
    void validate(const string& equation, const Kinetics& kin) override;

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

    //! Flag labeling reaction as chemically activated
    bool m_chemicallyActivated = false;
    //! Flag indicating whether negative A values are permitted
    bool m_negativeA_ok = false;

    double m_rc_low = NAN; //!< Evaluated reaction rate in the low-pressure limit
    double m_rc_high = NAN; //!< Evaluated reaction rate in the high-pressure limit
    vector<double> m_work; //!< Work vector
};


//! The Lindemann falloff parameterization.
/**
 * This class implements the trivial falloff function F = 1.0 @cite lindemann1922.
 *
 * @ingroup falloffGroup
 */
class LindemannRate final : public FalloffRate
{
public:
    LindemannRate() = default;

    LindemannRate(const AnyMap& node, const UnitStack& rate_units={});

    LindemannRate(const ArrheniusRate& low, const ArrheniusRate& high,
                  const vector<double>& c);

    unique_ptr<MultiRateBase> newMultiRate() const override{
        return make_unique<MultiRate<LindemannRate, FalloffData>>();
    }

    const string subType() const override {
        return "Lindemann";
    }
};


//! The 3- or 4-parameter Troe falloff parameterization.
/*!
 * The falloff function defines the value of @f$ F @f$ in the following
 * rate expression @cite gilbert1983
 *
 * @f[
 *  k = k_{\infty} \left( \frac{P_r}{1 + P_r} \right) F(T, P_r)
 * @f]
 * where
 * @f[
 *  P_r = \frac{k_0 [M]}{k_{\infty}}
 * @f]
 *
 * This parameterization is defined by
 *
 * @f[
 *  \log_{10} F(T, P_r) = \frac{\log_{10} F_{cent}(T)}{1 + f_1^2}
 * @f]
 * where
 * @f[
 *  F_{cent}(T) = (1 - A)\exp\left(\frac{-T}{T_3}\right)
 *      + A \exp\left(\frac{-T}{T_1}\right) + \exp\left(\frac{-T_2}{T}\right)
 * @f]
 *
 * @f[
 *  f_1 = \frac{\log_{10} P_r + C}{N - 0.14 (\log_{10} P_r + C)}
 * @f]
 *
 * @f[
 *  C = -0.4 - 0.67 \log_{10} F_{cent}
 * @f]
 *
 * @f[
 *  N = 0.75 - 1.27 \log_{10} F_{cent}
 * @f]
 *
 *  - If @f$ T_3 @f$ is zero, then the corresponding term is set to zero.
 *  - If @f$ T_1 @f$ is zero, then the corresponding term is set to zero.
 *  - If @f$ T_2 @f$ is zero, then the corresponding term is set to zero.
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

    TroeRate(const AnyMap& node, const UnitStack& rate_units={});
    TroeRate(const ArrheniusRate& low, const ArrheniusRate& high,
             const vector<double>& c);

    unique_ptr<MultiRateBase> newMultiRate() const override {
        return make_unique<MultiRate<TroeRate, FalloffData>>();
    }

    //! Set coefficients used by parameterization
    /*!
     * @param c Vector of three or four doubles: The doubles are the parameters,
     *          a, T_3, T_1, and (optionally) T_2 of the Troe parameterization
     */
    void setFalloffCoeffs(const vector<double>& c) override;

    void getFalloffCoeffs(vector<double>& c) const override;

    //! Update the temperature parameters in the representation
    /*!
     *   @param T         Temperature (Kelvin)
     *   @param work      Vector of working space, length 1, representing the
     *                    temperature-dependent part of the parameterization.
     */
    void updateTemp(double T, double* work) const override;

    double F(double pr, const double* work) const override;

    const string subType() const override {
        return "Troe";
    }

    size_t nParameters() const override {
        return 4;
    }

    void setParameters(const AnyMap& node, const UnitStack& rate_units) override;

    void getParameters(AnyMap& node) const override;

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
 * This falloff function is based on the one originally due to Stewart et al.
 * @cite stewart1989, which required three parameters @f$ a @f$, @f$ b @f$, and
 * @f$ c @f$. Kee et al. @cite kee1989 generalized this slightly by adding two more
 * parameters @f$ d @f$ and @f$ e @f$. (The original form corresponds to @f$ d = 1 @f$
 * and @f$ e = 0 @f$.) In keeping with the nomenclature of Kee et al. @cite kee1989,
 * the rate is referred to as the *SRI falloff function*.
 *
 * The falloff function defines the value of @f$ F @f$ in the following
 * rate expression
 * @f[
 *  k = k_{\infty} \left( \frac{P_r}{1 + P_r} \right) F
 * @f]
 * where
 * @f[
 *  P_r = \frac{k_0 [M]}{k_{\infty}}
 * @f]
 *
 * @f[
 *  F(T, P_r) = {\left[ a \; \exp\left(\frac{-b}{T}\right)
 *      + \exp\left(\frac{-T}{c}\right)\right]}^n \; d \; T^e
 * @f]
 * where
 * @f[
 *  n = \frac{1.0}{1.0 + (\log_{10} P_r)^2}
 * @f]
 *
 * @f$ c @f$ is required to be greater than or equal to zero. If it is zero, then the
 * corresponding term is set to zero. @f$ d @f$ is required to be greater than zero.
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

    SriRate(const AnyMap& node, const UnitStack& rate_units={});

    SriRate(const ArrheniusRate& low, const ArrheniusRate& high, const vector<double>& c)
        : SriRate()
    {
        m_lowRate = low;
        m_highRate = high;
        setFalloffCoeffs(c);
    }

    unique_ptr<MultiRateBase> newMultiRate() const override {
        return make_unique<MultiRate<SriRate, FalloffData>>();
    }

    //! Set coefficients used by parameterization
    /*!
     * @param c Vector of three or five doubles: The doubles are the parameters,
     *          a, b, c, d (optional; default 1.0), and e (optional; default
     *          0.0) of the SRI parameterization
     */
    void setFalloffCoeffs(const vector<double>& c) override;

    void getFalloffCoeffs(vector<double>& c) const override;

    //! Update the temperature parameters in the representation
    /*!
     *   @param T         Temperature (Kelvin)
     *   @param work      Vector of working space, length 2, representing the
     *                    temperature-dependent part of the parameterization.
     */
    void updateTemp(double T, double* work) const override;

    double F(double pr, const double* work) const override;

    const string subType() const override {
        return "SRI";
    }

    size_t nParameters() const override {
        return 5;
    }

    void setParameters(const AnyMap& node, const UnitStack& rate_units) override;
    void getParameters(AnyMap& node) const override;

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
 *  It provides a constant or linear in temperature value for @f$ F_{cent} @f$:
 *  @f[ F_{cent} = A + B*T @f]
 *
 *  The value of @f$ F_{cent} @f$ is then applied to Troe's model for the
 *  determination of the value of @f$ F(T, P_r) @f$:
 *  @f[ \log_{10} F(T, P_r) = \frac{\log_{10} F_{cent}(T)}{1 + f_1^2} @f]
 *  where
 *  @f[ f_1 = \frac{\log_{10} P_r + C}{N - 0.14 (\log_{10} P_r + C)} @f]
 *
 *  @f[ C = -0.4 - 0.67 \log_{10} F_{cent} @f]
 *
 *  @f[ N = 0.75 - 1.27 \log_{10} F_{cent} @f]
 *
 *  References:
 *  * Example of reaction database developed by Tsang utilizing this format
 *      @cite tsang1991
 *  * Example of Chemkin implementation of Tsang format (supplemental materials)
 *      @cite lucassen2011
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

    TsangRate(const AnyMap& node, const UnitStack& rate_units={});

    TsangRate(const ArrheniusRate& low, const ArrheniusRate& high, const vector<double>& c)
        : TsangRate()
    {
        m_lowRate = low;
        m_highRate = high;
        setFalloffCoeffs(c);
    }

    unique_ptr<MultiRateBase> newMultiRate() const override {
        return make_unique<MultiRate<TsangRate, FalloffData>>();
    }

    //! Set coefficients used by parameterization
    /*!
     * @param c Vector of one or two doubles: The doubles are the parameters,
     *          a and (optionally) b of the Tsang F_cent parameterization
     */
    void setFalloffCoeffs(const vector<double>& c) override;

    void getFalloffCoeffs(vector<double>& c) const override;

    //! Update the temperature parameters in the representation
    /*!
     *   @param T         Temperature (Kelvin)
     *   @param work      Vector of working space, length 1, representing the
     *                    temperature-dependent part of the parameterization.
     */
    void updateTemp(double T, double* work) const override;

    double F(double pr, const double* work) const override;

    const string subType() const override {
        return "Tsang";
    }

    size_t nParameters() const override {
        return 2;
    }

    void setParameters(const AnyMap& node, const UnitStack& rate_units) override;

    void getParameters(AnyMap& node) const override;

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
