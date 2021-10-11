// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_FALLOFF_H
#define CT_FALLOFF_H

#include "cantera/kinetics/Arrhenius.h"

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

/**
 * Base class for falloff function calculators. Each instance of a subclass of
 * Falloff computes one falloff function. This base class implements the
 * trivial falloff function F = 1.0.
 *
 * @todo remove Lindemann-specific portions to establish a 'true' base class
 */
class Falloff
{
public:
    Falloff()
        : allow_negative_pre_exponential_factor(false)
        , third_body_concentration(NAN)
        , rate_index(npos)
        , m_chemicallyActivated(false)
        , m_rc_low(NAN)
        , m_rc_high(NAN)
    {
        m_work.resize(workSize());
    }

    virtual ~Falloff() {}

    /**
     * Initialize. Must be called before any other method is invoked.
     *
     * @param c Vector of coefficients of the parameterization. The number and
     *     meaning of these coefficients is subclass-dependent.
     *
     * @todo  deprecate; superseded by setData
     */
    void init(const vector_fp& c);

    /**
     * Set coefficients of the falloff parametetrization.
     *
     * @param c Vector of coefficients of the parameterization. The number and
     *     meaning of these coefficients is subclass-dependent.
     */
    virtual void setData(const vector_fp& c);

    /**
     * Retrieve coefficients of the falloff parametetrization.
     *
     * @param c Vector of coefficients of the parameterization. The number and
     *     meaning of these coefficients is subclass-dependent.
     */
    virtual void getData(vector_fp& c) const;

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
        FalloffData data;
        data.update(T);
        update(data);
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

    //! Return a string representing the type of the Falloff parameterization.
    virtual const std::string type() const {
        return "Falloff";
    }

    //! Returns the number of parameters used by this parameterization. The
    //! values of these parameters can be obtained from getParameters().
    virtual size_t nParameters() const {
        return 0;
    }

    //! Perform object setup based on AnyMap node information
    //! @param node  AnyMap containing rate information
    //! @param rate_units  unit definitions specific to rate information
    virtual void setParameters(const AnyMap& node, const UnitsVector& rate_units);

    //! Get the values of the parameters for this object. *params* must be an
    //! array of at least nParameters() elements.
    /**
     * @todo  deprecate; superseded by getData
     */
    virtual void getParameters(double* params) const {}

    //! Store the falloff-related parameters needed to reconstruct an identical
    //! Reaction using the newReaction(AnyMap&, Kinetics&) function.
    virtual void getParameters(AnyMap& rateNode) const;

    //! Flag indicating that information specific to reaction rate is required
    const static bool usesUpdate() {
        return true;
    }

    //! Update information specific to reaction
    //! @param shared_data  data shared by all reactions of a given type
    virtual void update(const FalloffData& shared_data) {
        updateTemp(shared_data.temperature, m_work.data());
        m_rc_low = m_lowRate.eval(shared_data.logT, shared_data.recipT);
        m_rc_high = m_highRate.eval(shared_data.logT, shared_data.recipT);
        if (shared_data.finalized && rate_index != npos) {
            third_body_concentration = shared_data.conc_3b[rate_index];
        }
    }

    //! Evaluate reaction rate
    //! @param shared_data  data shared by all reactions of a given type
    virtual double eval(const FalloffData& shared_data) const {
        double pr = third_body_concentration * m_rc_low / (m_rc_high + SmallNumber);
        // AssertFinite(pr, "Falloff::eval", "pr is not finite.");

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

    //! Evaluate derivative of reaction rate with respect to temperature
    //! @param shared_data  data shared by all reactions of a given type
    virtual double ddT(const FalloffData& shared_data) const
    {
        return NAN; // @todo
    }

    //! Check the reaction rate expression
    void check(const std::string& equation, const AnyMap& node);

    //! Validate the reaction rate expression
    void validate(const std::string& equation) {}

    bool allow_negative_pre_exponential_factor; // Flag is directly accessible
    double third_body_concentration; //!< Buffered third-body concentration

    //! Get flag indicating whether reaction is chemically activated
    const bool chemicallyActivated() const {
        return m_chemicallyActivated;
    }

    //! Set flag indicating whether reaction is chemically activated
    void setChemicallyActivated(bool activated) {
        m_chemicallyActivated = activated;
    }

    //! Get reaction rate in the low-pressure limit
    ArrheniusBase& lowRate() {
        return m_lowRate;
    }

    //! Set reaction rate in the low-pressure limit
    void setLowRate(const ArrheniusBase& low);

    //! Get reaction rate in the high-pressure limit
    ArrheniusBase& highRate() {
        return m_highRate;
    }

    //! Set reaction rate in the high-pressure limit
    void setHighRate(const ArrheniusBase& high);

    size_t rate_index; //!< Reaction rate index within kinetics evaluator

protected:
    ArrheniusBase m_lowRate; //!< The reaction rate in the low-pressure limit
    ArrheniusBase m_highRate; //!< The reaction rate in the high-pressure limit

    bool m_chemicallyActivated; //!< Flag indicating whether reaction is chemically activated
    double m_rc_low; //!< Evaluated reaction rate in the low-pressure limit
    double m_rc_high; //!< Evaluated reaction rate in the high-pressure limit
    vector_fp m_work; //!< Work vector

    double m_conc_3b; //!< Third body concentration
};


//! The Lindemann falloff parameterization.
/**
 * This class implements the trivial falloff function F = 1.0.
 *
 * @ingroup falloffGroup
 */
class Lindemann : public Falloff
{
public:
    Lindemann() : Falloff() {
        m_work.resize(workSize());
    }

    Lindemann(const ArrheniusBase& low, const ArrheniusBase& high, const vector_fp& c)
        : Falloff()
    {
        m_lowRate = low;
        m_highRate = high;
        setData(c);
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
class Troe : public Falloff
{
public:
    //! Constructor
    Troe() : Falloff(), m_a(NAN), m_rt3(0.0), m_rt1(0.0), m_t2(0.0) {
        m_work.resize(workSize());
    }

    Troe(const ArrheniusBase& low, const ArrheniusBase& high, const vector_fp& c)
        : Troe()
    {
        m_lowRate = low;
        m_highRate = high;
        setData(c);
    }

    //! Set coefficients used by parameterization
    /*!
     * @param c Vector of three or four doubles: The doubles are the parameters,
     *          a, T_3, T_1, and (optionally) T_2 of the Troe parameterization
     */
    virtual void setData(const vector_fp& c);

    virtual void getData(vector_fp& c) const;

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

    virtual void setParameters(const AnyMap& node, const UnitsVector& rate_units);

    //! Sets params to contain, in order, \f[ (A, T_3, T_1, T_2) \f]
    /**
     * @todo  deprecate; superseded by getData
     */
    virtual void getParameters(double* params) const;

    virtual void getParameters(AnyMap& rateNode) const;

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
class SRI : public Falloff
{
public:
    //! Constructor
    SRI() : Falloff(), m_a(NAN), m_b(-1.0), m_c(-1.0), m_d(-1.0), m_e(-1.0) {
        m_work.resize(workSize());
    }

    SRI(const ArrheniusBase& low, const ArrheniusBase& high, const vector_fp& c)
        : SRI()
    {
        m_lowRate = low;
        m_highRate = high;
        setData(c);
    }

    //! Set coefficients used by parameterization
    /*!
     * @param c Vector of three or five doubles: The doubles are the parameters,
     *          a, b, c, d (optional; default 1.0), and e (optional; default
     *          0.0) of the SRI parameterization
     */
    virtual void setData(const vector_fp& c);

    virtual void getData(vector_fp& c) const;

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

    virtual void setParameters(const AnyMap& node, const UnitsVector& rate_units);

    //! Sets params to contain, in order, \f[ (a, b, c, d, e) \f]
    /**
     * @todo  deprecate; superseded by getData
     */
    virtual void getParameters(double* params) const;

    virtual void getParameters(AnyMap& rateNode) const;

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
class Tsang : public Falloff
{
public:
    //! Constructor
    Tsang() : Falloff(), m_a(NAN), m_b(0.0) {
        m_work.resize(workSize());
    }

    Tsang(const ArrheniusBase& low, const ArrheniusBase& high, const vector_fp& c)
        : Tsang()
    {
        m_lowRate = low;
        m_highRate = high;
        setData(c);
    }

    //! Set coefficients used by parameterization
    /*!
     * @param c Vector of one or two doubles: The doubles are the parameters,
     *          a and (optionally) b of the Tsang F_cent parameterization
     */
    virtual void setData(const vector_fp& c);

    virtual void getData(vector_fp& c) const;

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

    virtual void setParameters(const AnyMap& node, const UnitsVector& rate_units);

    //! Sets params to contain, in order, \f[ (A, B) \f]
    /**
     * @todo  deprecate; superseded by getData
     */
    virtual void getParameters(double* params) const;

    virtual void getParameters(AnyMap& rateNode) const;

protected:
    //! parameter a in the Tsang F_cent formulation. Dimensionless
    double m_a;

    //! parameter b in the Tsang F_cent formulation. [K^-1]
    double m_b;
};

}

#endif
