/**
 * @file ReactionRate.h
 *
 * @warning This file is an experimental part of the %Cantera API and
 *    may be changed or removed without notice.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_REACTIONRATE_H
#define CT_REACTIONRATE_H

#include "cantera/base/AnyMap.h"
#include "cantera/base/Units.h"
#include "cantera/kinetics/RxnRates.h"
#include "cantera/kinetics/ReactionData.h"
#include "cantera/base/ctexceptions.h"

namespace Cantera
{

class Func1;
class MultiRateBase;


//! Abstract base class for reaction rate definitions
/**
 * Because this class has no template parameters, derived objects can be
 * accessed via `shared_ptr<ReactionRateBase>`. For performance reasons
 * it is essential that derived classes use the keyword `final` to
 * de-virtualize `virtual` methods.
 *
 * Methods defined for the abstract base class are not aware of specialized
 * data handlers defined by the template class `ReactionRate<DataType>`
 * and thus can be exposed to the API.
 */
class ReactionRateBase
{
public:
    ReactionRateBase() : units(0.), m_index(npos) {}
    virtual ~ReactionRateBase() {}

    //! Identifier of reaction type
    virtual std::string type() const = 0;

    //! Update reaction rate data based on temperature
    //! @param T  temperature [K]
    virtual void update(double T) = 0;

    //! Update reaction rate data based on temperature and pressure
    //! @param T  temperature [K]
    //! @param P  pressure [Pa]
    //! @param concm  third-body concentration (if applicable)
    virtual void update(double T, double P, double concm=0.) = 0;

    //! Update reaction rate data based on bulk phase
    //! @param bulk  object representing bulk phase
    //! @param concm  third-body concentration (if applicable)
    virtual void update(const ThermoPhase& bulk, double concm=0.) = 0;

    //! Evaluate reaction rate based on temperature
    //! @param T  temperature [K]
    virtual double eval(double T) const = 0;

    //! Evaluate reaction rate based on temperature and pressure
    //! @param T  temperature [K]
    //! @param P  pressure [Pa]
    //! @param concm  third-body concentration (if applicable)
    virtual double eval(double T, double P, double concm=0.) const = 0;

    //! Evaluate reaction rate based on bulk phase
    //! @param bulk  object representing bulk phase
    //! @param concm  third-body concentration (if applicable)
    virtual double eval(const ThermoPhase& bulk, double concm=0.) const = 0;

    //! Evaluate reaction rate derivative based on temperature
    //! @param T  temperature [K]
    virtual double ddT(double T) const = 0;

    //! Evaluate reaction rate derivative based on temperature and pressure
    //! @param T  temperature [K]
    //! @param P  pressure [Pa]
    virtual double ddT(double T, double P) const = 0;

    //! Evaluate reaction rate derivative based on bulk phase
    //! @param bulk  object representing bulk phase
    //! @param concm  third-body concentration (if applicable)
    virtual double ddT(const ThermoPhase& bulk, double concm=0.) const = 0;

    //! Validate the reaction rate expression
    virtual void validate(const std::string& equation) = 0;

    //! Return the parameters such that an identical Reaction could be reconstructed
    //! using the newReaction() function. Behavior specific to derived classes is
    //! handled by the getParameters() method.
    //! @param rate_units  units used for rate parameters
    AnyMap parameters(const Units& rate_units) const;

    //! Return parameters using original unit system
    AnyMap parameters() const;

    //! Set parameters
    //! @param node  AnyMap object containing reaction rate specification
    //! @param rate_units  Description of units used for rate parameters
    virtual void setParameters(const AnyMap& node, const Units& rate_units);

    //! Indicate whether reaction is linked to a rate evaluator
    bool linked() { return bool(m_evaluator); }

    //! Link reaction to MultiRateBase evaluator
    void linkEvaluator(size_t index, const shared_ptr<MultiRateBase> evaluator);

    //! Release (unlink) reaction from MultiRateBase evaluator
    void releaseEvaluator();

    //! Return index of reaction within the Kinetics object owning the rate
    //! evaluator. Raises an exception if the reaction is not linked.
    size_t index();

protected:
    //! Get parameters
    //! Store the parameters of a ReactionRate needed to reconstruct an identical
    //! object. Does not include user-defined fields available in the #input map.
    virtual void getParameters(AnyMap& rateNode, const Units& rate_units) const {
        throw CanteraError("ReactionRate::getParameters",
                           "Not implemented by derived ReactionRate object.");
    }

    //! Input data used for specific models
    AnyMap input;

    //! The units of the rate constant. These are determined by the units of the
    //! standard concentration of the reactant species' phases of the phase
    //! where the reaction occurs.
    Units units;

    //! Evaluator handling the reaction rate
    std::shared_ptr<MultiRateBase> m_evaluator;

    //! Index of reaction within kinetics object (if applicable)
    size_t m_index;
};


//! Class template for reaction rate definitions with specialized DataType
/**
 * This class template ensures that derived objects are aware of specialized
 * data types, which are passed by MultiRateBase evaluators.
 */
template <class DataType>
class ReactionRate : public ReactionRateBase
{
public:
    ReactionRate() = default;

    //! Update information specific to reaction
    //! @param shared_data  data shared by all reactions of a given type
    virtual void update(const DataType& shared_data,
                        double concm=0.) {}

    virtual void update(double T) override {
        DataType data;
        data.update(T);
        update(data);
    }

    virtual void update(double T, double P, double concm=0.) override {
        DataType data;
        data.update(T, P);
        update(data, concm);
    }

    virtual void update(const ThermoPhase& bulk, double concm=0.) override {
        DataType data;
        data.update(bulk);
        update(data);
    }

    //! Evaluate reaction rate
    //! @param shared_data  data shared by all reactions of a given type
    //! @param concm  third-body concentration (if applicable)
    virtual double eval(const DataType& shared_data, double concm=0.) const = 0;

    virtual double eval(double T) const override {
        DataType data;
        data.update(T);
        return eval(data);
    }

    virtual double eval(double T, double P, double concm=0.) const override {
        DataType data;
        data.update(T, P);
        return eval(data, concm);
    }

    virtual double eval(const ThermoPhase& bulk, double concm=0.) const override {
        DataType data;
        data.update(bulk);
        return eval(data, concm);
    }

    //! Evaluate derivative of reaction rate with respect to temperature
    //! @param shared_data  data shared by all reactions of a given type
    virtual double ddT(const DataType& shared_data, double concm=0.) const {
        throw CanteraError("ReactionRate::ddT",
                           "Not implemented by derived ReactionRate object.");
    }

    virtual double ddT(double T) const override {
        DataType data;
        data.update(T);
        return ddT(data);
    }

    virtual double ddT(double T, double P) const override {
        DataType data;
        data.update(T, P);
        return ddT(data);
    }

    virtual double ddT(const ThermoPhase& bulk, double concm=0.) const override {
        DataType data;
        data.update(bulk);
        return ddT(data, concm);
    }
};


//! The ArrheniusRate reaction rate type depends only on temperature
/**
 *  A reaction rate coefficient of the following form.
 *
 *   \f[
 *        k_f =  A T^b \exp (-E/RT)
 *   \f]
 *
 * Note: ArrheniusRate acts as a wrapper for the Arrhenius class.
 */
class ArrheniusRate final : public ReactionRate<ArrheniusData>, public Arrhenius
{
public:
    //! Default constructor.
    ArrheniusRate();

    //! Constructor.
    //! @param A  pre-exponential. The unit system is
    //!     (kmol, m, s). The actual units depend on the reaction
    //!     order and the dimensionality (surface or bulk).
    //! @param b  Temperature exponent. Non-dimensional.
    //! @param E  Activation energy. J/kmol.
    ArrheniusRate(double A, double b, double E);

    //! Constructor using AnyMap content
    //! @param node  AnyMap containing rate information
    //! @param rate_units  unit definitions used for rate information
    ArrheniusRate(const AnyMap& node, const Units& rate_units);

    //! Constructor using AnyMap content
    //! @param node  AnyMap containing rate information
    ArrheniusRate(const AnyMap& node);

    //! Constructor based on Arrhenius object
    //! @param arr  Arrhenius object
    //! @param allow_negative_A  allow negative pre-exponential factor
    //!      (optional, default is false)
    ArrheniusRate(const Arrhenius& arr, bool allow_negative_A=false);

    virtual void setParameters(const AnyMap& node, const Units& rate_units) override;
    virtual void getParameters(AnyMap& rateNode,
                               const Units& rate_units) const override;

    virtual std::string type() const override { return "ArrheniusRate"; }

    //! Update information specific to reaction
    static bool usesUpdate() { return false; }

    virtual double eval(const ArrheniusData& shared_data,
                        double concm=0.) const override {
        return updateRC(shared_data.m_logT, shared_data.m_recipT);
    }

    virtual double ddT(const ArrheniusData& shared_data,
                       double concm=0.) const override {
        return updateRC(shared_data.m_logT, shared_data.m_recipT) *
            (m_b + m_E * shared_data.m_recipT) * shared_data.m_recipT;
    }

    //! Set the pre-exponential factor *A* (in m, kmol, s to powers depending
    //! on the reaction order)
    void setPreExponentialFactor(double A);

    //! Set the temperature exponent *b*
    void setTemperatureExponent(double b);

    //! Set the activation energy [J/kmol]
    void setActivationEnergy(double E);

    //! Return the activation energy [J/kmol]
    double activationEnergy() const {
        return m_E * GasConstant;
    }

    virtual void validate(const std::string& equation) override;

    bool allow_negative_pre_exponential_factor;
};


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
class PlogRate final : public ReactionRate<PlogData>, public Plog
{
public:
    PlogRate() = default;

    //! Constructor from Arrhenius rate expressions at a set of pressures
    explicit PlogRate(const std::vector<std::pair<double, Arrhenius>>& rates);

    //! Constructor using AnyMap content
    //! @param node  AnyMap containing rate information
    //! @param rate_units  unit definitions used for rate information
    PlogRate(const AnyMap& node, const Units& rate_units);

    //! Constructor using AnyMap content
    //! @param node  AnyMap containing rate information
    PlogRate(const AnyMap& node);

    virtual std::string type() const override { return "PlogRate"; }

    virtual void setParameters(const AnyMap& node, const Units& rate_units) override;
    virtual void getParameters(AnyMap& rateNode,
                               const Units& rate_units) const override;

    //! Update information specific to reaction
    static bool usesUpdate() { return true; }

    virtual void update(const PlogData& shared_data, double concm=0.) override {
        update_C(&shared_data.m_logP);
    }

    virtual double eval(const PlogData& shared_data,
                        double concm=0.) const override {
        return updateRC(shared_data.m_logT, shared_data.m_recipT);
    }

    //! Set pressures and Arrhenius expressions which comprise this reaction.
    void setRates(const std::vector<std::pair<double, Arrhenius>>& rates);

    virtual void validate(const std::string& equation) override {
        Plog::validate(equation);
    }
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
 * A Chebyshev rate expression is specified in terms of the coefficient matrix
 * \f$ \alpha \f$ and the temperature and pressure ranges. Note that the
 * Chebyshev polynomials are not defined outside the interval (-1,1), and
 * therefore extrapolation of rates outside the range of temperatures and
 * pressures for which they are defined is strongly discouraged.
 */
class ChebyshevRate3 final : public ReactionRate<ChebyshevData>, public Chebyshev
{
public:
    //! Default constructor.
    ChebyshevRate3() {}

    //! Constructor directly from coefficient array
    /*
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

    //! Constructor using AnyMap content
    //! @param node  AnyMap containing rate information
    //! @param rate_units  unit definitions used for rate information
    ChebyshevRate3(const AnyMap& node, const Units& rate_units);

    //! Constructor using AnyMap content
    //! @param node  AnyMap containing rate information
    ChebyshevRate3(const AnyMap& node);

    virtual std::string type() const override { return "ChebyshevRate"; }

    virtual void setParameters(const AnyMap& node, const Units& rate_units) override;
    virtual void getParameters(AnyMap& rateNode,
                               const Units& rate_units) const override;

    //! Update information specific to reaction
    static bool usesUpdate() { return true; }

    virtual void update(const ChebyshevData& shared_data, double concm=0.) override {
        update_C(&shared_data.m_log10P);
    }

    virtual double eval(const ChebyshevData& shared_data,
                        double concm=0.) const override {
        return updateRC(0., shared_data.m_recipT);
    }

    //! Access the Chebyshev coefficients as 2-dimensional array.
    const Array2D& coeffs() const;

    //! Set the Chebyshev coefficients as 2-dimensional array.
    void setCoeffs(const Array2D& coeffs);

    virtual void validate(const std::string& equation) override;
};


//! Custom reaction rate depending only on temperature
/**
 * The rate expression is provided by a `Func1` object taking a single
 * argument (temperature) and does not use a formalized parameterization.
 *
 * @warning This class is an experimental part of the %Cantera API and
 *    may be changed or removed without notice.
 */
class CustomFunc1Rate final : public ReactionRate<CustomFunc1Data>
{
public:
    CustomFunc1Rate();

    //! Constructor does nothing, as there is no formalized parameterization
    //! @param node  AnyMap object containing reaction rate specification
    //! @param rate_units  Description of units used for rate parameters
    CustomFunc1Rate(const AnyMap& rate, const Units& rate_units) {}

    virtual std::string type() const override { return "custom-function"; }

    virtual void setParameters(const AnyMap& node, const Units& rate_units) override {
        units = rate_units;
    }

    //! Update information specific to reaction
    static bool usesUpdate() { return false; }

    //! Set custom rate
    /**
     * The call to the Func1 object takes a single argument (temperature) and
     * does not depend on parameters handled in C++.
     */
    void setRateFunction(shared_ptr<Func1> f);

    virtual double eval(const CustomFunc1Data& shared_data,
                        double concm=0.) const override;

    virtual void validate(const std::string& equation) override {}

protected:
    shared_ptr<Func1> m_ratefunc;
};

}

#endif
