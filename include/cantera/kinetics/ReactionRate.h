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
#include "cantera/kinetics/MultiRate.h"
#include "cantera/kinetics/RxnRates.h"
#include "cantera/kinetics/Falloff.h"
#include "cantera/base/ctexceptions.h"

namespace Cantera
{

class Func1;
class Kinetics;


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
    ~ReactionRate() {}

    //! Update information specific to reaction
    //! @param shared_data  data shared by all reactions of a given type
    virtual void update(const DataType& shared_data) {}

    virtual void update(double T) override {
        DataType data;
        data.update(T);
        update(data);
    }

    virtual void update(double T, double P) override {
        DataType data;
        data.update(T, P);
        update(data);
    }

    virtual void update(const ThermoPhase& bulk, const Kinetics& kin) override {
        DataType data;
        data.update(bulk, kin);
        update(data);
    }

    //! Evaluate reaction rate
    //! @param shared_data  data shared by all reactions of a given type
    virtual double eval(const DataType& shared_data) const = 0;

    virtual double eval(double T) const override {
        DataType data;
        data.update(T);
        return eval(data);
    }

    virtual double eval(double T, double P) const override {
        DataType data;
        data.update(T, P);
        return eval(data);
    }

    virtual double eval(const ThermoPhase& bulk, const Kinetics& kin) const override {
        DataType data;
        data.update(bulk, kin);
        return eval(data);
    }

    //! Evaluate derivative of reaction rate with respect to temperature
    //! @param shared_data  data shared by all reactions of a given type
    virtual double ddT(const DataType& shared_data) const {
        throw NotImplementedError("ReactionRate::ddT",
                                  "Not implemented by '{}' object.", type());
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

    virtual double ddT(const ThermoPhase& bulk, const Kinetics& kin) const override {
        DataType data;
        data.update(bulk, kin);
        return ddT(data);
    }
};


//! Class template for reaction rate evaluators with specialized RateType and DataType
/**
 * This class template ensures that derived objects are aware of specialized
 * data types, which are passed by MultiRateBase evaluators.
 */
template <class RateType, class DataType>
class RateEvaluator : public ReactionRateBase, public RateType
{
public:
    RateEvaluator() = default;
    ~RateEvaluator() {}

    //! Constructor using AnyMap content
    //! @param node  AnyMap containing rate information
    //! @param rate_units  unit definitions used for rate information
    RateEvaluator(const AnyMap& node, const UnitsVector& units={}) {
        setParameters(node, units);
    }

    using RateType::RateType;
    using RateType::usesUpdate;
    using RateType::update;
    using RateType::eval;

    // using RateType::type;
    virtual const std::string type() const override {
        // need to be explicitly included due to multiple inheritance
        return RateType::type();
    }

    virtual void setParameters(const AnyMap& node, const UnitsVector& units) override {
        ReactionRateBase::setParameters(node, units);
        RateType::setParameters(node, units);
    }

    // using RateType::getParameters;
    virtual void getParameters(AnyMap& node) const override {
        RateType::getParameters(node);
    }

    virtual unique_ptr<MultiRateBase> newMultiRate() const override {
        return unique_ptr<MultiRateBase>(
            new MultiBulkRate<RateEvaluator<RateType, DataType>, DataType>);
    }

    virtual void update(double T) override {
        DataType data;
        data.update(T);
        update(data);
    }

    virtual void update(double T, double P) override {
        DataType data;
        data.update(T, P);
        update(data);
    }

    virtual void update(const ThermoPhase& bulk, const Kinetics& kin) override {
        DataType data;
        data.update(bulk, kin);
        update(data);
    }

    virtual double eval(double T) const override {
        DataType data;
        data.update(T);
        return eval(data);
    }

    virtual double eval(double T, double P) const override {
        DataType data;
        data.update(T, P);
        return eval(data);
    }

    virtual double eval(const ThermoPhase& bulk, const Kinetics& kin) const override {
        DataType data;
        data.update(bulk, kin);
        return eval(data);
    }

    //! Evaluate derivative of reaction rate with respect to temperature
    //! @param shared_data  data shared by all reactions of a given type
    virtual double ddT(const DataType& shared_data) const {
        throw NotImplementedError("RateEvaluator::ddT",
                                  "Not implemented by '{}' object.", type());
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

    virtual double ddT(const ThermoPhase& bulk, const Kinetics& kin) const override {
        DataType data;
        data.update(bulk, kin);
        return ddT(data);
    }

    virtual void check(const std::string& equation, const AnyMap& node) override {
        RateType::check(equation, node);
    }

    virtual void validate(const std::string& equation) override {
        RateType::validate(equation);
    }

    virtual void setRateIndex(size_t idx) {
        ReactionRateBase::setRateIndex(idx);
        RateType::rate_index = idx;
    }
};


//! Arrhenius reaction rate type; @see Arrhenius3
class ArrheniusRate final : public RateEvaluator<Arrhenius3, ArrheniusData>
{
public:
    ArrheniusRate() = default;
    ~ArrheniusRate() {}

    // inherit constructors
    using RateEvaluator<Arrhenius3, ArrheniusData>::RateEvaluator;
};


//! Blowers-Masel reaction rate type; @see BlowersMasel3
class BlowersMaselRate final : public RateEvaluator<BlowersMasel3, BlowersMaselData>
{
public:
    BlowersMaselRate() = default;
    ~BlowersMaselRate() {}

    // inherit constructors
    using RateEvaluator<BlowersMasel3, BlowersMaselData>::RateEvaluator;
};


//! Class template for falloff reaction rate definitions; @see Falloff
template <class RateType=Lindemann>
class FalloffRate final : public RateEvaluator<RateType, FalloffData>
{
public:
    FalloffRate() = default;
    ~FalloffRate() {}

    // inherit constructors
    using RateEvaluator<RateType, FalloffData>::RateEvaluator;
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
    ~PlogRate() {}

    //! Constructor from Arrhenius rate expressions at a set of pressures
    explicit PlogRate(const std::multimap<double, Arrhenius>& rates);

    //! Constructor using AnyMap content
    //! @param node  AnyMap containing rate information
    //! @param units  Description of units used for rate parameters
    PlogRate(const AnyMap& node, const UnitsVector& units={});

    virtual const std::string type() const override {
        return "pressure-dependent-Arrhenius";
    }

    virtual unique_ptr<MultiRateBase> newMultiRate() const override;

    virtual void setParameters(const AnyMap& node, const UnitsVector& units) override;
    virtual void getParameters(AnyMap& rateNode) const override;

    //! Update information specific to reaction
    static bool usesUpdate() { return true; }

    virtual void update(const PlogData& shared_data) override {
        update_C(&shared_data.logP);
    }

    virtual double eval(const PlogData& shared_data) const override {
        return updateRC(shared_data.logT, shared_data.recipT);
    }

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
 *
 * @TODO  rename to ChebyshevRate when the legacy ChebyshevRate class is removed
 *      from RxnRates.h after Cantera 2.6.
 */
class ChebyshevRate3 final : public ReactionRate<ChebyshevData>, public Chebyshev
{
public:
    //! Default constructor.
    ChebyshevRate3() = default;
    ~ChebyshevRate3() {}

    //! Constructor using coefficient array
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
    //! @param units  Description of units used for rate parameters
    ChebyshevRate3(const AnyMap& node, const UnitsVector& rate_units={});

    virtual const std::string type() const override { return "Chebyshev"; }

    virtual unique_ptr<MultiRateBase> newMultiRate() const override;

    virtual void setParameters(const AnyMap& node, const UnitsVector& units) override;
    virtual void getParameters(AnyMap& rateNode) const override;

    //! Update information specific to reaction
    static bool usesUpdate() { return true; }

    virtual void update(const ChebyshevData& shared_data) override {
        update_C(&shared_data.log10P);
    }

    virtual double eval(const ChebyshevData& shared_data) const override {
        return updateRC(0., shared_data.recipT);
    }

protected:
    Units m_rate_units; //!< Reaction rate units
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
    ~CustomFunc1Rate() {}

    //! Constructor does nothing, as there is no formalized parameterization
    //! @param node  AnyMap object containing reaction rate specification
    //! @param units  Description of units used for rate parameters
    CustomFunc1Rate(const AnyMap& rate, const UnitsVector& units) {}

    virtual const std::string type() const override { return "custom-rate-function"; }

    virtual unique_ptr<MultiRateBase> newMultiRate() const override;

    virtual void setParameters(const AnyMap& node, const UnitsVector& units) override {
        ReactionRateBase::setParameters(node, units);
    }

    //! Update information specific to reaction
    static bool usesUpdate() { return false; }

    //! Set custom rate
    /**
     * The call to the Func1 object takes a single argument (temperature) and
     * does not depend on parameters handled in C++.
     */
    void setRateFunction(shared_ptr<Func1> f);

    virtual double eval(const CustomFunc1Data& shared_data) const override;

protected:
    shared_ptr<Func1> m_ratefunc;
};

}

#endif
