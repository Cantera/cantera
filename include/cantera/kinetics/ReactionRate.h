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
//! between Arrhenius rate expressions at various pressures; @see Plog
class PlogRate final : public RateEvaluator<Plog, PlogData>
{
public:
    PlogRate() = default;
    ~PlogRate() {}

    // inherit constructors
    using RateEvaluator<Plog, PlogData>::RateEvaluator;
};


//! Pressure-dependent rate expression where the rate coefficient is expressed
//! as a bivariate Chebyshev polynomial in temperature and pressure; @see Chebyshev
class ChebyshevRate3 final : public RateEvaluator<Chebyshev, ChebyshevData>
{
public:
    ChebyshevRate3() = default;
    ~ChebyshevRate3() {}

    // inherit constructors
    using RateEvaluator<Chebyshev, ChebyshevData>::RateEvaluator;
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
