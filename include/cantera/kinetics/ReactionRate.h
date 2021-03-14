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

#include "cantera/kinetics/RxnRates.h"
#include "cantera/base/ctexceptions.h"

namespace Cantera
{

class Func1;
class ThermoPhase;
class AnyMap;


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
 *
 * @warning This class is an experimental part of the %Cantera API and
 *    may be changed or removed without notice.
 */
class ReactionRateBase
{
public:
    virtual ~ReactionRateBase() {}

    //! Identifier of reaction type
    virtual std::string type() const = 0;

    //! Evaluate reaction rate based on temperature
    virtual double eval(double T) const = 0;

    //! Evaluate reaction rate based on temperature and pressure
    virtual double eval(double T, double P) const = 0;

    //! Evaluate reaction rate based on bulk phase
    virtual double eval(const ThermoPhase& bulk) const = 0;

    //! Evaluate reaction rate derivative based on temperature
    virtual double ddT(double T) const = 0;

    //! Evaluate reaction rate derivative based on temperature and pressure
    virtual double ddT(double T, double P) const = 0;

    //! Evaluate reaction rate derivative based on bulk phase
    virtual double ddT(const ThermoPhase& bulk) const = 0;
};


//! Class template for reaction rate definitions with specialized DataType
/**
 * This class template ensures that derived objects are aware of specialized
 * data types, which are passed by `MultiRateBase` evaluators.
 *
 * @warning This class is an experimental part of the %Cantera API and
 *    may be changed or removed without notice.
 */
template <class DataType>
class ReactionRate : public ReactionRateBase
{
public:
    //! Constructor
    ReactionRate() = default;

    //! Update information specific to reaction
    virtual void update(const DataType& shared_data, const ThermoPhase& bulk) {}

    //! Evaluate reaction rate
    virtual double eval(const DataType& shared_data) const = 0;

    virtual double eval(double T) const override {
        return eval(DataType(T));
    }

    virtual double eval(double T, double P) const override {
        return eval(DataType(T, P));
    }

    virtual double eval(const ThermoPhase& bulk) const override {
        return eval(DataType(bulk));
    }

    //! Evaluate derivative of reaction rate with respect to temperature
    virtual double ddT(const DataType& shared_data) const {
        throw CanteraError("ReactionRate::ddT",
                           "Not implemented by derived ReactionRate object.");
    }

    virtual double ddT(double T) const override {
        return ddT(DataType(T));
    }

    virtual double ddT(double T, double P) const override {
        return ddT(DataType(T, P));
    }

    virtual double ddT(const ThermoPhase& bulk) const {
        return ddT(DataType(bulk));
    }

    // [...] signatures for derivative with respect to pressure and concentration
    // are not created (yet)
};


//! Data container holding shared data specific to `ArrheniusRate`
/**
 * The data container `ArrheniusData` holds precalculated data common to
 * all `ArrheniusRate` objects.
 *
 * @warning This class is an experimental part of the %Cantera API and
 *    may be changed or removed without notice.
 */
struct ArrheniusData {
    ArrheniusData() : m_temperature(1.), m_logT(0.), m_recipT(1.) {}

    //! Constructor based on temperature *T*
    ArrheniusData(double T) : m_temperature(T) {
        m_logT = std::log(T);
        m_recipT = 1./T;
    }

    //! Constructor based on temperature *T*
    ArrheniusData(double T, double P) : ArrheniusData(T) {}

    //! Constructor accessing bulk phase definitions
    ArrheniusData(const ThermoPhase& bulk) { update(bulk); }

    //! Update based on bulk phase definitions
    void update(const ThermoPhase& bulk);

    double m_temperature;
    double m_logT;
    double m_recipT;
};


//! Arrhenius reaction rate type depends only on temperature
/**
 * Wrapped Arrhenius rate.
 *
 * @warning This class is an experimental part of the %Cantera API and
 *    may be changed or removed without notice.
 */
class ArrheniusRate final : public ReactionRate<ArrheniusData>, public Arrhenius
{
public:
    //! Constructor.
    ArrheniusRate();

    ArrheniusRate(double A, double b, double E);

    ArrheniusRate(const AnyMap& node, const Units& rc_units);

    //! @internal access to reaction type is used by `MultiRateBase` evaluator
    static std::string staticType() { return "ArrheniusRate"; }

    virtual std::string type() const override { return staticType(); }

    //! Update information specific to reaction
    static bool uses_update() { return false; }

    virtual double eval(const ArrheniusData& shared_data) const override {
        return updateRC(shared_data.m_logT, shared_data.m_recipT);
    }

    virtual double ddT(const ArrheniusData& shared_data) const override {
        return updateRC(shared_data.m_logT, shared_data.m_recipT) *
            (m_b + m_E * shared_data.m_recipT) * shared_data.m_recipT;
    }
};


//! Data container holding shared data specific to `CustomFunc1Rate`
/**
 * @warning This class is an experimental part of the %Cantera API and
 *    may be changed or removed without notice.
 */
struct CustomFunc1Data {
    CustomFunc1Data() : m_temperature(1.) {}

    //! Constructor based on temperature *T*
    CustomFunc1Data(double T) : m_temperature(T) {}

    //! Constructor based on temperature *T*
    CustomFunc1Data(double T, double P) : CustomFunc1Data(T) {}

    //! Constructor accessing bulk phase definitions
    CustomFunc1Data(const ThermoPhase& bulk) { update(bulk); }

    //! Update based on reacting phase definitions
    void update(const ThermoPhase& bulk);

    double m_temperature;
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
    //! Constructor.
    CustomFunc1Rate();

    // Does nothing, as there is no formalized parameterization.
    CustomFunc1Rate(const AnyMap& rate, const Units& rc_units) {}

    //! @internal access to reaction type is used by `MultiRateBase` evaluator
    static std::string staticType() { return "custom-function"; }

    virtual std::string type() const override { return staticType(); }

    //! Update information specific to reaction
    static bool uses_update() { return false; }

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
