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
#include "cantera/kinetics/ReactionData.h"
#include "cantera/base/ctexceptions.h"

namespace Cantera
{

class Func1;
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
    //! @param T  temperature [K]
    virtual double eval(double T) const = 0;

    //! Evaluate reaction rate based on temperature and pressure
    //! @param T  temperature [K]
    //! @param P  pressure [Pa]
    virtual double eval(double T, double P) const = 0;

    //! Evaluate reaction rate based on bulk phase
    //! @param bulk  object representing bulk phase
    virtual double eval(const ThermoPhase& bulk) const = 0;

    //! Evaluate reaction rate derivative based on temperature
    //! @param T  temperature [K]
    virtual double ddT(double T) const = 0;

    //! Evaluate reaction rate derivative based on temperature and pressure
    //! @param T  temperature [K]
    //! @param P  pressure [Pa]
    virtual double ddT(double T, double P) const = 0;

    //! Evaluate reaction rate derivative based on bulk phase
    //! @param bulk  object representing bulk phase
    virtual double ddT(const ThermoPhase& bulk) const = 0;
};


//! Class template for reaction rate definitions with specialized DataType
/**
 * This class template ensures that derived objects are aware of specialized
 * data types, which are passed by MultiRateBase evaluators.
 *
 * @warning This class is an experimental part of the %Cantera API and
 *    may be changed or removed without notice.
 */
template <class DataType>
class ReactionRate : public ReactionRateBase
{
public:
    ReactionRate() = default;

    //! Update information specific to reaction
    //! @param shared_data  data shared by all reactions of a given type
    virtual void update(const DataType& shared_data, const ThermoPhase& bulk) {}

    //! Evaluate reaction rate
    //! @param shared_data  data shared by all reactions of a given type
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
    //! @param shared_data  data shared by all reactions of a given type
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

    virtual double ddT(const ThermoPhase& bulk) const override {
        return ddT(DataType(bulk));
    }
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
    ArrheniusRate();

    ArrheniusRate(double A, double b, double E);

    //! Constructor
    //! @param node  AnyMap object containing reaction rate specification
    //! @param rate_units  Description of units used for rate parameters
    ArrheniusRate(const AnyMap& node, const Units& rate_units);

    virtual std::string type() const override { return "ArrheniusRate"; }

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
