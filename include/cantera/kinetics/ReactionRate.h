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
    virtual ~ReactionRateBase() {}

    //! Set parameters
    //! @param node  AnyMap object containing reaction rate specification
    //! @param rate_units  Description of units used for rate parameters
    virtual bool setParameters(const AnyMap& node, const Units& rate_units) {
        throw CanteraError("ReactionRate::setParameters",
                           "Not implemented by derived ReactionRate object.");
    }

    //! Get parameters
    //! Store the parameters of a ReactionRate needed to reconstruct an identical
    //! object. Does not include user-defined fields available in the #input map.
    virtual void getParameters(AnyMap& rateNode, const Units& rate_units) const {
        throw CanteraError("ReactionRate::getParameters",
                           "Not implemented by derived ReactionRate object.");
    }

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
        update(DataType(T));
    }

    virtual void update(double T, double P, double concm=0.) override {
        update(DataType(T, P), concm);
    }

    virtual void update(const ThermoPhase& bulk, double concm=0.) override {
        update(DataType(bulk), concm);
    }


    //! Evaluate reaction rate
    //! @param shared_data  data shared by all reactions of a given type
    //! @param concm  third-body concentration (if applicable)
    virtual double eval(const DataType& shared_data, double concm=0.) const = 0;

    virtual double eval(double T) const override {
        return eval(DataType(T));
    }

    virtual double eval(double T, double P, double concm=0.) const override {
        return eval(DataType(T, P), concm);
    }

    virtual double eval(const ThermoPhase& bulk, double concm=0.) const override {
        return eval(DataType(bulk), concm);
    }

    //! Evaluate derivative of reaction rate with respect to temperature
    //! @param shared_data  data shared by all reactions of a given type
    virtual double ddT(const DataType& shared_data, double concm=0.) const {
        throw CanteraError("ReactionRate::ddT",
                           "Not implemented by derived ReactionRate object.");
    }

    virtual double ddT(double T) const override {
        return ddT(DataType(T));
    }

    virtual double ddT(double T, double P) const override {
        return ddT(DataType(T, P));
    }

    virtual double ddT(const ThermoPhase& bulk, double concm=0.) const override {
        return ddT(DataType(bulk), concm);
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

    /// Constructor.
    /// @param A  pre-exponential. The unit system is
    ///     (kmol, m, s). The actual units depend on the reaction
    ///     order and the dimensionality (surface or bulk).
    /// @param b  Temperature exponent. Non-dimensional.
    /// @param E  Activation energy. J/kmol.
    ArrheniusRate(double A, double b, double E);

    //! Constructor using AnyMap content
    /// @param node  AnyMap containing rate information
    /// @param rate_units  unit definitions used for rate information
    ArrheniusRate(const AnyMap& node, const Units& rate_units);

    virtual bool setParameters(const AnyMap& node, const Units& rate_units) override;
    virtual void getParameters(AnyMap& rateNode,
                               const Units& rate_units) const override;

    virtual std::string type() const override { return "ArrheniusRate"; }

    //! Update information specific to reaction
    static bool uses_update() { return false; }

    virtual double eval(const ArrheniusData& shared_data,
                        double concm=0.) const override {
        return updateRC(shared_data.m_logT, shared_data.m_recipT);
    }

    virtual double ddT(const ArrheniusData& shared_data,
                       double concm=0.) const override {
        return updateRC(shared_data.m_logT, shared_data.m_recipT) *
            (m_b + m_E * shared_data.m_recipT) * shared_data.m_recipT;
    }

    //! Return the activation energy [J/kmol]
    double activationEnergy() const {
        return m_E * GasConstant;
    }
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
    PlogRate();

    //! Constructor from Arrhenius rate expressions at a set of pressures
    explicit PlogRate(const std::multimap<double, Arrhenius>& rates);

    //! Constructor using AnyMap content
    /// @param node  AnyMap containing rate information
    /// @param rate_units  unit definitions used for rate information
    PlogRate(const AnyMap& node, const Units& rate_units);

    virtual bool setParameters(const AnyMap& node, const Units& rate_units) override;
    virtual void getParameters(AnyMap& rateNode,
                               const Units& rate_units) const override;

    //! Update information specific to reaction
    static bool uses_update() { return true; }

    virtual void update(const PlogData& shared_data, double concm=0.) override {
        update_C(shared_data.logP());
    }

    virtual std::string type() const override { return "PlogRate"; }

    virtual double eval(const PlogData& shared_data,
                        double concm=0.) const override {
        return updateRC(shared_data.m_logT, shared_data.m_recipT);
    }
};


//! Custom reaction rate depending only on temperature
/**
 * The rate expression is provided by a `Func1` object taking a single
 * argument (temperature) and does not use a formalized parameterization.
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

    virtual double eval(const CustomFunc1Data& shared_data,
                        double concm=0.) const override;

protected:
    shared_ptr<Func1> m_ratefunc;
};

}

#endif
