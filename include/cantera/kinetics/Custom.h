/**
 * @file Custom.h
 *
 * @warning This file is an experimental part of the %Cantera API and
 *    may be changed or removed without notice.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_CUSTOM_H
#define CT_CUSTOM_H

#include "cantera/base/ct_defs.h"
#include "cantera/base/Units.h"
#include "cantera/kinetics/Arrhenius.h"
#include "ReactionRate.h"
#include "MultiRate.h"

namespace Cantera
{

class AnyValue;
class AnyMap;
class Func1;


//! Custom reaction rate depending only on temperature
/**
 * The rate expression is provided by a Func1 object taking a single
 * argument (temperature) and does not use a formalized parameterization.
 *
 * @warning This class is an experimental part of the %Cantera API and
 *    may be changed or removed without notice.
 * @ingroup otherRateGroup
 */
class CustomFunc1Rate final : public ReactionRate
{
public:
    CustomFunc1Rate() = default;
    CustomFunc1Rate(const AnyMap& node, const UnitStack& rate_units);

    unique_ptr<MultiRateBase> newMultiRate() const override {
        return make_unique<MultiRate<CustomFunc1Rate, ArrheniusData>>();
    }

    const string type() const override { return "custom-rate-function"; }

    void getParameters(AnyMap& rateNode, const Units& rate_units=Units(0.)) const;
    using ReactionRate::getParameters;

    void validate(const string& equation, const Kinetics& kin) override;

    //! Update information specific to reaction
    /*!
     *  @param shared_data  data shared by all reactions of a given type
     */
    double evalFromStruct(const ArrheniusData& shared_data) const;

    //! Set custom rate
    /**
     * The call to the Func1 object takes a single argument (temperature) and
     * does not depend on parameters handled in C++.
     */
    void setRateFunction(shared_ptr<Func1> f);

protected:
    shared_ptr<Func1> m_ratefunc;
};


}

#endif
