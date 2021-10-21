// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_CUSTOM_H
#define CT_CUSTOM_H

#include "cantera/base/ct_defs.h"
#include "cantera/base/Units.h"
#include "cantera/kinetics/ReactionData.h"

namespace Cantera
{

class AnyValue;
class AnyMap;
class Func1;


//! Custom reaction rate depending only on temperature
/**
 * The rate expression is provided by a `Func1` object taking a single
 * argument (temperature) and does not use a formalized parameterization.
 *
 * @warning This class is an experimental part of the %Cantera API and
 *    may be changed or removed without notice.
 */
class CustomFunc1
{
public:
    CustomFunc1();
    ~CustomFunc1() {}

    const std::string type() const { return "custom-rate-function"; }

    void setParameters(const AnyMap& node, const UnitsVector& units) {}

    void getParameters(AnyMap& rateNode, const Units& rate_units=Units(0.)) const;

    //! Update information specific to reaction
    static bool usesUpdate() { return false; }

    //! Update information specific to reaction
    /*!
     *  @param shared_data  data shared by all reactions of a given type
     */
    void update(const CustomFunc1Data& shared_data) {}

    //! Update information specific to reaction
    /*!
     *  @param shared_data  data shared by all reactions of a given type
     */
    double eval(const CustomFunc1Data& shared_data) const;

    //! Check the reaction rate expression
    void check(const std::string& equation, const AnyMap& node) {}

    void validate(const std::string& equation) {}

    //! Set custom rate
    /**
     * The call to the Func1 object takes a single argument (temperature) and
     * does not depend on parameters handled in C++.
     */
    void setRateFunction(shared_ptr<Func1> f);

    size_t rate_index; //!< Reaction rate index within kinetics evaluator

protected:
    shared_ptr<Func1> m_ratefunc;
};

}

#endif
