/**
 * @file ReactionRate.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_REACTIONRATE_H
#define CT_REACTIONRATE_H

#include "MultiRateBase.h"
#include "cantera/base/AnyMap.h"
#include "cantera/base/Units.h"
#include "cantera/base/ctexceptions.h"

namespace Cantera
{

//! Abstract base class for reaction rate definitions; this base class is used by
//! user-facing APIs to access reaction rate objects
//!
//! In addition to the pure virtual methods declared in this class, complete derived
//! classes must implement the method `evalFromStruct(const DataType& shared_data)`,
//! where `DataType` is a container for parameters needed to evaluate reactions of that
//! type.
//!
//! Derived classes may also implement the method
//! `ddTScaledFromStruct(const DataType& shared_data)`.
class ReactionRate
{
public:
    ReactionRate() : m_rate_index(npos) {}

    // Copy constructor and assignment operator need to be defined because of the
    // #m_evaluator member that can't (and shouldn't) be copied.
    ReactionRate(const ReactionRate& other)
        : m_input(other.m_input)
        , m_rate_index(other.m_rate_index)
    {}

    ReactionRate& operator=(const ReactionRate& other) {
        if (this == &other) {
            return *this;
        }
        m_input = other.m_input;
        m_rate_index = other.m_rate_index;
        return *this;
    }

    virtual ~ReactionRate() = default;

    //! Create a rate evaluator for reactions of a particular derived type.
    //! Derived classes usually implement this as:
    //!
    //! ```.cpp
    //! unique_ptr<MultiRateBase> newMultiRate() const override {
    //!     return unique_ptr<MultiRateBase>(new MultiBulkRate<RateType, DataType>);
    //! ```
    //!
    //! where `RateType` is the derived class name and `DataType` is the corresponding
    //! container for parameters needed to evaluate reactions of that type.
    virtual unique_ptr<MultiRateBase> newMultiRate() const = 0;

    virtual const std::string type() const = 0;

    //! Set parameters
    //! @param node  AnyMap object containing reaction rate specification
    //! @param units  unit definitions specific to rate information
    virtual void setParameters(const AnyMap& node, const UnitStack& units) {
        m_input = node;
    }

    //! Return the parameters such that an identical Reaction could be reconstructed
    //! using the newReaction() function. Behavior specific to derived classes is
    //! handled by the getParameters() method.
    AnyMap parameters() const {
        AnyMap out;
        getParameters(out);
        return out;
    }

    //! Check basic syntax and settings of reaction rate expression
    virtual void check(const std::string& equation, const AnyMap& node) {}

    //! Validate the reaction rate expression
    virtual void validate(const std::string& equation, const Kinetics& kin) {}

    //! Validate the reaction rate expression (legacy call)
    //! @todo deprecate in favor of two-parameter version
    virtual void validate(const std::string& equation) {}

    //! Reaction rate index within kinetics evaluator
    size_t rateIndex() const {
        return m_rate_index;
    }

    //! Set reaction rate index within kinetics evaluator
    void setRateIndex(size_t idx) {
        m_rate_index = idx;
    }

    //! Evaluate reaction rate based on temperature
    //! @param T  temperature [K]
    //! Used in conjunction with MultiRateBase::evalSingle / ReactionRate::eval.
    //! This method allows for testing of a reaction rate expression outside of
    //! Kinetics reaction rate evaluators.
    double eval(double T) {
        _evaluator().update(T);
        return _evaluator().evalSingle(*this);
    }

    //! Evaluate reaction rate based on temperature and an extra parameter.
    //! Specific rate parameterizations may require an additional parameter, which
    //! is specific to the derived ReactionRate object.
    //! @param T  temperature [K]
    //! @param extra  extra parameter used by parameterization
    //! Used in conjunction with MultiRateBase::evalSingle / ReactionRate::eval.
    //! This method allows for testing of a reaction rate expression outside of
    //! Kinetics reaction rate evaluators.
    double eval(double T, double extra) {
        _evaluator().update(T, extra);
        return _evaluator().evalSingle(*this);
    }

protected:
    //! Get parameters
    //! @param node  AnyMap containing rate information
    //! Store the parameters of a ReactionRate needed to reconstruct an identical
    //! object. Does not include user-defined fields available in the #m_input map.
    virtual void getParameters(AnyMap& node) const {
        throw NotImplementedError("ReactionRate::getParameters",
                                  "Not implemented by '{}' object.", type());
    }

    //! Input data used for specific models
    AnyMap m_input;

    //! Index of reaction rate within kinetics evaluator
    size_t m_rate_index;

private:
    //! Return an object that be used to evaluate the rate by converting general input
    //! such as temperature and pressure into the `DataType` struct that is particular
    //! to the rate model.
    MultiRateBase& _evaluator() {
        if (!m_evaluator) {
            m_evaluator = newMultiRate();
        }
        return *m_evaluator;
    }

    unique_ptr<MultiRateBase> m_evaluator;
};

}

#endif
