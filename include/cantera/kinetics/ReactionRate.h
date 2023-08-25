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
#include "cantera/base/global.h"
#include "cantera/base/ctexceptions.h"

namespace Cantera
{


class Reaction;

//! Abstract base class for reaction rate definitions; this base class is used by
//! user-facing APIs to access reaction rate objects
//!
//! In addition to the pure virtual methods declared in this class, complete derived
//! classes must implement the method `evalFromStruct(const DataType& shared_data)`,
//! where `DataType` is a container for parameters needed to evaluate reactions of that
//! type. In addition, derived classes may also implement the method
//! `updateFromStruct(const DataType& shared_data)` to update buffered data that
//! is specific to a given reaction rate.
//!
//! The calculation of derivatives (or Jacobians) relies on the following methods:
//!  -  Derived classes may implement the method
//!     `ddTScaledFromStruct(const DataType& shared_data)` for an analytical derivative
//!     with respect to temperature.
//!  -  Associated `DataType` containers may overload the method
//!     `perturbTemperature(double deltaT)`, which is used for the calculation of
//!     numerical derivatives with respect to temperature if an analytic implementation
//!     is not available.
//!  -  For reaction rate constants that depend on pressure or third-body collision
//!     partners, associated `DataType` containers should implement the methods
//!     `perturbPressure(double deltaP)` and/or `perturbThirdBodies(double deltaM)`,
//!     which allow for the calculation of numerical derivatives.
//!  -  For additional information, refer to the @ref kinDerivs "Kinetics Derivatives"
//!     documentation.
//! @ingroup reactionGroup
class ReactionRate
{
public:
    ReactionRate() {}

    // Copy constructor and assignment operator need to be defined because of the
    // #m_evaluator member that can't (and shouldn't) be copied.
    ReactionRate(const ReactionRate& other)
        : m_input(other.m_input)
        , m_rate_index(other.m_rate_index)
        , m_valid(other.m_valid)
        , m_conversion_units(other.m_conversion_units)
    {}

    ReactionRate& operator=(const ReactionRate& other) {
        if (this == &other) {
            return *this;
        }
        m_input = other.m_input;
        m_rate_index = other.m_rate_index;
        m_valid = other.m_valid;
        m_conversion_units = other.m_conversion_units;
        return *this;
    }

    virtual ~ReactionRate() = default;

    //! Create a rate evaluator for reactions of a particular derived type.
    //! Derived classes usually implement this as:
    //!
    //! ```.cpp
    //! unique_ptr<MultiRateBase> newMultiRate() const override {
    //!     return make_unique<MultiRate<RateType, DataType>>();
    //! ```
    //!
    //! where `RateType` is the derived class name and `DataType` is the corresponding
    //! container for parameters needed to evaluate reactions of that type.
    virtual unique_ptr<MultiRateBase> newMultiRate() const {
        throw NotImplementedError("ReactionRate::newMultiRate",
            "Not implemented by '{}' object.", type());
    }

    //! String identifying reaction rate specialization
    virtual const string type() const = 0;

    //! String identifying sub-type of reaction rate specialization
    virtual const string subType() const {
        return "";
    }

    //! Set parameters
    //! @param node  AnyMap object containing reaction rate specification
    //! @param units  unit definitions specific to rate information
    virtual void setParameters(const AnyMap& node, const UnitStack& units) {
        setRateUnits(units);
        m_input = node;
    }

    //! Return the parameters such that an identical Reaction could be reconstructed
    //! using the newReaction() function. Behavior specific to derived classes is
    //! handled by the getParameters() method.
    AnyMap parameters() const {
        AnyMap out;
        out["type"] = type();
        getParameters(out);
        return out;
    }

    //! Get the units for converting the leading term in the reaction rate expression.
    //!
    //! These units are often the same as the units of the rate expression itself, but
    //! not always; sticking coefficients are a notable exception.
    //! @since New in %Cantera 3.0
    const Units& conversionUnits() const {
        return m_conversion_units;
    }

    //! Set the units of the reaction rate expression
    //!
    //! Used to determine the units that should be used for converting terms in the
    //! reaction rate expression, which often have the same units (for example, the
    //! Arrhenius pre-exponential) but may also be different (for example, sticking
    //! coefficients).
    //! @since New in %Cantera 3.0
    virtual void setRateUnits(const UnitStack& rate_units) {
        if (rate_units.size() > 1) {
            m_conversion_units = rate_units.product();
        } else {
            m_conversion_units = rate_units.standardUnits();
        }
    }

    //! Check basic syntax and settings of reaction rate expression
    virtual void check(const string& equation) {}

    //! Validate the reaction rate expression
    virtual void validate(const string& equation, const Kinetics& kin) {}

    //! Reaction rate index within kinetics evaluator
    size_t rateIndex() const {
        return m_rate_index;
    }

    //! Set reaction rate index within kinetics evaluator
    void setRateIndex(size_t idx) {
        m_rate_index = idx;
    }

    //! Set context of reaction rate evaluation
    //! @param rxn  Reaction object associated with rate
    //! @param kin  Kinetics object used for rate evaluation
    //! This method allows for passing of information specific to the associated
    //! reaction when a ReactionRate object is added a MultiRate reaction evaluator.
    virtual void setContext(const Reaction& rxn, const Kinetics& kin) {
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

    //! Evaluate reaction rate based on temperature and an extra vector parameter.
    //! Specific rate parameterizations may require additional parameters, which
    //! are specific to the derived ReactionRate object.
    //! @param T  temperature [K]
    //! @param extra  extra vector parameter used by parameterization
    //! Used in conjunction with MultiRateBase::evalSingle / ReactionRate::eval.
    //! This method allows for testing of a reaction rate expression outside of
    //! Kinetics reaction rate evaluators.
    //! @warning  This method is an experimental part of the %Cantera API and
    //!     may be changed or removed without notice.
    double eval(double T, const vector<double>& extra) {
        _evaluator().update(T, extra);
        return _evaluator().evalSingle(*this);
    }

    //! Get flag indicating whether reaction rate is set up correctly
    bool valid() const {
        return m_valid;
    }

    //! Boolean indicating whether rate has compositional dependence
    //! @since New in %Cantera 3.0
    bool compositionDependent() {
        return m_composition_dependent_rate;
    }

    //! Set rate compositional dependence
    //! @since New in %Cantera 3.0
    void setCompositionDependence(bool comp_dep) {
        m_composition_dependent_rate = comp_dep;
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
    size_t m_rate_index = npos;

    //! Flag indicating whether reaction rate is set up correctly
    bool m_valid = false;

    //! Flag indicating composition dependent rate
    bool m_composition_dependent_rate = false;

    //! Units of the leading term in the reaction rate expression
    Units m_conversion_units{0.};

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
