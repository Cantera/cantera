/**
 * @file ReactionRateBase.h
 *
 * @warning This file is an experimental part of the %Cantera API and
 *    may be changed or removed without notice.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/ctexceptions.h"
#include "cantera/base/AnyMap.h"
#include "cantera/base/Units.h"

#ifndef CT_REACTIONRATEBASE_H
#define CT_REACTIONRATEBASE_H

namespace Cantera
{

class MultiRateBase;
class ThermoPhase;
class Kinetics;

//! Header-only abstract base class for reaction rate definitions
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
    ReactionRateBase() : units(0.), m_idx(npos) {}
    virtual ~ReactionRateBase() {}

public:
    //! Identifier of reaction type
    virtual const std::string type() const = 0;

    //! Create multi-rate evaluator
    virtual unique_ptr<MultiRateBase> newMultiRate() const = 0;

    //! Update reaction rate data based on temperature
    //! @param T  temperature [K]
    virtual void update(double T) = 0;

    //! Update reaction rate data based on temperature and pressure
    //! @param T  temperature [K]
    //! @param P  pressure [Pa]
    virtual void update(double T, double P) = 0;

    //! Update reaction rate data based on bulk phase
    //! @param bulk  object representing bulk phase
    //! @param kin  object representing kinetics (not required for all rate types)
    virtual void update(const ThermoPhase& bulk, const Kinetics& kin) = 0;

    //! Evaluate reaction rate based on temperature
    //! @param T  temperature [K]
    virtual double eval(double T) const = 0;

    //! Evaluate reaction rate based on temperature and pressure
    //! @param T  temperature [K]
    //! @param P  pressure [Pa]
    virtual double eval(double T, double P) const = 0;

    //! Evaluate reaction rate based on bulk phase
    //! @param bulk  object representing bulk phase
    //! @param kin  object representing kinetics (not required for all rate types)
    virtual double eval(const ThermoPhase& bulk, const Kinetics& kin) const = 0;

    //! Evaluate reaction rate derivative based on temperature
    //! @param T  temperature [K]
    virtual double ddT(double T) const = 0;

    //! Evaluate reaction rate derivative based on temperature and pressure
    //! @param T  temperature [K]
    //! @param P  pressure [Pa]
    virtual double ddT(double T, double P) const = 0;

    //! Evaluate reaction rate derivative based on bulk phase
    //! @param bulk  object representing bulk phase
    //! @param kin  object representing kinetics (not required for all rate types)
    virtual double ddT(const ThermoPhase& bulk, const Kinetics& kin) const = 0;

    //! Validate the reaction rate expression
    virtual void validate(const std::string& equation) {};

    //! Return the parameters such that an identical Reaction could be reconstructed
    //! using the newReaction() function. Behavior specific to derived classes is
    //! handled by the getParameters() method.
    //! @param rate_units  units used for rate parameters
    AnyMap parameters(const Units& rate_units) const {
        AnyMap out;
        getParameters(out, rate_units);
        return out;
    }

    //! Return parameters using original unit system
    AnyMap parameters() const {
        return parameters(units);
    }

    //! Set parameters
    //! @param node  AnyMap object containing reaction rate specification
    //! @param rate_units  Description of units used for rate parameters
    virtual void setParameters(const AnyMap& node, const Units& rate_units) {
        units = rate_units;
        input = node;
    }

    //! Set rate units
    //! @param rate_units  Description of units used for rate parameters
    virtual void setUnits(const Units& rate_units) {
        units = rate_units;
    }

    //! Reaction rate index within kinetics evaluator
    size_t rateIndex() const {
        return m_idx;
    }

    //! Set reaction rate index within kinetics evaluator
    virtual void setRateIndex(size_t idx) {
        m_idx = idx;
    }

protected:
    //! Get parameters
    //! Store the parameters of a ReactionRate needed to reconstruct an identical
    //! object. Does not include user-defined fields available in the #input map.
    virtual void getParameters(AnyMap& rateNode, const Units& rate_units) const {
        throw NotImplementedError("ReactionRate::getParameters",
                                  "Not implemented by '{}' object.", type());
    }

    //! Input data used for specific models
    AnyMap input;

    //! The units of the rate constant. These are determined by the units of the
    //! standard concentration of the reactant species' phases of the phase
    //! where the reaction occurs.
    Units units;

    size_t m_idx; //!< Index of reaction rate within kinetics evaluator
};

}

#endif
