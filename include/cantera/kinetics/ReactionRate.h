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

#include "MultiRateBase.h"
#include "cantera/base/AnyMap.h"
#include "cantera/base/Units.h"
#include "cantera/base/ctexceptions.h"

namespace Cantera
{


//! Abstract base class for reaction rate definitions; this base class is used by
//! user-facing APIs to access reaction rate objects.
class ReactionRate
{
public:
    ReactionRate() : m_rate_index(npos) {}

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

    virtual unique_ptr<MultiRateBase> newMultiRate() const = 0;

    virtual const std::string type() const = 0;

    virtual bool usesUpdate() const { return true; };

    //! Set parameters
    //! @param node  AnyMap object containing reaction rate specification
    virtual void setParameters(const AnyMap& node, const UnitsVector& units) {
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

    virtual void check(const std::string& equation, const AnyMap& node) {}

    //! Validate the reaction rate expression
    //! @todo add Kinetics object
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
    double eval(double T) {
        _evaluator().update(T);
        return _evaluator().evalSingle(*this);
    }

    //! Evaluate reaction rate based on temperature and pressure
    //! @param T  temperature [K]
    //! @param P  pressure [Pa]
    double eval(double T, double P) {
        _evaluator().update(T, P);
        return _evaluator().evalSingle(*this);
    }

    //! Evaluate reaction rate based on bulk phase
    //! @param bulk  object representing bulk phase
    //! @param kin  object representing kinetics (not required for all rate types)
    double eval(const ThermoPhase& bulk, const Kinetics& kin) {
        _evaluator().update(bulk, kin);
        return _evaluator().evalSingle(*this);
    }

    //! Update reaction rate data based on temperature
    //! @param T  temperature [K]
    void update(double T) {
        _evaluator().update(T);
    }

    //! Update reaction rate data based on temperature and pressure
    //! @param T  temperature [K]
    //! @param P  pressure [Pa]
    void update(double T, double P) {
        _evaluator().update(T, P);
    }

    //! Update reaction rate data based on bulk phase
    //! @param bulk  object representing bulk phase
    //! @param kin  object representing kinetics (not required for all rate types)
    void update(const ThermoPhase& bulk, const Kinetics& kin) {
        _evaluator().update(bulk, kin);
    }

    //! Evaluate reaction rate derivative based on temperature
    //! @param T  temperature [K]
    double ddT(double T) {
        _evaluator().update(T);
        return _evaluator().ddTSingle(*this);
    }

    //! Evaluate reaction rate derivative based on temperature and pressure
    //! @param T  temperature [K]
    //! @param P  pressure [Pa]
    double ddT(double T, double P) {
        _evaluator().update(T, P);
        return _evaluator().ddTSingle(*this);
    }

    //! Evaluate reaction rate derivative based on bulk phase
    //! @param bulk  object representing bulk phase
    //! @param kin  object representing kinetics (not required for all rate types)
    double ddT(const ThermoPhase& bulk, const Kinetics& kin) {
        _evaluator().update(bulk, kin);
        return _evaluator().ddTSingle(*this);
    }

protected:
    //! Get parameters
    //! Store the parameters of a ReactionRate needed to reconstruct an identical
    //! object. Does not include user-defined fields available in the #input map.
    virtual void getParameters(AnyMap& rateNode) const {
        throw NotImplementedError("ReactionRate::getParameters",
                                  "Not implemented by '{}' object.", type());
    }

    //! Input data used for specific models
    AnyMap m_input;

    //! Index of reaction rate within kinetics evaluator
    size_t m_rate_index;

private:
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
