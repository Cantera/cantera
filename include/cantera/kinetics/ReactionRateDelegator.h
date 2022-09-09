//! @file ReactionRateDelegator.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_REACTION_RATE_DELEGATOR_H
#define CT_REACTION_RATE_DELEGATOR_H

#include "ReactionRate.h"
#include "cantera/base/Delegator.h"
#include "Arrhenius.h"

namespace Cantera
{

//! Delegate methods of the ReactionRate class to external functions
//!
//! @since New in Cantera 3.0
class ReactionRateDelegator : public Delegator, public ReactionRate
{
public:
    ReactionRateDelegator() {
        install("evalFromStruct", m_evalFromStruct,
            [](void*) {
                throw NotImplementedError("ReactionRateDelegator::evalFromStruct");
                return 0.0; // necessary to set lambda's function signature
            }
        );
        install("setParameters", m_setParameters,
            [this](const AnyMap& node, const UnitStack& units) {
                ReactionRate::setParameters(node, units); });
    }

    virtual unique_ptr<MultiRateBase> newMultiRate() const override {
        return unique_ptr<MultiRateBase>(
            new MultiRate<ReactionRateDelegator, ArrheniusData>);
    }

    //! Set the reaction type based on the user-provided reaction rate parameterization
    void setType(const std::string& type) {
        m_rateType = type;
    }

    virtual const std::string type() const override {
        return m_rateType;
    }

    // Delegatable methods

    //! Evaluate reaction rate
    //!
    //! @param shared_data  data shared by all reactions of a given type
    double evalFromStruct(const ArrheniusData& shared_data) {
        // @TODO: replace passing pointer to temperature with a language-specific
        //     wrapper of the ReactionData object
        double T = shared_data.temperature;
        return m_evalFromStruct(&T);
    }

    void setParameters(const AnyMap& node, const UnitStack& units) override {
        m_setParameters(node, units);
    }

private:
    std::string m_rateType;
    std::function<double(void*)> m_evalFromStruct;
    std::function<void(const AnyMap&, const UnitStack&)> m_setParameters;
};

}

#endif
