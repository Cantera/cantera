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
    }

    virtual unique_ptr<MultiRateBase> newMultiRate() const override {
        return unique_ptr<MultiRateBase>(
            new MultiRate<ReactionRateDelegator, ArrheniusData>);
    }

    virtual const std::string type() const override {
        return "ReactionRateDelegator";
    }

    // Delegatable methods

    double evalFromStruct(const ArrheniusData& shared_data) {
        // @TODO: replace passing pointer to temperature with a language-specific
        //     wrapper of the ReactionData object
        double T = shared_data.temperature;
        return m_evalFromStruct(&T);
    }

private:
    std::function<double(void*)> m_evalFromStruct;
};

}

#endif
