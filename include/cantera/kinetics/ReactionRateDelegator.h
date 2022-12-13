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

//! Delegate methods of the ReactionData class to external functions
//!
//! @since New in Cantera 3.0
class ReactionDataDelegator : public Delegator, public ReactionData
{
public:
    ReactionDataDelegator();

    bool update(const ThermoPhase& phase, const Kinetics& kin) override;

    void update(double T) override {
        throw NotImplementedError("ReactionDataDelegator",
            "Not implemented for delegated reaction rates");
    }

    using ReactionData::update;

    void setType(const std::string& name) {
        m_rateType = name;
    }

    shared_ptr<ExternalHandle> getWrapper() const {
        return m_wrappedData;
    }

    void setWrapper(shared_ptr<ExternalHandle> wrapper) {
        m_wrappedData = wrapper;
    }

    void setSolutionWrapperType(const std::string& type) {
        m_solutionWrapperType = type;
    }

protected:
    std::string m_rateType;
    std::string m_solutionWrapperType;
    shared_ptr<ExternalHandle> m_wrappedSolution;
    shared_ptr<ExternalHandle> m_wrappedData;

    std::function<double(void*)> m_update;
};

//! Delegate methods of the ReactionRate class to external functions
//!
//! @since New in Cantera 3.0
class ReactionRateDelegator : public Delegator, public ReactionRate
{
public:
    ReactionRateDelegator();

    virtual unique_ptr<MultiRateBase> newMultiRate() const override;

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
    double evalFromStruct(const ReactionDataDelegator& shared_data) {
        return m_evalFromStruct(shared_data.getWrapper()->get());
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
