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
//! @since New in %Cantera 3.0
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

    //! Set the type of the ReactionData class. This should match the corresponding
    //! ReactionRate class's type
    void setType(const string& name) {
        m_rateType = name;
    }

    //! Get the external language wrapper for this ReactionData object
    shared_ptr<ExternalHandle> getWrapper() const {
        return m_wrappedData;
    }

    //! Set the external language wrapper for this ReactionData object
    void setWrapper(shared_ptr<ExternalHandle> wrapper) {
        m_wrappedData = wrapper;
    }

protected:
    //! The reaction rate type
    string m_rateType;

    //! An external language's wrapper for the Solution object where this ReactionData
    //! object is being used
    shared_ptr<ExternalHandle> m_wrappedSolution;

    //! An external language's wrapper for this ReactionData object
    shared_ptr<ExternalHandle> m_wrappedData;

    //! Delegated `update` method taking the Solution wrapper as its argument
    function<double(void*)> m_update;
};

//! Delegate methods of the ReactionRate class to external functions
//!
//! @since New in %Cantera 3.0
//! @ingroup otherRateGroup
class ReactionRateDelegator : public Delegator, public ReactionRate
{
public:
    ReactionRateDelegator();

    unique_ptr<MultiRateBase> newMultiRate() const override;

    //! Set the reaction type based on the user-provided reaction rate parameterization
    void setType(const string& type) {
        m_rateType = type;
    }

    const string type() const override {
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

    void getParameters(AnyMap& node) const override {
        m_getParameters(node);
    }

    void validate(const string& equation, const Kinetics& kin) override;

private:
    //! The name of the reaction rate type
    string m_rateType;

    //! Delegated `evalFromStruct` method taking a pointer to the corresponding
    //! ReactionData wrapper object
    function<double(void*)> m_evalFromStruct;

    function<void(const string&, void*)> m_validate;
    function<void(const AnyMap&, const UnitStack&)> m_setParameters;
    function<void(AnyMap&)> m_getParameters;
};

}

#endif
