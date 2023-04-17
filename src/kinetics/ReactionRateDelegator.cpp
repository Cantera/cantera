//! @file ReactionRateDelegator.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/ReactionRateDelegator.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/base/Solution.h"

namespace Cantera
{

ReactionDataDelegator::ReactionDataDelegator()
{
    install("update", m_update,
        [](void*) {
            throw NotImplementedError("ReactionDataDelegator::update");
            return 0.0; // necessary to set lambda's function signature
        }
    );
}

bool ReactionDataDelegator::update(const ThermoPhase& phase, const Kinetics& kin)
{
    if (!m_wrappedSolution) {
        auto wrapperType = ExtensionManager::getSolutionWrapperType(m_rateType);
        auto soln = kin.root();
        if (!soln) {
            throw CanteraError("ReactionDataDelegator::update",
                "Phase must be instantiated as a Solution to use extensible "
                "reactions of type '{}'", m_rateType);
        }
        if (soln->getExternalHandle(wrapperType)) {
            m_wrappedSolution = soln->getExternalHandle(wrapperType);
        } else {
            m_wrappedSolution = ExtensionManager::wrapSolution(wrapperType, soln);
        }
    }
    double needsUpdate = m_update(m_wrappedSolution->get());
    return needsUpdate != 0.0;
}

ReactionRateDelegator::ReactionRateDelegator()
{
    install("evalFromStruct", m_evalFromStruct,
        [](void*) {
            throw NotImplementedError("ReactionRateDelegator::evalFromStruct");
            return 0.0; // necessary to set lambda's function signature
        }
    );
    install("setParameters", m_setParameters,
        [this](const AnyMap& node, const UnitStack& units) {
            ReactionRate::setParameters(node, units); });
    install("getParameters", m_getParameters,
        [this](AnyMap& node) { ReactionRate::getParameters(node); });
    install("validate", m_validate,
        [](const string& equation, void* soln) {
            throw NotImplementedError("ReactionRateDelegator::validate"); });
}

unique_ptr<MultiRateBase> ReactionRateDelegator::newMultiRate() const
{
    auto multirate = make_unique<MultiRate<ReactionRateDelegator,
                                           ReactionDataDelegator>>();
    multirate->sharedData().setType(m_rateType);
    ExtensionManager::wrapReactionData(m_rateType, multirate->sharedData());
    return multirate;
}

void ReactionRateDelegator::validate(const string& equation, const Kinetics& kin)
{
    auto soln = kin.root();
    if (!soln) {
        throw CanteraError("ReactionRateDelegator::validate",
            "Phase must be instantiated as a Solution to use extensible "
            "reactions of type '{}'", m_rateType);
    }
    auto wrapperType = ExtensionManager::getSolutionWrapperType(m_rateType);
    auto wrappedSoln = soln->getExternalHandle(wrapperType);
    if (!wrappedSoln) {
        wrappedSoln = ExtensionManager::wrapSolution(wrapperType, soln);
    }

    try {
        m_validate(equation, wrappedSoln->get());
    } catch (CanteraError& err) {
        throw InputFileError("'" + m_rateType + "' validate", m_input,
            err.getMessage());
    }
}

}
