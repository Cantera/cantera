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
        auto soln = kin.root();
        if (!soln) {
            throw CanteraError("ReactionDataDelegator::update",
                "Phase must be instantiated as a Solution to use extensible "
                "reactions of type '{}'", m_rateType);
        }
        if (soln->getExternalHandle(m_solutionWrapperType)) {
            m_wrappedSolution = soln->getExternalHandle(m_solutionWrapperType);
        } else {
            m_wrappedSolution = ExtensionManager::wrapSolution(m_solutionWrapperType,
                                                               soln);
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
}

unique_ptr<MultiRateBase> ReactionRateDelegator::newMultiRate() const
{
    auto multirate = make_unique<MultiRate<ReactionRateDelegator,
                                           ReactionDataDelegator>>();
    multirate->sharedData().setType(m_rateType);
    ExtensionManager::wrapReactionData(m_rateType, multirate->sharedData());
    return multirate;
}

}
