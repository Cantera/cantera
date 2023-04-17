//! @file ExtensionManager.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/ExtensionManager.h"
#include "cantera/base/global.h"

namespace Cantera
{

map<string, function<void(ReactionDataDelegator&)>> ExtensionManager::s_ReactionData_linkers = {};
map<string, function<shared_ptr<ExternalHandle>(shared_ptr<Solution>)>> ExtensionManager::s_Solution_linkers = {};
map<string, string> ExtensionManager::s_userTypeToWrapperType = {};

void ExtensionManager::wrapReactionData(const string& rateName,
                                        ReactionDataDelegator& data)
{
    if (s_ReactionData_linkers.count(rateName) != 0) {
        s_ReactionData_linkers.at(rateName)(data);
    } else {
        throw CanteraError("ExtensionManager::wrapReactionData",
                           "No ReactionData delegator for type {} registered",
                           rateName);
    }
}

void ExtensionManager::registerReactionDataLinker(const string& rateName,
    const string& wrapperName, function<void(ReactionDataDelegator& delegator)> link)
{
    s_ReactionData_linkers[rateName] = link;
    s_userTypeToWrapperType[rateName] = wrapperName;
}

shared_ptr<ExternalHandle> ExtensionManager::wrapSolution(
    const string& rateName, shared_ptr<Solution> soln)
{
    if (s_Solution_linkers.count(rateName) != 0) {
        return s_Solution_linkers.at(rateName)(soln);
    } else {
        throw CanteraError("ExtensionManager::wrapSolution",
                           "No Solution linker for type {} registered",
                           rateName);
    }
}

void ExtensionManager::registerSolutionLinker(const string& rateName,
    function<shared_ptr<ExternalHandle>(shared_ptr<Solution>)> link)
{
    s_Solution_linkers[rateName] = link;
}

string ExtensionManager::getSolutionWrapperType(const string& userType)
{
    if (s_userTypeToWrapperType.count(userType)) {
        return s_userTypeToWrapperType[userType];
    } else {
        throw CanteraError("ExtensionManager::getSolutionWrapperType",
                           "No Solution linker for type {} registered", userType);
    }
}

}
