//! @file ExtensionManager.h

#ifndef CT_EXTENSIONMANAGER_H
#define CT_EXTENSIONMANAGER_H

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/ctexceptions.h"
#include <functional>

namespace Cantera
{

class ReactionDataDelegator;
class Solution;

//! A base class for managing the lifetime of an external object, such as a Python
//! object used by a Delegator
class ExternalHandle
{
public:
    ExternalHandle() {}
    ExternalHandle(const ExternalHandle&) = delete;
    virtual ~ExternalHandle() = default;

    //! Get the underlying external object
    virtual void* get() {
        throw NotImplementedError("ExternalHandle::get");
    }
};

//! Base class for managing user-defined Cantera extensions written in other languages
//!
//! @since New in Cantera 3.0
class ExtensionManager
{
public:
    virtual ~ExtensionManager() = default;

    //! Register ReactionRate defined in a user extension with ReactionRateFactory
    //! @param extensionName
    virtual void registerRateBuilders(const std::string& extensionName) {
        throw NotImplementedError("ExtensionManager::registerRateBuilders");
    };

    static void wrapReactionData(const std::string& rateName,
                                 ReactionDataDelegator& data);

    static shared_ptr<ExternalHandle> wrapSolution(const std::string& rateName,
                                                   shared_ptr<Solution> soln);

    static void registerReactionDataLinker(const std::string& rateName,
        std::function<void(ReactionDataDelegator&)> link);

    static void registerSolutionLinker(const std::string& rateName,
        std::function<shared_ptr<ExternalHandle>(shared_ptr<Solution>)> link);

protected:
    static std::map<std::string,
        std::function<void(ReactionDataDelegator&)>> s_ReactionData_linkers;

    static std::map<std::string,
        std::function<shared_ptr<ExternalHandle>(shared_ptr<Solution>)>> s_Solution_linkers;

};

}

#endif
