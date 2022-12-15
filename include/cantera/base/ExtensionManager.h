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

    //! Create an object in an external language that wraps the specified ReactionData
    //! object
    //!
    //! @param rateName  The name of the reaction rate type, which corresponds to the
    //!     name used to register the wrapper generator using registerReactionDataLinker
    //! @param data  The ReactionData object to be wrapped
    static void wrapReactionData(const std::string& rateName,
                                 ReactionDataDelegator& data);

    //! Create an object in an external language that wraps the specified Solution
    //! object.
    //!
    //! @param wrapperType  A name specifying the wrapper type, which corresponds to
    //!     the name used to register the wrapper generator using registerSolutionLinker
    //! @param soln  The Solution object to be wrapped
    static shared_ptr<ExternalHandle> wrapSolution(const std::string& wrapperType,
                                                   shared_ptr<Solution> soln);

    //! Register a function that can be used to create wrappers for ReactionData objects
    //! in an external language and link them to the corresponding C++ object
    static void registerReactionDataLinker(const std::string& rateName,
        std::function<void(ReactionDataDelegator&)> link);

    //! Register a function that can be used to create wrappers for Solution objects in
    //! an external language and link it to the corresponding C++ objects
    static void registerSolutionLinker(const std::string& wrapperName,
        std::function<shared_ptr<ExternalHandle>(shared_ptr<Solution>)> link);

protected:
    //! Functions for wrapping and linking ReactionData objects
    static std::map<std::string,
        std::function<void(ReactionDataDelegator&)>> s_ReactionData_linkers;

    //! Functions for wrapping and linking Solution objects
    static std::map<std::string,
        std::function<shared_ptr<ExternalHandle>(shared_ptr<Solution>)>> s_Solution_linkers;

};

}

#endif
