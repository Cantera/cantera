//! @file ExtensionManager.h

#ifndef CT_EXTENSIONMANAGER_H
#define CT_EXTENSIONMANAGER_H

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/ctexceptions.h"

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

//! Base class for managing user-defined %Cantera extensions written in other languages
//!
//! @since New in %Cantera 3.0
class ExtensionManager
{
public:
    virtual ~ExtensionManager() = default;

    //! Register ReactionRate defined in a user extension with ReactionRateFactory
    //! @param extensionName
    virtual void registerRateBuilders(const string& extensionName) {
        throw NotImplementedError("ExtensionManager::registerRateBuilders");
    };

    //! Register a user-defined ReactionRate implementation with ReactionRateFactory
    //! @param extensionName The name of the library/module containing the user-defined
    //!     rate. For example, the module name for rates implemented in Python.
    //! @param className The name of the rate in the user's code. For example, the
    //!     Python class name
    //! @param rateName The name used to construct a rate of this type using
    //!     the newReactionRate() function or from a YAML input file
    virtual void registerRateBuilder(const string& extensionName,
        const string& className, const string& rateName)
    {
        throw NotImplementedError("ExtensionManager::registerRateBuilder");
    }

    //! Register a user-defined ReactionData implementation
    //! @param extensionName The name of the library/module containing the user-defined
    //!     type. For example, the module name for rates implemented in Python.
    //! @param className The name of the data object in the user's code. For example,
    //!     the Python class name
    //! @param rateName The name of the corresponding reaction rate type
    virtual void registerRateDataBuilder(const string& extensionName,
        const string& className, const string& rateName)
    {
        throw NotImplementedError("ExtensionManager::registerRateDataBuilder");
    }

    //! Create an object in an external language that wraps the specified ReactionData
    //! object
    //!
    //! @param rateName  The name of the reaction rate type, which corresponds to the
    //!     name used to register the wrapper generator using registerReactionDataLinker
    //! @param data  The ReactionData object to be wrapped
    static void wrapReactionData(const string& rateName, ReactionDataDelegator& data);

    //! Create an object in an external language that wraps the specified Solution
    //! object.
    //!
    //! @param wrapperType  A name specifying the wrapper type, which corresponds to
    //!     the name used to register the wrapper generator using registerSolutionLinker
    //! @param soln  The Solution object to be wrapped
    static shared_ptr<ExternalHandle> wrapSolution(const string& wrapperType,
                                                   shared_ptr<Solution> soln);

    //! Register a function that can be used to create wrappers for ReactionData objects
    //! in an external language and link them to the corresponding C++ object
    //!
    //! @param rateName  The name of the reaction rate type
    //! @param wrapperName  The name used for Solution wrappers to be used with this
    //!     object, corresponding to a type registered with registerSolutionLinker().
    //! @param link  Function that creates ReactionData wrapper and links it to the
    //!     provided C++ object
    static void registerReactionDataLinker(const string& rateName,
        const string& wrapperName, function<void(ReactionDataDelegator&)> link);

    //! Register a function that can be used to create wrappers for Solution objects in
    //! an external language and link it to the corresponding C++ objects
    static void registerSolutionLinker(const string& wrapperName,
        function<shared_ptr<ExternalHandle>(shared_ptr<Solution>)> link);

    //! Get the Solution wrapper type corresponding to the specified user-defined
    //! reaction rate type.
    static string getSolutionWrapperType(const string& userType);

protected:
    //! Functions for wrapping and linking ReactionData objects
    static map<string, function<void(ReactionDataDelegator&)>> s_ReactionData_linkers;

    //! Functions for wrapping and linking Solution objects
    static map<string,
        function<shared_ptr<ExternalHandle>(shared_ptr<Solution>)>> s_Solution_linkers;

    //! Mapping from user-defined rate types to Solution wrapper types
    static map<string, string> s_userTypeToWrapperType;
};

}

#endif
