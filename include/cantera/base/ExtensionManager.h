//! @file ExtensionManager.h

#ifndef CT_EXTENSIONMANAGER_H
#define CT_EXTENSIONMANAGER_H

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/ctexceptions.h"

namespace Cantera
{

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
};

//! A base class for managing the lifetime of an external object, such as a Python
//! object used by a Delegator
class ExternalHandle
{
public:
    virtual ~ExternalHandle() = default;
};

}

#endif
