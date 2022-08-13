//! @file ExtensionManager.h

#ifndef CT_EXTENSIONMANAGER_H
#define CT_EXTENSIONMANAGER_H

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/ctexceptions.h"

namespace Cantera
{

//! Base class for managing user-defined Cantera extensions written in other languages
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

}

#endif
