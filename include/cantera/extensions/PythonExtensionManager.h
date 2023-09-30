//! @file PythonExtensionManager.h

#ifndef CT_PYTHONEXTENSIONMANAGER_H
#define CT_PYTHONEXTENSIONMANAGER_H

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/ExtensionManager.h"

namespace Cantera
{

//! Class for managing user-defined %Cantera extensions written in Python
//!
//! Handles Python initialization if the main application is not the Python interpreter.
//!
//! Imports a user-specified module, which must be on the Python path and registers
//! user-defined classes that are marked with the `@extension` decorator. See the
//! documentation for [\@extension](../python/utilities.html#cantera.extension)
//! in the Python documentation for more information.
//!
//! @since New in %Cantera 3.0
class PythonExtensionManager : public ExtensionManager
{
public:
    void registerRateBuilders(const string& extensionName) override;

    void registerRateBuilder(const string& moduleName,
        const string& className, const string& rateName) override;

    static void registerSelf();

    void registerRateDataBuilder(const string& moduleName,
        const string& className, const string& rateName) override;

private:
    PythonExtensionManager() = default;
};

}

#endif
