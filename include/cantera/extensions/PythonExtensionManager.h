//! @file PythonExtensionManager.h

#ifndef CT_PYTHONEXTENSIONMANAGER_H
#define CT_PYTHONEXTENSIONMANAGER_H

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/ExtensionManager.h"

namespace Cantera
{

//! Class for managing user-defined Cantera extensions written in Python
//!
//! Handles Python initialization if the main application is not the Python interpreter.
class PythonExtensionManager : public ExtensionManager
{
public:
    PythonExtensionManager();
    virtual void registerRateBuilders(const std::string& extensionName) override;
};

}

#endif
