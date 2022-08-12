//! @file ExtensionManagerFactory.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/ExtensionManagerFactory.h"

#ifdef CT_HAS_PYTHON
#include "cantera/extensions/PythonExtensionManager.h"
#endif

using namespace std;

namespace Cantera
{

ExtensionManagerFactory* ExtensionManagerFactory::s_factory = 0;
mutex ExtensionManagerFactory::s_mutex;

ExtensionManagerFactory::ExtensionManagerFactory()
{
    #ifdef CT_HAS_PYTHON
    reg("python", []() { return new PythonExtensionManager(); });
    #endif
}

ExtensionManagerFactory& ExtensionManagerFactory::factory()
{
    unique_lock<mutex> lock(s_mutex);
    if (!s_factory) {
        s_factory = new ExtensionManagerFactory();
    }
    return *s_factory;
}

void ExtensionManagerFactory::deleteFactory()
{
    unique_lock<mutex> lock(s_mutex);
    delete s_factory;
    s_factory = 0;
}

}
