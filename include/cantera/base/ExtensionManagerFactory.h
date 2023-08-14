//! @file ExtensionManagerFactory.h

#ifndef CT_EXTENSIONMANAGERFACTORY_H
#define CT_EXTENSIONMANAGERFACTORY_H

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "FactoryBase.h"
#include "ExtensionManager.h"

namespace Cantera
{

//! A factory class for creating ExtensionManager objects
//!
//! @since New in %Cantera 3.0
class ExtensionManagerFactory : public Factory<ExtensionManager>
{
public:
    //! Create a new ExtensionManager
    static shared_ptr<ExtensionManager> build(const string& extensionType) {
        return shared_ptr<ExtensionManager>(factory().create(extensionType));
    }

    //! Delete the static instance of this factory
    void deleteFactory() override;

    //! Static function that returns the static instance of the factory, creating it
    //! if necessary.
    static ExtensionManagerFactory& factory();

private:
    //! static member of the single factory instance
    static ExtensionManagerFactory* s_factory;

    //! Private constructor prevents direct usage
    ExtensionManagerFactory() = default;

    //! Decl for locking mutex for thermo factory singleton
    static std::mutex s_mutex;
};

}

#endif
