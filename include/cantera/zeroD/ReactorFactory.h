//! @file ReactorFactory.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef REACTOR_FACTORY_H
#define REACTOR_FACTORY_H

#include "cantera/zeroD/ReactorBase.h"
#include "cantera/base/FactoryBase.h"

namespace Cantera
{

//! Factory class to create reactor objects
//!
//! This class is mainly used via the newReactor3() function, for example:
//!
//! ```cpp
//!     shared_ptr<ReactorBase> r1 = newReactor3("IdealGasReactor");
//! ```
//!
//! @ingroup ZeroD
class ReactorFactory : public Factory<ReactorBase>
{
public:
    static ReactorFactory* factory();

    virtual void deleteFactory();

    //! Create a new reactor by type name.
    /*!
     * @param reactorType the type to be created.
     * @deprecated  To be removed after Cantera 3.0; replaceable by newReactor3.
     */
    virtual ReactorBase* newReactor(const std::string& reactorType);

private:
    static ReactorFactory* s_factory;
    static std::mutex reactor_mutex;
    ReactorFactory();
};

//! Create a Reactor object of the specified type
//! @deprecated  To be changed after Cantera 3.0; for new behavior, see newReactor3.
//! @ingroup ZeroD
ReactorBase* newReactor(const string& model);

//! Create a Reactor object of the specified type
//! @ingroup ZeroD
//! @since New in Cantera 3.0.
//! @todo Transition back to newReactor after Cantera 3.0
shared_ptr<ReactorBase> newReactor3(const string& model);

}

#endif
