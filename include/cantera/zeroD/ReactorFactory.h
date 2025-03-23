//! @file ReactorFactory.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef REACTOR_FACTORY_H
#define REACTOR_FACTORY_H

#include "cantera/zeroD/ReactorNode.h"
#include "cantera/base/FactoryBase.h"

namespace Cantera
{

//! Factory class to create reactor objects
//!
//! This class is mainly used via the newReactor() function, for example:
//!
//! ```cpp
//!     shared_ptr<ReactorNode> r1 = newReactor("IdealGasReactor");
//! ```
class ReactorFactory : public Factory<ReactorNode, shared_ptr<Solution>, const string&>
{
public:
    static ReactorFactory* factory();

    void deleteFactory() override;

private:
    static ReactorFactory* s_factory;
    static std::mutex reactor_mutex;
    ReactorFactory();
};

//! @defgroup reactorGroup Reactors
//! Zero-dimensional objects representing stirred reactors.
//! Reactors simulate time-dependent behavior considering gas-phase chemistry.
//! Reactor objects should be instantiated via the newReactor() function, for example:
//!
//! ```cpp
//!     shared_ptr<ReactorNode> r1 = newReactor("IdealGasReactor");
//! ```
//! @ingroup zerodGroup
//! @{

//! Create a Reactor object of the specified type and contents
//! @since Starting in %Cantera 3.1, this method requires a valid Solution object and
//!     returns a `shared_ptr<ReactorNode>` instead of a `ReactorNode*`.
shared_ptr<ReactorNode> newReactor(
    const string& model, shared_ptr<Solution> contents, const string& name="(none)");

//! @}

}

#endif
