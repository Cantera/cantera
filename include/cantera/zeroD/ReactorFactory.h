//! @file ReactorFactory.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef REACTOR_FACTORY_H
#define REACTOR_FACTORY_H

#include "cantera/zeroD/ReactorBase.h"
#include "cantera/base/FactoryBase.h"

namespace Cantera
{

class Reactor;
class Reservoir;

//! Factory class to create reactor objects
//!
//! This class is mainly used via the newReactor() function, for example:
//!
//! ```cpp
//!     shared_ptr<ReactorBase> r1 = newReactor("IdealGasReactor");
//! ```
class ReactorFactory : public Factory<ReactorBase, shared_ptr<Solution>, bool, const string&>
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
//!     shared_ptr<ReactorBase> r1 = newReactor("IdealGasReactor");
//! ```
//! @ingroup zerodGroup
//! @{

//! Create a ReactorBase object of the specified type and contents.
//! @since  New in %Cantera 3.2.
shared_ptr<ReactorBase> newReactorBase(
    const string& model, shared_ptr<Solution> contents, bool clone=true,
    const string& name="(none)");

//! Create a Reactor object of the specified type and contents
//! @since Starting in %Cantera 3.1, this method requires a valid Solution object and
//!     returns a `shared_ptr<ReactorBase>` instead of a `ReactorBase*`.
//! @deprecated  Behavior changes after %Cantera 3.2, when a `shared_ptr<Reactor>` will
//!     be returned. For new behavior, see `newReactor4`.
shared_ptr<ReactorBase> newReactor(
    const string& model, shared_ptr<Solution> contents, const string& name="(none)");

//! Create a Reactor object of the specified type and contents
//! @since  New in %Cantera 3.2. Transitional method returning a `Reactor` object.
shared_ptr<Reactor> newReactor4(
    const string& model, shared_ptr<Solution> contents, bool clone=true,
    const string& name="(none)");

//! Create a Reservoir object with the specified contents
//! @since New in %Cantera 3.2.
shared_ptr<Reservoir> newReservoir(
    shared_ptr<Solution> contents, bool clone=true, const string& name="(none)");

//! Create a ReactorSurface object with the specified contents and adjacent reactors
//! participating in surface reactions.
//! @since  New in %Cantera 3.2.
shared_ptr<ReactorSurface> newReactorSurface(
    shared_ptr<Solution> soln, const vector<shared_ptr<ReactorBase>>& reactors,
    bool clone=true, const string& name="(none)");

//! @}

}

#endif
