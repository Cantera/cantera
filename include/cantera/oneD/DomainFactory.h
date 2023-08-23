//! @file DomainFactory.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef DOMAIN_FACTORY_H
#define DOMAIN_FACTORY_H

#include "cantera/oneD/Domain1D.h"
#include "cantera/base/FactoryBase.h"

namespace Cantera
{

//! Factory class to create domain objects
//!
//! This class is mainly used via the newDomain() function, for example:
//!
//! ```cpp
//!     shared_ptr<Domain1D> d1 = newDomain("Inlet", sol, "reactants");
//! ```
class DomainFactory : public Factory<Domain1D, shared_ptr<Solution>, const string&>
{
public:
    /**
     * Return a pointer to the factory. On the first call, a new instance is
     * created. Since there is no need to instantiate more than one factory,
     * on all subsequent calls, a pointer to the existing factory is returned.
     */
    static DomainFactory* factory();

    void deleteFactory() override;

private:
    //! Pointer to the single instance of the factory
    static DomainFactory* s_factory;

    //! default constructor, which is defined as private
    DomainFactory();

    //! Mutex for use when calling the factory
    static std::mutex domain_mutex;
};

//! Create a Domain object of the specified type. An optional template argument will
//! dynamically cast Domain1D to the desired specialization.
//! @param domainType  string identifying domain type.
//! @param solution  Solution holding ThermoPhase, Kinetics and Transport objects.
//! @param id  string identifier describing domain. If omitted, id defaults to the
//!     domain type identifier.
//! @ingroup onedGroup
template <class T=Domain1D>
shared_ptr<T> newDomain(
    const string& domainType, shared_ptr<Solution> solution, const string& id="")
{
    string id_ = id;
    if (id_ == "") {
        id_ = domainType;
    }
    auto ret = std::dynamic_pointer_cast<T>(
        shared_ptr<Domain1D>(
            DomainFactory::factory()->create(domainType, solution, id_)));
    if (!ret) {
        throw CanteraError("newDomain",
            "Invalid cast: unable to access 'Domain1D' as '{}'.", demangle(typeid(T)));
    }
    return ret;
}

}

#endif
