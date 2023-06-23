//! @file Func1Factory.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef FUNC1_FACTORY_H
#define FUNC1_FACTORY_H

#include "cantera/numerics/Func1.h"
#include "cantera/base/FactoryBase.h"

namespace Cantera
{

//! Factory class to create domain objects
//!
//! This class is mainly used via the newFunc1() function, for example:
//!
//! ```cpp
//!     shared_ptr<Func1> d1 = newFunc1("sin", 0, {1.0});
//! ```
class Func1Factory : public Factory<Func1, size_t, const vector<double>&>
{
public:
    /**
     * Return a pointer to the factory. On the first call, a new instance is
     * created. Since there is no need to instantiate more than one factory,
     * on all subsequent calls, a pointer to the existing factory is returned.
     */
    static Func1Factory* factory();

    virtual void deleteFactory();

private:
    //! Pointer to the single instance of the factory
    static Func1Factory* s_factory;

    //! default constructor, which is defined as private
    Func1Factory();

    //! Mutex for use when calling the factory
    static std::mutex domain_mutex;
};

//! Create a new Func1 functor object.
//! @param func1Type  string identifying function type.
//! @param n  Integer; definition depends on the function type.
//! @param params  Parameter vector; definition depends on the function type.
shared_ptr<Func1> newFunc1(const string& func1Type, size_t n=0,
                           const vector<double>& params={})
{
    return shared_ptr<Func1>(
        Func1Factory::factory()->create(func1Type, n, params));
}

}

#endif
