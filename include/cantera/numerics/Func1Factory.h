//! @file Func1Factory.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef FUNC1_FACTORY_H
#define FUNC1_FACTORY_H

#include "cantera/numerics/Func1.h"
#include "cantera/base/FactoryBase.h"

namespace Cantera
{

//! Factory class to create Func1 objects
//!
//! This class is mainly used via the newFunc1() function, for example:
//!
//! ```cpp
//!     shared_ptr<Func1> d1 = newFunc1("sin", {1.0});
//! ```
class Func1Factory : public Factory<Func1, const vector<double>&, size_t>
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


//! Factory class to create Func1 math objects - version A
//!
//! This class is mainly used via the newMath1(const string&, const shared_ptr<Func1>,
//! const shared_ptr<Func1>) function, for example:
//!
//! ```cpp
//!     shared_ptr<Func1> d1 = newMath1("sum", f1, f2);
//! ```
class Math1FactoryA
    : public Factory<Func1, const shared_ptr<Func1>, const shared_ptr<Func1>>
{
public:
    /**
     * Return a pointer to the factory. On the first call, a new instance is
     * created. Since there is no need to instantiate more than one factory,
     * on all subsequent calls, a pointer to the existing factory is returned.
     */
    static Math1FactoryA* factory();

    virtual void deleteFactory();

private:
    //! Pointer to the single instance of the factory
    static Math1FactoryA* s_factory;

    //! default constructor, which is defined as private
    Math1FactoryA();

    //! Mutex for use when calling the factory
    static std::mutex domain_mutex;
};


//! Factory class to create Func1 math objects - version B
//!
//! This class is mainly used via the newMath1(const string&, const shared_ptr<Func1>,
//! double) function, for example:
//!
//! ```cpp
//!     shared_ptr<Func1> d1 = newMath1("plus-constant", f, 1.);
//! ```
class Math1FactoryB : public Factory<Func1, const shared_ptr<Func1>, double>
{
public:
    /**
     * Return a pointer to the factory. On the first call, a new instance is
     * created. Since there is no need to instantiate more than one factory,
     * on all subsequent calls, a pointer to the existing factory is returned.
     */
    static Math1FactoryB* factory();

    virtual void deleteFactory();

private:
    //! Pointer to the single instance of the factory
    static Math1FactoryB* s_factory;

    //! default constructor, which is defined as private
    Math1FactoryB();

    //! Mutex for use when calling the factory
    static std::mutex domain_mutex;
};


//! Create a new Func1 functor object.
//! @param func1Type  string identifying function type.
//! @param coeff  Coefficient; definition depends on the function type.
shared_ptr<Func1> newFunc1(const string& func1Type, double coeff=1.);

//! Create a new Func1 functor object.
//! @param func1Type  string identifying function type.
//! @param params  Parameter vector; definition depends on the function type.
//! @param n  Integer; definition depends on function type and may or may not be used.
shared_ptr<Func1> newFunc1(const string& func1Type,
                           const vector<double>& params, size_t n=1);

//! Create a new Func1 functor object based on mathematical operations.
//! @param mathType  String identifying operation.
//! @param f1  First Func1 object.
//! @param f2  Second Func1 object.
shared_ptr<Func1> newMath1(const string& mathType, const shared_ptr<Func1> f1,
                           const shared_ptr<Func1> f2);

//! Create a new Func1 functor object based on mathematical operations.
//! @param mathType  String identifying operation.
//! @param f  Func1 object.
//! @param c  Coefficient.
shared_ptr<Func1> newMath1(const string& mathType, const shared_ptr<Func1> f, double c);

}

#endif
