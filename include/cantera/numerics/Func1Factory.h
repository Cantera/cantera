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
//! @since New in %Cantera 3.0
class Func1Factory : public Factory<Func1, const vector<double>&>
{
public:
    /**
     * Return a pointer to the factory. On the first call, a new instance is
     * created. Since there is no need to instantiate more than one factory,
     * on all subsequent calls, a pointer to the existing factory is returned.
     */
    static Func1Factory* factory();

    void deleteFactory() override;

private:
    //! Pointer to the single instance of the factory
    static Func1Factory* s_factory;

    //! default constructor, which is defined as private
    Func1Factory();

    //! Mutex for use when calling the factory
    static std::mutex s_mutex;
};


//! Factory class to create Func1 compound objects - version A
//!
//! This class is mainly used via the newFunc1(const string&, const shared_ptr<Func1>,
//! const shared_ptr<Func1>) function, for example:
//!
//! ```cpp
//!     shared_ptr<Func1> d1 = newFunc1("sum", f1, f2);
//! ```
//! @since New in %Cantera 3.0
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

    void deleteFactory() override;

private:
    //! Pointer to the single instance of the factory
    static Math1FactoryA* s_factory;

    //! default constructor, which is defined as private
    Math1FactoryA();

    //! Mutex for use when calling the factory
    static std::mutex s_mutex;
};


//! Factory class to create Func1 compound objects - version B
//!
//! This class is mainly used via the newFunc1(const string&, const shared_ptr<Func1>,
//! double) function, for example:
//!
//! ```cpp
//!     shared_ptr<Func1> d1 = newFunc1("plus-constant", f, 1.);
//! ```
//! @since New in %Cantera 3.0
class Math1FactoryB : public Factory<Func1, const shared_ptr<Func1>, double>
{
public:
    /**
     * Return a pointer to the factory. On the first call, a new instance is
     * created. Since there is no need to instantiate more than one factory,
     * on all subsequent calls, a pointer to the existing factory is returned.
     */
    static Math1FactoryB* factory();

    void deleteFactory() override;

private:
    //! Pointer to the single instance of the factory
    static Math1FactoryB* s_factory;

    //! default constructor, which is defined as private
    Math1FactoryB();

    //! Mutex for use when calling the factory
    static std::mutex s_mutex;
};


//! Create a new basic functor object (see @ref func1basic).
//! @param func1Type  String identifying functor type.
//! @param coeff  Coefficient; definition depends on functor type.
//! @ingroup func1basic
//! @since New in %Cantera 3.0
shared_ptr<Func1> newFunc1(const string& func1Type, double coeff=1.);

//! Create a new advanced functor object (see @ref func1advanced).
//! @param func1Type  String identifying functor type.
//! @param params  Parameter vector; definition depends on functor type.
//! @ingroup func1advanced
//! @since New in %Cantera 3.0
shared_ptr<Func1> newFunc1(const string& func1Type, const vector<double>& params);

//! Create a new compound functor object (see @ref func1compound).
//! @param func1Type  String identifying functor type.
//! @param f1  First Func1 object.
//! @param f2  Second Func1 object.
//! @ingroup func1compound
//! @since New in %Cantera 3.0
shared_ptr<Func1> newFunc1(const string& func1Type,
                           const shared_ptr<Func1> f1, const shared_ptr<Func1> f2);

//! Create a new modified functor object (see @ref func1modified).
//! @param func1Type  String identifying functor type.
//! @param f  Func1 object.
//! @param coeff  Coefficient; definition depends on functor type.
//! @ingroup func1modified
//! @since New in %Cantera 3.0
shared_ptr<Func1> newFunc1(const string& func1Type,
                           const shared_ptr<Func1> f, double coeff);

//! Check definition of functor object.
//! @param func1Type  String identifying functor type.
//! @return string indicating functor type: @c "undefined" if not defined;
//!     @c "standard" if @ref func1basic or @ref func1advanced (defined in
//!     Func1Factory); @c "compound" if @ref func1compound (defined in Math1FactoryA);
//!     or @c "modified" if @ref func1modified (defined in Math1FactoryB).
//! @internal Not intended for use in external API's (Python, MATLAB, etc.).
//! @ingroup func1helper
//! @since New in %Cantera 3.1
string checkFunc1(const string& func1Type);

}

#endif
