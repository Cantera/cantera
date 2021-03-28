//! @file Delegator.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_DELEGATOR_H
#define CT_DELEGATOR_H

#include "cantera/base/global.h"
#include "cantera/base/ctexceptions.h"
#include <array>

namespace Cantera
{

//! Delegate member functions of a C++ class to externally-specified functions
/*
 * This base class provides functions for setting delegates for the member
 * functions of a C++ class at runtime. The purpose of this capability is to
 * allow the class to be extended using functions defined in any programming
 * language that provides a C API for calling functions in that language.
 *
 * Delegates are specified as std::function objects that are responsible for
 * encapsulating the data specific to the target language and calling the
 * appropriate function in the target language. The std::function has a
 * modified function signature compared to the method that it is replacing or
 * augmenting:
 * - Methods with no return value and scalar arguments are treated the same
 * - Methods with a return value have that value as the first reference argument
 *   of their delegate function, and return an int. The delegate should return
 *   zero if it does not set the arguments value, and a non-zero value if it
 *   does.
 * - Methods with pointers to arrays as arguments have an additional argument
 *   introduced to indicate the length of each array argument. This argument
 *   occurs either first, or after the return value reference, and is a
 *   std::array<size_t, N> where N is the number of array arguments.
 *
 * Delegated methods can be specified to either "replace" the original class's
 * method, or to run "before" or "after" the original method, using the `when`
 * parameter of the `setDelegate` method.
 */
class Delegator
{
public:
    //! A method overridden in derived classes to set delegates for member
    //! functions with the signature `void()`.
    virtual void setDelegate(const std::string& name,
                             const std::function<void()>& func,
                             const std::string& when)
    {
        throw NotImplementedError("Delegator::setDelegate for void(double)");
    }

    //! A method overridden in derived classes to set delegates for member
    //! functions with the signature `void(bool)`
    virtual void setDelegate(const std::string& name,
                             const std::function<void(bool)>& func,
                             const std::string& when)
    {
        throw NotImplementedError("Delegator::setDelegate for void(bool)");
    }

    //! A method overridden in derived classes to set delegates for member
    //! functions with the signature `void(double)`
    virtual void setDelegate(const std::string& name,
                             const std::function<void(double)>& func,
                             const std::string& when)
    {
        throw NotImplementedError("Delegator::setDelegate for void(double)");
    }

    //! A method overridden in derived classes to set delegates for member
    //! functions with the signature `void(double*)`
    virtual void setDelegate(
        const std::string& name,
        const std::function<void(std::array<size_t, 1>, double*)>& func,
        const std::string& when)
    {
        throw NotImplementedError("Delegator::setDelegate"
            " for void(array<size_t, 1>, double*)");
    }

    //! A method overridden in derived classes to set delegates for member
    //! functions with the signature `void(double, double*)`
    virtual void setDelegate(
        const std::string& name,
        const std::function<void(std::array<size_t, 1>, double, double*)>& func,
        const std::string& when)
    {
        throw NotImplementedError("Delegator::setDelegate"
            " for void(array<size_t, 1>, double, double*)");
    }

    //! A method overridden in derived classes to set delegates for member
    //! functions with the signature `double(double, double*)`
    virtual void setDelegate(
        const std::string& name,
        const std::function<int(double&, std::array<size_t, 1>, double, double*)>& func,
        const std::string& when)
    {
        throw NotImplementedError("Delegator::setDelegate"
            " for void(double&, array<size_t, 1>, double, double*)");
    }

    //! A method overridden in derived classes to set delegates for member
    //! functions with the signature `string(size_t)`
    virtual void setDelegate(
        const std::string& name,
        const std::function<int(std::string&, size_t)>& func,
        const std::string& when)
    {
        throw NotImplementedError("Delegator::setDelegate"
            " for void(string&, size_t)");
    }

    //! A method overridden in derived classes to set delegates for member
    //! functions with the signature `size_t(string)`
    virtual void setDelegate(
        const std::string& name,
        const std::function<int(size_t&, const std::string&)>& func,
        const std::string& when)
    {
        throw NotImplementedError("Delegator::setDelegate"
            " for void(size_t&, string)");
    }

protected:
    //! Create a delegate for a function with no return value and no array
    //! arguments
    template <typename BaseFunc, class ... Args>
    std::function<void(Args ...)> makeDelegate(
        const std::function<void(Args ...)>& func,
        const std::string& when,
        BaseFunc base)
    {
        if (when == "before") {
            return [base, func](Args ... args) {
                func(args ...);
                base(args ...);
            };
        } else if (when == "after") {
            return [base, func](Args ... args) {
                base(args ...);
                func(args ...);
            };
        } else if (when == "replace") {
            return func;
        } else {
            throw CanteraError("Delegator::makeDelegate",
                "'when' must be one of 'before', 'after', or 'replace';"
                " not '{}", when);
        }
    }

    //! Create a delegate for a function with array arguments and no return
    //! value
    template <int nArrays, typename BaseFunc, class ... Args>
    std::function<void(Args ...)> makeDelegate(
        const std::function<void(std::array<size_t, nArrays>, Args ...)>& func,
        const std::function<std::array<size_t, nArrays>()>& sizeGetter,
        const std::string& when,
        BaseFunc base)
    {
        if (when == "before") {
            return [base, func, sizeGetter](Args ... args) {
                func(sizeGetter(), args ...);
                base(args ...);
            };
        } else if (when == "after") {
            return [base, func, sizeGetter](Args ... args) {
                base(args ...);
                func(sizeGetter(), args ...);
            };
        } else if (when == "replace") {
            return [func, sizeGetter](Args ... args) {
                func(sizeGetter(), args ...);
            };
        } else {
            throw CanteraError("Delegator::makeDelegate",
                "'when' must be one of 'before', 'after', or 'replace';"
                " not '{}", when);
        }
    }

    //! Create a delegate for a function with a return value and no array
    //! arguments
    template <typename ReturnType, typename BaseFunc, class ... Args>
    std::function<ReturnType(Args ...)> makeDelegate(
        const std::function<int(ReturnType&, Args ...)>& func,
        const std::string& when,
        const BaseFunc& base)
    {
        if (when == "before") {
            return [base, func](Args ... args) {
                // Call the provided delegate first. If it sets the return
                // value, return that, otherwise return the value from the
                // original method
                ReturnType ret;
                int done = func(ret, args ...);
                if (done) {
                    return ret;
                } else {
                    return base(args ...);
                }
            };
        } else if (when == "after") {
            return [base, func](Args ... args) {
                // Add the value returned by the original method and the
                // provided delegate
                ReturnType ret1 = base(args ...);
                ReturnType ret2;
                int done = func(ret2, args ...);
                if (done) {
                    return ret1 + ret2;
                } else {
                    return ret1;
                }
            };
        } else if (when == "replace") {
            return [func](Args ... args) {
                ReturnType ret;
                func(ret, args ...);
                return ret;
            };
        } else {
            throw CanteraError("Delegator::makeDelegate",
                "'when' must be one of 'before', 'after', or 'replace';"
                " not '{}", when);
        }
    }

    //! Create a delegate for a function with a return value and array arguments
    template <int nArrays, typename ReturnType, typename BaseFunc, class ... Args>
    std::function<ReturnType(Args ...)> makeDelegate(
        const std::function<int(ReturnType&, std::array<size_t, nArrays>, Args ...)>& func,
        const std::function<std::array<size_t, nArrays>()>& sizeGetter,
        const std::string& when,
        BaseFunc base)
    {
        if (when == "before") {
            return [base, func, sizeGetter](Args ... args) {
                // Call the provided delegate first. If it sets the return
                // value, return that, otherwise return the value from the
                // original method.
                ReturnType ret;
                int done = func(ret, sizeGetter(), args ...);
                if (done) {
                    return ret;
                } else {
                    return base(args ...);
                }
            };
        } else if (when == "after") {
            return [base, func, sizeGetter](Args ... args) {
                // Add the value returned by the original method and the
                // provided delegate.
                ReturnType ret1 = base(args ...);
                ReturnType ret2;
                int done = func(ret2, sizeGetter(), args ...);
                if (done) {
                    return ret1 + ret2;
                } else {
                    return ret1;
                }
            };
        } else if (when == "replace") {
            return [func, sizeGetter](Args ... args) {
                ReturnType ret;
                func(ret, sizeGetter(), args ...);
                return ret;
            };
        } else {
            throw CanteraError("Delegator::makeDelegate",
                "'when' must be one of 'before', 'after', or 'replace';"
                " not '{}", when);
        }
    }
};

}

#endif
