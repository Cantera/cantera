//! @file Delegator.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_DELEGATOR_H
#define CT_DELEGATOR_H

#include "cantera/base/global.h"
#include "cantera/base/Units.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/ExtensionManager.h"
#include <array>
#include <list>

namespace Cantera
{

//! Delegate member functions of a C++ class to externally-specified functions
/*!
 * This base class provides functions for setting delegates for the member
 * functions of a C++ class at runtime. The purpose of this capability is to
 * allow the class to be extended using functions defined in any programming
 * language that provides a C API for calling functions in that language.
 *
 * Delegates are specified as `std::function` objects that are responsible for
 * encapsulating the data specific to the target language and calling the
 * appropriate function in the target language. The `std::function` has a
 * modified function signature compared to the method that it is replacing or
 * augmenting:
 * - Methods with no return value and scalar arguments are treated the same
 * - Methods with a return value have that value as the first reference argument
 *   of their delegate function, and return an `int`. The delegate should return
 *   zero if it does not set the arguments value, and a non-zero value if it
 *   does.
 * - Methods with pointers to arrays as arguments have an additional argument
 *   introduced to indicate the length of each array argument. This argument
 *   occurs either first, or after the return value reference, and is a
 *   `std::array<size_t, N>` where N is the number of array arguments.
 *
 * Delegated methods can be specified to either "replace" the original class's
 * method, or to run "before" or "after" the original method, using the `when`
 * parameter of the `setDelegate` method. There are two special cases for delegates of
 * methods with return values:
 * - If the "before" delegate specifies a value for the return parameter (and then
 *   returns with a non-zero status value), this value will be returned and the original
 *   method will not be called. Otherwise, the original method will be called and its
 *   value will be returned.
 * - If an "after" delegate specifies a return value (and returns with a non-zero
 *   status value), this value will be added to the value returned by the original
 *   method, and this combined result will be then be returned. The meaning of "added"
 *   is determined by the `+` operator for the return type, for example addition for
 *   numeric types or concatenation for strings.
 *
 * ## Implementation for each delegated function type
 *
 * Several functions and member variables are defined in this base class for each
 * distinct function signature that can be delegated by a derived class such as
 * ReactorDelegator. These are:
 * - The `install` function stores the address of a function that will be called by
 *   a derived delegator class in the corresponding `m_funcs_...` map, and sets the
 *   default implementation of this method (which should be to call the base class
 *   implementation). For delegates with a return type, a copy of this default
 *   implementation is also stored in the corresponding `m_base_...` map
 * - The `setDelegate` function wraps the user-provided delegate in a function that
 *   handles whether the delegate is to be called before, after, or instead of the base
 *   class's implementation. This function is then stored at the address specified in
 *   the corresponding `m_funcs_...` map.
 * - `m_funcs_...` is a mapping between member function names and the addresses
 *   of the functions that will be called to implement these functions in the derived
 *   delegator class.
 * - `m_base_...` is a mapping between member function names and the default
 *   implementations of these member functions, for delegates with return values
 *
 * Additional implementation for each function type is specific to the programming
 * language that the delegate is written in. For the Python delegates, see additional
 * documentation in `delegator.pyx`.
 *
 * ## Implementation for specific delegated functions
 *
 * Beyond the implementation of particular function signatures, there are no elements
 * of the Delegator class that are specific to individual delegated functions, which
 * are handled by derived classes such as ReactorDelegator, which will also inherit from
 * another base class such as Reactor.
 *
 * Delegation of a member function (for example, `Reactor::eval`) is handled by several
 * elements:
 * - A `std::function` member variable to hold a delegate function, or its default
 *   implementation. This function type includes the additional `std::array` argument
 *   of array sizes for functions with array arguments.
 * - An override of the member function whose implementation calls the stored delegate.
 *   For delegates that need an array of array sizes, this function first calculates the
 *   necessary values and passes them as an additional argument to the delegate.
 * - A call to `install` from the constructor of the derived delegator class, which
 *   takes the member function name, a reference to the `std::function` member variable
 *   described above, and a lambda function that implements the default behavior, that
 *   is, calling the equivalent base class method.
 *
 * Additional implementation for each function is specific to the programming language
 * that the delegate is written in. For Python delegates, see additional documentation
 * in `delegator.pyx`.
 */
class Delegator
{
public:
    //! Get the name of the user-defined class in the extension language
    string delegatorName() const {
        return m_delegatorName;
    }

    //! Set the name of the user-defined class in the extension language
    void setDelegatorName(const string& delegatorName) {
        m_delegatorName = delegatorName;
    }

    //! Set delegates for member functions with the signature `void()`.
    void setDelegate(const string& name, const function<void()>& func,
                     const string& when)
    {
        if (!m_funcs_v.count(name)) {
            throw NotImplementedError("Delegator::setDelegate",
                "for function named '{}' with signature 'void()'.", name);
        }
        *m_funcs_v[name] = makeDelegate(func, when, *m_funcs_v[name]);
    }

    //! set delegates for member functions with the signature `void(bool)`
    void setDelegate(const string& name, const function<void(bool)>& func,
                     const string& when)
    {
        if (!m_funcs_v_b.count(name)) {
            throw NotImplementedError("Delegator::setDelegate",
                "for function named '{}' with signature 'void(bool)'.", name);
        }
        *m_funcs_v_b[name] = makeDelegate(func, when, *m_funcs_v_b[name]);
    }

    //! set delegates for member functions with the signature `void(double)`
    void setDelegate(const string& name, const function<void(double)>& func,
                     const string& when)
    {
        if (!m_funcs_v_d.count(name)) {
            throw NotImplementedError("Delegator::setDelegate",
                "for function named '{}' with signature 'void(double)'.", name);
        }
        *m_funcs_v_d[name] = makeDelegate(func, when, *m_funcs_v_d[name]);
    }

    //! set delegates for member functions with the signature `void(AnyMap&)`
    void setDelegate(const string& name, const function<void(AnyMap&)>& func,
                     const string& when)
    {
        if (!m_funcs_v_AMr.count(name)) {
            throw NotImplementedError("Delegator::setDelegate",
                "for function named '{}' with signature 'void(AnyMap&)'.", name);
        }
        *m_funcs_v_AMr[name] = makeDelegate(func, when, *m_funcs_v_AMr[name]);
    }

    //! set delegates for member functions with the signature
    //! `void(AnyMap&, UnitStack&)`
    void setDelegate(const string& name,
                     const function<void(const AnyMap&, const UnitStack&)>& func,
                     const string& when)
    {
        if (!m_funcs_v_cAMr_cUSr.count(name)) {
            throw NotImplementedError("Delegator::setDelegate",
                "for function named '{}' with signature "
                "'void(const AnyMap&, const UnitStack&)'.",
                name);
        }
        *m_funcs_v_cAMr_cUSr[name] = makeDelegate(func, when, *m_funcs_v_cAMr_cUSr[name]);
    }

    //! set delegates for member functions with the signature
    //! `void(const string&, void*)`
    void setDelegate(const string& name,
                     const function<void(const string&, void*)>& func,
                     const string& when)
    {
        if (!m_funcs_v_csr_vp.count(name)) {
            throw NotImplementedError("Delegator::setDelegate",
                "for function named '{}' with signature 'void(const string&, void*)'.");
        }
        *m_funcs_v_csr_vp[name] = makeDelegate(func, when, *m_funcs_v_csr_vp[name]);
    }

    //! Set delegates for member functions with the signature `void(double*)`
    void setDelegate(const string& name,
                     const function<void(std::array<size_t, 1>, double*)>& func,
                     const string& when)
    {
        if (!m_funcs_v_dp.count(name)) {
            throw NotImplementedError("Delegator::setDelegate",
                "for function named '{}' with signature 'void(double*)'.", name);
        }
        *m_funcs_v_dp[name] = makeDelegate(func, when, *m_funcs_v_dp[name]);
    }

    //! Set delegates for member functions with the signature `void(double, double*)`
    void setDelegate(
        const string& name,
        const function<void(std::array<size_t, 1>, double, double*)>& func,
        const string& when)
    {
        if (!m_funcs_v_d_dp.count(name)) {
            throw NotImplementedError("Delegator::setDelegate",
                "for function named '{}' with signature 'void(double, double*)'.",
                name);
        }
        *m_funcs_v_d_dp[name] = makeDelegate(func, when, *m_funcs_v_d_dp[name]);
    }

    //! Set delegates for member functions with the signature
    //! `void(double, double*, double*)`
    void setDelegate(
        const string& name,
        const function<void(std::array <size_t, 2>, double, double*, double*)>& func,
        const string& when)
    {
        if (!m_funcs_v_d_dp_dp.count(name)) {
            throw NotImplementedError("Delegator::setDelegate",
                "for function named '{}' with signature "
                "'void(double, double*, double*)'.", name);
        }
        *m_funcs_v_d_dp_dp[name] = makeDelegate(func, when, *m_funcs_v_d_dp_dp[name]);
    }

    //! Set delegates for member functions with the signature
    //! `void(double*, double*, double*)`
    void setDelegate(
        const string& name,
        const function<void(std::array<size_t, 3>, double*, double*, double*)>& func,
        const string& when)
    {
        if (!m_funcs_v_dp_dp_dp.count(name)) {
            throw NotImplementedError("Delegator::setDelegate",
                "for function named '{}' with signature "
                "'void(double*, double*, double*)'.", name);
        }
        *m_funcs_v_dp_dp_dp[name] = makeDelegate(func, when, *m_funcs_v_dp_dp_dp[name]);
    }

    //! set delegates for member functions with the signature `double(void*)`
    void setDelegate(const string& name,
                     const function<int(double&, void*)>& func,
                     const string& when)
    {
        if (!m_funcs_d_vp.count(name)) {
            throw NotImplementedError("Delegator::setDelegate",
                "for function named '{}' with signature 'double(void*)'.", name);
        }
        *m_funcs_d_vp[name] = makeDelegate(name, func, when, m_base_d_vp[name]);
    }

    //! Set delegates for member functions with the signature `string(size_t)`
    void setDelegate(const string& name,
                     const function<int(string&, size_t)>& func,
                     const string& when)
    {
        if (!m_funcs_s_sz.count(name)) {
            throw NotImplementedError("Delegator::setDelegate",
                "for function named '{}' with signature "
                "'string(size_t)'.", name);
        }
        *m_funcs_s_sz[name] = makeDelegate(name, func, when, m_base_s_sz[name]);
    }

    //! Set delegates for member functions with the signature `size_t(string)`
    void setDelegate(const string& name,
                     const function<int(size_t&, const string&)>& func,
                     const string& when)
    {
        if (!m_funcs_sz_csr.count(name)) {
            throw NotImplementedError("Delegator::setDelegate",
                "for function '{}' with signature "
                "'size_t(const string&)'.", name);
        }
        *m_funcs_sz_csr[name] = makeDelegate(name, func, when, m_base_sz_csr[name]);
    }

    //! Store a handle to a wrapper for the delegate from an external language interface
    void holdExternalHandle(const string& name,
                            const shared_ptr<ExternalHandle>& handle) {
        m_handles[name] = handle;
    }

    //! Get the handle for a wrapper for the delegate from the external language
    //! interface specified by *name*.
    //! Returns a null pointer if the requested handle does not exist.
    shared_ptr<ExternalHandle> getExternalHandle(const string& name) const {
        if (m_handles.count(name)) {
            return m_handles.at(name);
        } else {
            return shared_ptr<ExternalHandle>();
        }
    }

protected:
    //! Install a function with the signature `void()` as being delegatable
    void install(const string& name, function<void()>& target,
                 const function<void()>& func)
    {
        target = func;
        m_funcs_v[name] = &target;
    }

    //! Install a function with the signature `void(bool)` as being delegatable
    void install(const string& name, function<void(bool)>& target,
                 const function<void(bool)>& func)
    {
        target = func;
        m_funcs_v_b[name] = &target;
    }

    //! Install a function with the signature `void(double)` as being delegatable
    void install(const string& name, function<void(double)>& target,
                 const function<void(double)>& func)
    {
        target = func;
        m_funcs_v_d[name] = &target;
    }

    //! Install a function with the signature `void(AnyMap&)` as being delegatable
    void install(const string& name, function<void(AnyMap&)>& target,
                 const function<void(AnyMap&)>& func)
    {
        target = func;
        m_funcs_v_AMr[name] = &target;
    }

    //! Install a function with the signature `void(const AnyMap&, const UnitStack&)`
    //! as being delegatable
    void install(const string& name,
                 function<void(const AnyMap&, const UnitStack&)>& target,
                 const function<void(const AnyMap&, const UnitStack&)>& func)
    {
        target = func;
        m_funcs_v_cAMr_cUSr[name] = &target;
    }

    //! Install a function with the signature `void(const string&, void*) as being
    //! delegatable
    void install(const string& name, function<void(const string&, void*)>& target,
                 const function<void(const string&, void*)>& func)
    {
        target = func;
        m_funcs_v_csr_vp[name] = &target;
    }

    //! Install a function with the signature `void(double*)` as being delegatable
    void install(const string& name,
                 function<void(std::array<size_t, 1>, double*)>& target,
                 const function<void(std::array<size_t, 1>, double*)>& func)
    {
        target = func;
        m_funcs_v_dp[name] = &target;
    }

    //! Install a function with the signature `void(double, double*)` as being delegatable
    void install(const string& name,
                 function<void(std::array<size_t, 1>, double, double*)>& target,
                 const function<void(std::array<size_t, 1>, double, double*)>& func)
    {
        target = func;
        m_funcs_v_d_dp[name] = &target;
    }

    //! Install a function with the signature `void(double, double*, double*)` as being
    //! delegatable
    void install(const string& name,
                 function<void(std::array<size_t, 2>, double, double*, double*)>& target,
                 const function<void(std::array<size_t, 2>, double, double*, double*)>& func)
    {
        target = func;
        m_funcs_v_d_dp_dp[name] = &target;
    }

    //! Install a function with the signature
    //! `void(double*, double*, double*)` as being delegatable
    void install(const string& name,
                 function<void(std::array<size_t, 3>, double*, double*, double*)>& target,
                 const function<void(std::array<size_t, 3>, double*, double*, double*)>& base)
    {
        target = base;
        m_funcs_v_dp_dp_dp[name] = &target;
    }

    //! Install a function with the signature `double(void*)` as being delegatable
    void install(const string& name, function<double(void*)>& target,
                 const function<double(void*)>& func)
    {
        target = func;
        m_funcs_d_vp[name] = &target;
    }

    //! Install a function with the signature `string(size_t)` as being delegatable
    void install(const string& name,
                 function<string(size_t)>& target,
                 const function<string(size_t)>& base)
    {
        target = base;
        m_funcs_s_sz[name] = &target;
        m_base_s_sz[name] = base;
    }

    //! Install a function with the signature `size_t(string)` as being delegatable
    void install(const string& name,
                 function<size_t(const string&)>& target,
                 const function<size_t(const string&)>& base)
    {
        target = base;
        m_funcs_sz_csr[name] = &target;
        m_base_sz_csr[name] = base;
    }

    //! Create a delegate for a function with no return value
    template <typename BaseFunc, class ... Args>
    function<void(Args ...)> makeDelegate(
        const function<void(Args ...)>& func,
        const string& when,
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
            return [func](Args ... args) {
                func(args ...);
            };
        } else {
            throw CanteraError("Delegator::makeDelegate",
                "'when' must be one of 'before', 'after', or 'replace';"
                " not '{}", when);
        }
    }

    //! Create a delegate for a function with a return value
    template <typename ReturnType, class ... Args>
    function<ReturnType(Args ...)> makeDelegate(
        const string& name,
        const function<int(ReturnType&, Args ...)>& func,
        const string& when,
        const function<ReturnType(Args ...)>& base)
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
            return [base, name, func, this](Args ... args) {
                ReturnType ret;
                int has_ret = func(ret, args ...);
                if (!has_ret) {
                    throw CanteraError("Lambda generated by Delegator::makeDelegate",
                        "Method '{}' of class '{}' did not return a value of type '{}'.",
                        name, delegatorName(), demangle(typeid(ret)));
                }
                return ret;
            };
        } else {
            throw CanteraError("Delegator::makeDelegate",
                "For function named '{}':\n"
                "'when' must be one of 'before', 'after', or 'replace';"
                " not '{}'", name, when);
        }
    }

    //! @name Containers for delegates and base class implementations.
    //!
    //! Maps named `m_funcs_*` hold the current delegate. For functions with a return
    //! value, `m_base_*` holds a pointer to the base class's implementation of the
    //! function.
    //!
    //! These containers use a naming scheme based on the signature of the corresponding
    //! member functions. Following the prefix `_m_funcs_` is a notation of the return
    //! type, followed by an underscore, then the notations for each argument, separated
    //! by underscores. The following shorthand is used for different return / argument
    //! types:
    //!
    //! - `v` for `void`
    //! - `b` for `bool`
    //! - `d` for `double`
    //! - `s` for `string`
    //! - `sz` for `size_t`
    //! - `AM` for `AnyMap`
    //! - `US` for `UnitStack`
    //! - prefix `c` for `const` arguments
    //! - suffix `r` for reference arguments
    //! - suffix `p` for pointer arguments
    //! @{

    // Delegates with no return value
    map<string, function<void()>*> m_funcs_v;
    map<string, function<void(bool)>*> m_funcs_v_b;
    map<string, function<void(double)>*> m_funcs_v_d;
    map<string, function<void(AnyMap&)>*> m_funcs_v_AMr;
    map<string, function<void(const AnyMap&, const UnitStack&)>*> m_funcs_v_cAMr_cUSr;
    map<string, function<void(const string&, void*)>*> m_funcs_v_csr_vp;
    map<string, function<void(std::array<size_t, 1>, double*)>*> m_funcs_v_dp;
    map<string, function<void(std::array<size_t, 1>, double, double*)>*> m_funcs_v_d_dp;
    map<string,
        function<void(std::array<size_t, 2>, double, double*, double*)>*> m_funcs_v_d_dp_dp;
    map<string,
        function<void(std::array<size_t, 3>, double*, double*, double*)>*> m_funcs_v_dp_dp_dp;

    // Delegates with a return value
    map<string, function<double(void*)>> m_base_d_vp;
    map<string, function<double(void*)>*> m_funcs_d_vp;

    map<string, function<string(size_t)>> m_base_s_sz;
    map<string, function<string(size_t)>*> m_funcs_s_sz;

    map<string, function<size_t(const string&)>> m_base_sz_csr;
    map<string, function<size_t(const string&)>*> m_funcs_sz_csr;
    //! @}

    //! Handles to wrappers for the delegated object in external language interfaces.
    //! Used to provide access to these wrappers as well as managing cleanup functions
    //! called from the destructor of the derived ExternalHandle class.
    map<string, shared_ptr<ExternalHandle>> m_handles;

    //! Name of the class in the extension language
    string m_delegatorName;
};

}

#endif
