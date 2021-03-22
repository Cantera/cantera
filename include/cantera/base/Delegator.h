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

class Delegator
{
public:
    virtual void setDelegate(const std::string& name,
                             const std::function<void()>& func,
                             const std::string& when)
    {
        throw NotImplementedError("Delegator::setDelegate for void(double)");
    }

    virtual void setDelegate(const std::string& name,
                             const std::function<void(bool)>& func,
                             const std::string& when)
    {
        throw NotImplementedError("Delegator::setDelegate for void(bool)");
    }

    virtual void setDelegate(const std::string& name,
                             const std::function<void(double)>& func,
                             const std::string& when)
    {
        throw NotImplementedError("Delegator::setDelegate for void(double)");
    }

    // For functions with the signature void(double*)
    virtual void setDelegate(
        const std::string& name,
        const std::function<void(std::array<size_t, 1>, double*)>& func,
        const std::string& when)
    {
        throw NotImplementedError("Delegator::setDelegate"
            " for void(array<size_t, 1>, double*)");
    }

    // For functions with the signature void(double, double*, double*, double*)
    virtual void setDelegate(
        const std::string& name,
        const std::function<void(std::array<size_t, 3>, double, double*, double*, double*)>& func,
        const std::string& when)
    {
        throw NotImplementedError("Delegator::setDelegate"
            " for void(array<size_t, 3>, double, double*, double*, double*)");
    }

    // For functions with the signature double(double, double*)
    virtual void setDelegate(
        const std::string& name,
        const std::function<int(double&, std::array<size_t, 1>, double, double*)>& func,
        const std::string& when)
    {
        throw NotImplementedError("Delegator::setDelegate"
            " for void(double&, array<size_t, 1>, double, double*)");
    }

    // For functions with the signature string(size_t)
    virtual void setDelegate(
        const std::string& name,
        const std::function<int(std::string&, size_t)>& func,
        const std::string& when)
    {
        throw NotImplementedError("Delegator::setDelegate"
            " for void(string&, size_t)");
    }

    // For functions with the signature size_t(string)
    virtual void setDelegate(
        const std::string& name,
        const std::function<int(size_t&, const std::string&)>& func,
        const std::string& when)
    {
        throw NotImplementedError("Delegator::setDelegate"
            " for void(size_t&, string)");
    }

protected:
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

    template <typename ReturnType, typename BaseFunc, class ... Args>
    std::function<ReturnType(Args ...)> makeDelegate(
        const std::function<int(ReturnType&, Args ...)>& func,
        const std::string& when,
        const BaseFunc& base)
    {
        if (when == "before") {
            return [base, func](Args ... args) {
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

    template <int nArrays, typename ReturnType, typename BaseFunc, class ... Args>
    std::function<ReturnType(Args ...)> makeDelegate(
        const std::function<int(ReturnType&, std::array<size_t, nArrays>, Args ...)>& func,
        const std::function<std::array<size_t, nArrays>()>& sizeGetter,
        const std::string& when,
        BaseFunc base)
    {
        if (when == "before") {
            return [base, func, sizeGetter](Args ... args) {
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
