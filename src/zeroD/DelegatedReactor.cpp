//! @file DelegatedReactor.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/DelegatedReactor.h"

namespace Cantera
{

DelegatedReactor::DelegatedReactor()
{
    m_initialize = [this](double t0) { Reactor::initialize(t0); };
    m_evalEqs = [this](double t, double* y, double* ydot, double* params) {
        Reactor::evalEqs(t, y, ydot, params);
    };
}

void DelegatedReactor::setInitialize(
    const std::function<void(double)>& func,
    const std::string& when)
{
    setOverride(
        &m_initialize, func, when,
        [this](double t0) { Reactor::initialize(t0); }
    );
}

void DelegatedReactor::setEvalEqs(
    const std::function<void(std::array<size_t, 3>, double, double*, double*, double*)>& func,
    const std::string& when)
{
    setOverride<3>(
        &m_evalEqs, func,
        [this]() {
            return std::array<size_t, 3>{neq(), neq(), nSensParams()};
        },
        when,
        [this](double t, double* y, double* ydot, double* params) {
            Reactor::evalEqs(t, y, ydot, params);
        }
    );
}

template <typename BaseFunc, class ... Args>
void DelegatedReactor::setOverride(
    std::function<void(Args ...)>* target,
    const std::function<void(Args ...)>& func,
    const std::string& when,
    BaseFunc base)
{
    if (when == "before") {
        *target = [base, func](Args ... args) {
            func(args ...);
            base(args ...);
        };
    } else if (when == "after") {
        *target = [base, func](Args ... args) {
            base(args ...);
            func(args ...);
        };
    } else if (when == "replace") {
        *target = func;
    } else {
        throw CanteraError("DelegatedReactor::setOverride",
            "'when' must be one of 'before', 'after', or 'replace';"
            " not '{}", when);
    }
}

template <int nArrays, typename BaseFunc, class ... Args>
void DelegatedReactor::setOverride(
    std::function<void(Args ...)>* target,
    const std::function<void(std::array<size_t, nArrays>, Args ...)>& func,
    const std::function<std::array<size_t, nArrays>()>& sizeGetter,
    const std::string& when,
    BaseFunc base)
{
    if (when == "before") {
        *target = [base, func, sizeGetter](Args ... args) {
            func(sizeGetter(), args ...);
            base(args ...);
        };
    } else if (when == "after") {
        *target = [base, func, sizeGetter](Args ... args) {
            base(args ...);
            func(sizeGetter(), args ...);
        };
    } else if (when == "replace") {
        *target = [func, sizeGetter](Args ... args) {
            func(sizeGetter(), args ...);
        };
    } else {
        throw CanteraError("DelegatedReactor::setOverride",
            "'when' must be one of 'before', 'after', or 'replace';"
            " not '{}", when);
    }
}


} // end namespace Cantera