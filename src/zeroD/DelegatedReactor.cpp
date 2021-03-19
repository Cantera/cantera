//! @file DelegatedReactor.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/DelegatedReactor.h"

namespace Cantera
{

DelegatedReactor::DelegatedReactor()
{
    m_initialize = [this](double t0) { Reactor::initialize(t0); };
    m_syncState = [this]() { Reactor::syncState(); };
    m_getState = [this](double* y) { Reactor::getState(y); };
    m_updateState = [this](double* y) { Reactor::updateState(y); };
    m_updateSurfaceState = [this](double* y) { Reactor::updateSurfaceState(y); };
    m_getSurfaceInitialConditions = [this](double* y) {
        Reactor::getSurfaceInitialConditions(y);
    };
    m_updateConnected = [this](bool updatePressure) {
        Reactor::updateConnected(updatePressure);
    };
    m_evalEqs = [this](double t, double* y, double* ydot, double* params) {
        Reactor::evalEqs(t, y, ydot, params);
    };
    m_evalWalls = [this](double t) { Reactor::evalWalls(t); };
    m_evalSurfaces = [this](double t, double* ydot) {
        return Reactor::evalSurfaces(t, ydot);
    };
    m_componentName = [this](size_t k) { return Reactor::componentName(k); };
    m_componentIndex = [this](const std::string& nm) {
        return Reactor::componentIndex(nm);
    };
    m_speciesIndex = [this](const std::string& nm) {
        return Reactor::speciesIndex(nm);
    };
}

void DelegatedReactor::setInitialize(
    const std::function<void(double)>& func,
    const std::string& when)
{
    setDelegate(
        &m_initialize, func, when,
        [this](double t0) { Reactor::initialize(t0); }
    );
}

void DelegatedReactor::setSyncState(
    const std::function<void()>& func,
    const std::string& when)
{
    setDelegate(
        &m_syncState, func, when,
        [this]() { Reactor::syncState(); }
    );
}

void DelegatedReactor::setGetState(
    const std::function<void(std::array<size_t, 1>, double*)>& func,
    const std::string& when)
{
    setDelegate<1>(
        &m_getState, func,
        [this]() {
            return std::array<size_t, 1>{neq()};
        },
        when,
        [this](double* y) {
            Reactor::getState(y);
        }
    );
}

void DelegatedReactor::setUpdateState(
    const std::function<void(std::array<size_t, 1>, double*)>& func,
    const std::string& when)
{
    setDelegate<1>(
        &m_updateState, func,
        [this]() {
            return std::array<size_t, 1>{neq()};
        },
        when,
        [this](double* y) {
            Reactor::updateState(y);
        }
    );
}

void DelegatedReactor::setUpdateSurfaceState(
    const std::function<void(std::array<size_t, 1>, double*)>& func,
    const std::string& when)
{
    setDelegate<1>(
        &m_updateSurfaceState, func,
        [this]() {
            return std::array<size_t, 1>{m_nv_surf};
        },
        when,
        [this](double* y) {
            Reactor::updateSurfaceState(y);
        }
    );
}

void DelegatedReactor::setGetSurfaceInitialConditions(
    const std::function<void(std::array<size_t, 1>, double*)>& func,
    const std::string& when)
{
    setDelegate<1>(
        &m_getSurfaceInitialConditions, func,
        [this]() {
            return std::array<size_t, 1>{m_nv_surf};
        },
        when,
        [this](double* y) {
            Reactor::getSurfaceInitialConditions(y);
        }
    );
}

void DelegatedReactor::setUpdateConnected(
    const std::function<void(bool)>& func,
    const std::string& when)
{
    setDelegate(
        &m_updateConnected, func, when,
        [this](bool updatePressure) {
            Reactor::updateConnected(updatePressure);
        }
    );
}

void DelegatedReactor::setEvalEqs(
    const std::function<void(std::array<size_t, 3>, double, double*, double*, double*)>& func,
    const std::string& when)
{
    setDelegate<3>(
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

void DelegatedReactor::setEvalWalls(
    const std::function<void(double)>& func,
    const std::string& when)
{
    setDelegate(
        &m_evalWalls, func, when,
        [this](double t) {
            Reactor::evalWalls(t);
        }
    );
}

void DelegatedReactor::setEvalSurfaces(
    const std::function<int(double&, std::array<size_t, 1>, double, double*)>& func,
    const std::string& when)
{
    setDelegate<1>(
        &m_evalSurfaces, func,
        [this]() {
            return std::array<size_t, 1>{m_nv_surf};
        },
        when,
        [this](double t, double* ydot) {
            return Reactor::evalSurfaces(t, ydot);
        }
    );
}

void DelegatedReactor::setComponentName(
    const std::function<int(std::string&, size_t)>& func,
    const std::string& when)
{
    setDelegate(
        &m_componentName, func, when,
        [this](size_t k) {
            return Reactor::componentName(k);
        }
    );
}

void DelegatedReactor::setComponentIndex(
    const std::function<int(size_t&, const std::string& nm)>& func,
    const std::string& when)
{
    setDelegate(
        &m_componentIndex, func, when,
        [this](const std::string& nm) {
            return Reactor::componentIndex(nm);
        }
    );
}

void DelegatedReactor::setSpeciesIndex(
    const std::function<int(size_t&, const std::string& nm)>& func,
    const std::string& when)
{
    setDelegate(
        &m_speciesIndex, func, when,
        [this](const std::string& nm) {
            return Reactor::speciesIndex(nm);
        }
    );
}

template <typename BaseFunc, class ... Args>
void DelegatedReactor::setDelegate(
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
        throw CanteraError("DelegatedReactor::setDelegate",
            "'when' must be one of 'before', 'after', or 'replace';"
            " not '{}", when);
    }
}

template <int nArrays, typename BaseFunc, class ... Args>
void DelegatedReactor::setDelegate(
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
        throw CanteraError("DelegatedReactor::setDelegate",
            "'when' must be one of 'before', 'after', or 'replace';"
            " not '{}", when);
    }
}

template <typename ReturnType, typename BaseFunc, class ... Args>
void DelegatedReactor::setDelegate(
    std::function<ReturnType(Args ...)>* target,
    const std::function<int(ReturnType&, Args ...)>& func,
    const std::string& when,
    const BaseFunc& base)
{
    if (when == "before") {
        *target = [base, func](Args ... args) {
            ReturnType ret;
            int done = func(ret, args ...);
            if (done) {
                return ret;
            } else {
                return base(args ...);
            }
        };
    } else if (when == "after") {
        *target = [base, func](Args ... args) {
            ReturnType ret = base(args ...);
            func(ret, args ...);
            return ret;
        };
    } else if (when == "replace") {
        *target = [func](Args ... args) {
            ReturnType ret;
            func(ret, args ...);
            return ret;
        };
    } else {
        throw CanteraError("DelegatedReactor::setDelegate",
            "'when' must be one of 'before', 'after', or 'replace';"
            " not '{}", when);
    }
}

template <int nArrays, typename ReturnType, typename BaseFunc, class ... Args>
void DelegatedReactor::setDelegate(
    std::function<ReturnType(Args ...)>* target,
    const std::function<int(ReturnType&, std::array<size_t, nArrays>, Args ...)>& func,
    const std::function<std::array<size_t, nArrays>()>& sizeGetter,
    const std::string& when,
    BaseFunc base)
{
    if (when == "before") {
        *target = [base, func, sizeGetter](Args ... args) {
            ReturnType ret;
            int done = func(ret, sizeGetter(), args ...);
            if (done) {
                return ret;
            } else {
                return base(args ...);
            }
        };
    } else if (when == "after") {
        *target = [base, func, sizeGetter](Args ... args) {
            ReturnType ret = base(args ...);
            func(ret, sizeGetter(), args ...);
            return ret;
        };
    } else if (when == "replace") {
        *target = [func, sizeGetter](Args ... args) {
            ReturnType ret;
            func(ret, sizeGetter(), args ...);
            return ret;
        };
    } else {
        throw CanteraError("DelegatedReactor::setDelegate",
            "'when' must be one of 'before', 'after', or 'replace';"
            " not '{}", when);
    }
}


} // end namespace Cantera