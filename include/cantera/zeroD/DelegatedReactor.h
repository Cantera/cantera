//! @file DelegatedReactor.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_DELEGATEDREACTOR_H
#define CT_DELEGATEDREACTOR_H

#include "Reactor.h"
#include "cantera/base/global.h"
#include <array>

namespace Cantera
{

class DelegatedReactor : public Reactor
{
public:
    DelegatedReactor();

    void setInitialize(const std::function<void(double)>& func,
                       const std::string& when);
    virtual void initialize(double t0) {
        m_initialize(t0);
    }

    void setSyncState(const std::function<void()>& func,
                      const std::string& when);
    virtual void syncState() {
        m_syncState();
    }

    void setGetState(const std::function<void(std::array<size_t, 1>, double*)>& func,
                     const std::string& when);
    virtual void getState(double* y) {
        m_getState(y);
    }

    void setUpdateState(const std::function<void(std::array<size_t, 1>, double*)>& func,
                        const std::string& when);
    virtual void updateState(double* y) {
        m_updateState(y);
    }

    void setUpdateSurfaceState(
        const std::function<void(std::array<size_t, 1>, double*)>& func,
        const std::string& when);
    virtual void updateSurfaceState(double* y) {
        m_updateSurfaceState(y);
    }

    void setGetSurfaceInitialConditions(
        const std::function<void(std::array<size_t, 1>, double*)>& func,
        const std::string& when);
    virtual void getSurfaceInitialConditions(double* y) {
        m_getSurfaceInitialConditions(y);
    }

    void setUpdateConnected(const std::function<void(bool)>& func,
                            const std::string& when);
    virtual void updateConnected(bool updatePressure) {
        m_updateConnected(updatePressure);
    }

    void setEvalEqs(
        const std::function<void(std::array<size_t, 3>, double, double*, double*, double*)>& func,
        const std::string& when);
    virtual void evalEqs(double t, double* y, double* ydot, double* params) {
        m_evalEqs(t, y, ydot, params);
    }

    void setEvalWalls(const std::function<void(double)>& func,
                      const std::string& when);
    virtual void evalWalls(double t) {
        m_evalWalls(t);
    }

    void setEvalSurfaces(
        const std::function<int(double&, std::array<size_t, 1>, double, double*)>& func,
        const std::string& when);
    virtual double evalSurfaces(double t, double* ydot) {
        return m_evalSurfaces(t, ydot);
    }

    void setComponentName(const std::function<int(std::string&, size_t)>& func,
                          const std::string& when);
    virtual std::string componentName(size_t k) {
        return m_componentName(k);
    }

    void setComponentIndex(const std::function<int(size_t&, const std::string& nm)>& func,
                           const std::string& when);
    virtual size_t componentIndex(const std::string& nm) const {
        return m_componentIndex(nm);
    }

    void setSpeciesIndex(const std::function<int(size_t&, const std::string& nm)>& func,
                         const std::string& when);
    virtual size_t speciesIndex(const std::string& nm) const {
        return m_speciesIndex(nm);
    }


private:
    template <typename BaseFunc, class ... Args>
    void setDelegate(
        std::function<void(Args ...)>* target,
        const std::function<void(Args ...)>& func,
        const std::string& when,
        BaseFunc base);

    template <int nArrays, typename BaseFunc, class ... Args>
    void setDelegate(
        std::function<void(Args ...)>* target,
        const std::function<void(std::array<size_t, nArrays>, Args ...)>& func,
        const std::function<std::array<size_t, nArrays>()>& sizeGetter,
        const std::string& when,
        BaseFunc base);

    template <typename ReturnType, typename BaseFunc, class ... Args>
    void setDelegate(
        std::function<ReturnType(Args ...)>* target,
        const std::function<int(ReturnType&, Args ...)>& func,
        const std::string& when,
        const BaseFunc& base);

    template <int nArrays, typename ReturnType, typename BaseFunc, class ... Args>
    void setDelegate(
        std::function<ReturnType(Args ...)>* target,
        const std::function<int(ReturnType&, std::array<size_t, nArrays>, Args ...)>& func,
        const std::function<std::array<size_t, nArrays>()>& sizeGetter,
        const std::string& when,
        BaseFunc base);

    std::function<void(double)> m_initialize;
    std::function<void()> m_syncState;
    std::function<void(double*)> m_getState;
    std::function<void(double*)> m_updateState;
    std::function<void(double*)> m_updateSurfaceState;
    std::function<void(double*)> m_getSurfaceInitialConditions;
    std::function<void(bool)> m_updateConnected;
    std::function<void(double, double*, double*, double*)> m_evalEqs;
    std::function<void(double)> m_evalWalls;
    std::function<double(double, double*)> m_evalSurfaces;
    std::function<std::string(size_t)> m_componentName;
    std::function<size_t(const std::string&)> m_componentIndex;
    std::function<size_t(const std::string&)> m_speciesIndex;
};

}
#endif
