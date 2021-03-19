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

    void setEvalEqs(const std::function<void(std::array<size_t, 3>, double, double*, double*, double*)>& func,
                    const std::string& when);
    virtual void evalEqs(double t, double* y, double* ydot, double* params) {
        m_evalEqs(t, y, ydot, params);
    }

    void setComponentName(const std::function<int(std::string&, size_t)>& func,
                          const std::string& when);
    virtual std::string componentName(size_t k) {
        return m_componentName(k);
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

    std::function<void(double)> m_initialize;
    std::function<void(double, double*, double*, double*)> m_evalEqs;
    std::function<std::string(size_t)> m_componentName;
};

}
#endif
