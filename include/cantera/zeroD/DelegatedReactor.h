//! @file DelegatedReactor.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_DELEGATEDREACTOR_H
#define CT_DELEGATEDREACTOR_H

#include "Reactor.h"
#include "cantera/base/global.h"

namespace Cantera
{

class DelegatedReactor : public Reactor
{
public:
    DelegatedReactor() {
        m_initialize = [this](double t0) { Reactor::initialize(t0); };
    }

    void setInitialize(const std::function<void(double)>& func,
                       const std::string& when) {
        setOverride(&m_initialize, func, when,
                    &DelegatedReactor::baseInitialize);
    }
    virtual void initialize(double t0) {
        m_initialize(t0);
    }
    void baseInitialize(double t0) { Reactor::initialize(t0); }

private:
    template <class ... Args>
    void setOverride(
        std::function<void(Args ... args)>* target,
        const std::function<void(Args ... args)>& func,
        const std::string& when,
        void(DelegatedReactor::* base)(Args ... args))
    {

        if (when == "before") {
            *target = [this, base, func](Args ... args) {
                func(args ...);
                (this->*base)(args ...);
            };
        } else if (when == "after") {
            *target = [this, base, func](Args ... args) {
                (this->*base)(args ...);
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

    std::function<void(double)> m_initialize;
};

}
#endif
