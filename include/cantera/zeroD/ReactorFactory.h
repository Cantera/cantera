//! @file ReactorFactory.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef REACTOR_FACTORY_H
#define REACTOR_FACTORY_H

#include "ReactorBase.h"
#include "cantera/base/FactoryBase.h"

namespace Cantera
{

class ReactorFactory : public Factory<ReactorBase>
{
public:
    static ReactorFactory* factory() {
        std::unique_lock<std::mutex> lock(reactor_mutex);
        if (!s_factory) {
            s_factory = new ReactorFactory;
        }
        return s_factory;
    }

    virtual void deleteFactory() {
        std::unique_lock<std::mutex> lock(reactor_mutex);
        delete s_factory;
        s_factory = 0;
    }

    /**
     * Create a new reactor.
     * @param n the type to be created.
     */
    virtual ReactorBase* newReactor(int n);
    virtual ReactorBase* newReactor(const std::string& reactorType);

private:
    static ReactorFactory* s_factory;
    static std::mutex reactor_mutex;
    ReactorFactory();
};

//! Create a Reactor object of the specified type
inline ReactorBase* newReactor(const std::string& model)
{
    return ReactorFactory::factory()->newReactor(model);
}

}

#endif
