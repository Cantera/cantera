/**
 *  @file ReactorFactory.h
 */
// Copyright 2001  California Institute of Technology


#ifndef REACTOR_FACTORY_H
#define REACTOR_FACTORY_H

#include "ReactorBase.h"
#include "cantera/base/FactoryBase.h"
#include "cantera/base/ct_thread.h"

namespace Cantera
{

class ReactorFactory : Cantera::FactoryBase
{

public:

    static ReactorFactory* factory() {
        ScopedLock lock(reactor_mutex);
        if (!s_factory) {
            s_factory = new ReactorFactory;
        }
        return s_factory;
    }

    virtual void deleteFactory() {
        ScopedLock lock(reactor_mutex);
        if (s_factory) {
            delete s_factory;
            s_factory = 0;
        }
    }

    /**
     * Create a new reactor.
     * @param n the type to be created.
     */
    virtual ReactorBase* newReactor(int n);
    virtual ReactorBase* newReactor(const std::string& reactorType);

private:
    static ReactorFactory* s_factory;
    static mutex_t reactor_mutex;
    ReactorFactory() {}
};

inline ReactorBase* newReactor(const std::string& model,
                               ReactorFactory* f=0)
{
    if (f == 0) {
        f = ReactorFactory::factory();
    }
    return f->newReactor(model);
}

}

#endif



