/**
 *  @file ReactorFactory.h
 */

/*
 * $Author: dggoodwin $
 * $Revision: 1.3 $
 * $Date: 2007/05/10 03:28:32 $
 */

// Copyright 2001  California Institute of Technology


#ifndef REACTOR_FACTORY_H
#define REACTOR_FACTORY_H

#include "ReactorBase.h"
#include "FactoryBase.h"

#if defined(THREAD_SAFE_CANTERA)
#include <boost/thread/mutex.hpp>
#endif

namespace CanteraZeroD {


    class ReactorFactory : FactoryBase {

    public:

        static ReactorFactory* factory() {
            #if defined(THREAD_SAFE_CANTERA)
               boost::mutex::scoped_lock   lock(reactor_mutex) ;
            #endif
            if (!s_factory) s_factory = new ReactorFactory;
            return s_factory;
        }

        virtual void deleteFactory() {
       #if defined(THREAD_SAFE_CANTERA)
         boost::mutex::scoped_lock   lock(reactor_mutex) ;
       #endif
        if (s_factory) {
          delete s_factory;
          s_factory = 0;
        }
    }

    /**
         * Destructor doesn't do anything. 
         */
        virtual ~ReactorFactory() {}

        /**
         * Create a new reactor.
         * @param n the type to be created.
         */
        virtual ReactorBase* newReactor(int n);
        virtual ReactorBase* newReactor(std::string reactorType);

    private:

        static ReactorFactory* s_factory;
         #if defined(THREAD_SAFE_CANTERA)
            static boost::mutex reactor_mutex ;
         #endif
        ReactorFactory(){}
    };

    inline ReactorBase* newReactor(std::string model,  
        ReactorFactory* f=0) {
        if (f == 0) {
            f = ReactorFactory::factory();
        }
        return f->newReactor(model);
    }

}

#endif



