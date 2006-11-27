/**
 *  @file ReactorFactory.h
 */

/*
 * $Author$
 * $Revision$
 * $Date$
 */

// Copyright 2001  California Institute of Technology


#ifndef REACTOR_FACTORY_H
#define REACTOR_FACTORY_H

#include "ReactorBase.h"

namespace CanteraZeroD {


    class ReactorFactory {

    public:

        static ReactorFactory* factory() {
            if (!s_factory) s_factory = new ReactorFactory;
            return s_factory;
        }

	static void deleteFactory() {
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


