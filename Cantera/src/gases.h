/**
 *  @file gases.h
 *
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_GASES_H
#define CT_GASES_H

#include "IdealGasMix.h"
#include "TransportFactory.h"

namespace Cantera {

    template<class T>
    class GasWithTransport : public IdealGasMix, public T {
    public:
        GasWithTransport(map<string, string>& params, 
            TransportFactory* f = 0) {
            if (f == 0) 
                f = TransportFactory::factory();
            f->initTransport(*self, dbase);
        }
    };
}

#endif
