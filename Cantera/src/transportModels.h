#ifndef CT_TRANSPORT_MODELS_H
#define CT_TRANSPORT_MODELS_H


#include "TransportFactory.h"

namespace Cantera {

    inline Transport* MultiTransport(mixture_t& mix, string file, 
        int loglevel=0) {
        TransportFactory* f = TransportFactory::factory();
        Transport* t = f->newTransport(Multicomponent, file, mix, loglevel);
        mix->setTransport(t);
    }

}
