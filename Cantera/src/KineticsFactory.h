/**
 *  @file KineticsFactory.h
 */

/*
 * $Author$
 * $Revision$
 * $Date$
 */

// Copyright 2001  California Institute of Technology


#ifndef KINETICS_FACTORY_H
#define KINETICS_FACTORY_H

#include "Kinetics.h"
#include "xml.h"

namespace Cantera {

    /**
     * Factory for kinetics managers.
     */
    class KineticsFactory {

    public:

        static KineticsFactory* factory() {
            if (!__factory) __factory = new KineticsFactory;
            return __factory;
        }

        virtual ~KineticsFactory() {
            delete __factory;
            __factory = 0;
        }

        /**
         * Create a new kinetics manager.
         */ 
        virtual Kinetics* newKinetics(XML_Node& phase,
            vector<ThermoPhase*> th);

        virtual Kinetics* newKinetics(string model);

    private:

        static KineticsFactory* __factory;
        KineticsFactory(){}
    };


    /**
     *  Create a new kinetics manager.
     */
    inline Kinetics* newKineticsMgr(XML_Node& phase,  
        vector<ThermoPhase*> th, KineticsFactory* f=0) {
        if (f == 0) {
            f = KineticsFactory::factory();
        }
        Kinetics* kin = f->newKinetics(phase, th);
        return kin;
    }

    /**
     *  Create a new kinetics manager.
     */
    inline Kinetics* newKineticsMgr(string model, KineticsFactory* f=0) {
        if (f == 0) {
            f = KineticsFactory::factory();
        }
        Kinetics* kin = f->newKinetics(model);
        return kin;
    }
}

#endif


