/**
 *  @file ThermoFactory.h
 */

/*
 * $Author$
 * $Revision$
 * $Date$
 */

// Copyright 2001  California Institute of Technology


#ifndef THERMO_FACTORY_H
#define THERMO_FACTORY_H

#include "ThermoPhase.h"
#include "xml.h"

namespace Cantera {

    /**
     * Factory for thermodynamic property managers.
     */
    class ThermoFactory {

    public:

        static ThermoFactory* factory() {
            if (!__factory) __factory = new ThermoFactory;
            return __factory;
        }

        virtual ~ThermoFactory() {
            delete __factory;
            __factory = 0;
        }

        /**
         * Create a new thermodynamic property manager.
         * @param type the type to be created.
         */ 
        virtual ThermoPhase* newThermo(XML_Node& node, string id);
        virtual ThermoPhase* newThermoPhase(string model);

    private:

        static ThermoFactory* __factory;
        ThermoFactory(){}
    };


    /**
     *  Create a new thermo manager instance.
     */
    inline ThermoPhase* newThermoMgr(XML_Node& root, string id,  
        ThermoFactory* f=0) {
        if (f == 0) {
            f = ThermoFactory::factory();
        }
        ThermoPhase* therm = f->newThermo(root, id);
        return therm;
    }

    /**
     *  Create a new thermo manager instance.
     */
    inline ThermoPhase* newThermoPhase(string model,  
        ThermoFactory* f=0) {
        if (f == 0) {
            f = ThermoFactory::factory();
        }
        return f->newThermoPhase(model);
    }

}

#endif


