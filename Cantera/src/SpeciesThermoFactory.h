/**
 *  @file SpeciesThermoFactory.h
 */

/*
 * $Author$
 * $Revision$
 * $Date$
 */

// Copyright 2001  California Institute of Technology


#ifndef SPECIESTHERMO_FACTORY_H
#define SPECIESTHERMO_FACTORY_H

#include "SpeciesThermo.h"
//#include "xml.h"

namespace Cantera {

    class XML_Node;

    /**
     * Factory to build instances of classes that manage the
     * standard-state thermodynamic properties of a set of species.
     */
    class SpeciesThermoFactory {

    public:

        static SpeciesThermoFactory* factory() {
            if (!__factory) __factory = new SpeciesThermoFactory;
            return __factory;
        }

        virtual ~SpeciesThermoFactory() {
            delete __factory;
            __factory = 0;
        }

        /**
         * Create a new species property manager.
         * @param type the type to be created.
         */ 
        virtual SpeciesThermo* newSpeciesThermo(int type);
        virtual SpeciesThermo* newSpeciesThermo(XML_Node* node);
        virtual SpeciesThermo* newSpeciesThermo(vector<XML_Node*> nodes);

    private:
        static SpeciesThermoFactory* __factory;
        SpeciesThermoFactory(){}
    };


    /**
     *  Create a new species thermo manager instance.
     */
    inline SpeciesThermo* newSpeciesThermoMgr(int type, 
        SpeciesThermoFactory* f=0) {
        if (f == 0) {
            f = SpeciesThermoFactory::factory();
        }
        SpeciesThermo* sptherm = f->newSpeciesThermo(type);
        return sptherm;
    }

    /**
     *  Create a new species thermo manager instance.
     */
    inline SpeciesThermo* newSpeciesThermoMgr(XML_Node* node, 
        SpeciesThermoFactory* f=0) {
        if (f == 0) {
            f = SpeciesThermoFactory::factory();
        }
        SpeciesThermo* sptherm = f->newSpeciesThermo(node);
        return sptherm;
    }

    inline SpeciesThermo* newSpeciesThermoMgr(vector<XML_Node*> nodes, 
        SpeciesThermoFactory* f=0) {
        if (f == 0) {
            f = SpeciesThermoFactory::factory();
        }
        SpeciesThermo* sptherm = f->newSpeciesThermo(nodes);
        return sptherm;
    }

}

#endif


