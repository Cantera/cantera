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
#include "ctexceptions.h"

namespace Cantera {

    class XML_Node;

    /**
     * Throw a named error for an unknown or missing species thermo
     * model. 
     */
    class UnknownSpeciesThermoModel: public CanteraError {
    public:
	UnknownSpeciesThermoModel(string proc, string spName,
				  string speciesThermoModel) :
	    CanteraError(proc, "species :" + spName + 
			 ": Specified speciesThermoPhase model "   
			 + speciesThermoModel + 
			 " does not match any known type.") {}
	virtual ~UnknownSpeciesThermoModel() {}
    };

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
	static void deleteFactory() {
	    if (__factory) {
	      delete __factory;
	      __factory = 0;
	    }
	}
	
	/**
         * Destructor doesn't do anything. We do not delete statically
	 * created single instance of this class here, because it would
	 * create an infinite loop if destructor is called for that
	 * single instance.
         */
        virtual ~SpeciesThermoFactory() {
        }

        /**
         * Create a new species property manager.
         * @param type the type to be created.
         */ 
        virtual SpeciesThermo* newSpeciesThermo(int type);
        virtual SpeciesThermo* newSpeciesThermo(XML_Node* node);
        virtual SpeciesThermo* newSpeciesThermo(vector<XML_Node*> nodes);
        virtual SpeciesThermo* newSpeciesThermoOpt(vector<XML_Node*> nodes);

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
        SpeciesThermoFactory* f=0, bool opt=false) {
        if (f == 0) {
            f = SpeciesThermoFactory::factory();
        }
	SpeciesThermo* sptherm;
	if (opt) {
	  sptherm = f->newSpeciesThermoOpt(nodes);
	} else { 
	  sptherm = f->newSpeciesThermo(nodes);
	}
        return sptherm;
    }

}

#endif


