/**
 *  @file SpeciesThermoFactory.h
 */

/*
 * $Author: dggoodwin $
 * $Revision: 1.7 $
 * $Date: 2006/05/03 19:46:40 $
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
	    CanteraError(proc, "species " + spName + 
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

        /**
         * This class is implemented as a singleton -- one in which
         * only one instance is needed.  The recommended way to access
         * the factory is to call this static method, which
         * instantiates the class if it is the first call, but
         * otherwise simply returns the pointer to the existing
         * instance.
         */
        static SpeciesThermoFactory* factory() {
            if (!s_factory) s_factory = new SpeciesThermoFactory;
            return s_factory;
        }

        /**
         * If it is necessary to explicitly delete the factory before
         * the process terminates (for example, when checking for
         * memory leaks) then this method can be called to delete it.
         */
	static void deleteFactory() {
	    if (s_factory) {
	      delete s_factory;
	      s_factory = 0;
	    }
	}
	

	/**
         * Destructor. Doesn't do anything. We do not delete statically
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


        virtual void installThermoForSpecies(int k, const XML_Node& s, 
            SpeciesThermo& spthermo);

    private:

        /// pointer to the sole instance of this class
        static SpeciesThermoFactory* s_factory;

        /// Constructor. This is made private, so that only the static
        /// method factory() can instantiate the class.
        SpeciesThermoFactory(){}
    };


    ////////////////////// Convenience functions ////////////////////
    //
    //  These functions allow using a different factory class that
    //  derives from SpeciesThermoFactory.
    //
    //////////////////////////////////////////////////////////////////


    /**
     *  Create a new species thermo manager instance, by specifying
     * the type and (optionally) a pointer to the factory to use to
     * create it.
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


