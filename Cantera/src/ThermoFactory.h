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


    class UnknownThermoPhaseModel : public CanteraError {
    public:
	UnknownThermoPhaseModel(string proc, string thermoModel) :
	    CanteraError(proc, "Specified ThermoPhase model "   
			 + thermoModel + 
			 " does not match any known type.") {}
	virtual ~UnknownThermoPhaseModel() {}
    };

    /**
     * Factory for thermodynamic property managers.
     */
    class ThermoFactory {

    public:

        static ThermoFactory* factory() {
            if (!__factory) __factory = new ThermoFactory;
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
        virtual ~ThermoFactory() {
        }

        /**
         * Create a new thermodynamic property manager.
         * @param type the type to be created.
         */ 
        //virtual ThermoPhase* newThermo(XML_Node& node, string id);
        virtual ThermoPhase* newThermoPhase(string model);

    private:

        static ThermoFactory* __factory;
        ThermoFactory(){}
    };


    /**
     *  Create a new thermo manager instance.
     */
//     inline ThermoPhase* newThermoMgr(XML_Node& root, string id,  
//         ThermoFactory* f=0) {
//         if (f == 0) {
//             f = ThermoFactory::factory();
//         }
//         ThermoPhase* therm = f->newThermo(root, id);
//         return therm;
//     }

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


