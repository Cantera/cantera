/**
 *  @file ThermoFactory.h
 */

/*
 * $Author: dggoodwin $
 * $Revision: 1.5 $
 * $Date: 2005/11/22 17:59:04 $
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
            if (!s_factory) s_factory = new ThermoFactory;
            return s_factory;
        }

	static void deleteFactory() {
	    if (s_factory) {
	      delete s_factory;
	      s_factory = 0;
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

        static ThermoFactory* s_factory;
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


