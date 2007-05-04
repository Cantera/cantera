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


    class UnknownKineticsModel : public CanteraError {
    public:
	UnknownKineticsModel(std::string proc, std::string kineticsModel) :
	    CanteraError(proc, "Specified Kinetics model "   
			 + kineticsModel + 
			 " does not match any known type.") {}
	virtual ~UnknownKineticsModel() {}
    };


    /**
     * Factory for kinetics managers.
     */
    class KineticsFactory {

    public:

        static KineticsFactory* factory() {
            if (!s_factory) s_factory = new KineticsFactory;
            return s_factory;
        }

        virtual ~KineticsFactory() {
            delete s_factory;
            s_factory = 0;
        }

        /**
         * Create a new kinetics manager.
         */ 
        virtual Kinetics* newKinetics(XML_Node& phase,
            std::vector<ThermoPhase*> th);

        virtual Kinetics* newKinetics(std::string model);

    private:

        static KineticsFactory* s_factory;
        KineticsFactory(){}
    };


    /**
     *  Create a new kinetics manager.
     */
    inline Kinetics* newKineticsMgr(XML_Node& phase,  
        std::vector<ThermoPhase*> th, KineticsFactory* f=0) {
        if (f == 0) {
            f = KineticsFactory::factory();
        }
        Kinetics* kin = f->newKinetics(phase, th);
        return kin;
    }

    /**
     *  Create a new kinetics manager.
     */
    inline Kinetics* newKineticsMgr(std::string model, KineticsFactory* f=0) {
        if (f == 0) {
            f = KineticsFactory::factory();
        }
        Kinetics* kin = f->newKinetics(model);
        return kin;
    }
}

#endif


