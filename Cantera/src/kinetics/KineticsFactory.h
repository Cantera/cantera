/**
 *  @file KineticsFactory.h
 */

/*
 * $Author: hkmoffa $
 * $Revision: 1.3 $
 * $Date: 2009/02/11 01:50:58 $
 */

// Copyright 2001  California Institute of Technology


#ifndef KINETICS_FACTORY_H
#define KINETICS_FACTORY_H

#include "Kinetics.h"
#include "xml.h"
#include "FactoryBase.h"

#if defined(THREAD_SAFE_CANTERA)
#include <boost/thread/mutex.hpp>
#endif

namespace Cantera {


    class UnknownKineticsModel : public CanteraError {
    public:
    UnknownKineticsModel(std::string proc, std::string kineticsModel) :
        CanteraError(proc, "Specified Kinetics model "   
             + kineticsModel + 
             " does not match any known type.") {}
    virtual ~UnknownKineticsModel() throw() {}
    };


    /**
     * Factory for kinetics managers.
     */
    class KineticsFactory : public FactoryBase {

    public:

        static KineticsFactory* factory() {
            #if defined(THREAD_SAFE_CANTERA)
               boost::mutex::scoped_lock   lock(kinetics_mutex) ;
            #endif
            if (!s_factory) s_factory = new KineticsFactory;
            return s_factory;
        }

        virtual ~KineticsFactory() {
            //delete s_factory;
            //s_factory = 0;
        }

        virtual void deleteFactory() {
             #if defined(THREAD_SAFE_CANTERA)
               boost::mutex::scoped_lock   lock(kinetics_mutex) ;
            #endif
          if ( s_factory ) {
               delete s_factory ;
               s_factory = 0 ;
           }
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
      #if defined(THREAD_SAFE_CANTERA)
        static boost::mutex kinetics_mutex ;
      #endif
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



