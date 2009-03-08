/**
 *
 *  @file TransportFactory.h
 *
 *  Header file defining class TransportFactory
 */

/*
 *  $Author: dggoodwin $
 *  $Date: 2005/11/22 18:50:19 $
 *  $Revision: 1.6 $
 *
 *  Copyright 2001 California Institute of Technology
 *
 */

#ifndef CT_TRANFACTORY_H
#define CT_TRANFACTORY_H


// turn off warnings under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

// STL includes
#include <vector>
#include <string>
#include <iostream>
#include <new>

using namespace std;

// Cantera includes
#include "../ct_defs.h"
#include "TransportBase.h"

namespace Cantera {

    /**
     * Struct to hold data read from a transport property database file.
     */
    struct GasTransportData {
        GasTransportData() : speciesName("-"), 
			     geometry(-1), wellDepth(-1.0),
			     diameter(-1.0), 
			     dipoleMoment(-1.0), 
			     polarizability(-1.0),
			     rotRelaxNumber(-1.0) {}
        
        string speciesName;
        int geometry;
        doublereal wellDepth;
        doublereal diameter;
        doublereal dipoleMoment;
        doublereal polarizability;
        doublereal rotRelaxNumber;
    };

    // forward references
    class MMCollisionInt;
    class TransportParams;
    class XML_Node;

    /**
     * The purpose of TransportFactory is to create new instances of
     * 'transport managers', which are classes that provide transport
     * properties and are derived from base class
     * Transport. TransportFactory handles all initialization
     * required, including evaluation of collision integrals and
     * generating polynomial fits.  Transport managers can also be
     * created in other ways. @ingroup transportgroup
     * @ingroup transportProps
     */
    class TransportFactory {

    public:

        /**
         * Return a pointer to a TransportFactory
         * instance. TransportFactory is implemented as a 'singleton',
         * which means that at most one instance may be created. The
         * constructor is private. When a TransportFactory instance is
         * required, call static method factory() to return a pointer
         * to the TransportFactory instance.
         *
         * @code
         * TransportFactory* f;
         * f = TransportFactory::factory();
         * @endcode
         */
        static TransportFactory* factory() {
	  if (!s_factory) {
	      s_factory = new TransportFactory();
	    }
            return s_factory;
        }

	/**
	 * Deletes the statically malloced instance.
	 */
	static void deleteTransportFactory();

        /**
         * Destructor 
	 *
	 * We do not delete statically
	 * created single instance of this class here, because it would
	 * create an infinite loop if destructor is called for that
	 * single instance. 
         */
        virtual ~TransportFactory();
	

        /// Build a new transport manager
        virtual Transport*
        newTransport(string model="", thermo_t* thermo=0, int log_level=0);

        /// Initialize an existing transport manager
        virtual void initTransport(Transport* tr,  
            thermo_t* thermo=0, int mode=0, int log_level=0);


    private:

        static TransportFactory* s_factory;

        // The constructor is private; use static method factory() to
        // get a pointer to a factory instance
        TransportFactory();

        /// Read in transport parameters from a database
        //void readTransportDatabase(ostream& logfile, 
        //    XML_Node* db, 
        //    const vector<string>& names, 
        //    TransportParams& tr);

        void getTransportData(const XML_Node* db,  
            XML_Node& log, const vector<string>& names, 
            TransportParams& tr);

        /** Generate polynomial fits to viscosity, conductivity, and
         *  binary diffusion coefficients */
        void fitProperties(TransportParams& tr, ostream& logfile=cout);

        /// Generate polynomial fits to collision integrals
        void fitCollisionIntegrals(ostream& logfile, 
            TransportParams& tr);

        MMCollisionInt* m_integrals;

        void setupMM(ostream& flog,  const XML_Node* transport_database, 
            thermo_t* thermo, int mode, int log_level, 
            TransportParams& tr);

        /// construct a new power-law transport manager
        //Transport* newPowerTransport(const string& transport_database, 
        //    phase_t* mix);

        /// construct a new multicomponent transport manager
        //        Transport* newMultiTransport(const string& fname, 
        //    thermo_t* thermo, int mode = 0, int log_level = 0);

        /// construct a new mixture-averaged transport server
        //Transport* newMixTransport(const string& fname, 
        //    thermo_t* thermo, int mode = 0, int log_level = 0);


        /// Second-order correction to the binary diffusion coefficients
        void getBinDiffCorrection(doublereal t, 
            const TransportParams& tr, int k, int j, doublereal xk, doublereal xj, 
            doublereal& fkj, doublereal& fjk);

        /// Corrections for polar-nonpolar binary diffusion coefficients
        void makePolarCorrections(int i, int j, 
            const TransportParams& tr, doublereal& f_eps, doublereal& f_sigma);

        map<string, int> m_models;
    };


    /**
     *  Create a new transport manager instance.
     * @ingroup transportProps
     */
    inline Transport* newTransportMgr(string transportModel="", 
        thermo_t* thermo=0, int loglevel=0, TransportFactory* f=0) {
        if (f == 0) {
            f = TransportFactory::factory();
        }
        Transport* ptr = f->newTransport(transportModel, thermo, loglevel);
	/*
	 * Note: We delete the static s_factory instance here, instead of in
	 *       appdelete() in misc.cpp, to avoid linking problems involving
	 *       the need for multiple cantera and transport library statements
	 *       for applications that don't have transport in them.
	 */
	TransportFactory::deleteTransportFactory();
	return ptr;
    }


}
#endif






