/**
 *  @file ThermoFactory.cpp
 */

/*
 * $Author$
 * $Revision$
 * $Date$
 */

// Copyright 2001  California Institute of Technology


#ifdef WIN32
#pragma warning(disable:4786)
#endif

#include "ThermoFactory.h"

#include "SpeciesThermoFactory.h"
#include "IdealGasPhase.h"
#include "ConstDensityThermo.h"
#include "SurfPhase.h"
#include "MetalPhase.h"
#include "SolidCompound.h"
#include "importCTML.h"

namespace Cantera {

    ThermoFactory* ThermoFactory::__factory = 0;

    static int ntypes = 5;
    static string _types[] = {"IdealGas", "Incompressible", 
                              "Surface", "Metal", "SolidCompound"};
    static int _itypes[]   = {cIdealGas, cIncompressible, 
                              cSurf, cMetal, cSolidCompound};


    ThermoPhase* ThermoFactory::newThermoPhase(string model) {

        int ieos=-1;
        for (int n = 0; n < ntypes; n++) {
            if (model == _types[n]) ieos = _itypes[n];
        }

        ThermoPhase* th=0;
        map<string, double> d;
        switch (ieos) {

        case cIdealGas:
            th = new IdealGasPhase;
            break;

        case cIncompressible:
            th = new ConstDensityThermo;
            break;

        case cSurf:
            th = new SurfPhase;
            break;

        case cMetal:
            th = new MetalPhase;
            break;

        case cSolidCompound:
            th = new SolidCompound;
            break;

        default:
            throw CanteraError("newThermo",
                "newThermo: unknown equation of state: "+model);
        }
        return th;
    }

    /**
     * Return a thermo manager to handle the parameterizations
     * specified in a CTML phase specification.
     */
    ThermoPhase* ThermoFactory::newThermo(XML_Node& root, string id) {

        // Find the node with the specified id, check that it is
        // a 'phase' node, and set the phase id to 'id'.
        XML_Node* ph;
        ph = root.findID(id);
        if (ph == 0) return 0; // false;   // id not found
        XML_Node& node = *ph;

        if (node.name() != "phase") 
            throw CanteraError("newThermo","node with id = "+id
                +" is not a phase object.");

        //Phase* p = new Phase;
        //p->setID(id);                // set the phase id

        // get equaton of state type
        XML_Node& eos = node.child("thermo");
        string eostype = eos["model"];
        int ieos=-1;
        for (int n = 0; n < ntypes; n++) {
            if (eostype == _types[n]) ieos = _itypes[n];
        }

        // build species thermo manager
        SpeciesThermo* spthermo = newSpeciesThermoMgr(&node);

        ThermoPhase* th=0;
        //        doublereal dens;
        map<string, double> d;
        switch (ieos) {

        case cIdealGas:
            th = new IdealGasPhase;
            break;

        case cIncompressible:
            th = new ConstDensityThermo;
            break;

        case cSurf:
            th = new SurfPhase;
            break;

        case cMetal:
            th = new MetalPhase;
            break;

        case cSolidCompound:
            th = new SolidCompound;
            break;

        default:
            throw CanteraError("newThermo",
                "newThermo: unknown equation of state: "+eostype);
        }
        th->setSpeciesThermo(spthermo);

        // import the phase specification
        importPhase(node, th);

        return th;
    }

}
