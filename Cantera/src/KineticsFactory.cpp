/**
 *  @file KineticsFactory.cpp
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

#include "KineticsFactory.h"

#include "GasKinetics.h"
#include "GRI_30_Kinetics.h"
#include "InterfaceKinetics.h"

#include "importCTML.h"

namespace Cantera {

    KineticsFactory* KineticsFactory::__factory = 0;

    static int ntypes = 3;
    static string _types[] = {"GasKinetics", "GRI30", "Interface"};
    static int _itypes[]   = {cGasKinetics, cGRI30, cInterfaceKinetics};


    /**
     * Return a new kinetics manager that implements a reaction
     * mechanism specified in a CTML file.
     */
    Kinetics* KineticsFactory::newKinetics(XML_Node& phase, 
        vector<ThermoPhase*> th) {

        string kintype = phase.child("kinetics")["model"];
        int ikin=-1;
        int n;
        for (n = 0; n < ntypes; n++) {
            if (kintype == _types[n]) ikin = _itypes[n];
        }
        Kinetics* k=0;
        switch (ikin) {

        case cGasKinetics:
            k = new GasKinetics;
            break;

        case cGRI30:
            k = new GRI_30_Kinetics;
            break;

        case cInterfaceKinetics:
            k = new InterfaceKinetics;
            break;

        default:
            throw CanteraError("newKinetics",
                "unknown kinetics manager: "+kintype);
        }

        // import the reaction mechanism
        importKinetics(phase, th, k);
        return k;
    }


    /**
     * Return a new, empty kinetics manager.
     */
    Kinetics* KineticsFactory::newKinetics(string model) {

        int ikin = -1;
        int n;
        for (n = 0; n < ntypes; n++) {
            if (model == _types[n]) ikin = _itypes[n];
        }
        Kinetics* k=0;
        switch (ikin) {

        case cGasKinetics:
            k = new GasKinetics;
            break;

        case cGRI30:
            k = new GRI_30_Kinetics;
            break;

        case cInterfaceKinetics:
            k = new InterfaceKinetics;
            break;

        default:
            throw CanteraError("newKinetics",
                "unknown kinetics manager: "+model);
        }
        return k;
    }

}
