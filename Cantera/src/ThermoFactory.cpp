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
	    throw UnknownThermoPhaseModel("ThermoFactory::newThermoPhase",
					  model);
        }
        return th;
    }
}
