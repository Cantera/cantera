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
#include "PureFluidPhase.h"
#include "ConstDensityThermo.h"
#include "SurfPhase.h"
#include "EdgePhase.h"
#include "MetalPhase.h"
#include "SolidCompound.h"
#include "importCTML.h"

namespace Cantera {

    ThermoFactory* ThermoFactory::__factory = 0;

    static int ntypes = 7;
    static string _types[] = {"IdealGas", "Incompressible", 
                              "Surface", "Edge", "Metal", "SolidCompound",
                              "PureFluid"};

    static int _itypes[]   = {cIdealGas, cIncompressible, 
                              cSurf, cEdge, cMetal, cSolidCompound,
                              cPureFluid};

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

        case cEdge:
            th = new EdgePhase;
            break;

        case cMetal:
            th = new MetalPhase;
            break;

        case cSolidCompound:
            th = new SolidCompound;
            break;

#ifdef INCL_PURE_FLUIDS
        case cPureFluid:
            th = new PureFluid;
            break;
#endif

        default:
	    throw UnknownThermoPhaseModel("ThermoFactory::newThermoPhase",
					  model);
        }
        return th;
    }
}
