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
//#include "SolidCompound.h"
#include "StoichSubstance.h"
#include "importCTML.h"

namespace Cantera {

    ThermoFactory* ThermoFactory::__factory = 0;

    static int ntypes = 7;
    static string _types[] = {"IdealGas", "Incompressible", 
                              "Surface", "Edge", "Metal", "StoichSubstance",
                              "PureFluid"};

    static int _itypes[]   = {cIdealGas, cIncompressible, 
                              cSurf, cEdge, cMetal, cStoichSubstance,
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

        case cStoichSubstance:
            th = new StoichSubstance;
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


//     void setEOSParameters(const XML_Node& xmlphase, ThermoPhase* th) {

//         // if no thermo model is specified for the phase, simply
//         // return
//         if (!phase.hasChild("thermo")) return;

//         const XML_Node& eos = phase.child("thermo");

//         // set the parameters for the particular equation of state type,
//         // and 
//         if (eos["model"] == "Incompressible") {
//             if (th->eosType() == cIncompressible) {
//                 doublereal rho = getFloat(eos, "density", "-");
//                 th->setParameters(1, &rho);
//             }
//             else {
//                 eoserror = true;
//             }
//         }
//             else if (eos["model"] == "StoichSubstance") {
//                 if (th->eosType() == cStoichSubstance) {
//                     doublereal rho = getFloat(eos, "density", "-");
//                     th->setDensity(rho);
//                 }
//                 else {
//                     eoserror = true;
//                 }
//             }
//             else if (eos["model"] == "Surface") {
//                 if (th->eosType() == cSurf) {
//                     doublereal n = getFloat(eos, "site_density", "-");
//                     if (n <= 0.0) 
//                         throw CanteraError("importCTML",
//                             "missing or negative site density");
//                     th->setParameters(1, &n);
//                 }
//                 else {
//                     eoserror = true;
//                 }
//             }
//             else if (eos["model"] == "Edge") {
//                 if (th->eosType() == cEdge) {
//                     doublereal n = getFloat(eos, "site_density", "-");
//                     if (n <= 0.0) 
//                         throw CanteraError("importCTML",
//                             "missing or negative site density");
//                     th->setParameters(1, &n);
//                 }
//                 else {
//                     eoserror = true;
//                 }
//             }
// #ifdef INCL_PURE_FLUIDS
//             else if (eos["model"] == "PureFluid") {
//                 if (th->eosType() == cPureFluid) {
//                     subflag = atoi(eos["fluid_type"].c_str());
//                     if (subflag < 0) 
//                         throw CanteraError("importCTML",
//                             "missing fluid type flag");
//                 }
//                 else {
//                     eoserror = true;
//                 }
//             }
// #endif
//             if (eoserror) {
//                 string msg = "Wrong equation of state type for phase "+phase["id"]+"\n";
//                 msg += eos["model"]+" is not consistent with eos type "+int2str(th->eosType());
//                 throw CanteraError("importCTML",msg);
//             }

}
