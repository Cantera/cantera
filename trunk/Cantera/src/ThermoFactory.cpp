/**
 *  @file ThermoFactory.cpp
 */

/*
 * $Author: dggoodwin $
 * $Revision: 1.13 $
 * $Date: 2006/05/30 16:15:48 $
 */

// Copyright 2001  California Institute of Technology


#ifdef WIN32
#pragma warning(disable:4786)
#endif

#include "ThermoFactory.h"

#include "SpeciesThermoFactory.h"
#include "IdealGasPhase.h"

#ifdef WITH_PURE_FLUIDS
#include "PureFluidPhase.h"
#endif

#include "ConstDensityThermo.h"
#include "SurfPhase.h"
#include "EdgePhase.h"

#ifdef WITH_METAL
#include "MetalPhase.h"
#endif

#ifdef WITH_STOICH_SUBSTANCE
#include "StoichSubstance.h"
#endif

#include "importCTML.h"

#ifdef WITH_LATTICE_SOLID
#include "LatticeSolidPhase.h"
#include "LatticePhase.h"
#endif

namespace Cantera {

    ThermoFactory* ThermoFactory::s_factory = 0;

    static int ntypes = 9;
    static string _types[] = {"IdealGas", "Incompressible", 
                              "Surface", "Edge", "Metal", "StoichSubstance",
                              "PureFluid", "LatticeSolid", "Lattice"};

    static int _itypes[]   = {cIdealGas, cIncompressible, 
                              cSurf, cEdge, cMetal, cStoichSubstance,
                              cPureFluid, cLatticeSolid, cLattice};

    /**
     * This method returns a new instance of a subclass of ThermoPhase
     */ 
    ThermoPhase* ThermoFactory::newThermoPhase(string model) {

        int ieos=-1;

        for (int n = 0; n < ntypes; n++) {
            if (model == _types[n]) ieos = _itypes[n];
        }

        ThermoPhase* th=0;
        //        map<string, double> d;
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

#ifdef WITH_METAL
        case cMetal:
            th = new MetalPhase;
            break;
#endif

#ifdef WITH_STOICH_SUBSTANCE
        case cStoichSubstance:
            th = new StoichSubstance;
            break;
#endif

#ifdef WITH_LATTICE_SOLID
        case cLatticeSolid:
            th = new LatticeSolidPhase;
            break;

        case cLattice:
            th = new LatticePhase;
            break;
#endif

#ifdef WITH_PURE_FLUIDS
        case cPureFluid:
            th = new PureFluidPhase;
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
