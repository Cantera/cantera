/**
 *  @file SpeciesThermoFactory.cpp
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

#include "SpeciesThermoFactory.h"

#include "SpeciesThermo.h"
#include "NasaThermo.h"
#include "ShomateThermo.h"
#include "PolyThermoMgr.h"
#include "SimpleThermo.h"

#include "SpeciesThermoMgr.h"
#include "speciesThermoTypes.h"

#include "xml.h"

namespace Cantera {

    SpeciesThermoFactory* SpeciesThermoFactory::__factory = 0;


    static void getSpeciesThermoTypes(XML_Node* node, 
        int& has_nasa, int& has_shomate, int& has_simple) {
        XML_Node& sparray = *node;
        vector<XML_Node*> sp;
        sparray.getChildren("species",sp);
        int ns = sp.size();
        for (int n = 0; n < ns; n++) {
            XML_Node& th = sp[n]->child("thermo");
            if (th.hasChild("NASA")) has_nasa = 1;
            if (th.hasChild("Shomate")) has_shomate = 1;
            if (th.hasChild("const_cp")) has_simple = 1;
            if (th.hasChild("poly")) {
                if (th.child("poly")["order"] == "1") has_simple = 1;
                else throw CanteraError("newSpeciesThermo",
                    "poly with order > 1 not yet supported");
            }
        }
    }


    /**
     * Return a species thermo manager to handle the parameterizations
     * specified in a CTML phase specification.
     */
    SpeciesThermo* SpeciesThermoFactory::newSpeciesThermo(XML_Node* node) {
        int inasa = 0, ishomate = 0, isimple = 0;
        getSpeciesThermoTypes(node, inasa, ishomate, isimple);
        return newSpeciesThermo(NASA*inasa
            + SHOMATE*ishomate + SIMPLE*isimple);
    }

    SpeciesThermo* SpeciesThermoFactory::newSpeciesThermo(vector<XML_Node*> nodes) {
        int n = nodes.size();
        int inasa = 0, ishomate = 0, isimple = 0;
        for (int j = 0; j < n; j++) {
            getSpeciesThermoTypes(nodes[j], inasa, ishomate, isimple);
        }
        return newSpeciesThermo(NASA*inasa
            + SHOMATE*ishomate + SIMPLE*isimple);
    }


    SpeciesThermo* SpeciesThermoFactory::newSpeciesThermo(int type) {
        
        switch (type) {

        case NASA:
            return new NasaThermo;
        case SHOMATE:
            return new ShomateThermo;
        case SIMPLE:
            return new SimpleThermo;
        case NASA + SHOMATE:
            return new SpeciesThermoDuo<NasaThermo, ShomateThermo>;
        case NASA + SIMPLE:
            return new SpeciesThermoDuo<NasaThermo, SimpleThermo>;
        default:
            throw UnknownSpeciesThermo(
                "SpeciesThermoFactory::newSpeciesThermo",type);
            return 0; 
        }
    }

}
