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
        const XML_Node& sparray = *node;
        vector<XML_Node*> sp;
        sparray.getChildren("species",sp);
        int ns = sp.size();
        for (int n = 0; n < ns; n++) {
	  XML_Node* spNode = sp[n];
	  if (spNode->hasChild("thermo")) {
            const XML_Node& th = sp[n]->child("thermo");
            if (th.hasChild("NASA")) has_nasa = 1;
            if (th.hasChild("Shomate")) has_shomate = 1;
            if (th.hasChild("const_cp")) has_simple = 1;
            if (th.hasChild("poly")) {
                if (th.child("poly")["order"] == "1") has_simple = 1;
                else throw CanteraError("newSpeciesThermo",
                    "poly with order > 1 not yet supported");
            }
	  } else {
	    throw UnknownSpeciesThermoModel("getSpeciesThermoTypes:",
					    spNode->attrib("name"), "missing");
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


    SpeciesThermo* SpeciesThermoFactory::newSpeciesThermoOpt(vector<XML_Node*> nodes) {
        int n = nodes.size();
        int inasa = 0, ishomate = 0, isimple = 0;
        for (int j = 0; j < n; j++) {
	  try {
            getSpeciesThermoTypes(nodes[j], inasa, ishomate, isimple);
	  } catch (UnknownSpeciesThermoModel) {

	  }
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


        /// Check the continuity of properties at the midpoint
        /// temperature, and adjust the high-T coefficients to
        /// make the properties exactly continuous at Tmid.
    void NasaThermo::checkContinuity(double tmid, const doublereal* clow,
            doublereal* chigh) {

            // heat capacity
            doublereal cplow = poly4(tmid, clow+2);
            doublereal cphigh = poly4(tmid, chigh+2);
            doublereal delta = cplow - cphigh;
            if (fabs(delta/cplow) > 0.001) {
                writelog("\n**** WARNING ****\nDiscontinuity in cp/R detected at Tmid = "
                    +fp2str(tmid)+"\n");
                writelog("\tValue computed using low-temperature polynomial:  "+fp2str(cplow)+".\n");
                writelog("\tValue computed using high-temperature polynomial: "+fp2str(cphigh)+".\n");
            }
            //chigh[2] += cplow - cphigh;

            // enthalpy
            doublereal hrtlow = enthalpy_RT(tmid, clow);
            doublereal hrthigh = enthalpy_RT(tmid, chigh);
            delta = hrtlow - hrthigh;
            //            chigh[0] += tmid*(hrtlow - hrthigh);
            if (fabs(delta/hrtlow) > 0.001) {
                writelog("\n**** WARNING ****\nDiscontinuity in h/RT detected at Tmid = "
                    +fp2str(tmid)+"\n");
                writelog("\tValue computed using low-temperature polynomial:  "+fp2str(hrtlow)+".\n");
                writelog("\tValue computed using high-temperature polynomial: "+fp2str(hrthigh)+".\n");
            }

            // entropy
            doublereal srlow = entropy_R(tmid, clow);
            doublereal srhigh = entropy_R(tmid, chigh);
            //chigh[1] += srlow - srhigh;
            delta = srlow - srhigh;
            if (fabs(delta/srlow) > 0.001) {
                writelog("\n**** WARNING ****\nDiscontinuity in s/R detected at Tmid = "
                    +fp2str(tmid)+"\n");
                writelog("\tValue computed using low-temperature polynomial:  "+fp2str(srlow)+".\n");
                writelog("\tValue computed using high-temperature polynomial: "+fp2str(srhigh)+".\n");
            }
        }

}
