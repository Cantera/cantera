/**
 *  @file SpeciesThermoFactory.cpp
 *    Definitions for factory to build instances of classes that manage the
 *    standard-state thermodynamic properties of a set of species 
 *    (see \ref spthermo and class \link Cantera::SpeciesThermoFactory SpeciesThermoFactory\endlink);
 */
/*
 * $Id$
 */
// Copyright 2001  California Institute of Technology

#ifdef WIN32
#pragma warning(disable:4786)
#endif


#include "SpeciesThermoFactory.h"
using namespace std;

#include "SpeciesThermo.h"
#include "NasaThermo.h"
#include "ShomateThermo.h"
#include "SimpleThermo.h"
#include "GeneralSpeciesThermo.h"
#include "Mu0Poly.h"
#include "Nasa9PolyMultiTempRegion.h"
#include "Nasa9Poly1.h"

#include "SpeciesThermoMgr.h"
#include "speciesThermoTypes.h"

#include "xml.h"
#include "ctml.h"

using namespace ctml;


namespace Cantera {

  SpeciesThermoFactory* SpeciesThermoFactory::s_factory = 0;

#if defined(THREAD_SAFE_CANTERA)
  boost::mutex SpeciesThermoFactory::species_thermo_mutex ;
#endif
 

  /**
   * Examine the types of species thermo parameterizations,
   * and return a flag indicating the type of parameterization
   * needed by the species.
   *
   *  @param spData_node Species Data XML node. This node contains a list
   *                     of species XML nodes underneath it.
   * 
   * @todo Make sure that spDadta_node is species Data XML node by checking its name is speciesData
   */
  static void getSpeciesThermoTypes(XML_Node* spData_node, 
                    int& has_nasa, int& has_shomate, int& has_simple,
                    int &has_other) {
    const XML_Node& sparray = *spData_node;
    std::vector<XML_Node*> sp;

    // get all of the species nodes
    sparray.getChildren("species",sp);
    size_t n, ns = sp.size();
    for (n = 0; n < ns; n++) {
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
	if (th.hasChild("Mu0")) has_other = 1;
	if (th.hasChild("NASA9")) has_other = 1;
	if (th.hasChild("NASA9MULTITEMP")) has_other = 1;
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
  SpeciesThermo* SpeciesThermoFactory::newSpeciesThermo(XML_Node* spData_node) {
    int inasa = 0, ishomate = 0, isimple = 0, iother = 0;
    try {
      getSpeciesThermoTypes(spData_node, inasa, ishomate, isimple, iother);
    } catch (UnknownSpeciesThermoModel) {
      iother = 1;
      popError();
    }
    if (iother) {
      writelog("returning new GeneralSpeciesThermo");
      return new GeneralSpeciesThermo();
    }
    return newSpeciesThermo(NASA*inasa
                + SHOMATE*ishomate + SIMPLE*isimple);
  }
    
  SpeciesThermo* SpeciesThermoFactory::
  newSpeciesThermo(std::vector<XML_Node*> spData_nodes) {
    int n = static_cast<int>(spData_nodes.size());
    int inasa = 0, ishomate = 0, isimple = 0, iother = 0;
    for (int j = 0; j < n; j++) {
      try {
    getSpeciesThermoTypes(spData_nodes[j], inasa, ishomate, isimple, iother);
      } catch (UnknownSpeciesThermoModel) {
    iother = 1;
    popError();
      }
    }
    if (iother) {
      return new GeneralSpeciesThermo();
    }
    return newSpeciesThermo(NASA*inasa
                + SHOMATE*ishomate + SIMPLE*isimple);
  }


  /*
   * @todo is this used? 
   */
  SpeciesThermo* SpeciesThermoFactory::
  newSpeciesThermoOpt(std::vector<XML_Node*> nodes) {
    int n = static_cast<int>(nodes.size());
    int inasa = 0, ishomate = 0, isimple = 0, iother = 0;
    for (int j = 0; j < n; j++) {
      try {
    getSpeciesThermoTypes(nodes[j], inasa, ishomate, isimple, iother);
      } catch (UnknownSpeciesThermoModel) {
    iother = 1;
    popError();
      }
    }
    if (iother) {
      return new GeneralSpeciesThermo();
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
    case SHOMATE + SIMPLE:
      return new SpeciesThermoDuo<ShomateThermo, SimpleThermo>;
    default:
      throw UnknownSpeciesThermo(
                 "SpeciesThermoFactory::newSpeciesThermo",type);
      return 0; 
    }
  }


  /*
   * Check the continuity of properties at the midpoint
   * temperature.
   */
  void NasaThermo::checkContinuity(std::string name, double tmid, const doublereal* clow,
                   doublereal* chigh) {

        // heat capacity
        doublereal cplow = poly4(tmid, clow);
        doublereal cphigh = poly4(tmid, chigh);
        doublereal delta = cplow - cphigh;
        if (fabs(delta/cplow) > 0.001) {
            writelog("\n\n**** WARNING ****\nFor species "+name+
                ", discontinuity in cp/R detected at Tmid = "
                +fp2str(tmid)+"\n");
            writelog("\tValue computed using low-temperature polynomial:  "
                +fp2str(cplow)+".\n");
            writelog("\tValue computed using high-temperature polynomial: "
                +fp2str(cphigh)+".\n");
        }

        // enthalpy
        doublereal hrtlow = enthalpy_RT(tmid, clow);
        doublereal hrthigh = enthalpy_RT(tmid, chigh);
        delta = hrtlow - hrthigh;
        if (fabs(delta/hrtlow) > 0.001) {
            writelog("\n\n**** WARNING ****\nFor species "+name+
                ", discontinuity in h/RT detected at Tmid = "
                +fp2str(tmid)+"\n");
            writelog("\tValue computed using low-temperature polynomial:  "
                +fp2str(hrtlow)+".\n");
            writelog("\tValue computed using high-temperature polynomial: "
                +fp2str(hrthigh)+".\n");
        }

        // entropy
        doublereal srlow = entropy_R(tmid, clow);
        doublereal srhigh = entropy_R(tmid, chigh);
        delta = srlow - srhigh;
        if (fabs(delta/srlow) > 0.001) {
            writelog("\n\n**** WARNING ****\nFor species "+name+
                ", discontinuity in s/R detected at Tmid = "
                +fp2str(tmid)+"\n");
            writelog("\tValue computed using low-temperature polynomial:  "
                +fp2str(srlow)+".\n");
            writelog("\tValue computed using high-temperature polynomial: "
                +fp2str(srhigh)+".\n");
        }
    }


    /** 
     * Install a NASA polynomial thermodynamic property
     * parameterization for species k into a SpeciesThermo instance.
     * This is called by method installThermoForSpecies if a NASA
     * block is found in the XML input.
     */
    static void installNasaThermoFromXML(std::string speciesName,
        SpeciesThermo& sp, int k, 
        const XML_Node* f0ptr, const XML_Node* f1ptr) {
        doublereal tmin0, tmax0, tmin1, tmax1, tmin, tmid, tmax;

        const XML_Node& f0 = *f0ptr;

        // default to a single temperature range
        bool dualRange = false;

        // but if f1ptr is suppled, then it is a two-range
        // parameterization
        if (f1ptr) {dualRange = true;}

        tmin0 = fpValue(f0["Tmin"]);
        tmax0 = fpValue(f0["Tmax"]);
        tmin1 = tmax0;
        tmax1 = tmin1 + 0.0001;
        if (dualRange) {
            tmin1 = fpValue((*f1ptr)["Tmin"]);
            tmax1 = fpValue((*f1ptr)["Tmax"]);
        }

        vector_fp c0, c1;
        if (fabs(tmax0 - tmin1) < 0.01) {
            // f0 has the lower T data, and f1 the higher T data
            tmin = tmin0;
            tmid = tmax0;
            tmax = tmax1;
            getFloatArray(f0.child("floatArray"), c0, false);
            if (dualRange)
                getFloatArray(f1ptr->child("floatArray"), c1, false);
            else {
                // if there is no higher range data, then copy c0 to c1.
                c1.resize(7,0.0);
                copy(c0.begin(), c0.end(), c1.begin());
            }
        }
        else if (fabs(tmax1 - tmin0) < 0.01) {
            // f1 has the lower T data, and f0 the higher T data
            tmin = tmin1;
            tmid = tmax1;
            tmax = tmax0;
            getFloatArray(f1ptr->child("floatArray"), c0, false);
            getFloatArray(f0.child("floatArray"), c1, false);
        }
        else {
            throw CanteraError("installNasaThermo",
                   "non-continuous temperature ranges.");
        }

        // The NasaThermo species property manager expects the
        // coefficients in a different order, so rearrange them.
        array_fp c(15);
        c[0] = tmid;
        doublereal p0 = OneAtm;
        c[1] = c0[5];
        c[2] = c0[6];
        copy(c0.begin(), c0.begin()+5, c.begin() + 3);
        c[8] = c1[5];
        c[9] = c1[6];
        copy(c1.begin(), c1.begin()+5, c.begin() + 10);
        sp.install(speciesName, k, NASA, &c[0], tmin, tmax, p0);
    }

#ifdef INCL_NASA96

    /** 
     * Install a NASA96 polynomial thermodynamic property
     * parameterization for species k into a SpeciesThermo instance.
     */
    static void installNasa96ThermoFromXML(std::string speciesName,
        SpeciesThermo& sp, int k, 
        const XML_Node* f0ptr, const XML_Node* f1ptr) {
        doublereal tmin0, tmax0, tmin1, tmax1, tmin, tmid, tmax;

        const XML_Node& f0 = *f0ptr;
        bool dualRange = false;
        if (f1ptr) {dualRange = true;}
        tmin0 = fpValue(f0["Tmin"]);
        tmax0 = fpValue(f0["Tmax"]);
        tmin1 = tmax0;
        tmax1 = tmin1 + 0.0001;
        if (dualRange) {
            tmin1 = fpValue((*f1ptr)["Tmin"]);
            tmax1 = fpValue((*f1ptr)["Tmax"]);
        }

        vector_fp c0, c1;
        if (fabs(tmax0 - tmin1) < 0.01) {
            tmin = tmin0;
            tmid = tmax0;
            tmax = tmax1;
            getFloatArray(f0.child("floatArray"), c0, false);
            if (dualRange)
                getFloatArray(f1ptr->child("floatArray"), c1, false);
            else {
                c1.resize(7,0.0);
                copy(c0.begin(), c0.end(), c1.begin());
            }
        }
        else if (fabs(tmax1 - tmin0) < 0.01) {
            tmin = tmin1;
            tmid = tmax1;
            tmax = tmax0;
            getFloatArray(f1ptr->child("floatArray"), c0, false);
            getFloatArray(f0.child("floatArray"), c1, false);
        }
        else {
            throw CanteraError("installNasaThermo",
                   "non-continuous temperature ranges.");
        }
        array_fp c(15);
        c[0] = tmid;
        doublereal p0 = OneAtm;
        c[1] = c0[5];
        c[2] = c0[6];
        copy(c0.begin(), c0.begin()+5, c.begin() + 3);
        c[8] = c1[5];
        c[9] = c1[6];
        copy(c1.begin(), c1.begin()+5, c.begin() + 10);
        sp.install(speciesName, k, NASA, &c[0], tmin, tmax, p0);
    }

#endif


    /** 
     * Install a Shomate polynomial thermodynamic property
     * parameterization for species k.
     */
    static void installShomateThermoFromXML(std::string speciesName, 
        SpeciesThermo& sp, int k, 
        const XML_Node* f0ptr, const XML_Node* f1ptr) {
        doublereal tmin0, tmax0, tmin1, tmax1, tmin, tmid, tmax;

        const XML_Node& f0 = *f0ptr;
        bool dualRange = false;
        if (f1ptr) {dualRange = true;}
        tmin0 = fpValue(f0["Tmin"]);
        tmax0 = fpValue(f0["Tmax"]);
        tmin1 = tmax0;
        tmax1 = tmin1 + 0.0001;
        if (dualRange) {
            tmin1 = fpValue((*f1ptr)["Tmin"]);
            tmax1 = fpValue((*f1ptr)["Tmax"]);
        }

        vector_fp c0, c1;
        if (fabs(tmax0 - tmin1) < 0.01) {
            tmin = tmin0;
            tmid = tmax0;
            tmax = tmax1;
            getFloatArray(f0.child("floatArray"), c0, false);
            if (dualRange)
                getFloatArray(f1ptr->child("floatArray"), c1, false);
            else {
          c1.resize(7,0.0);
          copy(c0.begin(), c0.begin()+7, c1.begin());
        }
        }
        else if (fabs(tmax1 - tmin0) < 0.01) {
            tmin = tmin1;
            tmid = tmax1;
            tmax = tmax0;
            getFloatArray(f1ptr->child("floatArray"), c0, false);
            getFloatArray(f0.child("floatArray"), c1, false);
        }
        else {
            throw CanteraError("installShomateThermoFromXML",
                   "non-continuous temperature ranges.");
        }
        array_fp c(15);
        c[0] = tmid;
        doublereal p0 = OneAtm;
        copy(c0.begin(), c0.begin()+7, c.begin() + 1);
        copy(c1.begin(), c1.begin()+7, c.begin() + 8);
        sp.install(speciesName, k, SHOMATE, &c[0], tmin, tmax, p0);
    }



    /** 
     * Install a constant-cp thermodynamic property
     * parameterization for species k.
     */
    static void installSimpleThermoFromXML(std::string speciesName, 
        SpeciesThermo& sp, int k, 
        const XML_Node& f) {
        doublereal tmin, tmax;
        tmin = fpValue(f["Tmin"]);
        tmax = fpValue(f["Tmax"]);
        if (tmax == 0.0) tmax = 1.0e30;

        vector_fp c(4);
        c[0] = getFloat(f, "t0", "-");
        c[1] = getFloat(f, "h0", "-");
        c[2] = getFloat(f, "s0", "-");
        c[3] = getFloat(f, "cp0", "-");
        doublereal p0 = OneAtm;
        sp.install(speciesName, k, SIMPLE, &c[0], tmin, tmax, p0);
    }

   /** 
     * Install a NASA9 polynomial thermodynamic property
     * parameterization for species k into a SpeciesThermo instance.
     * This is called by method installThermoForSpecies if a NASA9
     * block is found in the XML input.
     */
    static void installNasa9ThermoFromXML(std::string speciesName,
					  SpeciesThermo& sp, int k, 
					  const std::vector<XML_Node*>& tp)
    { 				
      const XML_Node * fptr = tp[0];
      int nRegTmp = tp.size();
      int nRegions = 0;
      vector_fp cPoly;
      Nasa9Poly1 *np_ptr = 0; 
      std::vector<Nasa9Poly1 *> regionPtrs;
      doublereal tmin, tmax, pref;
      // Loop over all of the possible temperature regions
      for (int i = 0; i < nRegTmp; i++) {
	fptr = tp[i];
	if (fptr) {
	  if (fptr->name() == "NASA9") {
	    if (fptr->hasChild("floatArray")) {

	      tmin = fpValue((*fptr)["Tmin"]);
	      tmax = fpValue((*fptr)["Tmax"]);
	      pref = fpValue((*fptr)["P0"]);

	      getFloatArray(fptr->child("floatArray"), cPoly, false);
	      if (cPoly.size() != 9) {
		throw CanteraError("installNasa9ThermoFromXML",
				   "Expected 9 coeff polynomial");
	      }
	      np_ptr = new Nasa9Poly1(k, tmin, tmax, pref,
				      DATA_PTR(cPoly));
	      regionPtrs.push_back(np_ptr);
	      nRegions++;
	    } 
	  }
	}
      }
      if (nRegions == 0) {
	throw UnknownSpeciesThermoModel("installThermoForSpecies", 
					speciesName, "  " );
      } else if (nRegions == 1)  {
	sp.install_STIT(np_ptr);
      } else {
	Nasa9PolyMultiTempRegion*  npMulti_ptr = new  Nasa9PolyMultiTempRegion(regionPtrs);
	sp.install_STIT(npMulti_ptr);
      }
    }

    /**
     * Install a species thermodynamic property parameterization
     * for one species into a species thermo manager.
     * @param k species number
     * @param s XML node specifying species
     * @param spthermo species thermo manager
     */
  void SpeciesThermoFactory::
  installThermoForSpecies(int k, const XML_Node& s, 
			  SpeciesThermo& spthermo) {
    /*
     * Check to see that the species block has a thermo block
     * before processing. Throw an error if not there.
     */
    if (!(s.hasChild("thermo"))) {
      throw UnknownSpeciesThermoModel("installThermoForSpecies", 
				      s["name"], "<nonexistent>");
    }
    const XML_Node& thermo = s.child("thermo");
    const std::vector<XML_Node*>& tp = thermo.children();
    int nc = static_cast<int>(tp.size());
    if (nc == 1) {
      const XML_Node* f = tp[0];
      if (f->name() == "Shomate") {
	installShomateThermoFromXML(s["name"], spthermo, k, f, 0);
      }
      else if (f->name() == "const_cp") {
	installSimpleThermoFromXML(s["name"], spthermo, k, *f);
      }
      else if (f->name() == "NASA") {
	installNasaThermoFromXML(s["name"], spthermo, k, f, 0);
      }
      else if (f->name() == "Mu0") {
	installMu0ThermoFromXML(s["name"], spthermo, k, f);
      }
      else if (f->name() == "NASA9") {
	installNasa9ThermoFromXML(s["name"], spthermo, k, tp);
      } 
      else {
	throw UnknownSpeciesThermoModel("installThermoForSpecies", 
					s["name"], f->name());
      }
    }
    else if (nc == 2) {
      const XML_Node* f0 = tp[0];
      const XML_Node* f1 = tp[1];
      if (f0->name() == "NASA" && f1->name() == "NASA") {
	installNasaThermoFromXML(s["name"], spthermo, k, f0, f1);
      } 
      else if (f0->name() == "Shomate" && f1->name() == "Shomate") {
	installShomateThermoFromXML(s["name"], spthermo, k, f0, f1);
      } 
      else if (f0->name() == "NASA9" && f1->name() == "NASA9") {
	installNasa9ThermoFromXML(s["name"], spthermo, k, tp);
      } else {
	throw UnknownSpeciesThermoModel("installThermoForSpecies", s["name"], 
					f0->name() + " and "
					+ f1->name());
      }
    }
    else if (nc >= 2) {
      const XML_Node* f0 = tp[0];
      if (f0->name() == "NASA9") {
	installNasa9ThermoFromXML(s["name"], spthermo, k, tp);
      } else {
	throw UnknownSpeciesThermoModel("installThermoForSpecies", s["name"], 
					"multiple");
      }
    } else {
	throw UnknownSpeciesThermoModel("installThermoForSpecies", s["name"], 
					"multiple");
      }
  }

}
