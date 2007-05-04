/**
 *  @file ThermoFactory.cpp
 *     Definitions for the factory class that can create known %ThermoPhase objects
 *     (see \ref thermoprops and class \link Cantera::ThermoFactory ThermoFactory\endlink).
 *
 
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

#include "speciesThermoTypes.h"
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

#undef USE_SSTP
#ifdef WITH_STOICH_SUBSTANCE
#ifdef USE_SSTP
#include "StoichSubstanceSSTP.h"
#else
#include "StoichSubstance.h"
#endif
#endif

//#include "importCTML.h"

#ifdef WITH_LATTICE_SOLID
#include "LatticeSolidPhase.h"
#include "LatticePhase.h"
#endif

using namespace std;

namespace Cantera {

    ThermoFactory* ThermoFactory::s_factory = 0;

    static int ntypes = 9;
    static string _types[] = {"IdealGas", "Incompressible", 
                              "Surface", "Edge", "Metal", "StoichSubstance",
                              "PureFluid", "LatticeSolid", "Lattice"};

    static int _itypes[]   = {cIdealGas, cIncompressible, 
                              cSurf, cEdge, cMetal, cStoichSubstance,
                              cPureFluid, cLatticeSolid, cLattice};

    /*
     * This method returns a new instance of a subclass of ThermoPhase
     */ 
    ThermoPhase* ThermoFactory::newThermoPhase(std::string model) {

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
#ifdef USE_SSTP
            th = new StoichSubstanceSSTP;
#else
            th = new StoichSubstance;
#endif
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



  /*
   * Create a new ThermoPhase object and initializes it according to
   * the XML tree database.  This routine first looks up the
   * identity of the model for the solution thermodynamics in the
   * model attribute of the thermo child of the xml phase
   * node. Then, it does a string lookup on the model to figure out
   * what ThermoPhase derived class is assigned. It creates a new
   * instance of that class, and then calls importPhase() to
   * populate that class with the correct parameters from the XML
   * tree.
   */
  ThermoPhase* newPhase(XML_Node& xmlphase) {
    const XML_Node& th = xmlphase.child("thermo");
    string model = th["model"];
    ThermoPhase* t = newThermoPhase(model);
    importPhase(xmlphase, t);
    return t;
  }

  ThermoPhase* newPhase(std::string infile, std::string id) {
    XML_Node* root = get_XML_File(infile); 
    if (id == "-") id = "";
    XML_Node* x = get_XML_Node(string("#")+id, root);
    if (x) 
      return newPhase(*x);
    else
      return 0;
  }


        
  /*
   * Import a phase specification.
   *   Here we read an XML description of the phase.
   *   We import descriptions of the elements that make up the
   *   species in a phase.
   *   We import information about the species, including their
   *   reference state thermodynamic polynomials. We then freeze
   *   the state of the species, and finally call initThermoXML(phase, id)
   *   a member function of the ThermoPhase object to "finish"
   *   the description.
   *
   *
   * @param phase This object must be the phase node of a
   *             complete XML tree
   *             description of the phase, including all of the
   *             species data. In other words while "phase" must
   *             point to an XML phase object, it must have
   *             sibling nodes "speciesData" that describe
   *             the species in the phase.
   * @param th   Pointer to the ThermoPhase object which will
   *             handle the thermodynamics for this phase.
   *             We initialize part of the Thermophase object
   *             here, especially for those objects which are
   *             part of the Cantera Kernel.
   */
  bool importPhase(XML_Node& phase, ThermoPhase* th, 
		   SpeciesThermoFactory* spfactory) {

    // Check the the supplied XML node in fact represents a 
    // phase.
    if (phase.name() != "phase") 
      throw CanteraError("importPhase",
			 "Current const XML_Node is not a phase element.");

    // if no species thermo factory was supplied,
    // use the default one. 
    if (!spfactory) 
      spfactory = SpeciesThermoFactory::factory();

    // set the id attribute of the phase to the 'id' attribute 
    // in the XML tree.
    th->setID(phase.id());
    th->setName(phase.id());

    // Number of spatial dimensions. Defaults to 3 (bulk phase)
    if (phase.hasAttrib("dim")) {
      int idim = intValue(phase["dim"]);
      if (idim < 1 || idim > 3)
	throw CanteraError("importPhase",
			   "unphysical number of dimensions: "+phase["dim"]);
      th->setNDim(idim);
    }
    else
      th->setNDim(3);     // default


	
    // Set equation of state parameters. The parameters are
    // specific to each subclass of ThermoPhase, so this is done
    // by method setParametersFromXML in each subclass.
    if (phase.hasChild("thermo")) {
      const XML_Node& eos = phase.child("thermo");
      th->setParametersFromXML(eos);
    }


    /***************************************************************
     * Add the elements.
     ***************************************************************/
    th->addElementsFromXML(phase);


    /***************************************************************
     * Add the species. 
     *
     * Species definitions may be imported from multiple
     * sources. For each one, a speciesArray element must be
     * present.
     ***************************************************************/
    XML_Node* db = 0;
    vector<XML_Node*> sparrays;
    phase.getChildren("speciesArray", sparrays);
    int jsp, nspa = static_cast<int>(sparrays.size());
    vector<XML_Node*> dbases;
    vector_int sprule(nspa,0);

    // loop over the speciesArray elements
    for (jsp = 0; jsp < nspa; jsp++) {

      const XML_Node& species = *sparrays[jsp];

      // If the speciesArray element has a child element
      //   <skip element="undeclared"> 
      // then set sprule[jsp] to 1, so
      // that any species with an undeclared element will be
      // quietly skipped when importing species.
      if (species.hasChild("skip")) {
	const XML_Node& sk = species.child("skip");
	string eskip = sk["element"];
	if (eskip == "undeclared") {
	  sprule[jsp] = 1;
	}
	string dskip = sk["species"];
	if (dskip == "duplicate") {
	  sprule[jsp] += 10;
	}
      }

      string fname, idstr;

      // get a pointer to the node containing the species
      // definitions for the species declared in this 
      // speciesArray element. This may be in the local file
      // containing the phase element, or may be in another
      // file.            
      db = get_XML_Node(species["datasrc"], &phase.root());

      // add this node to the list of species database nodes.
      dbases.push_back(db);
    }


    // if the phase has a species thermo manager already installed,
    // delete it since we are adding new species.
    delete &th->speciesThermo();

    // create a new species thermo manager.  Function
    // 'newSpeciesThermoMgr' looks at the species in the database
    // to see what thermodynamic property parameterizations are
    // used, and selects a class that can handle the
    // parameterizations found.
    SpeciesThermo* spth = newSpeciesThermoMgr(dbases);

    // install it in the phase object
    th->setSpeciesThermo(spth);
    SpeciesThermo& spthermo = th->speciesThermo();

    // used to check that each species is declared only once
    map<string,bool> declared;

    int i, k = 0;

    // loop over the species arrays
    for (jsp = 0; jsp < nspa; jsp++) {

      const XML_Node& species = *sparrays[jsp]; 
      db = dbases[jsp];

      // Get the array of species name strings.
      vector<string> spnames;
      getStringArray(species, spnames);
      int nsp = static_cast<int>(spnames.size());

      // if 'all' is specified, then add all species 
      // defined in this database to the phase
      if (nsp == 1 && spnames[0] == "all") {
	vector<XML_Node*> allsp;
	db->getChildren("species",allsp);
	nsp = static_cast<int>(allsp.size());
	spnames.resize(nsp);
	for (int nn = 0; nn < nsp; nn++) {
	  spnames[nn] = (*allsp[nn])["name"];
	}
      }
      else if (nsp == 1 && spnames[0] == "unique") {
	vector<XML_Node*> uniquesp;
	db->getChildren("species",uniquesp);
	nsp = static_cast<int>(uniquesp.size());
	spnames.clear();
	spnames.resize(nsp);
	string spnm;
	for (int nn = 0; nn < nsp; nn++) {
	  spnm = (*uniquesp[nn])["name"];
	  if (!declared[spnm]) spnames[nn] = spnm;
	}
      }

      string name;
      bool skip;
      for (i = 0; i < nsp; i++) {
	name = spnames[i];
	skip = false;
	if (name == "") skip = true;
	// Check that every species is only declared once
	if (declared[name]) {
	  if (sprule[jsp] >= 10) 
	    skip = true;
	  else
	    throw CanteraError("importPhase",
			       "duplicate species: "+name);
	}
	if (!skip) {
	  declared[name] = true;

	  // Find the species in the database by name.
	  XML_Node* s = db->findByAttr("name",spnames[i]);
	  if (s) {
	    if (installSpecies(k, *s, *th, spthermo, sprule[jsp], 
			       spfactory)) 
	      ++k;
	  }
	  else {
	    throw CanteraError("importPhase","no data for species "
			       +name);
	  }
	}
      }
    }

    // done adding species. 
    th->freezeSpecies();

    th->saveSpeciesData(db);

    // Perform any required subclass-specific initialization.
    string id = "";
    th->initThermoXML(phase, id);

    return true;
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



  /*
   * Install a species into a ThermoPhase object, which defines
   * the phase thermodynamics and speciation.
   *
   *  This routine first gathers the information from the Species XML
   *  tree and calls addUniqueSpecies() to add it to the
   *  ThermoPhase object, p.
   *  This information consists of:
   *         ecomp[] = element composition of species.
   *         chgr    = electric charge of species
   *         name    = string name of species
   *         sz      = size of the species 
   *                 (option double used a lot in thermo)
   *
   *  Then, the routine processes the "thermo" XML element and
   *  calls underlying utility routines to read the XML elements
   *  containing the thermodynamic information for the reference
   *  state of the species. Failures or lack of information trigger
   *  an "UnknownSpeciesThermoModel" exception being thrown.
   */
  bool installSpecies(int k, const XML_Node& s, thermo_t& p, 
		      SpeciesThermo& spthermo, int rule, 
		      SpeciesThermoFactory* factory) {

    // get the composition of the species
    const XML_Node& a = s.child("atomArray");
    map<string,string> comp;
    getMap(a, comp);

    // check that all elements in the species
    // exist in 'p'. If rule != 0, quietly skip 
    // this species if some elements are undeclared;
    // otherwise, throw an exception
    map<string,string>::const_iterator _b = comp.begin();
    for (; _b != comp.end(); ++_b) {
      if (p.elementIndex(_b->first) < 0) {
	if (rule == 0) {
	  throw CanteraError("installSpecies", 
			     "Species " + s["name"] + 
			     " contains undeclared element " + _b->first);
	}
	else
	  return false;
      }
    }

    // construct a vector of atom numbers for each 
    // element in phase p. Elements not declared in the
    // species (i.e., not in map comp) will have zero
    // entries in the vector.
    int m, nel = p.nElements();
    vector_fp ecomp(nel, 0.0);            
    for (m = 0; m < nel; m++) {
      ecomp[m] = atoi(comp[p.elementName(m)].c_str());
    }


    // get the species charge, if any. Note that the charge need
    // not be explicitly specified if special element 'E'
    // (electron) is one of the elements.
    doublereal chrg = 0.0;
    if (s.hasChild("charge")) chrg = getFloat(s, "charge");

    // get the species size, if any. (This is used by surface
    // phases to represent how many sites a species occupies.)
    doublereal sz = 1.0;
    if (s.hasChild("size")) sz = getFloat(s, "size");

    // add the species to phase p.
    p.addUniqueSpecies(s["name"], &ecomp[0], chrg, sz);

    // install the thermo parameterization for this species into
    // the species thermo manager for phase p.
    factory->installThermoForSpecies(k, s, spthermo);
        
    return true;
  }

  /*
   * Search an XML tree for species data.
   *
   *   This utility routine will search the XML tree for the species
   *   named by the string, kname. It will return the XML_Node
   *   pointer.
   *   Failures of any kind return the null pointer.
   */
  const XML_Node *speciesXML_Node(std::string kname,
				  const XML_Node *phaseSpeciesData) {
    /*
     * First look at the species database.
     *  -> Look for the subelement "stoichIsMods"
     *     in each of the species SS databases.
     */
    if (!phaseSpeciesData) return ((const XML_Node *) 0);
    string jname;
    vector<XML_Node*> xspecies;
    phaseSpeciesData->getChildren("species", xspecies);
    int jj = xspecies.size();
    for (int j = 0; j < jj; j++) {
      const XML_Node& sp = *xspecies[j];
      jname = sp["name"];
      if (jname == kname) {
	return &sp;
      }
    }
    return ((const XML_Node *) 0);
  }

}
