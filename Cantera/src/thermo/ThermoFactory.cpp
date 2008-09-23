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
#include "VPSSMgr.h"
#include "VPSSMgrFactory.h"

#ifdef WITH_IDEAL_SOLUTIONS
#include "IdealSolidSolnPhase.h"
#endif

#ifdef WITH_PURE_FLUIDS
#include "PureFluidPhase.h"
#endif

#include "ConstDensityThermo.h"
#include "SurfPhase.h"
#include "EdgePhase.h"

#ifdef WITH_METAL
#include "MetalPhase.h"
#endif

#ifdef WITH_SEMICONDUCTOR
#include "SemiconductorPhase.h"
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

#ifdef WITH_ELECTROLYTES
#include "HMWSoln.h"
#include "DebyeHuckel.h"
#include "IdealMolalSoln.h"
#endif

#include "IdealSolnGasVPSS.h"

using namespace std;

namespace Cantera {

    ThermoFactory* ThermoFactory::s_factory = 0;
#if defined(THREAD_SAFE_CANTERA)
    boost::mutex ThermoFactory::thermo_mutex;
#endif

    static int ntypes = 14;
    static string _types[] = {"IdealGas", "Incompressible", 
                              "Surface", "Edge", "Metal", "StoichSubstance",
                              "PureFluid", "LatticeSolid", "Lattice",
                              "HMW", "IdealSolidSolution", "DebyeHuckel", 
                              "IdealMolalSolution", "IdealGasVPSS"
    };

    static int _itypes[]   = {cIdealGas, cIncompressible, 
                              cSurf, cEdge, cMetal, cStoichSubstance,
                              cPureFluid, cLatticeSolid, cLattice,
                              cHMW, cIdealSolidSolnPhase, cDebyeHuckel,
                              cIdealMolalSoln, cVPSS_IdealGas
    };

    /*
     * This method returns a new instance of a subclass of ThermoPhase
     */ 
    ThermoPhase* ThermoFactory::newThermoPhase(std::string model) {

        int ieos=-1;

        for (int n = 0; n < ntypes; n++) {
            if (model == _types[n]) ieos = _itypes[n];
        }

        ThermoPhase* th=0;
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

#ifdef WITH_IDEAL_SOLUTIONS
        case cIdealSolidSolnPhase:
            th = new IdealSolidSolnPhase();
            break;
#endif

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
#ifdef WITH_ELECTROLYTES
        case cHMW:
            th = new HMWSoln;
            break;

        case cDebyeHuckel:
            th = new DebyeHuckel;
            break;

        case cIdealMolalSoln:
            th = new IdealMolalSoln;
            break;
#endif

	case cVPSS_IdealGas:
	  th = new IdealSolnGasVPSS;
	  break;

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
#ifdef WITH_ELECTROLYTES
    if (model == "HMW") {
      HMWSoln* p = (HMWSoln*)t;
      p->constructPhaseXML(xmlphase,"");
    }
    else
#endif
      importPhase(xmlphase, t);
    return t;
  }

  ThermoPhase* newPhase(std::string infile, std::string id) {
    XML_Node* root = get_XML_File(infile); 
    if (id == "-") id = "";
    XML_Node* xphase = get_XML_NameID("phase", std::string("#")+id, root);
    if (!xphase) {
      throw CanteraError("newPhase",
			  "Couldn't find phase named \"" + id + "\" in file, " + infile);
    }
    if (xphase) 
      return newPhase(*xphase);
    else
      return (ThermoPhase *) 0;
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
    if (phase.name() != "phase") {
      throw CanteraError("importPhase",
			 "Current const XML_Node named, " + phase.name() + ", is not a phase element.");
    }

    // set the id attribute of the phase to the 'id' attribute 
    // in the XML tree.
    th->setID(phase.id());
    th->setName(phase.id());

    // Number of spatial dimensions. Defaults to 3 (bulk phase)
    if (phase.hasAttrib("dim")) {
      int idim = intValue(phase["dim"]);
      if (idim < 1 || idim > 3)
	throw CanteraError("importPhase",
			   "phase, " + th->id() + ", has unphysical number of dimensions: " + phase["dim"]);
      th->setNDim(idim);
    }
    else {
      th->setNDim(3);     // default
    }
	
    // Set equation of state parameters. The parameters are
    // specific to each subclass of ThermoPhase, so this is done
    // by method setParametersFromXML in each subclass.
    if (phase.hasChild("thermo")) {
      const XML_Node& eos = phase.child("thermo");
      th->setParametersFromXML(eos);
    } else {
      throw CanteraError("importPhase", 
			 " phase, " + th->id() + ", XML_Node does not have a \"thermo\" XML_Node");
    }

    VPStandardStateTP *vpss_ptr = 0;
    int ssConvention = th->standardStateConvention();
    if (ssConvention == cSS_CONVENTION_VPSS) {
      vpss_ptr = dynamic_cast <VPStandardStateTP *>(th);
      if (vpss_ptr == 0) {
	throw CanteraError("importPhase",
			   "phase, " + th->id() + ", was VPSS, but dynamic cast failed");
      }
    } 

    // if no species thermo factory was supplied,
    // use the default one. 
    if (!spfactory) {
      spfactory = SpeciesThermoFactory::factory();
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
    if (nspa == 0) {
      throw CanteraError("importPhase",
			 "phase, " + th->id() + ", has zero \"speciesArray\" XML nodes.\n"
			 + " There must be at least one speciesArray nodes with one or more species");
    }
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

    // Decide whether the the phase has a variable pressure ss or not
    SpeciesThermo* spth = 0;
    VPSSMgr* vp_spth = 0;
    if (ssConvention == cSS_CONVENTION_TEMPERATURE) {
      // Create a new species thermo manager.  Function
      // 'newSpeciesThermoMgr' looks at the species in the database
      // to see what thermodynamic property parameterizations are
      // used, and selects a class that can handle the
      // parameterizations found.
      spth = newSpeciesThermoMgr(dbases);
      
      // install it in the phase object
      th->setSpeciesThermo(spth);
      // SpeciesThermo& spthermo = th->speciesThermo();
    } else {
      vp_spth = newVPSSMgr(vpss_ptr, &phase, dbases);
      vpss_ptr->setVPSSMgr(vp_spth);
      spth = vp_spth->SpeciesThermoMgr();
      th->setSpeciesThermo(spth);
    }


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
			       "duplicate species: \"" + name + "\"");
	}
	if (!skip) {
	  declared[name] = true;

	  // Find the species in the database by name.
	  XML_Node* s = db->findByAttr("name",spnames[i]);
	  if (s) {
	    if (installSpecies(k, *s, *th, spth, sprule[jsp], 
                               &phase, vp_spth, spfactory)) 
	      ++k;
	  }
	  else {
	    throw CanteraError("importPhase","no data for species, \""
			       + name + "\"");
	  }
	}
      }
    }

    // done adding species. 
    th->freezeSpecies();

    th->saveSpeciesData(db);

    // Perform any required subclass-specific initialization.
    th->initThermo();
    string id = "";
    th->initThermoXML(phase, id);

    return true;
  }

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
   * *
   * @param k     Species Index in the phase
   * @param s     XML_Node containing the species data for this species.
   * @param p     Reference to the ThermoPhase object.
   * @param spthermo Reference to the SpeciesThermo object, where
   *              the standard state thermo properties for this
   *              species will be installed.
   * @param rule  Parameter that handles what to do with species
   *              who have elements that aren't declared.
   *              Check that all elements in the species
   *              exist in 'p'. If rule != 0, quietly skip
   *              this species if some elements are undeclared;
   *              otherwise, throw an exception
   * @param phaseNode_ptr Pointer to the XML_Node for this phase
   *              (defaults to 0)
   * @param factory Pointer to the SpeciesThermoFactory .
   *              (defaults to 0)
   *
   * @return
   *  Returns true if everything is ok, false otherwise.
   */
  bool installSpecies(int k, const XML_Node& s, thermo_t& p, 
		      SpeciesThermo *spthermo_ptr, int rule, 
		      XML_Node *phaseNode_ptr,
		      VPSSMgr *vpss_ptr,
		      SpeciesThermoFactory* factory) {

    std::string xname = s.name();
    if (xname != "species") {
      throw CanteraError("installSpecies",
			 "Unexpected XML name of species XML_Node: " + xname);
    }
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

    if (vpss_ptr) {
      VPStandardStateTP *vp_ptr = dynamic_cast<VPStandardStateTP *>(&p);
      factory->installVPThermoForSpecies(k, s, vp_ptr, vpss_ptr, spthermo_ptr,
					 phaseNode_ptr);
    } else {
      // install the thermo parameterization for this species into
      // the species thermo manager for phase p.
      factory->installThermoForSpecies(k, s, *spthermo_ptr, phaseNode_ptr);
    }
    
    
    return true;
  }


  //  Search an XML tree for species data.
  /*
   *   This utility routine will search the XML tree for the species
   *   named by the string, kname. It will return the XML_Node
   *   pointer to the species data for that species.
   *   Failures of any kind return the null pointer.
   *
   * @param kname String containing the name of the species.
   * @param phaseSpeciesData   Pointer to the XML speciesData element
   *              containing the species data for that phase.
   *
   */
  const XML_Node *speciesXML_Node(std::string kname,
				  const XML_Node *phaseSpeciesData) {
    if (!phaseSpeciesData) return ((const XML_Node *) 0);
    string jname = phaseSpeciesData->name();
    if (jname != "speciesData") {
      throw CanteraError("speciesXML_Node()",
			 "Unexpected phaseSpeciesData name: " + jname);
    }
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
