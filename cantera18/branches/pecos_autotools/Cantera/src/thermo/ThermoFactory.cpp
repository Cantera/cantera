/**
 *  @file ThermoFactory.cpp
 *     Definitions for the factory class that can create known %ThermoPhase objects
 *     (see \ref thermoprops and class \link Cantera::ThermoFactory ThermoFactory\endlink).
 *
 */

/*
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
#include "PerfectGasPhase.h"
#include "VPSSMgr.h"
#include "VPSSMgrFactory.h"

#ifdef WITH_IDEAL_SOLUTIONS
#include "IdealSolidSolnPhase.h"
#include "MargulesVPSSTP.h"
#include "IonsFromNeutralVPSSTP.h"
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

#ifdef WITH_STOICH_SUBSTANCE
#include "MineralEQ3.h"
#include "MetalSHEelectrons.h"
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

#include <cstdlib>

using namespace std;

namespace Cantera {

    ThermoFactory* ThermoFactory::s_factory = 0;
#if defined(THREAD_SAFE_CANTERA)
    boost::mutex ThermoFactory::thermo_mutex;
#endif

    static int ntypes = 18;
    static string _types[] = {"IdealGas", "Incompressible", 
                              "Surface", "Edge", "Metal", "StoichSubstance",
                              "PureFluid", "LatticeSolid", "Lattice",
                              "HMW", "IdealSolidSolution", "DebyeHuckel", 
                              "IdealMolalSolution", "IdealGasVPSS",
			      "MineralEQ3", "MetalSHEelectrons", "Margules",
                              "IonsFromNeutralMolecule","PerfectGas"
    };

    static int _itypes[]   = {cIdealGas, cIncompressible, 
                              cSurf, cEdge, cMetal, cStoichSubstance,
                              cPureFluid, cLatticeSolid, cLattice,
                              cHMW, cIdealSolidSolnPhase, cDebyeHuckel,
                              cIdealMolalSoln, cVPSS_IdealGas,
			      cMineralEQ3, cMetalSHEelectrons,
			      cMargulesVPSSTP, cIonsFromNeutral, 
			      cPerfectGas
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

      // adding perfect gas
    case cPerfectGas:
      th = new PerfectGasPhase;
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

    case cMargulesVPSSTP:
      th = new MargulesVPSSTP();
      break;

    case cIonsFromNeutral:
      th = new IonsFromNeutralVPSSTP();
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

#ifdef WITH_STOICH_SUBSTANCE
    case cMineralEQ3:
      th = new MineralEQ3();
      break;
#endif

#ifdef WITH_STOICH_SUBSTANCE
    case cMetalSHEelectrons:
      th = new MetalSHEelectrons();
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

  // Translate the eosType id into a string
  /*
   *  Returns a string representation of the eosType id for a phase.
   *  @param ieos  eosType id of the phase. This is unique for the phase
   *  @param length maximum length of the return string. Defaults to 100
   *
   *  @return returns a string representation.
   */
  std::string eosTypeString(int ieos, int length)
  {
    std::string ss = "UnknownPhaseType";
    // bool found = false;
    for (int n = 0; n <  ntypes; n++) {
      if (_itypes[n] == ieos) {
	ss = _types[n];
	//found = true;
      }  
    }
    return ss;
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
    if (model == "singing cows") {
      throw CanteraError(" newPhase", "Cows don't sing");
    } 
#ifdef WITH_ELECTROLYTES
    else if (model == "HMW") {
      HMWSoln* p = dynamic_cast<HMWSoln*>(t);
      p->constructPhaseXML(xmlphase,"");
    }
#endif
#ifdef WITH_IDEAL_SOLUTIONS
    else if (model == "IonsFromNeutralMolecule") {
      IonsFromNeutralVPSSTP* p = dynamic_cast<IonsFromNeutralVPSSTP*>(t);
      p->constructPhaseXML(xmlphase,"");
    }
#endif
    else {
      importPhase(xmlphase, t);
    }
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


  static void formSpeciesXMLNodeList(std::vector<XML_Node *> &spDataNodeList,
				     std::vector<std::string> &spNamesList,
				     std::vector<int> &spRuleList,
				     const std::vector<XML_Node *> spArray_names,
				     const std::vector<XML_Node *> spArray_dbases,
				     const vector_int sprule) {
    
    // used to check that each species is declared only once
    std::map<std::string, bool> declared;
    
    int nspa = spArray_dbases.size();
    int nSpecies = 0;
    bool skip;

    for (int jsp = 0; jsp < nspa; jsp++) {
      const XML_Node& speciesArray = *spArray_names[jsp]; 
      
      // Get the top XML for the database
      const XML_Node *db = spArray_dbases[jsp];

      // Get the array of species name strings and the count them
      std::vector<std::string> spnames;
      getStringArray(speciesArray, spnames);
      int nsp = static_cast<int>(spnames.size());

      // if 'all' is specified as the one and only species in the
      // spArray_names field, then add all species 
      // defined in the corresponding database to the phase
      if (nsp == 1 && spnames[0] == "all") {
	std::vector<XML_Node *> allsp;
	db->getChildren("species", allsp);
	nsp = static_cast<int>(allsp.size());
	spnames.resize(nsp);
	for (int nn = 0; nn < nsp; nn++) {
	  string stemp = (*allsp[nn])["name"];
	  bool skip = false;
	  if (declared[stemp]) {
	    if (sprule[jsp] >= 10) {
	      skip = true;
	    } else {
	      throw CanteraError("ThermoFactory::formSpeciesXMLNodeList()",
				 "duplicate species: \"" + stemp + "\"");
	    }
	  }
	  if (!skip) {
	    declared[stemp] = true;
	    nSpecies++;
	    spNamesList.resize(nSpecies);
	    spDataNodeList.resize(nSpecies, 0);
	    spRuleList.resize(nSpecies, 0);
	    spNamesList[nSpecies-1] = stemp;
	    spDataNodeList[nSpecies-1] = allsp[nn];
	    spRuleList[nSpecies-1] = sprule[jsp];
	  }
	}
      }
      else if (nsp == 1 && spnames[0] == "unique") {
	std::vector<XML_Node *> allsp;
	db->getChildren("species", allsp);
	nsp = static_cast<int>(allsp.size());
	spnames.resize(nsp);
	for (int nn = 0; nn < nsp; nn++) {
	  string stemp = (*allsp[nn])["name"];
	  bool skip = false;
	  if (declared[stemp]) {
	      skip = true;
	  }
	  if (!skip) {
	    declared[stemp] = true;
	    nSpecies++;
	    spNamesList.resize(nSpecies);
	    spDataNodeList.resize(nSpecies, 0);
	    spRuleList.resize(nSpecies, 0);
	    spNamesList[nSpecies-1] = stemp;
	    spDataNodeList[nSpecies-1] = allsp[nn];
	    spRuleList[nSpecies-1] = sprule[jsp];
	  }
	}
      } else {
	for (int k = 0; k < nsp; k++) {
	  string stemp = spnames[k];
	  skip = false;
	  if (declared[stemp]) {
	    if (sprule[jsp] >= 10) {
	      skip = true;
	    } else {
	      throw CanteraError("ThermoFactory::formSpeciesXMLNodeList()",
				 "duplicate species: \"" + stemp + "\"");
	    }
	  }
	  if (!skip) {
	    declared[stemp] = true;
	    // Find the species in the database by name.
	    XML_Node* s = db->findByAttr("name", stemp); 
	    if (!s) {
	      throw CanteraError("importPhase","no data for species, \""
				 + stemp + "\"");
	    }
	    nSpecies++;
	    spNamesList.resize(nSpecies);
	    spDataNodeList.resize(nSpecies, 0);
	    spRuleList.resize(nSpecies, 0);
	    spNamesList[nSpecies-1] = stemp;
	    spDataNodeList[nSpecies-1] = s;
	    spRuleList[nSpecies-1] = sprule[jsp];
	  }
	}
      }
    }
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
			 "Current const XML_Node named, " + phase.name() + 
			 ", is not a phase element.");
    }

    /*
     * In this section of code, we get the reference to the 
     * phase xml tree within the ThermoPhase object. Then,
     * we clear it and fill it with the current information that
     * we are about to use to construct the object. We will then
     * be able to resurrect the information later by calling xml().
     */
    XML_Node &phaseNode_XML = th->xml();
    phaseNode_XML.clear();
    phase.copy(&phaseNode_XML);

    // set the id attribute of the phase to the 'id' attribute 
    // in the XML tree.
    th->setID(phase.id());
    th->setName(phase.id());

    // Number of spatial dimensions. Defaults to 3 (bulk phase)
    if (phase.hasAttrib("dim")) {
      int idim = intValue(phase["dim"]);
      if (idim < 1 || idim > 3)
	throw CanteraError("importPhase",
			   "phase, " + th->id() + 
			   ", has unphysical number of dimensions: " + phase["dim"]);
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
			 " phase, " + th->id() + 
			 ", XML_Node does not have a \"thermo\" XML_Node");
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
			 + " There must be at least one speciesArray nodes "
			 "with one or more species");
    }
    vector<XML_Node*> dbases;
    vector_int sprule(nspa,0);

    // loop over the speciesArray elements
    for (jsp = 0; jsp < nspa; jsp++) {

      const XML_Node& speciesArray = *sparrays[jsp];

      // If the speciesArray element has a child element
      //
      //   <skip element="undeclared"> 
      //
      // then set sprule[jsp] to 1, so
      // that any species with an undeclared element will be
      // quietly skipped when importing species.
      // Additionally, if the skip node has the following attribute:
      //
      // <skip species="duplicate">
      //
      // then duplicate species names will not cause Cantera to
      // throw an exception. Instead, the duplicate entry will 
      // be discarded.
      if (speciesArray.hasChild("skip")) {
	const XML_Node& sk = speciesArray.child("skip");
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

      // Get a pointer to the node containing the species
      // definitions for the species declared in this 
      // speciesArray element. This may be in the local file
      // containing the phase element, or may be in another
      // file.            
      db = get_XML_Node(speciesArray["datasrc"], &phase.root());
      if (db == 0) {
	throw CanteraError("importPhase",
			   " Can not find XML node for species database: " 
			   + speciesArray["datasrc"]);
      }

      // add this node to the list of species database nodes.
      dbases.push_back(db);
    }

    // Now, collect all the species names and all the XML_Node * pointers
    // for those species in a single vector. This is where we decide what
    // species are to be included in the phase.
    // The logic is complicated enough that we put it in a separate routine.
    std::vector<XML_Node *>  spDataNodeList;
    std::vector<std::string> spNamesList;
    std::vector<int> spRuleList;
    formSpeciesXMLNodeList(spDataNodeList, spNamesList, spRuleList,
			   sparrays, dbases, sprule);

    // If the phase has a species thermo manager already installed,
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
      spth = newSpeciesThermoMgr(spDataNodeList);
      
      // install it in the phase object
      th->setSpeciesThermo(spth);
    } else {
      vp_spth = newVPSSMgr(vpss_ptr, &phase, spDataNodeList);
      vpss_ptr->setVPSSMgr(vp_spth);
      spth = vp_spth->SpeciesThermoMgr();
      th->setSpeciesThermo(spth);
    }


    int k = 0;

    int nsp = spDataNodeList.size();
    for (int i = 0; i < nsp; i++) {
      XML_Node *s = spDataNodeList[i];
      AssertTrace(s != 0);
      bool ok = installSpecies(k, *s, *th, spth, spRuleList[i], 
			       &phase, vp_spth, spfactory);
      if (ok) {
	th->saveSpeciesData(k, s);
	++k;
      }
    }

    // done adding species. 
    th->freezeSpecies();

    // Perform any required subclass-specific initialization.
    th->initThermo();

    // Perform any required subclass-specific initialization
    // that requires the XML phase object
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
  bool installSpecies(int k, const XML_Node& s, thermo_t& th, 
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
      if (th.elementIndex(_b->first) < 0) {
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
    // element in phase th. Elements not declared in the
    // species (i.e., not in map comp) will have zero
    // entries in the vector.
    int m, nel = th.nElements();
    vector_fp ecomp(nel, 0.0);            
    for (m = 0; m < nel; m++) {
      const char *es = comp[th.elementName(m)].c_str();
      if (strlen(es) > 0) {
        ecomp[m] = atofCheck(es);
      }
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

    // add the species to phase th 
    th.addUniqueSpecies(s["name"], &ecomp[0], chrg, sz);

    if (vpss_ptr) {
      VPStandardStateTP *vp_ptr = dynamic_cast<VPStandardStateTP *>(&th);
      factory->installVPThermoForSpecies(k, s, vp_ptr, vpss_ptr, spthermo_ptr,
					 phaseNode_ptr);
    } else {
      // install the thermo parameterization for this species into
      // the species thermo manager for phase th
      factory->installThermoForSpecies(k, s, &th, *spthermo_ptr, phaseNode_ptr);
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
