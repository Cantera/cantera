/**
 *  @file ThermoFactory.cpp
 *     Definitions for the factory class that can create known %ThermoPhase objects
 *     (see \ref thermoprops and class \link Cantera::ThermoFactory ThermoFactory\endlink).
 *
 */
// Copyright 2001  California Institute of Technology

#include "cantera/thermo/ThermoFactory.h"

#include "cantera/thermo/speciesThermoTypes.h"
#include "cantera/thermo/SpeciesThermoFactory.h"
#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/thermo/VPSSMgr.h"
#include "VPSSMgrFactory.h"

#include "cantera/thermo/IdealSolidSolnPhase.h"
#include "cantera/thermo/MargulesVPSSTP.h"
#include "cantera/thermo/RedlichKisterVPSSTP.h"
#include "cantera/thermo/IonsFromNeutralVPSSTP.h"
#include "cantera/thermo/PhaseCombo_Interaction.h"

#include "cantera/thermo/PureFluidPhase.h"
#include "cantera/thermo/RedlichKwongMFTP.h"

#include "cantera/thermo/ConstDensityThermo.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/thermo/EdgePhase.h"

#include "cantera/thermo/MetalPhase.h"
#include "cantera/thermo/SemiconductorPhase.h"

#undef USE_SSTP
#ifdef USE_SSTP
#include "cantera/thermo/StoichSubstanceSSTP.h"
#else
#include "cantera/thermo/StoichSubstance.h"
#endif

#include "cantera/thermo/MineralEQ3.h"
#include "cantera/thermo/MetalSHEelectrons.h"
#include "cantera/thermo/FixedChemPotSSTP.h"

#include "cantera/thermo/LatticeSolidPhase.h"
#include "cantera/thermo/LatticePhase.h"

#include "cantera/thermo/HMWSoln.h"
#include "cantera/thermo/DebyeHuckel.h"
#include "cantera/thermo/IdealMolalSoln.h"
#include "cantera/thermo/MolarityIonicVPSSTP.h"
#include "cantera/thermo/MixedSolventElectrolyte.h"
#include "cantera/thermo/IdealSolnGasVPSS.h"

#include "cantera/base/stringUtils.h"

using namespace std;
using namespace ctml;

namespace Cantera
{

ThermoFactory* ThermoFactory::s_factory = 0;
mutex_t ThermoFactory::thermo_mutex;

//! Define the number of %ThermoPhase types for use in this factory routine
static int ntypes = 26;

//! Define the string name of the %ThermoPhase types that are handled by this factory routine
static string _types[] = {"IdealGas", "Incompressible",
                          "Surface", "Edge", "Metal", "StoichSubstance",
                          "PureFluid", "LatticeSolid", "Lattice",
                          "HMW", "IdealSolidSolution", "DebyeHuckel",
                          "IdealMolalSolution", "IdealGasVPSS", "IdealSolnVPSS",
                          "MineralEQ3", "MetalSHEelectrons", "Margules", "PhaseCombo_Interaction",
                          "IonsFromNeutralMolecule", "FixedChemPot", "MolarityIonicVPSSTP",
                          "MixedSolventElectrolyte", "Redlich-Kister", "RedlichKwong",
                          "RedlichKwongMFTP"
                         };

//! Define the integer id of the %ThermoPhase types that are handled by this factory routine
static int _itypes[]   = {cIdealGas, cIncompressible,
                          cSurf, cEdge, cMetal, cStoichSubstance,
                          cPureFluid, cLatticeSolid, cLattice,
                          cHMW, cIdealSolidSolnPhase, cDebyeHuckel,
                          cIdealMolalSoln, cVPSS_IdealGas, cIdealSolnGasVPSS_iscv,
                          cMineralEQ3, cMetalSHEelectrons,
                          cMargulesVPSSTP,  cPhaseCombo_Interaction, cIonsFromNeutral, cFixedChemPot,
                          cMolarityIonicVPSSTP, cMixedSolventElectrolyte, cRedlichKisterVPSSTP,
                          cRedlichKwongMFTP, cRedlichKwongMFTP
                         };

ThermoPhase* ThermoFactory::newThermoPhase(const std::string& model)
{

    int ieos=-1;

    for (int n = 0; n < ntypes; n++) {
        if (model == _types[n]) {
            ieos = _itypes[n];
        }
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

    case cIdealSolidSolnPhase:
        th = new IdealSolidSolnPhase();
        break;

    case cMargulesVPSSTP:
        th = new MargulesVPSSTP();
        break;

    case cRedlichKisterVPSSTP:
        th = new RedlichKisterVPSSTP();
        break;

    case cMolarityIonicVPSSTP:
        th = new MolarityIonicVPSSTP();
        break;

    case cPhaseCombo_Interaction:
        th = new PhaseCombo_Interaction();
        break;

    case cIonsFromNeutral:
        th = new IonsFromNeutralVPSSTP();
        break;

    case cMetal:
        th = new MetalPhase;
        break;

    case cStoichSubstance:
#ifdef USE_SSTP
        th = new StoichSubstanceSSTP;
#else
        th = new StoichSubstance;
#endif
        break;

    case cFixedChemPot:
        th = new FixedChemPotSSTP;
        break;

    case cMineralEQ3:
        th = new MineralEQ3();
        break;

    case cMetalSHEelectrons:
        th = new MetalSHEelectrons();
        break;

    case cLatticeSolid:
        th = new LatticeSolidPhase;
        break;

    case cLattice:
        th = new LatticePhase;
        break;

    case cPureFluid:
        th = new PureFluidPhase;
        break;

    case cRedlichKwongMFTP:
        th = new RedlichKwongMFTP;
        break;

    case cHMW:
        th = new HMWSoln;
        break;

    case cDebyeHuckel:
        th = new DebyeHuckel;
        break;

    case cIdealMolalSoln:
        th = new IdealMolalSoln;
        break;

    case cVPSS_IdealGas:
        th = new IdealSolnGasVPSS;
        break;

    case cIdealSolnGasVPSS_iscv:
        th = new IdealSolnGasVPSS;
        break;

    default:
        throw UnknownThermoPhaseModel("ThermoFactory::newThermoPhase",
                                      model);
    }
    return th;
}

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

ThermoPhase* newPhase(XML_Node& xmlphase)
{
    const XML_Node& th = xmlphase.child("thermo");
    string model = th["model"];
    ThermoPhase* t = newThermoPhase(model);
    if (model == "singing cows") {
        throw CanteraError("ThermoPhase::newPhase", "Cows don't sing");
    } else if (model == "HMW") {
        HMWSoln* p = dynamic_cast<HMWSoln*>(t);
        p->constructPhaseXML(xmlphase,"");
    } else if (model == "IonsFromNeutralMolecule") {
        IonsFromNeutralVPSSTP* p = dynamic_cast<IonsFromNeutralVPSSTP*>(t);
        p->constructPhaseXML(xmlphase,"");
    } else {
        importPhase(xmlphase, t);
    }
    //return t;
    //importPhase(xmlphase, t);
    return t;
}

ThermoPhase* newPhase(const std::string& infile, std::string id)
{
    XML_Node* root = get_XML_File(infile);
    if (id == "-") {
        id = "";
    }
    XML_Node* xphase = get_XML_NameID("phase", std::string("#")+id, root);
    if (!xphase) {
        throw CanteraError("newPhase",
                           "Couldn't find phase named \"" + id + "\" in file, " + infile);
    }
    if (xphase) {
        return newPhase(*xphase);
    } else {
        return (ThermoPhase*) 0;
    }
}

//====================================================================================================================
//!  Gather a vector of pointers to XML_Nodes for a phase
/*!
 *   @param spDataNodeList   Output vector of pointer to XML_Nodes which contain the species XML_Nodes for the
 *                           species in the current phase.
 *   @param spNamesList      Output Vector of strings, which contain the names of the species in the phase
 *   @param spRuleList       Output Vector of ints, which contain the value of sprule for each species in the phase
 *   @param spArray_names    Vector of pointers to the XML_Nodes which contains the names of the
 *                           species in the phase
 *   @param spArray_dbases   Input vector of pointers to species data bases.
 *                           We search each data base for the required species names
 *   @param  sprule          Input vector of sprule values
 */
static void formSpeciesXMLNodeList(std::vector<XML_Node*> &spDataNodeList,
                                   std::vector<std::string> &spNamesList,
                                   std::vector<int> &spRuleList,
                                   const std::vector<XML_Node*> spArray_names,
                                   const std::vector<XML_Node*> spArray_dbases,
                                   const vector_int sprule)
{

    // used to check that each species is declared only once
    std::map<std::string, bool> declared;

    size_t nSpecies = 0;
    bool skip;

    for (size_t jsp = 0; jsp < spArray_dbases.size(); jsp++) {
        const XML_Node& speciesArray = *spArray_names[jsp];

        // Get the top XML for the database
        const XML_Node* db = spArray_dbases[jsp];

        // Get the array of species name strings and then count them
        std::vector<std::string> spnames;
        getStringArray(speciesArray, spnames);
        size_t nsp = spnames.size();

        // if 'all' is specified as the one and only species in the
        // spArray_names field, then add all species
        // defined in the corresponding database to the phase
        if (nsp == 1 && spnames[0] == "all") {
            std::vector<XML_Node*> allsp;
            db->getChildren("species", allsp);
            nsp = allsp.size();
            spnames.resize(nsp);
            for (size_t nn = 0; nn < nsp; nn++) {
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
        } else if (nsp == 1 && spnames[0] == "unique") {
            std::vector<XML_Node*> allsp;
            db->getChildren("species", allsp);
            nsp = allsp.size();
            spnames.resize(nsp);
            for (size_t nn = 0; nn < nsp; nn++) {
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
            std::map<std::string, XML_Node*> speciesNodes;
            for (size_t k = 0; k < db->nChildren(); k++) {
                XML_Node& child = db->child(k);
                speciesNodes[child["name"]] = &child;
            }
            for (size_t k = 0; k < nsp; k++) {
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
                    std::map<std::string, XML_Node*>::iterator iter = speciesNodes.find(stemp);
                    if (iter == speciesNodes.end()) {
                        throw CanteraError("importPhase","no data for species, \""
                                           + stemp + "\"");
                    }
                    XML_Node* s = iter->second;
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

bool importPhase(XML_Node& phase, ThermoPhase* th,
                 SpeciesThermoFactory* spfactory)
{

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
    XML_Node& phaseNode_XML = th->xml();
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
    } else {
        th->setNDim(3);     // default
    }

    // Set equation of state parameters. The parameters are
    // specific to each subclass of ThermoPhase, so this is done
    // by method setParametersFromXML in each subclass.
    const XML_Node& eos = phase.child("thermo");
    if (phase.hasChild("thermo")) {
        th->setParametersFromXML(eos);
    } else {
        throw CanteraError("importPhase",
                           " phase, " + th->id() +
                           ", XML_Node does not have a \"thermo\" XML_Node");
    }

    VPStandardStateTP* vpss_ptr = 0;
    int ssConvention = th->standardStateConvention();
    if (ssConvention == cSS_CONVENTION_VPSS) {
        vpss_ptr = dynamic_cast <VPStandardStateTP*>(th);
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
    if (ssConvention != cSS_CONVENTION_SLAVE) {
        th->addElementsFromXML(phase);
    }

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
    if (ssConvention != cSS_CONVENTION_SLAVE) {
        if (nspa == 0) {
            throw CanteraError("importPhase",
                               "phase, " + th->id() + ", has zero \"speciesArray\" XML nodes.\n"
                               + " There must be at least one speciesArray nodes "
                               "with one or more species");
        }
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
            throw CanteraError("importPhase()",
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
    std::vector<XML_Node*>  spDataNodeList;
    std::vector<std::string> spNamesList;
    std::vector<int> spRuleList;
    formSpeciesXMLNodeList(spDataNodeList, spNamesList, spRuleList,
                           sparrays, dbases, sprule);

    // If the phase has a species thermo manager already installed,
    // delete it since we are adding new species.
    //delete &th->speciesThermo();

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
        if (eos["allow_discontinuities"] == "true") {
            std::cout << "ALLOWING DISCONTINUOUS THERMO!" << std::endl;
            spth->m_allow_discontinuities = true;
        }

        // install it in the phase object
        th->setSpeciesThermo(spth);
    } else if (ssConvention == cSS_CONVENTION_SLAVE) {
        /*
         * No species thermo manager for this type
         */
    } else if (ssConvention == cSS_CONVENTION_VPSS) {
        vp_spth = newVPSSMgr(vpss_ptr, &phase, spDataNodeList);
        vpss_ptr->setVPSSMgr(vp_spth);
        spth = vp_spth->SpeciesThermoMgr();
        th->setSpeciesThermo(spth);
    } else {
        throw CanteraError("importPhase()", "unknown convention");
    }


    size_t k = 0;

    size_t nsp = spDataNodeList.size();
    if (ssConvention == cSS_CONVENTION_SLAVE) {
        if (nsp > 0) {
            throw CanteraError("importPhase()", "For Slave standard states, number of species must be zero: "
                               + int2str(nsp));
        }
    }
    for (size_t i = 0; i < nsp; i++) {
        XML_Node* s = spDataNodeList[i];
        AssertTrace(s != 0);
        bool ok = installSpecies(k, *s, *th, spth, spRuleList[i],
                                 &phase, vp_spth, spfactory);
        if (ok) {
            th->saveSpeciesData(k, s);
            ++k;
        }
    }

    if (ssConvention == cSS_CONVENTION_SLAVE) {
        th->installSlavePhases(&phase);
    }

    // done adding species.
    th->freezeSpecies();

    // Perform any required subclass-specific initialization.
    th->initThermo();

    // Perform any required subclass-specific initialization
    // that requires the XML phase object
    std::string id = "";
    th->initThermoXML(phase, id);

    return true;
}

bool installSpecies(size_t k, const XML_Node& s, thermo_t& th,
                    SpeciesThermo* spthermo_ptr, int rule,
                    XML_Node* phaseNode_ptr,
                    VPSSMgr* vpss_ptr,
                    SpeciesThermoFactory* factory)
{
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
        if (th.elementIndex(_b->first) == npos) {
            if (rule == 0) {
                throw CanteraError("installSpecies",
                                   "Species " + s["name"] +
                                   " contains undeclared element " + _b->first);
            } else {
                return false;
            }
        }
    }

    // construct a vector of atom numbers for each
    // element in phase th. Elements not declared in the
    // species (i.e., not in map comp) will have zero
    // entries in the vector.
    size_t nel = th.nElements();
    vector_fp ecomp(nel, 0.0);
    for (size_t m = 0; m < nel; m++) {
        std::string& es = comp[th.elementName(m)];
        if (!es.empty()) {
            ecomp[m] = fpValueCheck(es);
        }
    }


    // get the species charge, if any. Note that the charge need
    // not be explicitly specified if special element 'E'
    // (electron) is one of the elements.
    doublereal chrg = 0.0;
    if (s.hasChild("charge")) {
        chrg = getFloat(s, "charge");
    }

    // get the species size, if any. (This is used by surface
    // phases to represent how many sites a species occupies.)
    doublereal sz = 1.0;
    if (s.hasChild("size")) {
        sz = getFloat(s, "size");
    }

    // add the species to phase th
    th.addUniqueSpecies(s["name"], &ecomp[0], chrg, sz);

    if (vpss_ptr) {
        VPStandardStateTP* vp_ptr = dynamic_cast<VPStandardStateTP*>(&th);
        factory->installVPThermoForSpecies(k, s, vp_ptr, vpss_ptr, spthermo_ptr,
                                           phaseNode_ptr);
    } else {
        // install the thermo parameterization for this species into
        // the species thermo manager for phase th
        factory->installThermoForSpecies(k, s, &th, *spthermo_ptr, phaseNode_ptr);
    }

    return true;
}

const XML_Node* speciesXML_Node(const std::string& kname,
                                const XML_Node* phaseSpeciesData)
{
    if (!phaseSpeciesData) {
        return 0;
    }
    string jname = phaseSpeciesData->name();
    if (jname != "speciesData") {
        throw CanteraError("speciesXML_Node()",
                           "Unexpected phaseSpeciesData name: " + jname);
    }
    vector<XML_Node*> xspecies;
    phaseSpeciesData->getChildren("species", xspecies);
    for (size_t j = 0; j < xspecies.size(); j++) {
        const XML_Node& sp = *xspecies[j];
        jname = sp["name"];
        if (jname == kname) {
            return &sp;
        }
    }
    return 0;
}

}
