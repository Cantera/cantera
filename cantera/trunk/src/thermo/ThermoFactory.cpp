/**
 *  @file ThermoFactory.cpp
 *     Definitions for the factory class that can create known ThermoPhase objects
 *     (see \ref thermoprops and class \link Cantera::ThermoFactory ThermoFactory\endlink).
 */
// Copyright 2001  California Institute of Technology

#include "cantera/thermo/ThermoFactory.h"

#include "cantera/thermo/speciesThermoTypes.h"
#include "cantera/thermo/SpeciesThermoFactory.h"
#include "cantera/thermo/GeneralSpeciesThermo.h"
#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/thermo/VPSSMgr.h"
#include "VPSSMgrFactory.h"

#include "cantera/thermo/IdealSolidSolnPhase.h"
#include "cantera/thermo/MaskellSolidSolnPhase.h"
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

#include "cantera/transport/TransportData.h"

#include "cantera/base/stringUtils.h"

using namespace std;
using namespace ctml;

namespace Cantera
{

ThermoFactory* ThermoFactory::s_factory = 0;
mutex_t ThermoFactory::thermo_mutex;

//! Define the number of ThermoPhase types for use in this factory routine
static int ntypes = 27;

//! Define the string name of the ThermoPhase types that are handled by this factory routine
static string _types[] = {"IdealGas", "Incompressible",
                          "Surface", "Edge", "Metal", "StoichSubstance",
                          "PureFluid", "LatticeSolid", "Lattice",
                          "HMW", "IdealSolidSolution", "DebyeHuckel",
                          "IdealMolalSolution", "IdealGasVPSS", "IdealSolnVPSS",
                          "MineralEQ3", "MetalSHEelectrons", "Margules", "PhaseCombo_Interaction",
                          "IonsFromNeutralMolecule", "FixedChemPot", "MolarityIonicVPSSTP",
                          "MixedSolventElectrolyte", "Redlich-Kister", "RedlichKwong",
                          "RedlichKwongMFTP", "MaskellSolidSolnPhase"
                         };

//! Define the integer id of the ThermoPhase types that are handled by this factory routine
static int _itypes[]   = {cIdealGas, cIncompressible,
                          cSurf, cEdge, cMetal, cStoichSubstance,
                          cPureFluid, cLatticeSolid, cLattice,
                          cHMW, cIdealSolidSolnPhase, cDebyeHuckel,
                          cIdealMolalSoln, cVPSS_IdealGas, cIdealSolnGasVPSS_iscv,
                          cMineralEQ3, cMetalSHEelectrons,
                          cMargulesVPSSTP,  cPhaseCombo_Interaction, cIonsFromNeutral, cFixedChemPot,
                          cMolarityIonicVPSSTP, cMixedSolventElectrolyte, cRedlichKisterVPSSTP,
                          cRedlichKwongMFTP, cRedlichKwongMFTP, cMaskellSolidSolnPhase
                         };

ThermoPhase* ThermoFactory::newThermoPhase(const std::string& model)
{
    int ieos=-1;

    for (int n = 0; n < ntypes; n++) {
        if (model == _types[n]) {
            ieos = _itypes[n];
            break;
        }
    }

    switch (ieos) {

    case cIdealGas:
        return new IdealGasPhase;
    case cIncompressible:
        return new ConstDensityThermo;
    case cSurf:
        return new SurfPhase;
    case cEdge:
        return new EdgePhase;
    case cIdealSolidSolnPhase:
        return new IdealSolidSolnPhase();
    case cMargulesVPSSTP:
        return new MargulesVPSSTP();
    case cRedlichKisterVPSSTP:
        return new RedlichKisterVPSSTP();
    case cMolarityIonicVPSSTP:
        return new MolarityIonicVPSSTP();
    case cPhaseCombo_Interaction:
        return new PhaseCombo_Interaction();
    case cIonsFromNeutral:
        return new IonsFromNeutralVPSSTP();
    case cMetal:
        return new MetalPhase;
    case cStoichSubstance:
#ifdef USE_SSTP
        return new StoichSubstanceSSTP;
#else
        return new StoichSubstance;
#endif
    case cFixedChemPot:
        return new FixedChemPotSSTP;
    case cMineralEQ3:
        return new MineralEQ3();
    case cMetalSHEelectrons:
        return new MetalSHEelectrons();
    case cLatticeSolid:
        return new LatticeSolidPhase;
    case cLattice:
        return new LatticePhase;
    case cPureFluid:
        return new PureFluidPhase;
    case cRedlichKwongMFTP:
        return new RedlichKwongMFTP;
    case cHMW:
        return new HMWSoln;
    case cDebyeHuckel:
        return new DebyeHuckel;
    case cIdealMolalSoln:
        return new IdealMolalSoln;
    case cVPSS_IdealGas:
        return new IdealSolnGasVPSS;
    case cIdealSolnGasVPSS_iscv:
        return new IdealSolnGasVPSS;
    case cMaskellSolidSolnPhase:
        return new MaskellSolidSolnPhase;
    default:
        throw UnknownThermoPhaseModel("ThermoFactory::newThermoPhase", model);
    }
}

std::string eosTypeString(int ieos, int length)
{
    for (int n = 0; n < ntypes; n++) {
        if (_itypes[n] == ieos) {
            return _types[n];
        }
    }
    return "UnknownPhaseType";
}

ThermoPhase* newPhase(XML_Node& xmlphase)
{
    string model = xmlphase.child("thermo")["model"];
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
    return t;
}

ThermoPhase* newPhase(const std::string& infile, std::string id)
{
    XML_Node* root = get_XML_File(infile);
    if (id == "-") {
        id = "";
    }
    XML_Node* xphase = get_XML_NameID("phase", "#"+id, root);
    if (!xphase) {
        throw CanteraError("newPhase",
                           "Couldn't find phase named \"" + id + "\" in file, " + infile);
    }
    return newPhase(*xphase);
}

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
            std::vector<XML_Node*> allsp = db->getChildren("species");
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
                    spNamesList.push_back(stemp);
                    spDataNodeList.push_back(allsp[nn]);
                    spRuleList.push_back(sprule[jsp]);
                }
            }
        } else if (nsp == 1 && spnames[0] == "unique") {
            std::vector<XML_Node*> allsp = db->getChildren("species");
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
                    spNamesList.push_back(stemp);
                    spDataNodeList.push_back(allsp[nn]);
                    spRuleList.push_back(sprule[jsp]);
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
                    // Find the species in the database by name.
                    std::map<std::string, XML_Node*>::iterator iter = speciesNodes.find(stemp);
                    if (iter == speciesNodes.end()) {
                        throw CanteraError("importPhase","no data for species, \""
                                           + stemp + "\"");
                    }
                    spNamesList.push_back(stemp);
                    spDataNodeList.push_back(iter->second);
                    spRuleList.push_back(sprule[jsp]);
                }
            }
        }
    }
}

bool importPhase(XML_Node& phase, ThermoPhase* th,
                 SpeciesThermoFactory* spfactory)
{
    // Check the the supplied XML node in fact represents a phase.
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
    th->setXMLdata(phase);

    // set the id attribute of the phase to the 'id' attribute in the XML tree.
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

    /***************************************************************
     * Add the elements.
     ***************************************************************/
    if (ssConvention != cSS_CONVENTION_SLAVE) {
        installElements(*th, phase);
    }

    /***************************************************************
     * Add the species.
     *
     * Species definitions may be imported from multiple
     * sources. For each one, a speciesArray element must be
     * present.
     ***************************************************************/
    vector<XML_Node*> sparrays = phase.getChildren("speciesArray");
    if (ssConvention != cSS_CONVENTION_SLAVE) {
        if (sparrays.empty()) {
            throw CanteraError("importPhase",
                               "phase, " + th->id() + ", has zero \"speciesArray\" XML nodes.\n"
                               + " There must be at least one speciesArray nodes "
                               "with one or more species");
        }
    }
    vector<XML_Node*> dbases;
    vector_int sprule(sparrays.size(),0);

    // loop over the speciesArray elements
    for (size_t jsp = 0; jsp < sparrays.size(); jsp++) {

        const XML_Node& speciesArray = *sparrays[jsp];

        // If the speciesArray element has a child element
        //
        //   <skip element="undeclared">
        //
        // then set sprule[jsp] to 1, so that any species with an undeclared
        // element will be quietly skipped when importing species. Additionally,
        // if the skip node has the following attribute:
        //
        // <skip species="duplicate">
        //
        // then duplicate species names will not cause Cantera to throw an
        // exception. Instead, the duplicate entry will be discarded.
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

        // Get a pointer to the node containing the species
        // definitions for the species declared in this
        // speciesArray element. This may be in the local file
        // containing the phase element, or may be in another
        // file.
        XML_Node* db = get_XML_Node(speciesArray["datasrc"], &phase.root());
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

    // Decide whether the the phase has a variable pressure ss or not
    SpeciesThermo* spth = &th->speciesThermo();
    VPSSMgr* vp_spth = 0;
    if (ssConvention == cSS_CONVENTION_VPSS) {
        vp_spth = newVPSSMgr(vpss_ptr, &phase, spDataNodeList);
        vpss_ptr->setVPSSMgr(vp_spth);
        spth = vp_spth->SpeciesThermoMgr();
        th->setSpeciesThermo(spth);
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

    // Done adding species. Perform any required subclass-specific
    // initialization.
    th->initThermo();

    // Perform any required subclass-specific initialization
    // that requires the XML phase object
    std::string id = "";
    th->initThermoXML(phase, id);

    return true;
}

void installElements(Phase& th, const XML_Node& phaseNode)
{
    // get the declared element names
    if (!phaseNode.hasChild("elementArray")) {
        throw CanteraError("installElements",
                           "phase xml node doesn't have \"elementArray\" XML Node");
    }
    XML_Node& elements = phaseNode.child("elementArray");
    vector<string> enames;
    getStringArray(elements, enames);

    // // element database defaults to elements.xml
    string element_database = "elements.xml";
    if (elements.hasAttrib("datasrc")) {
        element_database = elements["datasrc"];
    }

    XML_Node* doc = get_XML_File(element_database);
    XML_Node* dbe = &doc->child("ctml/elementData");

    XML_Node& root = phaseNode.root();
    XML_Node* local_db = 0;
    if (root.hasChild("ctml")) {
        if (root.child("ctml").hasChild("elementData")) {
            local_db = &root.child("ctml/elementData");
        }
    }

    for (size_t i = 0; i < enames.size(); i++) {
        // Find the element data
        XML_Node* e = 0;
        if (local_db) {
            e = local_db->findByAttr("name",enames[i]);
        }
        if (!e) {
            e = dbe->findByAttr("name",enames[i]);
        }
        if (!e) {
            throw CanteraError("addElementsFromXML","no data for element "
                               +enames[i]);
        }

        // Add the element
        doublereal weight = 0.0;
        if (e->hasAttrib("atomicWt")) {
            weight = fpValue(e->attrib("atomicWt"));
        }
        int anum = 0;
        if (e->hasAttrib("atomicNumber")) {
            anum = intValue(e->attrib("atomicNumber"));
        }
        string symbol = e->attrib("name");
        doublereal entropy298 = ENTROPY298_UNKNOWN;
        if (e->hasChild("entropy298")) {
            XML_Node& e298Node = e->child("entropy298");
            if (e298Node.hasAttrib("value")) {
                entropy298 = fpValueCheck(e298Node["value"]);
            }
        }
        if (weight != 0.0) {
            th.addElement(symbol, weight, anum, entropy298);
        } else {
            th.addElement(symbol);
        }
    }
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

    if (rule) {
        th.ignoreUndefinedElements();
    }

    // get the composition of the species
    const XML_Node& a = s.child("atomArray");
    map<string,string> comp;
    getMap(a, comp);

    // construct a vector of atom numbers for each element in phase th. Elements
    // not declared in the species (i.e., not in map comp) will have zero
    // entries in the vector.
    size_t nel = th.nElements();
    vector_fp ecomp(nel, 0.0);
    compositionMap comp_map = parseCompString(a.value());
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

    if (vpss_ptr) {
        th.addUniqueSpecies(s["name"], &ecomp[0], chrg, sz);
        VPStandardStateTP* vp_ptr = dynamic_cast<VPStandardStateTP*>(&th);
        vp_ptr->createInstallPDSS(k, s, phaseNode_ptr);
    } else {
        SpeciesThermoInterpType* st = newSpeciesThermoInterpType(s);
        Species sp(s["name"], comp_map, st, chrg, sz);

        // Read gas-phase transport data, if provided
        if (s.hasChild("transport") &&
                s.child("transport")["model"] == "gas_transport") {
            XML_Node& tr = s.child("transport");

            string geometry, dummy;
            getString(tr, "geometry", geometry, dummy);

            double diam = getFloat(tr, "LJ_diameter");
            double welldepth = getFloat(tr, "LJ_welldepth");

            double dipole = 0.0;
            getOptionalFloat(tr, "dipoleMoment", dipole);

            double polar = 0.0;
            getOptionalFloat(tr, "polarizability", polar);

            double rot = 0.0;
            getOptionalFloat(tr, "rotRelax", rot);
            double acentric = 0.0;
            getOptionalFloat(tr, "acentric_factor", acentric);

            GasTransportData* gastran = new GasTransportData;
            gastran->setCustomaryUnits(sp.name, geometry, diam, welldepth,
                                   dipole, polar, rot, acentric);
            sp.transport.reset(gastran);
            gastran->validate(sp);
        }
        th.addSpecies(sp);
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
    vector<XML_Node*> xspecies = phaseSpeciesData->getChildren("species");
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
