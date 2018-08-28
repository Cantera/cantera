/**
 *  @file ThermoFactory.cpp
 *     Definitions for the factory class that can create known ThermoPhase objects
 *     (see \ref thermoprops and class \link Cantera::ThermoFactory ThermoFactory\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/ThermoFactory.h"

#include "cantera/thermo/Species.h"
#include "cantera/thermo/speciesThermoTypes.h"
#include "cantera/thermo/SpeciesThermoFactory.h"
#include "cantera/thermo/PDSSFactory.h"
#include "cantera/thermo/MultiSpeciesThermo.h"
#include "cantera/thermo/IdealGasPhase.h"

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
#include "cantera/thermo/StoichSubstance.h"
#include "cantera/thermo/MineralEQ3.h"
#include "cantera/thermo/MetalSHEelectrons.h"
#include "cantera/thermo/FixedChemPotSSTP.h"
#include "cantera/thermo/LatticeSolidPhase.h"
#include "cantera/thermo/LatticePhase.h"
#include "cantera/thermo/HMWSoln.h"
#include "cantera/thermo/DebyeHuckel.h"
#include "cantera/thermo/IdealMolalSoln.h"
#include "cantera/thermo/MolarityIonicVPSSTP.h"
#include "cantera/thermo/IdealSolnGasVPSS.h"
#include "cantera/thermo/WaterSSTP.h"
#include "cantera/base/stringUtils.h"

using namespace std;

namespace Cantera
{

ThermoFactory* ThermoFactory::s_factory = 0;
std::mutex ThermoFactory::thermo_mutex;

ThermoFactory::ThermoFactory()
{
    reg("IdealGas", []() { return new IdealGasPhase(); });
    reg("Incompressible", []() { return new ConstDensityThermo(); });
    reg("Surface", []() { return new SurfPhase(); });
    reg("Edge", []() { return new EdgePhase(); });
    reg("Metal", []() { return new MetalPhase(); });
    reg("StoichSubstance", []() { return new StoichSubstance(); });
    reg("PureFluid", []() { return new PureFluidPhase(); });
    reg("LatticeSolid", []() { return new LatticeSolidPhase(); });
    reg("Lattice", []() { return new LatticePhase(); });
    reg("HMW", []() { return new HMWSoln(); });
    reg("IdealSolidSolution", []() { return new IdealSolidSolnPhase(); });
    reg("DebyeHuckel", []() { return new DebyeHuckel(); });
    reg("IdealMolalSolution", []() { return new IdealMolalSoln(); });
    reg("IdealGasVPSS", []() { return new IdealSolnGasVPSS(); });
    m_synonyms["IdealGasVPSS"] = "IdealSolnVPSS";
    reg("MineralEQ3", []() { return new MineralEQ3(); });
    reg("MetalSHEelectrons", []() { return new MetalSHEelectrons(); });
    reg("Margules", []() { return new MargulesVPSSTP(); });
    reg("PhaseCombo_Interaction", []() { return new PhaseCombo_Interaction(); });
    reg("IonsFromNeutralMolecule", []() { return new IonsFromNeutralVPSSTP(); });
    reg("FixedChemPot", []() { return new FixedChemPotSSTP(); });
    reg("MolarityIonicVPSSTP", []() { return new MolarityIonicVPSSTP(); });
    reg("Redlich-Kister", []() { return new RedlichKisterVPSSTP(); });
    reg("RedlichKwong", []() { return new RedlichKwongMFTP(); });
    m_synonyms["RedlichKwongMFTP"] = "RedlichKwong";
    reg("MaskellSolidSolnPhase", []() { return new MaskellSolidSolnPhase(); });
    reg("PureLiquidWater", []() { return new WaterSSTP(); });
}

ThermoPhase* ThermoFactory::newThermoPhase(const std::string& model)
{
    return create(model);
}

ThermoPhase* newPhase(XML_Node& xmlphase)
{
    string model = xmlphase.child("thermo")["model"];
    unique_ptr<ThermoPhase> t(newThermoPhase(model));
    importPhase(xmlphase, t.get());
    return t.release();
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
 *   @param spDataNodeList   Output vector of pointer to XML_Nodes which contain
 *       the species XML_Nodes for the species in the current phase.
 *   @param spNamesList      Output Vector of strings, which contain the names
 *       of the species in the phase
 *   @param spRuleList       Output Vector of ints, which contain the value of
 *       sprule for each species in the phase
 *   @param spArray_names    Vector of pointers to the XML_Nodes which contains
 *                           the names of the species in the phase
 *   @param spArray_dbases   Input vector of pointers to species data bases. We
 *                           search each data base for the required species
 *                           names
 *   @param  sprule          Input vector of sprule values
 */
static void formSpeciesXMLNodeList(std::vector<XML_Node*> &spDataNodeList,
                                   std::vector<std::string> &spNamesList,
                                   vector_int &spRuleList,
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
        // spArray_names field, then add all species defined in the
        // corresponding database to the phase
        if (nsp == 1 && spnames[0] == "all") {
            std::vector<XML_Node*> allsp = db->getChildren("species");
            nsp = allsp.size();
            spnames.resize(nsp);
            for (size_t nn = 0; nn < nsp; nn++) {
                string stemp = (*allsp[nn])["name"];
                if (!declared[stemp] || sprule[jsp] < 10) {
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
                if (!declared[stemp]) {
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
                if (!declared[stemp] || sprule[jsp] < 10) {
                    declared[stemp] = true;
                    // Find the species in the database by name.
                    auto iter = speciesNodes.find(stemp);
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

void importPhase(XML_Node& phase, ThermoPhase* th)
{
    // Check the the supplied XML node in fact represents a phase.
    if (phase.name() != "phase") {
        throw CanteraError("importPhase",
                           "Current const XML_Node named, " + phase.name() +
                           ", is not a phase element.");
    }

    // In this section of code, we get the reference to the phase XML tree
    // within the ThermoPhase object. Then, we clear it and fill it with the
    // current information that we are about to use to construct the object. We
    // will then be able to resurrect the information later by calling xml().
    th->setXMLdata(phase);

    // set the id attribute of the phase to the 'id' attribute in the XML tree.
    th->setID(phase.id());
    th->setName(phase.id());

    // Number of spatial dimensions. Defaults to 3 (bulk phase)
    if (phase.hasAttrib("dim")) {
        int idim = intValue(phase["dim"]);
        if (idim < 1 || idim > 3) {
            throw CanteraError("importPhase",
                               "phase, " + th->id() +
                               ", has unphysical number of dimensions: " + phase["dim"]);
        }
        th->setNDim(idim);
    } else {
        th->setNDim(3); // default
    }

    // Set equation of state parameters. The parameters are specific to each
    // subclass of ThermoPhase, so this is done by method setParametersFromXML
    // in each subclass.
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

    // Add the elements.
    if (ssConvention != cSS_CONVENTION_SLAVE) {
        installElements(*th, phase);
    }

    // Add the species.
    //
    // Species definitions may be imported from multiple sources. For each one,
    // a speciesArray element must be present.
    vector<XML_Node*> sparrays = phase.getChildren("speciesArray");
    if (ssConvention != cSS_CONVENTION_SLAVE && sparrays.empty()) {
        throw CanteraError("importPhase",
                           "phase, " + th->id() + ", has zero \"speciesArray\" XML nodes.\n"
                           + " There must be at least one speciesArray nodes "
                           "with one or more species");
    }
    vector<XML_Node*> dbases;
    vector_int sprule(sparrays.size(),0);

    // Default behavior when importing from CTI/XML is for undefined elements to
    // be treated as an error
    th->throwUndefinedElements();

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

        // Get a pointer to the node containing the species definitions for the
        // species declared in this speciesArray element. This may be in the
        // local file containing the phase element, or may be in another file.
        XML_Node* db = get_XML_Node(speciesArray["datasrc"], &phase.root());
        if (db == 0) {
            throw CanteraError("importPhase()",
                               " Can not find XML node for species database: "
                               + speciesArray["datasrc"]);
        }

        // add this node to the list of species database nodes.
        dbases.push_back(db);
    }

    // Now, collect all the species names and all the XML_Node * pointers for
    // those species in a single vector. This is where we decide what species
    // are to be included in the phase. The logic is complicated enough that we
    // put it in a separate routine.
    std::vector<XML_Node*> spDataNodeList;
    std::vector<std::string> spNamesList;
    vector_int spRuleList;
    formSpeciesXMLNodeList(spDataNodeList, spNamesList, spRuleList,
                           sparrays, dbases, sprule);

    size_t nsp = spDataNodeList.size();
    if (ssConvention == cSS_CONVENTION_SLAVE && nsp > 0) {
        throw CanteraError("importPhase()", "For Slave standard states, "
            "number of species must be zero: {}", nsp);
    }
    for (size_t k = 0; k < nsp; k++) {
        XML_Node* s = spDataNodeList[k];
        AssertTrace(s != 0);
        if (spRuleList[k]) {
           th->ignoreUndefinedElements();
        }
        th->addSpecies(newSpecies(*s));
        if (vpss_ptr) {
            const XML_Node* const ss = s->findByName("standardState");
            std::string ss_model = (ss) ? ss->attrib("model") : "ideal-gas";
            unique_ptr<PDSS> kPDSS(newPDSS(ss_model));
            kPDSS->setParametersFromXML(*s);
            vpss_ptr->installPDSS(k, std::move(kPDSS));
        }
        th->saveSpeciesData(k, s);
    }

    // Done adding species. Perform any required subclass-specific
    // initialization.
    th->initThermo();

    // Perform any required subclass-specific initialization that requires the
    // XML phase object
    std::string id = "";
    th->initThermoXML(phase, id);
}

void installElements(Phase& th, const XML_Node& phaseNode)
{
    // get the declared element names
    if (!phaseNode.hasChild("elementArray")) {
        throw CanteraError("installElements",
                           "phase XML node doesn't have \"elementArray\" XML Node");
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
    XML_Node* dbe = &doc->child("elementData");

    XML_Node& root = phaseNode.root();
    XML_Node* local_db = 0;
    if (root.hasChild("elementData")) {
        local_db = &root.child("elementData");
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
        th.addElement(symbol, weight, anum, entropy298);
    }
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
