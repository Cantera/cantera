/**
 *  @file ThermoFactory.cpp
 *     Definitions for the factory class that can create known ThermoPhase objects
 *     (see \ref thermoprops and class \link Cantera::ThermoFactory ThermoFactory\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

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
#include "cantera/thermo/PureFluidPhase.h"
#include "cantera/thermo/RedlichKwongMFTP.h"
#include "cantera/thermo/ConstDensityThermo.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/thermo/EdgePhase.h"
#include "cantera/thermo/MetalPhase.h"
#include "cantera/thermo/StoichSubstance.h"
#include "cantera/thermo/FixedChemPotSSTP.h"
#include "cantera/thermo/LatticeSolidPhase.h"
#include "cantera/thermo/LatticePhase.h"
#include "cantera/thermo/HMWSoln.h"
#include "cantera/thermo/DebyeHuckel.h"
#include "cantera/thermo/IdealMolalSoln.h"
#include "cantera/thermo/IdealSolnGasVPSS.h"
#include "cantera/thermo/WaterSSTP.h"
#include "cantera/thermo/BinarySolutionTabulatedThermo.h"
#include "cantera/base/stringUtils.h"

using namespace std;

namespace Cantera
{

ThermoFactory* ThermoFactory::s_factory = 0;
std::mutex ThermoFactory::thermo_mutex;

ThermoFactory::ThermoFactory()
{
    reg("ideal-gas", []() { return new IdealGasPhase(); });
    addAlias("ideal-gas", "IdealGas");
    reg("constant-density", []() { return new ConstDensityThermo(); });
    addAlias("constant-density", "Incompressible");
    reg("ideal-surface", []() { return new SurfPhase(); });
    addAlias("ideal-surface", "Surface");
    reg("edge", []() { return new EdgePhase(); });
    addAlias("edge", "Edge");
    reg("electron-cloud", []() { return new MetalPhase(); });
    addAlias("electron-cloud", "Metal");
    reg("fixed-stoichiometry", []() { return new StoichSubstance(); });
    addAlias("fixed-stoichiometry", "StoichSubstance");
    reg("pure-fluid", []() { return new PureFluidPhase(); });
    addAlias("pure-fluid", "PureFluid");
    reg("compound-lattice", []() { return new LatticeSolidPhase(); });
    addAlias("compound-lattice", "LatticeSolid");
    reg("lattice", []() { return new LatticePhase(); });
    addAlias("lattice", "Lattice");
    reg("HMW-electrolyte", []() { return new HMWSoln(); });
    addAlias("HMW-electrolyte", "HMW");
    reg("ideal-condensed", []() { return new IdealSolidSolnPhase(); });
    addAlias("ideal-condensed", "IdealSolidSolution");
    reg("Debye-Huckel", []() { return new DebyeHuckel(); });
    addAlias("Debye-Huckel", "DebyeHuckel");
    reg("ideal-molal-solution", []() { return new IdealMolalSoln(); });
    addAlias("ideal-molal-solution", "IdealMolalSolution");
    reg("ideal-solution-VPSS", []() { return new IdealSolnGasVPSS(); });
    reg("ideal-gas-VPSS", []() { return new IdealSolnGasVPSS(); });
    addAlias("ideal-solution-VPSS", "IdealSolnVPSS");
    addAlias("ideal-gas-VPSS", "IdealGasVPSS");
    reg("Margules", []() { return new MargulesVPSSTP(); });
    reg("ions-from-neutral-molecule", []() { return new IonsFromNeutralVPSSTP(); });
    addAlias("ions-from-neutral-molecule", "IonsFromNeutralMolecule");
    reg("fixed-chemical-potential", []() { return new FixedChemPotSSTP(); });
    addAlias("fixed-chemical-potential", "FixedChemPot");
    reg("Redlich-Kister", []() { return new RedlichKisterVPSSTP(); });
    reg("Redlich-Kwong", []() { return new RedlichKwongMFTP(); });
    addAlias("Redlich-Kwong", "RedlichKwongMFTP");
    addAlias("Redlich-Kwong", "RedlichKwong");
    reg("Maskell-solid-solution", []() { return new MaskellSolidSolnPhase(); });
    addAlias("Maskell-solid-solution", "MaskellSolidSolnPhase");
    reg("liquid-water-IAPWS95", []() { return new WaterSSTP(); });
    addAlias("liquid-water-IAPWS95", "PureLiquidWater");
    reg("binary-solution-tabulated", []() { return new BinarySolutionTabulatedThermo(); });
    addAlias("binary-solution-tabulated", "BinarySolutionTabulatedThermo");
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

unique_ptr<ThermoPhase> newPhase(AnyMap& phaseNode, const AnyMap& rootNode)
{
    unique_ptr<ThermoPhase> t(newThermoPhase(phaseNode["thermo"].asString()));
    setupPhase(*t, phaseNode, rootNode);
    return t;
}

ThermoPhase* newPhase(const std::string& infile, std::string id)
{
    size_t dot = infile.find_last_of(".");
    string extension;
    if (dot != npos) {
        extension = toLowerCopy(infile.substr(dot+1));
    }
    if (id == "-") {
        id = "";
    }

    if (extension == "yml" || extension == "yaml") {
        AnyMap root = AnyMap::fromYamlFile(infile);
        AnyMap& phase = root["phases"].getMapWhere("name", id);
        unique_ptr<ThermoPhase> t(newThermoPhase(phase["thermo"].asString()));
        setupPhase(*t, phase, root);
        return t.release();
    } else {
        XML_Node* root = get_XML_File(infile);
        XML_Node* xphase = get_XML_NameID("phase", "#"+id, root);
        if (!xphase) {
            throw CanteraError("newPhase",
                               "Couldn't find phase named \"" + id + "\" in file, " + infile);
        }
        return newPhase(*xphase);
    }
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
                        throw CanteraError("formSpeciesXMLNodeList",
                            "no data for species, \"{}\"", stemp);
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
    th->setName(phase.id());

    // Number of spatial dimensions. Defaults to 3 (bulk phase)
    if (phase.hasAttrib("dim")) {
        int idim = intValue(phase["dim"]);
        if (idim < 1 || idim > 3) {
            throw CanteraError("importPhase",
                               "phase, " + th->name() +
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
                           " phase, " + th->name() +
                           ", XML_Node does not have a \"thermo\" XML_Node");
    }

    VPStandardStateTP* vpss_ptr = 0;
    int ssConvention = th->standardStateConvention();
    if (ssConvention == cSS_CONVENTION_VPSS) {
        vpss_ptr = dynamic_cast <VPStandardStateTP*>(th);
        if (vpss_ptr == 0) {
            throw CanteraError("importPhase",
                               "phase, " + th->name() + ", was VPSS, but dynamic cast failed");
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
                           "phase, " + th->name() + ", has zero \"speciesArray\" XML nodes.\n"
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
            throw CanteraError("importPhase",
                               "Can not find XML node for species database: {}",
                               speciesArray["datasrc"]);
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
        throw CanteraError("importPhase", "For Slave standard states, "
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

void addDefaultElements(ThermoPhase& thermo, const vector<string>& element_names) {
    for (const auto& symbol : element_names) {
        thermo.addElement(symbol);
    }
}

void addElements(ThermoPhase& thermo, const vector<string>& element_names,
                 const AnyValue& elements, bool allow_default)
{
    const auto& local_elements = elements.asMap("symbol");
    for (const auto& symbol : element_names) {
        if (local_elements.count(symbol)) {
            auto& element = *local_elements.at(symbol);
            double weight = element["atomic-weight"].asDouble();
            long int number = element.getInt("atomic-number", 0);
            double e298 = element.getDouble("entropy298", ENTROPY298_UNKNOWN);
            thermo.addElement(symbol, weight, number, e298);
        } else if (allow_default) {
            thermo.addElement(symbol);
        } else {
            throw InputFileError("addElements", elements,
                                 "Element '{}' not found", symbol);
        }
    }
}

void addSpecies(ThermoPhase& thermo, const AnyValue& names, const AnyValue& species)
{
    if (names.is<vector<string>>()) {
        // 'names' is a list of species names which should be found in 'species'
        const auto& species_nodes = species.asMap("name");
        for (const auto& name : names.asVector<string>()) {
            if (species_nodes.count(name)) {
                thermo.addSpecies(newSpecies(*species_nodes.at(name)));
            } else {
                throw InputFileError("addSpecies", names, species,
                    "Could not find a species named '{}'.", name);
            }
        }
    } else if (names == "all") {
        // The keyword 'all' means to add all species from this source
        for (const auto& item : species.asVector<AnyMap>()) {
            thermo.addSpecies(newSpecies(item));
        }
    } else {
        throw InputFileError("addSpecies", names,
            "Could not parse species declaration of type '{}'", names.type_str());
    }
}

void setupPhase(ThermoPhase& thermo, AnyMap& phaseNode, const AnyMap& rootNode)
{
    thermo.setName(phaseNode["name"].asString());
    if (rootNode.hasKey("__file__")) {
        phaseNode["__file__"] = rootNode["__file__"];
    }

    if (phaseNode.hasKey("deprecated")) {
        string msg = phaseNode["deprecated"].asString();
        string filename = phaseNode.getString("__file__", "unknown file");
        string method = fmt::format("{}/{}", filename, phaseNode["name"].asString());
        warn_deprecated(method, msg);
    }

    // Add elements
    if (phaseNode.hasKey("elements")) {
        if (phaseNode.getBool("skip-undeclared-elements", false)) {
            thermo.ignoreUndefinedElements();
        } else {
            thermo.throwUndefinedElements();
        }

        if (phaseNode["elements"].is<vector<string>>()) {
            // 'elements' is a list of element symbols
            if (rootNode.hasKey("elements")) {
                addElements(thermo, phaseNode["elements"].asVector<string>(),
                            rootNode["elements"], true);
            } else {
                addDefaultElements(thermo, phaseNode["elements"].asVector<string>());
            }
        } else if (phaseNode["elements"].is<vector<AnyMap>>()) {
            // Each item in 'elements' is a map with one item, where the key is
            // a section in this file or another YAML file, and the value is a
            // list of element symbols to read from that section
            for (const auto& elemNode : phaseNode["elements"].asVector<AnyMap>()) {
                const string& source = elemNode.begin()->first;
                const auto& names = elemNode.begin()->second.asVector<string>();
                const auto& slash = boost::ifind_last(source, "/");
                if (slash) {
                    std::string fileName(source.begin(), slash.begin());
                    std::string node(slash.end(), source.end());
                    const AnyMap elements = AnyMap::fromYamlFile(fileName,
                        rootNode.getString("__file__", ""));
                    addElements(thermo, names, elements.at(node), false);
                } else if (rootNode.hasKey(source)) {
                    addElements(thermo, names, rootNode.at(source), false);
                } else if (source == "default") {
                    addDefaultElements(thermo, names);
                } else {
                    throw InputFileError("setupPhase", elemNode,
                        "Could not find elements section named '{}'", source);
                }
            }
        } else {
            throw InputFileError("setupPhase", phaseNode["elements"],
                "Could not parse elements declaration of type '{}'",
                phaseNode["elements"].type_str());
        }
    } else {
        // If no elements list is provided, just add elements as-needed from the
        // default list.
        thermo.addUndefinedElements();
    }

    // Add species
    if (phaseNode.hasKey("species")) {
        if (phaseNode["species"].is<vector<string>>()) {
            // 'species' is a list of species names to be added from the current
            // file's 'species' section
            addSpecies(thermo, phaseNode["species"], rootNode["species"]);
        } else if (phaseNode["species"].is<string>()) {
            // 'species' is a keyword applicable to the current file's 'species'
            // section
            addSpecies(thermo, phaseNode["species"], rootNode["species"]);
        } else if (phaseNode["species"].is<vector<AnyMap>>()) {
            // Each item in 'species' is a map with one item, where the key is
            // a section in this file or another YAML file, and the value is a
            // list of species names to read from that section
            for (const auto& speciesNode : phaseNode["species"].asVector<AnyMap>()) {
                const string& source = speciesNode.begin()->first;
                const auto& names = speciesNode.begin()->second;
                const auto& slash = boost::ifind_last(source, "/");
                if (slash) {
                    // source is a different input file
                    std::string fileName(source.begin(), slash.begin());
                    std::string node(slash.end(), source.end());
                    AnyMap species = AnyMap::fromYamlFile(fileName,
                        rootNode.getString("__file__", ""));
                    addSpecies(thermo, names, species[node]);
                } else if (rootNode.hasKey(source)) {
                    // source is in the current file
                    addSpecies(thermo, names, rootNode[source]);
                } else {
                    throw InputFileError("setupPhase", speciesNode,
                        "Could not find species section named '{}'", source);
                }
            }
        } else {
            throw InputFileError("setupPhase", phaseNode["species"],
                "Could not parse species declaration of type '{}'",
                phaseNode["species"].type_str());
        }
    } else if (rootNode.hasKey("species")) {
        // By default, add all species from the 'species' section
        addSpecies(thermo, AnyValue("all"), rootNode["species"]);
    }

    auto* vpssThermo = dynamic_cast<VPStandardStateTP*>(&thermo);
    if (vpssThermo) {
        for (size_t k = 0; k < thermo.nSpecies(); k++) {
            unique_ptr<PDSS> pdss;
            if (thermo.species(k)->input.hasKey("equation-of-state")) {
                // Use the first node which specifies a valid PDSS model
                auto& eos = thermo.species(k)->input["equation-of-state"];
                bool found = false;
                for (auto& node : eos.asVector<AnyMap>()) {
                    string model = node["model"].asString();
                    if (PDSSFactory::factory()->exists(model)) {
                        pdss.reset(newPDSS(model));
                        pdss->setParameters(node);
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    throw InputFileError("setupPhase", eos,
                        "Could not find an equation-of-state specification "
                        "which defines a known PDSS model.");
                }
            } else {
                pdss.reset(newPDSS("ideal-gas"));
            }
            vpssThermo->installPDSS(k, std::move(pdss));
        }
    }

    thermo.setParameters(phaseNode, rootNode);
    thermo.initThermo();

    if (phaseNode.hasKey("state")) {
        auto node = phaseNode["state"].as<AnyMap>();
        thermo.setState(node);
    } else {
        thermo.setState_TP(298.15, OneAtm);
    }
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
            throw CanteraError("installElements", "no data for element '{}'",
                               enames[i]);
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
        throw CanteraError("speciesXML_Node",
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
